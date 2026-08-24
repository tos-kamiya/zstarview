import math
from functools import lru_cache

import numpy as np
from astropy.time import Time
from PySide6.QtGui import QImage

from ..night_lights import sky_disc_ambient_sun_altitude_factor
from ..types import ScreenGeometry, ViewProjection
from .atmosphere import AEROSOL_REFERENCE_AOD550, atmospheric_sky_samples
from .qt_image import np_rgba_to_qimage

FLAT_SKY_DISC_RGB_U8 = np.array([10, 10, 10], dtype=np.uint8)
SKY_AMBIENT_RGB_U8 = np.array([0.5, 1.0, 2.5], dtype=np.float32)
_SKY_DISC_RENDER_CACHE_SIZE = 2
SKY_DISC_OVERSCAN_DEG = 0.75
SKY_DISC_RENDER_SCALE = 0.125
SKY_DISC_REFERENCE_WIDTH_PX = 1920
SKY_DISC_TWILIGHT_UPDATE_INTERVAL_SECONDS = 15
SKY_DISC_DEFAULT_UPDATE_INTERVAL_SECONDS = 60
SKY_DISC_TWILIGHT_MIN_SUN_ALT_DEG = -15.0
SKY_DISC_TWILIGHT_MAX_SUN_ALT_DEG = 15.0
SOLAR_HORIZON_COLOR_MIN_ALT_DEG = 0.0
SOLAR_HORIZON_COLOR_MAX_ALT_DEG = 10.0
SOLAR_HORIZON_COLOR_SAMPLES = 6
SKY_INTENSITY_MAX_LUMINANCE_BLEND = 0.65
SKY_INTENSITY_LUMINANCE_BLEND_GAMMA = 1.5
SKY_LUMINANCE_WEIGHTS = np.array([0.2126, 0.7152, 0.0722], dtype=np.float32)


def sky_disc_update_interval(sun_alt_deg: float | None) -> int:
    """Return the fixed sky-disc refresh interval for the current Sun altitude."""
    if (
        sun_alt_deg is not None
        and math.isfinite(float(sun_alt_deg))
        and SKY_DISC_TWILIGHT_MIN_SUN_ALT_DEG
        <= float(sun_alt_deg)
        <= SKY_DISC_TWILIGHT_MAX_SUN_ALT_DEG
    ):
        return SKY_DISC_TWILIGHT_UPDATE_INTERVAL_SECONDS
    return SKY_DISC_DEFAULT_UPDATE_INTERVAL_SECONDS


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """Perform smooth Hermite interpolation between 0 and 1."""
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return float(t * t * (3.0 - 2.0 * t))


def _apply_sky_intensity(colors: np.ndarray, intensity: float) -> np.ndarray:
    """Raise sky luminance faster than chroma at stronger intensities."""
    strength = float(np.clip(intensity, 0.0, 1.0))
    scaled = np.asarray(colors, dtype=np.float32) * strength
    luminance_blend = (
        SKY_INTENSITY_MAX_LUMINANCE_BLEND
        * strength**SKY_INTENSITY_LUMINANCE_BLEND_GAMMA
    )
    luminance = scaled @ SKY_LUMINANCE_WEIGHTS
    return scaled + luminance_blend * (luminance[:, np.newaxis] - scaled)


def _inverse_project_disc(
    width_px: int,
    height_px: int,
    geometry: ScreenGeometry,
    projection: ViewProjection,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project pixels up to the requested content FOV."""
    radius_px = max(1.0, float(geometry.radius))
    ys = (np.arange(height_px, dtype=np.float32) - float(geometry.center[1])) / radius_px
    xs = (np.arange(width_px, dtype=np.float32) - float(geometry.center[0])) / radius_px
    nx, ny = np.meshgrid(xs, ys)

    rr2 = nx * nx + ny * ny
    edge_fov = max(1.0e-6, float(projection.edge_fov_deg))
    max_r = max(0.0, float(projection.content_fov_deg) / edge_fov)
    inside = rr2 <= (max_r * max_r)
    if not np.any(inside):
        return np.array([], dtype=np.float32), np.array([], dtype=np.float32), inside

    r = np.sqrt(rr2[inside]).astype(np.float32)
    theta = np.radians(r * edge_fov)
    psi = np.arctan2(nx[inside], -ny[inside])

    alt_c, az_c = projection.view_center
    eps = 1e-3
    phi1 = np.float32(math.radians(np.clip(alt_c, -90.0 + eps, 90.0 - eps)))
    lam1 = np.float32(math.radians(az_c))
    sin_phi1 = np.sin(phi1)
    cos_phi1 = np.cos(phi1)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    sin_phi2 = sin_phi1 * cos_theta + cos_phi1 * sin_theta * np.cos(psi)
    sin_phi2 = np.clip(sin_phi2, -1.0, 1.0)
    phi2 = np.arcsin(sin_phi2)
    y = np.sin(psi) * sin_theta * cos_phi1
    x = cos_theta - sin_phi1 * sin_phi2
    lam2 = lam1 + np.arctan2(y, x)

    alt = np.degrees(phi2).astype(np.float32)
    az = (np.degrees(lam2) + 360.0) % 360.0
    return alt, az.astype(np.float32), inside


def sky_color_samples(
    view_alt_deg: np.ndarray,
    view_az_deg: np.ndarray,
    sun_altaz: tuple[float, float],
    *,
    exposure: float = 1.0,
    alpha: float = 1.0,
    eclipse_factor: float = 1.0,
    observer_height_m: float = 0.0,
    ambient_scale: float = 1.0,
    aerosol_optical_depth: float | None = None,
) -> np.ndarray:
    """Return Mie/Rayleigh sky RGB samples for local alt/az directions."""
    alt = np.asarray(view_alt_deg, dtype=np.float32)
    az = np.asarray(view_az_deg, dtype=np.float32)
    if alt.shape != az.shape:
        raise ValueError("view_alt_deg and view_az_deg must have the same shape")
    if alt.size == 0:
        return np.zeros((0, 3), dtype=np.float32)

    colors = atmospheric_sky_samples(
        alt.reshape(-1),
        az.reshape(-1),
        sun_altaz,
        observer_height_m=observer_height_m,
        exposure=exposure,
        aerosol_optical_depth=(
            AEROSOL_REFERENCE_AOD550
            if aerosol_optical_depth is None
            else aerosol_optical_depth
        ),
    )
    colors = _apply_sky_intensity(colors, alpha)
    colors *= max(0.0, float(eclipse_factor))
    colors += (
        SKY_AMBIENT_RGB_U8.astype(np.float32)
        * max(0.0, float(ambient_scale))
        / 255.0
    )
    return np.clip(colors, 0.0, 1.0).astype(np.float32)


def sky_color_at_direction(
    view_alt_deg: float,
    view_az_deg: float,
    sun_altaz: tuple[float, float],
    *,
    alpha: float = 1.0,
    exposure: float = 1.0,
    eclipse_factor: float = 1.0,
    observer_height_m: float = 0.0,
) -> tuple[int, int, int, int]:
    """Return the Mie/Rayleigh sky color at one direction."""
    colors = sky_color_samples(
        np.asarray([view_alt_deg], dtype=np.float32),
        np.asarray([view_az_deg], dtype=np.float32),
        sun_altaz,
        alpha=alpha,
        exposure=exposure,
        eclipse_factor=eclipse_factor,
        observer_height_m=observer_height_m,
    )
    red, green, blue = np.clip(np.round(colors[0] * 255.0), 0, 255).astype(int)
    return int(red), int(green), int(blue), 255


def sky_color_near_solar_horizon(
    sun_altaz: tuple[float, float],
    *,
    alpha: float = 1.0,
    exposure: float = 1.0,
    eclipse_factor: float = 1.0,
    observer_height_m: float = 0.0,
    aerosol_optical_depth: float | None = None,
) -> tuple[int, int, int, int]:
    """Return the mean sky color from 0 to 10 degrees at solar azimuth."""
    altitudes = np.linspace(
        SOLAR_HORIZON_COLOR_MIN_ALT_DEG,
        SOLAR_HORIZON_COLOR_MAX_ALT_DEG,
        SOLAR_HORIZON_COLOR_SAMPLES,
        dtype=np.float32,
    )
    colors = sky_color_samples(
        altitudes,
        np.full_like(altitudes, float(sun_altaz[1])),
        sun_altaz,
        alpha=alpha,
        exposure=exposure,
        eclipse_factor=eclipse_factor,
        observer_height_m=observer_height_m,
        aerosol_optical_depth=aerosol_optical_depth,
    )
    red, green, blue = np.clip(np.round(colors.mean(axis=0) * 255.0), 0, 255).astype(int)
    return int(red), int(green), int(blue), 255


@lru_cache(maxsize=_SKY_DISC_RENDER_CACHE_SIZE)
def _render_sky_color_disc_cached(
    width: int,
    height: int,
    center_x: int,
    center_y: int,
    radius: int,
    projection: ViewProjection,
    sun_alt_deg: float,
    sun_az_deg: float,
    exposure: float,
    alpha: float,
    disc_opacity: float,
    eclipse_factor: float,
    observer_height_m: float,
    ambient_scale: float,
    aerosol_optical_depth: float,
) -> QImage:
    local_geometry = ScreenGeometry(center=(center_x, center_y), radius=radius)
    alt, az, inside = _inverse_project_disc(width, height, local_geometry, projection)
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    if alt.size == 0:
        return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)

    colors = sky_color_samples(
        alt,
        az,
        (sun_alt_deg, sun_az_deg),
        exposure=exposure,
        alpha=alpha,
        eclipse_factor=eclipse_factor,
        observer_height_m=observer_height_m,
        ambient_scale=ambient_scale,
        aerosol_optical_depth=aerosol_optical_depth,
    )
    rgb_u8 = np.clip(np.round(colors * 255.0), 0, 255).astype(np.uint8)
    alpha_u8 = round(max(0.0, min(1.0, float(disc_opacity))) * 255.0)
    rgba[..., 3][inside] = alpha_u8
    rgba[..., :3][inside] = rgb_u8
    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)


def draw_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    sun_altaz: tuple[float, float],
    *,
    exposure: float = 1.0,
    alpha: float = 1.0,
    disc_opacity: float = 1.0,
    eclipse_factor: float = 1.0,
    observer_height_m: float = 0.0,
    time_obj: Time | None = None,
    timezone_name: str = "UTC",
    image_size: tuple[int, int] | None = None,
    aerosol_optical_depth: float | None = None,
) -> QImage:
    """Draw a Mie/Rayleigh sky disc using one NumPy inverse projection."""
    radius = int(geometry.radius)
    if image_size is None:
        width = height = max(2, radius * 2)
        local_geometry = ScreenGeometry(center=(radius, radius), radius=radius)
    else:
        width = max(2, int(image_size[0]))
        height = max(2, int(image_size[1]))
        local_geometry = geometry
    if radius < 1:
        return QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
    projection = ViewProjection(
        view_center=(float(view_center[0]), float(view_center[1])),
        edge_fov_deg=float(edge_fov_deg),
        content_fov_deg=float(content_fov_deg),
    )
    ambient_scale = 1.0 + (
        1.0 - sky_disc_ambient_sun_altitude_factor(float(sun_altaz[0]))
    )
    return _render_sky_color_disc_cached(
        width,
        height,
        int(local_geometry.center[0]),
        int(local_geometry.center[1]),
        int(local_geometry.radius),
        projection,
        float(sun_altaz[0]),
        float(sun_altaz[1]),
        float(exposure),
        float(alpha),
        float(disc_opacity),
        float(eclipse_factor),
        float(observer_height_m),
        ambient_scale,
        AEROSOL_REFERENCE_AOD550
        if aerosol_optical_depth is None
        else float(aerosol_optical_depth),
    )


def draw_uniform_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    *,
    image_size: tuple[int, int] | None = None,
    disc_opacity: float = 1.0,
) -> QImage:
    """Draw a flat disc used when sky-color shading is disabled."""
    radius = int(geometry.radius)
    if image_size is None:
        width = height = max(2, radius * 2)
        local_geometry = ScreenGeometry(center=(radius, radius), radius=radius)
    else:
        width = max(2, int(image_size[0]))
        height = max(2, int(image_size[1]))
        local_geometry = geometry
    if radius < 1:
        return QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)

    _, _, inside = _inverse_project_disc(
        width,
        height,
        local_geometry,
        ViewProjection(
            view_center=(float(view_center[0]), float(view_center[1])),
            edge_fov_deg=float(edge_fov_deg),
            content_fov_deg=float(content_fov_deg),
        ),
    )
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    alpha_u8 = round(max(0.0, min(1.0, float(disc_opacity))) * 255.0)
    rgba[..., 3][inside] = alpha_u8
    rgba[..., :3][inside] = FLAT_SKY_DISC_RGB_U8
    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
