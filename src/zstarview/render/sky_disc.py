import math
from functools import lru_cache
from typing import Tuple

import numpy as np
from PySide6.QtGui import QImage

from ..paths import PALETTE_NEVER_RISES_RGB
from ..types import ScreenGeometry
from .qt_image import np_rgba_to_qimage

TURBIDITY = 5  # 2 (clear blue sky) to 10 (hazy white sky)
FLAT_SKY_DISC_RGB_U8 = np.array([10, 10, 10], dtype=np.uint8)

NIGHT_SKY_RGB = np.array([0.01, 0.02, 0.05], dtype=np.float32)
HORIZON_DAY_RGB = np.array([0.96, 0.73, 0.50], dtype=np.float32)
ZENITH_DAY_RGB = np.array([0.24, 0.48, 0.86], dtype=np.float32)
HAZE_RGB = np.array([0.99, 0.96, 0.92], dtype=np.float32)
RAYLEIGH_BLUE_RGB = np.array([0.18, 0.34, 0.82], dtype=np.float32)
SUN_GLOW_RGB = np.array([1.00, 0.89, 0.70], dtype=np.float32)
SUNSET_RGB = np.array([1.00, 0.54, 0.20], dtype=np.float32)
ANTI_SOLAR_RGB = np.array([0.14, 0.18, 0.34], dtype=np.float32)

SUNSET_START_ALT_DEG = 4.0
SUNSET_END_ALT_DEG = 18.0
SUN_GLOW_EXPONENT_BASE = 1.75
ANTI_SOLAR_EXPONENT = 2.6
HAZE_STRENGTH_MIN = 0.16
HAZE_STRENGTH_MAX = 0.68
RAYLEIGH_STRENGTH = 0.58
SUN_GLOW_STRENGTH = 0.42
SUNSET_STRENGTH = 0.264
ANTI_SOLAR_STRENGTH = 0.16
SATURATION_CHROMA_SCALE = 0.35
# Periodic sky updates include continuously changing sun coordinates, so this
# cache mostly protects immediate duplicate requests. Keep it small because
# each entry can be a full-window QImage backed by Qt memory.
_SKY_DISC_RENDER_CACHE_SIZE = 2


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """Performs a smooth Hermite interpolation between 0 and 1."""
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return float(t * t * (3.0 - 2.0 * t))


def _inverse_project_disc(
    width_px: int,
    height_px: int,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    *,
    edge_fov_deg: float,
    content_fov_deg: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project pixels up to the requested content FOV."""
    radius_px = max(1.0, float(geometry.radius))
    ys = (np.arange(height_px, dtype=np.float32) - float(geometry.center[1])) / radius_px
    xs = (np.arange(width_px, dtype=np.float32) - float(geometry.center[0])) / radius_px
    nx, ny = np.meshgrid(xs, ys)

    rr2 = nx * nx + ny * ny
    edge_fov = max(1.0e-6, float(edge_fov_deg))
    max_r = max(0.0, float(content_fov_deg) / edge_fov)
    inside = rr2 <= (max_r * max_r)
    if not np.any(inside):
        return np.array([], dtype=np.float32), np.array([], dtype=np.float32), inside

    r = np.sqrt(rr2[inside]).astype(np.float32)
    theta = np.radians(r * edge_fov)

    # Bearing from local north (clockwise): north=(0,-1), east=(1,0).
    psi = np.arctan2(nx[inside], -ny[inside])

    alt_c, az_c = view_center
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


def _get_sky_color_vectorized(
    view_alt_deg: np.ndarray,
    view_az_deg: np.ndarray,
    sun_altaz: Tuple[float, float],
    *,
    saturation: float,
) -> np.ndarray:
    """Vectorized sky color model for many sky directions at once."""
    sun_alt_deg, sun_az_deg = sun_altaz
    n = view_alt_deg.shape[0]
    if n == 0:
        return np.zeros((0, 3), dtype=np.float32)
    if sun_alt_deg <= -10.0:
        return np.repeat(NIGHT_SKY_RGB[None, :], n, axis=0)

    tau = float(np.clip((TURBIDITY - 2.0) / 8.0, 0.0, 1.0))
    colorfulness = float(np.clip(1.0 + (float(saturation) - 1.0) * SATURATION_CHROMA_SCALE, 0.75, 1.25))
    t_alt = np.clip(view_alt_deg / 90.0, 0.0, 1.0)
    sun_up = _smoothstep(-8.0, 6.0, sun_alt_deg)
    twilight = _smoothstep(-10.0, 0.0, sun_alt_deg)
    sunset = 1.0 - _smoothstep(SUNSET_START_ALT_DEG, SUNSET_END_ALT_DEG, sun_alt_deg)
    low_altitude = 1.0 - t_alt
    high_altitude = t_alt

    horizon_day = HORIZON_DAY_RGB[None, :]
    zenith_day = ZENITH_DAY_RGB[None, :]
    base = horizon_day + (zenith_day - horizon_day) * t_alt[:, None]

    haze_strength = (HAZE_STRENGTH_MIN + HAZE_STRENGTH_MAX * low_altitude) * (0.72 + 0.48 * tau)
    haze_strength = np.clip(haze_strength * (0.88 + 0.12 * sunset), 0.0, 1.0)
    haze_color = HAZE_RGB[None, :]
    base = base + (haze_color - base) * haze_strength[:, None]

    a1 = np.radians(view_alt_deg)
    z1 = np.radians(view_az_deg)
    a2 = math.radians(sun_alt_deg)
    z2 = math.radians(sun_az_deg)
    cos_g = np.sin(a1) * math.sin(a2) + np.cos(a1) * math.cos(a2) * np.cos(z2 - z1)
    cos_g = np.clip(cos_g, -1.0, 1.0)

    forward = np.maximum(0.0, cos_g)
    back = np.maximum(0.0, -cos_g)

    rayleigh_angle = 1.0 - cos_g * cos_g
    rayleigh_spread = 0.55 + 0.75 * sunset
    rayleigh_amount = rayleigh_angle * rayleigh_spread * sun_up * (0.34 + 0.66 * high_altitude)
    rayleigh_amount *= 0.78 + 0.22 * tau
    rayleigh_strength = np.clip(RAYLEIGH_STRENGTH * rayleigh_amount * colorfulness, 0.0, 1.0)
    color = base + (RAYLEIGH_BLUE_RGB[None, :] - base) * rayleigh_strength[:, None]

    anti_amount = (back**ANTI_SOLAR_EXPONENT) * sun_up * (0.25 + 0.75 * high_altitude)
    anti_strength = np.clip(ANTI_SOLAR_STRENGTH * anti_amount * (0.85 + 0.15 * colorfulness), 0.0, 1.0)
    color = color + (ANTI_SOLAR_RGB[None, :] - color) * anti_strength[:, None]

    glow_exponent = SUN_GLOW_EXPONENT_BASE - 0.25 * tau
    sun_glow_amount = (forward**glow_exponent) * sun_up
    sun_glow_strength = np.clip(SUN_GLOW_STRENGTH * sun_glow_amount * (0.92 + 0.08 * colorfulness), 0.0, 1.0)
    color = color + (SUN_GLOW_RGB[None, :] - color) * sun_glow_strength[:, None]

    sunset_amount = sunset * low_altitude * (forward ** (1.20 - 0.20 * tau))
    sunset_amount *= 0.70 + 0.30 * sun_up
    sunset_strength = np.clip(SUNSET_STRENGTH * sunset_amount * (0.92 + 0.08 * colorfulness), 0.0, 1.0)
    color = color + (SUNSET_RGB[None, :] - color) * sunset_strength[:, None]

    color = NIGHT_SKY_RGB[None, :] + (color - NIGHT_SKY_RGB[None, :]) * twilight

    return np.clip(color, 0.0, 1.0).astype(np.float32)


def sky_color_samples(
    view_alt_deg: np.ndarray,
    view_az_deg: np.ndarray,
    sun_altaz: Tuple[float, float],
    *,
    exposure: float = 1.14,
    saturation: float = 1.35,
    alpha: float = 1.0,
    eclipse_factor: float = 1.0,
) -> np.ndarray:
    """Return sky RGB samples for (alt, az) arrays."""
    alt = np.asarray(view_alt_deg, dtype=np.float32)
    az = np.asarray(view_az_deg, dtype=np.float32)
    if alt.shape != az.shape:
        raise ValueError("view_alt_deg and view_az_deg must have the same shape")
    if alt.size == 0:
        return np.zeros((0, 3), dtype=np.float32)

    colors = _get_sky_color_vectorized(alt.reshape(-1), az.reshape(-1), sun_altaz, saturation=saturation)
    colors *= max(0.0, float(exposure))
    colors *= max(0.0, float(alpha)) * max(0.0, float(eclipse_factor))
    return np.clip(colors, 0.0, 1.0).astype(np.float32)


@lru_cache(maxsize=_SKY_DISC_RENDER_CACHE_SIZE)
def _render_sky_color_disc_cached(
    width: int,
    height: int,
    center_x: int,
    center_y: int,
    radius: int,
    view_alt_deg: float,
    view_az_deg: float,
    sun_alt_deg: float,
    sun_az_deg: float,
    exposure: float,
    saturation: float,
    alpha: float,
    disc_opacity: float,
    eclipse_factor: float,
    edge_fov_deg: float,
    content_fov_deg: float,
) -> QImage:
    local_geometry = ScreenGeometry(center=(center_x, center_y), radius=radius)
    alt, az, inside = _inverse_project_disc(
        width,
        height,
        local_geometry,
        (view_alt_deg, view_az_deg),
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )

    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    if alt.size == 0:
        return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)

    colors = sky_color_samples(
        alt,
        az,
        (sun_alt_deg, sun_az_deg),
        exposure=exposure,
        saturation=saturation,
        alpha=alpha,
        eclipse_factor=eclipse_factor,
    )
    colors = np.clip(colors, 0.0, 1.0)

    rgb_u8 = np.clip(np.round(colors * 255.0), 0, 255).astype(np.uint8)
    alpha_u8 = int(round(max(0.0, min(1.0, float(disc_opacity))) * 255.0))
    rgba[..., 3][inside] = alpha_u8
    rgba[..., :3][inside] = rgb_u8
    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)


def draw_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    sun_altaz: Tuple[float, float],
    *,
    exposure: float = 1.14,
    saturation: float = 1.35,
    alpha: float = 1.0,
    disc_opacity: float = 1.0,
    eclipse_factor: float = 1.0,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float,
    image_size: Tuple[int, int] | None = None,
) -> QImage:
    """
    Draw sky color disc using one-pass NumPy inverse projection.

    The disc is always rendered as a full circle (no horizon clipping inside the disc).
    """
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
    return _render_sky_color_disc_cached(
        width,
        height,
        int(local_geometry.center[0]),
        int(local_geometry.center[1]),
        int(local_geometry.radius),
        float(view_center[0]),
        float(view_center[1]),
        float(sun_altaz[0]),
        float(sun_altaz[1]),
        float(exposure),
        float(saturation),
        float(alpha),
        float(disc_opacity),
        float(eclipse_factor),
        float(edge_fov_deg),
        float(content_fov_deg),
    )


def draw_uniform_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    *,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float,
    image_size: Tuple[int, int] | None = None,
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
        view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    alpha_u8 = int(round(max(0.0, min(1.0, float(disc_opacity))) * 255.0))
    rgba[..., 3][inside] = alpha_u8
    rgba[..., :3][inside] = FLAT_SKY_DISC_RGB_U8
    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
