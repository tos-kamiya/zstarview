import math
from functools import lru_cache
from typing import Tuple

import numpy as np
from PySide6.QtGui import QImage

from ..types import ScreenGeometry, ViewProjection
from .qt_image import np_rgba_to_qimage

TURBIDITY = 6  # 2 (clear blue sky) to 10 (hazy white sky)
FLAT_SKY_DISC_RGB_U8 = np.array([10, 10, 10], dtype=np.uint8)

NIGHT_SKY_RGB = np.array([0.012, 0.024, 0.06], dtype=np.float32)
HORIZON_DAY_RGB = np.array([103, 171, 245], dtype=np.float32) / 255.0
# Daytime zenith color, specified as RGB(32, 136, 232).
ZENITH_DAY_RGB = np.array([32, 136, 232], dtype=np.float32) / 255.0
# Broad daytime blue-dome target color: RGB(132, 162, 219).
BLUE_DOME_RGB = np.array([132, 162, 219], dtype=np.float32) / 255.0
LOW_ALTITUDE_SKY_RGB = np.array([0.68, 0.75, 0.78], dtype=np.float32)
# A very weak, sun-independent atmospheric horizon tint.  This represents the
# warmer component that can remain at low altitude outside of sunset colors.
HORIZON_ATMOSPHERIC_WARM_RGB = np.array([0.78, 0.68, 0.62], dtype=np.float32)
LOW_HORIZON_WARM_RGB = np.array([0.80, 0.62, 0.50], dtype=np.float32)
SUN_GLOW_RGB = np.array([0.99, 0.98, 0.96], dtype=np.float32)
# A restrained yellow tint for the solar glow once the Sun is clearly above
# the horizon.  This is applied only through the angular solar-glow mask.
SUN_HIGH_ALT_GLOW_RGB = np.array([1.00, 0.94, 0.78], dtype=np.float32)
SUNSET_RGB = np.array([1.00, 0.40, 0.10], dtype=np.float32)
ANTI_SOLAR_RGB = np.array([0.14, 0.18, 0.34], dtype=np.float32)

SUNSET_START_ALT_DEG = -1.0
SUNSET_END_ALT_DEG = 4.0
SUNLIGHT_FLOOR_ALT_DEG = -12.0
SUN_ALT_BLUE_START_DEG = 0.0
SUN_ALT_BLUE_END_DEG = 45.0
SUN_GLOW_EXPONENT_BASE = 1.75
ANTI_SOLAR_EXPONENT = 2.6
LOW_ALTITUDE_WHITENING_STRENGTH = 0.55
LOW_ALTITUDE_WHITENING_EXPONENT = 2.0
HORIZON_ATMOSPHERIC_WARM_ALT_DEG = 8.0
HORIZON_ATMOSPHERIC_WARM_STRENGTH = 0.03
LOW_HORIZON_WARM_ALT_DEG = 4.0
LOW_HORIZON_WARM_ALT_EXPANSION_DEG = 2.0
LOW_HORIZON_WARM_STRENGTH = 0.10
LOW_HORIZON_WARM_MAX_STRENGTH_SCALE = 1.50
SUN_ALT_BLUE_STRENGTH = 0.05
# Temporary visual tuning value for the Rayleigh-disabled comparison.
SUN_GLOW_STRENGTH = 0.30
# The existing sunset layer remains separate; this only warms the solar glow
# itself as the Sun approaches the horizon.
SUN_GLOW_SUNSET_COLOR_MIX = 0.15
SUN_GLOW_HIGH_ALT_COLOR_MIX = 0.10
SUN_GLOW_HIGH_ALT_START_DEG = 4.0
SUN_GLOW_HIGH_ALT_END_DEG = 12.0
# Keep most of the solar glow additive while allowing part of the target color
# to replace the pre-glow color when a strong sky opacity makes that color too
# dominant over sunset hues.
SUN_GLOW_ADDITIVE_MIX = 0.70
SUN_GLOW_REPLACE_MIX = 0.30
SUNSET_STRENGTH = 0.55
# Keep the sunset tint at full strength in the exact solar direction.
SUNSET_SOLAR_GLARE_FACTOR = 1.0
ANTI_SOLAR_STRENGTH = 0.064
SATURATION_CHROMA_SCALE = 0.35
# Periodic sky updates include continuously changing sun coordinates, so this
# cache mostly protects immediate duplicate requests. Keep it small because
# each entry can be a full-window QImage backed by Qt memory.
_SKY_DISC_RENDER_CACHE_SIZE = 2
# Render masks and overlays through the same small margin as the sky disc so
# the ground reset cannot expose sky colors in the overscan strip.
SKY_DISC_OVERSCAN_DEG = 0.75
# The sky disc is a smooth background layer, so render it at one quarter of the
# viewport width and height and let the compositor scale it to the final
# surface.
SKY_DISC_RENDER_SCALE = 0.25


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """Performs a smooth Hermite interpolation between 0 and 1."""
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return float(t * t * (3.0 - 2.0 * t))


def _low_horizon_warm_amount(
    view_alt_deg: np.ndarray,
    sun_alt_deg: float,
) -> np.ndarray:
    """Return the low-horizon warm-haze amount before twilight fading."""
    sunset = 1.0 - _smoothstep(SUNSET_START_ALT_DEG, SUNSET_END_ALT_DEG, sun_alt_deg)
    warm_alt_deg = LOW_HORIZON_WARM_ALT_DEG + (
        LOW_HORIZON_WARM_ALT_EXPANSION_DEG * sunset
    )
    warm_strength = LOW_HORIZON_WARM_STRENGTH * (
        1.0 + ((LOW_HORIZON_WARM_MAX_STRENGTH_SCALE - 1.0) * sunset)
    )
    return (
        warm_strength
        * (1.0 - np.clip(view_alt_deg / warm_alt_deg, 0.0, 1.0))
        * sunset
    )


def _atmospheric_horizon_warm_amount(view_alt_deg: np.ndarray) -> np.ndarray:
    """Return a subtle sun-independent warm tint near the horizon."""
    return HORIZON_ATMOSPHERIC_WARM_STRENGTH * (
        1.0
        - np.clip(
            view_alt_deg / HORIZON_ATMOSPHERIC_WARM_ALT_DEG,
            0.0,
            1.0,
        )
    )


def _inverse_project_disc(
    width_px: int,
    height_px: int,
    geometry: ScreenGeometry,
    projection: ViewProjection,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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

    # Bearing from local north (clockwise): north=(0,-1), east=(1,0).
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
    if sun_alt_deg <= SUNLIGHT_FLOOR_ALT_DEG:
        return np.repeat(NIGHT_SKY_RGB[None, :], n, axis=0)

    tau = float(np.clip((TURBIDITY - 2.0) / 8.0, 0.0, 1.0))
    colorfulness = float(np.clip(1.0 + (float(saturation) - 1.0) * SATURATION_CHROMA_SCALE, 0.75, 1.25))
    t_alt = np.clip(view_alt_deg / 90.0 * 1.15, 0.0, 1.0)
    sun_up = _smoothstep(SUNLIGHT_FLOOR_ALT_DEG, 6.0, sun_alt_deg)
    twilight = _smoothstep(-11.0, 0.0, sun_alt_deg)
    sunset = 1.0 - _smoothstep(SUNSET_START_ALT_DEG, SUNSET_END_ALT_DEG, sun_alt_deg)
    low_altitude = 1.0 - t_alt
    high_altitude = t_alt

    horizon_day = HORIZON_DAY_RGB[None, :]
    zenith_day = ZENITH_DAY_RGB[None, :]
    base = horizon_day + (zenith_day - horizon_day) * t_alt[:, None]

    low_altitude_whitening = np.clip(
        LOW_ALTITUDE_WHITENING_STRENGTH
        * np.power(low_altitude, LOW_ALTITUDE_WHITENING_EXPONENT),
        0.0,
        1.0,
    )
    base = base + (
        LOW_ALTITUDE_SKY_RGB[None, :] - base
    ) * low_altitude_whitening[:, None]

    atmospheric_horizon_warm_amount = _atmospheric_horizon_warm_amount(
        view_alt_deg
    )
    base = base + (
        HORIZON_ATMOSPHERIC_WARM_RGB[None, :] - base
    ) * atmospheric_horizon_warm_amount[:, None]

    low_horizon_warm_amount = _low_horizon_warm_amount(
        view_alt_deg,
        sun_alt_deg,
    )
    low_horizon_warm_rgb = LOW_HORIZON_WARM_RGB + (
        SUNSET_RGB - LOW_HORIZON_WARM_RGB
    ) * sunset
    color = base + (
        low_horizon_warm_rgb[None, :] - base
    ) * low_horizon_warm_amount[:, None]

    a1 = np.radians(view_alt_deg)
    z1 = np.radians(view_az_deg)
    a2 = math.radians(sun_alt_deg)
    z2 = math.radians(sun_az_deg)
    cos_g = np.sin(a1) * math.sin(a2) + np.cos(a1) * math.cos(a2) * np.cos(z2 - z1)
    cos_g = np.clip(cos_g, -1.0, 1.0)

    forward = np.maximum(0.0, cos_g)
    back = np.maximum(0.0, -cos_g)

    # Add a broad daytime blue-dome contribution as the Sun rises.
    sun_alt_blue = _smoothstep(
        SUN_ALT_BLUE_START_DEG,
        SUN_ALT_BLUE_END_DEG,
        sun_alt_deg,
    )
    blue_dome_amount = sun_alt_blue * (0.35 + 0.65 * high_altitude)
    blue_dome_strength = np.clip(
        SUN_ALT_BLUE_STRENGTH * blue_dome_amount * colorfulness,
        0.0,
        1.0,
    )
    color = color + (BLUE_DOME_RGB[None, :] - color) * blue_dome_strength[:, None]

    anti_amount = (back**ANTI_SOLAR_EXPONENT) * sun_up * (0.25 + 0.75 * high_altitude)
    anti_strength = np.clip(ANTI_SOLAR_STRENGTH * anti_amount * (0.85 + 0.15 * colorfulness), 0.0, 1.0)
    color = color + (ANTI_SOLAR_RGB[None, :] - color) * anti_strength[:, None]

    glow_exponent = SUN_GLOW_EXPONENT_BASE - 0.25 * tau
    sun_glow_amount = (forward**glow_exponent) * sun_up
    sun_glow_strength = np.clip(SUN_GLOW_STRENGTH * sun_glow_amount * (0.92 + 0.08 * colorfulness), 0.0, 1.0)
    sun_glow_rgb = SUN_GLOW_RGB + (
        SUNSET_RGB - SUN_GLOW_RGB
    ) * (SUN_GLOW_SUNSET_COLOR_MIX * sunset)
    high_alt_glow = _smoothstep(
        SUN_GLOW_HIGH_ALT_START_DEG,
        SUN_GLOW_HIGH_ALT_END_DEG,
        sun_alt_deg,
    )
    sun_glow_rgb = sun_glow_rgb + (
        SUN_HIGH_ALT_GLOW_RGB - sun_glow_rgb
    ) * (SUN_GLOW_HIGH_ALT_COLOR_MIX * high_alt_glow)
    additive_glow_color = color + sun_glow_rgb[None, :] * sun_glow_strength[:, None]
    replace_glow_color = color + (
        sun_glow_rgb[None, :] - color
    ) * sun_glow_strength[:, None]
    color = (
        additive_glow_color * SUN_GLOW_ADDITIVE_MIX
        + replace_glow_color * SUN_GLOW_REPLACE_MIX
    )

    sunset_amount = sunset * low_altitude * (forward ** (1.20 - 0.20 * tau))
    sunset_amount *= 0.70 + 0.30 * sun_up
    solar_glare_sunset_factor = 1.0 - (
        1.0 - SUNSET_SOLAR_GLARE_FACTOR
    ) * forward
    sunset_amount *= solar_glare_sunset_factor
    sunset_strength = np.clip(SUNSET_STRENGTH * sunset_amount * (0.92 + 0.08 * colorfulness), 0.0, 1.0)
    color = color + (SUNSET_RGB[None, :] - color) * sunset_strength[:, None]

    # Keep a faint blue night-sky floor at every solar altitude, then add the
    # daylight contribution as it emerges through twilight.
    color = NIGHT_SKY_RGB[None, :] + color * twilight

    return np.clip(color, 0.0, 1.0).astype(np.float32)


def sky_color_samples(
    view_alt_deg: np.ndarray,
    view_az_deg: np.ndarray,
    sun_altaz: Tuple[float, float],
    *,
    exposure: float = 1.3,
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


def sky_color_at_direction(
    view_alt_deg: float,
    view_az_deg: float,
    sun_altaz: Tuple[float, float],
    *,
    alpha: float = 1.0,
    exposure: float = 1.3,
    eclipse_factor: float = 1.0,
) -> tuple[int, int, int, int]:
    """Return the rendered sky color at one horizontal-sky direction."""
    colors = sky_color_samples(
        np.asarray([view_alt_deg], dtype=np.float32),
        np.asarray([view_az_deg], dtype=np.float32),
        sun_altaz,
        alpha=alpha,
        exposure=exposure,
        eclipse_factor=eclipse_factor,
    )
    red, green, blue = np.clip(np.round(colors[0] * 255.0), 0, 255).astype(int)
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
    saturation: float,
    alpha: float,
    disc_opacity: float,
    eclipse_factor: float,
) -> QImage:
    local_geometry = ScreenGeometry(center=(center_x, center_y), radius=radius)
    alt, az, inside = _inverse_project_disc(
        width,
        height,
        local_geometry,
        projection,
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
    edge_fov_deg: float,
    content_fov_deg: float,
    sun_altaz: Tuple[float, float],
    *,
    exposure: float = 1.3,
    saturation: float = 1.35,
    alpha: float = 1.0,
    disc_opacity: float = 1.0,
    eclipse_factor: float = 1.0,
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
    projection = ViewProjection(
        view_center=(float(view_center[0]), float(view_center[1])),
        edge_fov_deg=float(edge_fov_deg),
        content_fov_deg=float(content_fov_deg),
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
        float(saturation),
        float(alpha),
        float(disc_opacity),
        float(eclipse_factor),
    )


def draw_uniform_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    *,
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
        ViewProjection(
            view_center=(float(view_center[0]), float(view_center[1])),
            edge_fov_deg=float(edge_fov_deg),
            content_fov_deg=float(content_fov_deg),
        ),
    )
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    alpha_u8 = int(round(max(0.0, min(1.0, float(disc_opacity))) * 255.0))
    rgba[..., 3][inside] = alpha_u8
    rgba[..., :3][inside] = FLAT_SKY_DISC_RGB_U8
    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
