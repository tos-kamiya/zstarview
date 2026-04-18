import math
from typing import Tuple

import numpy as np
from PySide6.QtGui import QImage

from ..paths import PALETTE_NEVER_RISES_RGB, TERRAIN_HORIZON_LINE_COLOR
from ..types import ScreenGeometry
from .qt_image import np_rgba_to_qimage


TURBIDITY = 5  # 2 (clear blue sky) to 10 (hazy white sky)
GROUND_TINT_RGB = np.array(TERRAIN_HORIZON_LINE_COLOR, dtype=np.float32) / 255.0
NEVER_RISES_TINT_RGB = np.array(PALETTE_NEVER_RISES_RGB, dtype=np.float32) / 255.0
NEVER_RISES_TINT_STRENGTH = 0.12
FLAT_SKY_DISC_RGB_U8 = np.array([10, 10, 10], dtype=np.uint8)


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
    content_fov_deg: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project pixels up to the requested content FOV."""
    radius_px = max(1.0, float(geometry.radius))
    ys = (np.arange(height_px, dtype=np.float32) - float(geometry.center[1])) / radius_px
    xs = (np.arange(width_px, dtype=np.float32) - float(geometry.center[0])) / radius_px
    nx, ny = np.meshgrid(xs, ys)

    rr2 = nx * nx + ny * ny
    max_r = max(0.0, float(content_fov_deg) / 90.0)
    inside = rr2 <= (max_r * max_r)
    if not np.any(inside):
        return np.array([], dtype=np.float32), np.array([], dtype=np.float32), inside

    r = np.sqrt(rr2[inside]).astype(np.float32)
    theta = r * (np.pi / 2.0)

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
) -> np.ndarray:
    """Vectorized sky color model for many sky directions at once."""
    sun_alt_deg, sun_az_deg = sun_altaz
    n = view_alt_deg.shape[0]
    if n == 0:
        return np.zeros((0, 3), dtype=np.float32)
    if sun_alt_deg <= -10.0:
        return np.zeros((n, 3), dtype=np.float32)

    tau = float(np.clip((TURBIDITY - 2.0) / 8.0, 0.0, 1.0))
    t_alt = np.clip(view_alt_deg / 90.0, 0.0, 1.0)
    sun_up = _smoothstep(-8.0, 6.0, sun_alt_deg)
    twilight = _smoothstep(-10.0, 0.0, sun_alt_deg)

    horizon_day = np.array([0.98, 0.70, 0.45], dtype=np.float32)
    zenith_day = np.array([0.18, 0.42, 0.93], dtype=np.float32)
    base = horizon_day[None, :] + (zenith_day - horizon_day)[None, :] * t_alt[:, None]

    haze = (0.22 + 0.48 * (1.0 - t_alt)) * (0.65 + 0.55 * tau)
    base = base + (1.0 - base) * haze[:, None]

    a1 = np.radians(view_alt_deg)
    z1 = np.radians(view_az_deg)
    a2 = math.radians(sun_alt_deg)
    z2 = math.radians(sun_az_deg)
    cos_g = np.sin(a1) * math.sin(a2) + np.cos(a1) * math.cos(a2) * np.cos(z2 - z1)
    cos_g = np.clip(cos_g, -1.0, 1.0)
    sun_angle = np.arccos(cos_g)

    f = np.maximum(0.0, np.cos(sun_angle))
    glow = f ** (1.9 - 0.55 * tau)
    sun_tint = np.array([1.0, 0.82, 0.58], dtype=np.float32) + (
        np.array([1.0, 0.92, 0.78], dtype=np.float32) - np.array([1.0, 0.82, 0.58], dtype=np.float32)
    ) * tau
    color = base + sun_tint[None, :] * (0.52 * glow * sun_up)[:, None]

    anti = np.maximum(0.0, np.cos(np.pi - sun_angle))
    anti_boost = anti**2.2
    color += np.array([0.06, 0.10, 0.20], dtype=np.float32)[None, :] * (0.24 * anti_boost * sun_up)[:, None]

    night = np.array([0.01, 0.02, 0.05], dtype=np.float32)
    color = night[None, :] + (color - night[None, :]) * twilight

    max_luma = 0.28 + 0.20 * sun_up
    luma = np.sum(color * np.array([0.299, 0.587, 0.114], dtype=np.float32)[None, :], axis=1)
    scale = np.ones_like(luma)
    clip_mask = luma > max_luma
    scale[clip_mask] = max_luma / np.maximum(luma[clip_mask], 1e-6)
    color = color * scale[:, None]

    return np.clip(color, 0.0, 1.0).astype(np.float32)


def _grade_color(
    color: np.ndarray,
    *,
    saturation: float = 1.0,
    exposure: float = 1.0,
    gamma: float = 1.0,
) -> np.ndarray:
    """
    Sky-disc-only color grading:
      1) luma-based saturation (BT.601/709)
      2) exposure scaling
      3) simple power-law gamma (on current space)
    Assumes color is a numpy array in [0,1] and *not* premultiplied.
    """
    luma = np.sum(color * np.array([0.299, 0.587, 0.114], dtype=np.float32)[None, :], axis=1, keepdims=True)

    color = luma + (color - luma) * saturation
    color *= exposure

    if gamma > 0.0 and gamma != 1.0:
        inv = 1.0 / gamma
        color = np.power(np.maximum(0.0, color), inv)

    color = color / (1.0 + 0.35 * color)

    return np.clip(color, 0.0, 1.0)


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
    """Return graded sky RGB samples for (alt, az) arrays."""
    alt = np.asarray(view_alt_deg, dtype=np.float32)
    az = np.asarray(view_az_deg, dtype=np.float32)
    if alt.shape != az.shape:
        raise ValueError("view_alt_deg and view_az_deg must have the same shape")
    if alt.size == 0:
        return np.zeros((0, 3), dtype=np.float32)

    colors = _get_sky_color_vectorized(alt.reshape(-1), az.reshape(-1), sun_altaz)
    gamma = (1.0 - alpha) * 0.2 + 1.0 if alpha < 1.0 else 1.0
    colors = _grade_color(colors, saturation=saturation, exposure=exposure, gamma=gamma)

    sky_scale = max(0.0, float(alpha))
    eclipse_scale = max(0.0, float(eclipse_factor))
    colors *= sky_scale * eclipse_scale
    return np.clip(colors, 0.0, 1.0).astype(np.float32)


def draw_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    sun_altaz: Tuple[float, float],
    observer_lat_deg: float | None = None,
    *,
    exposure: float = 1.14,
    saturation: float = 1.35,
    alpha: float = 1.0,
    disc_opacity: float = 1.0,
    eclipse_factor: float = 1.0,
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

    alt, az, inside = _inverse_project_disc(
        width,
        height,
        local_geometry,
        view_center,
        content_fov_deg=content_fov_deg,
    )

    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    if alt.size == 0:
        return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)

    colors = sky_color_samples(
        alt,
        az,
        sun_altaz,
        exposure=exposure,
        saturation=saturation,
        alpha=alpha,
        eclipse_factor=eclipse_factor,
    )
    # Mark declinations that never rise at the observer latitude.
    never_rises = np.zeros_like(alt, dtype=bool)
    if observer_lat_deg is not None:
        lat = float(np.clip(observer_lat_deg, -90.0, 90.0))
        lat_rad = math.radians(lat)
        alt_rad = np.radians(alt)
        az_rad = np.radians(az)
        # Convert (alt, az, lat) to declination:
        # sin(dec) = sin(alt)*sin(lat) + cos(alt)*cos(lat)*cos(az)
        # (az: 0=N, 90=E)
        sin_dec = np.sin(alt_rad) * math.sin(lat_rad) + np.cos(alt_rad) * math.cos(lat_rad) * np.cos(az_rad)
        sin_dec = np.clip(sin_dec, -1.0, 1.0)
        dec = np.degrees(np.arcsin(sin_dec))
        if lat > 0.0:
            threshold = lat - 90.0
            never_rises = dec <= threshold
        elif lat < 0.0:
            threshold = lat + 90.0
            never_rises = dec >= threshold
        else:
            never_rises = np.zeros_like(dec, dtype=bool)
    # Apply never-rises tint at the end so it remains visible after grading.
    if np.any(never_rises):
        colors[never_rises] = np.clip(
            colors[never_rises] + NEVER_RISES_TINT_RGB[None, :] * np.float32(NEVER_RISES_TINT_STRENGTH),
            0.0,
            1.0,
        )
    colors = np.clip(colors, 0.0, 1.0)

    rgb_u8 = np.clip(np.round(colors * 255.0), 0, 255).astype(np.uint8)
    alpha_u8 = int(round(max(0.0, min(1.0, float(disc_opacity))) * 255.0))
    rgba[..., 3][inside] = alpha_u8
    rgba[..., :3][inside] = rgb_u8

    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)


def draw_uniform_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    *,
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
        content_fov_deg=content_fov_deg,
    )
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    alpha_u8 = int(round(max(0.0, min(1.0, float(disc_opacity))) * 255.0))
    rgba[..., 3][inside] = alpha_u8
    rgba[..., :3][inside] = FLAT_SKY_DISC_RGB_U8
    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
