import math
from typing import Tuple

import numpy as np
from PySide6.QtGui import QImage

from ..types import ScreenGeometry
from ..utils.qt import np_rgba_to_qimage


TURBIDITY = 5  # 2 (clear blue sky) to 10 (hazy white sky)
GROUND_TINT_RGB = np.array([0.12, 0.19, 0.27], dtype=np.float32)
GROUND_TINT_STRENGTH = 0.2
GROUND_BASE_DARKNESS = 0.03


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """Performs a smooth Hermite interpolation between 0 and 1."""
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return float(t * t * (3.0 - 2.0 * t))


def _angle_between(alt1_deg: float, az1_deg: float, alt2_deg: float, az2_deg: float) -> float:
    """Calculates angular distance between two sky directions."""
    a1, z1 = math.radians(alt1_deg), math.radians(az1_deg)
    a2, z2 = math.radians(alt2_deg), math.radians(az2_deg)
    cos_g = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(z2 - z1)
    cos_g = max(-1.0, min(1.0, cos_g))
    return math.acos(cos_g)


def _lerp_color(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    """Linearly interpolates between two colors."""
    return a + (b - a) * t


def _luma_bt601(color: np.ndarray) -> float:
    """Return BT.601 luma from RGB in [0,1]."""
    return float(np.dot(color, np.array([0.299, 0.587, 0.114])))


def _soft_clip_luma(color: np.ndarray, max_luma: float) -> np.ndarray:
    """Keep hue/chroma while compressing highlights above `max_luma`."""
    y = _luma_bt601(color)
    if y <= 1e-6 or y <= max_luma:
        return color
    return np.clip(color * (max_luma / y), 0.0, 1.0)


def get_sun_color(sun_alt_deg: float) -> np.ndarray:
    """Determine sunlight color based on sun altitude."""
    zenith_color = np.array([0.3, 0.55, 0.98])
    horizon_color = np.array([0.95, 0.50, 0.30])
    night_color = np.array([0.01, 0.02, 0.05])

    t = float(np.clip((sun_alt_deg + 7.0) / 97.0, 0.0, 1.0))
    t = math.pow(t, 0.35)

    day_color = _lerp_color(horizon_color, zenith_color, t)
    fade = float(np.clip((-sun_alt_deg + 1.0) / 6.0, 0.0, 1.0))
    return _lerp_color(day_color, night_color, fade)


def get_sky_color(view_altaz: Tuple[float, float], sun_altaz: Tuple[float, float]) -> np.ndarray:
    """Calculate the sky color for one viewing direction."""
    sun_alt_deg, sun_az_deg = sun_altaz
    view_alt_deg, view_az_deg = view_altaz

    if sun_alt_deg <= -10.0:
        return np.array([0.0, 0.0, 0.0])

    tau = np.clip((TURBIDITY - 2.0) / 8.0, 0.0, 1.0)
    t_alt = np.clip(view_alt_deg / 90.0, 0.0, 1.0)
    sun_up = _smoothstep(-8.0, 6.0, sun_alt_deg)
    twilight = _smoothstep(-10.0, 0.0, sun_alt_deg)

    horizon_day = np.array([0.98, 0.70, 0.45])
    zenith_day = np.array([0.18, 0.42, 0.93])
    base = _lerp_color(horizon_day, zenith_day, t_alt)

    haze = (0.22 + 0.48 * (1.0 - t_alt)) * (0.65 + 0.55 * tau)
    base = _lerp_color(base, np.array([1.0, 1.0, 1.0]), haze)

    sun_angle = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)
    f = max(0.0, math.cos(sun_angle))
    glow = f ** (1.9 - 0.55 * tau)
    sun_tint = _lerp_color(np.array([1.0, 0.82, 0.58]), np.array([1.0, 0.92, 0.78]), tau)
    color = base + sun_tint * (0.52 * glow * sun_up)

    anti = max(0.0, math.cos(math.pi - sun_angle))
    anti_boost = anti ** 2.2
    color += np.array([0.06, 0.10, 0.20]) * (0.24 * anti_boost * sun_up)

    night = np.array([0.01, 0.02, 0.05])
    color = _lerp_color(night, color, twilight)

    max_luma = 0.28 + 0.20 * sun_up
    color = _soft_clip_luma(color, max_luma=max_luma)

    return np.clip(color, 0.0, 1.0)


def _inverse_project_disc(
    radius_px: int,
    view_center: Tuple[float, float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project all pixels inside the unit disc to (alt, az)."""
    size = radius_px * 2
    ys = (np.arange(size, dtype=np.float32) - radius_px) / float(radius_px)
    xs = (np.arange(size, dtype=np.float32) - radius_px) / float(radius_px)
    nx, ny = np.meshgrid(xs, ys)

    rr2 = nx * nx + ny * ny
    inside = rr2 <= 1.0
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


def grade_color(
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


def draw_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    sun_altaz: Tuple[float, float],
    *,
    exposure: float = 1.14,
    saturation: float = 1.35,
    alpha: float = 1.0,
    eclipse_factor: float = 1.0,
    mask_fov_deg: float = 90.0,
    sample_step_px: int = 10,
    min_ang_samples: int = 8,
    deriv_probe_deg: float = 0.25,
    min_theta_step_deg: float = 0.2,
    max_theta_step_deg: float = 6.0,
) -> QImage:
    """
    Draw sky color disc using one-pass NumPy inverse projection.

    The disc is always rendered as a full circle (no horizon clipping inside the disc).
    Unused legacy sampling parameters are kept for call-site compatibility.
    """
    _ = (mask_fov_deg, sample_step_px, min_ang_samples, deriv_probe_deg, min_theta_step_deg, max_theta_step_deg)

    radius = int(geometry.radius)
    size = max(2, radius * 2)
    if radius < 1:
        return QImage(size, size, QImage.Format.Format_ARGB32_Premultiplied)

    alt, az, inside = _inverse_project_disc(radius, view_center)

    rgba = np.zeros((size, size, 4), dtype=np.uint8)
    if alt.size == 0:
        return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)

    colors = _get_sky_color_vectorized(alt, az, sun_altaz)
    gamma = (1.0 - alpha) * 0.2 + 1.0 if alpha < 1.0 else 1.0
    colors = grade_color(colors, saturation=saturation, exposure=exposure, gamma=gamma)
    below_horizon = alt < 0.0
    if np.any(below_horizon):
        s = np.float32(GROUND_TINT_STRENGTH)
        base = colors[below_horizon] * np.float32(GROUND_BASE_DARKNESS)
        colors[below_horizon] = base * (1.0 - s) + GROUND_TINT_RGB[None, :] * s
    # Keep sky opacity for the upper sky only.
    # Ground tint should stay visible regardless of sky_opacity.
    sky_scale = max(0.0, float(alpha))
    eclipse_scale = max(0.0, float(eclipse_factor))
    if np.any(below_horizon):
        colors[~below_horizon] *= sky_scale * eclipse_scale
        colors[below_horizon] *= eclipse_scale
    else:
        colors *= sky_scale * eclipse_scale
    colors = np.clip(colors, 0.0, 1.0)

    rgb_u8 = np.clip(np.round(colors * 255.0), 0, 255).astype(np.uint8)
    rgba[..., 3][inside] = 255
    rgba[..., :3][inside] = rgb_u8

    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
