import math
from typing import Tuple

import numpy as np

from ..types import ScreenGeometry


def shortest_azimuth_delta_deg(az_deg: float, center_az_deg: float) -> float:
    """Return the signed shortest azimuth delta in degrees."""
    return ((float(az_deg) - float(center_az_deg) + 180.0) % 360.0) - 180.0


def altaz_to_normalized_xy_vectorized(
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    view_center_altaz_deg: Tuple[float, float],
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorized conversion of altitude/azimuth to normalized screen coordinates."""
    center_alt, center_az = view_center_altaz_deg
    alt1, az1 = np.radians(center_alt), np.radians(center_az)
    alt2, az2 = np.radians(alt_deg), np.radians(az_deg)

    cos_theta = (
        np.sin(alt1) * np.sin(alt2)
        + np.cos(alt1) * np.cos(alt2) * np.cos(az2 - az1)
    )
    theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))

    r = theta / (math.pi / 2)
    dx = np.cos(alt2) * np.sin(az2 - az1)
    dy = np.cos(alt1) * np.sin(alt2) - np.sin(alt1) * np.cos(alt2) * np.cos(az2 - az1)
    length = np.hypot(dx, dy)
    length[length == 0] = 1.0
    nx = r * dx / length
    ny = -r * dy / length
    return (nx, ny)


def normalized_to_screen_xy_vectorized(
    nx: np.ndarray, ny: np.ndarray, geometry: ScreenGeometry
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorized conversion of normalized coordinates to screen coordinates."""
    return (
        geometry.center[0] + nx * geometry.radius,
        geometry.center[1] + ny * geometry.radius,
    )


def get_screen_geometry(
    width_px: int, height_px: int, view_alt_deg: float
) -> ScreenGeometry:
    """Calculate circular viewport geometry."""
    margin_x = 10
    margin_y = 10
    avail_w = max(2, int(width_px) - margin_x * 2)
    avail_h = max(2, int(height_px) - margin_y * 2)
    alt = max(0.0, min(90.0, float(view_alt_deg)))

    if width_px >= height_px:
        r_height = int(avail_h / (1.0 + alt / 90.0))
        r_width = avail_w // 2
        radius_px = max(1, min(r_width, r_height))
        center = (int(width_px) // 2, margin_y + radius_px)
    else:
        radius_px = max(1, min(avail_w // 2, avail_h // 2))
        center = (int(width_px) // 2, int(height_px) // 2)
    return ScreenGeometry(center, radius_px)


def normalized_to_screen_xy(
    nx: float, ny: float, geometry: ScreenGeometry
) -> Tuple[float, float]:
    """Convert normalized coordinates to screen coordinates."""
    return geometry.center[0] + nx * geometry.radius, geometry.center[1] + ny * geometry.radius


def altaz_to_cylindrical_normalized_xy(
    alt_deg: float,
    az_deg: float,
    view_center_altaz_deg: Tuple[float, float],
    *,
    fov_deg: float = 180.0,
) -> Tuple[float, float]:
    """Convert altitude/azimuth to a cylindrical-style normalized screen position.

    This maps azimuth linearly across the horizontal field and altitude linearly
    across the vertical field, which makes near urban skylines read more like a
    silhouette on a surrounding cylinder than a line stuck to the sky dome.
    """
    center_alt, center_az = view_center_altaz_deg
    half_fov = max(1e-6, float(fov_deg) * 0.5)
    delta_az = shortest_azimuth_delta_deg(float(az_deg), float(center_az))
    delta_alt = float(alt_deg) - float(center_alt)
    nx = delta_az / half_fov
    ny = -(delta_alt / half_fov)
    return (nx, ny)


def altaz_to_urban_skyline_normalized_xy(
    alt_deg: float,
    az_deg: float,
    view_center_altaz_deg: Tuple[float, float],
    *,
    cylindrical_until_alt_deg: float = 0.0,
    spherical_from_alt_deg: float = 5.0,
) -> Tuple[float, float]:
    """Project urban skyline points with a cylindrical-to-spherical blend.

    Near the horizon, cylindrical mapping reads more naturally for city
    silhouettes. Higher above the horizon, the regular sky-dome projection is
    less visually surprising. This helper blends between the two.
    """
    cylindrical = altaz_to_cylindrical_normalized_xy(alt_deg, az_deg, view_center_altaz_deg)
    spherical_x, spherical_y = altaz_to_normalized_xy_vectorized(
        np.array([alt_deg], dtype=np.float64),
        np.array([az_deg], dtype=np.float64),
        view_center_altaz_deg,
    )
    spherical = (float(spherical_x[0]), float(spherical_y[0]))

    low = float(cylindrical_until_alt_deg)
    high = max(low, float(spherical_from_alt_deg))
    if float(alt_deg) <= low:
        return cylindrical
    if float(alt_deg) >= high or math.isclose(high, low):
        return spherical

    t = (float(alt_deg) - low) / (high - low)
    nx = cylindrical[0] * (1.0 - t) + spherical[0] * t
    ny = cylindrical[1] * (1.0 - t) + spherical[1] * t
    return (nx, ny)
