import math
from typing import Tuple

import numpy as np

from ..types import ScreenGeometry
from ..paths import OBSERVER_MIN_ALT_DEG


def _altaz_to_normalized_xy_vectorized(
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


def _normalized_to_screen_xy_vectorized(
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
    margin_x = 0
    margin_y = 0
    avail_w = max(2, int(width_px) - margin_x * 2)
    avail_h = max(2, int(height_px) - margin_y * 2)
    alt = max(OBSERVER_MIN_ALT_DEG, min(90.0, float(view_alt_deg)))

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
