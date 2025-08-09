import math
import numpy as np
from typing import Tuple, List

from PyQt5.QtCore import QPointF, QRectF
from PyQt5.QtGui import QPainter, QColor, QRadialGradient
from PyQt5.QtCore import Qt

from ..paths import FIELD_OF_VIEW_DEG
from ..types import ScreenGeometry, SkyData, ViewerData
from .draw import bv_to_qcolor


def draw_stars_fully_vectorized(painter: QPainter, geometry: ScreenGeometry, sky_data: SkyData, viewer_data: ViewerData, star_base_radius: float):
    """Fully vectorized star rendering function (including field-of-view check)."""

    if not sky_data.stars:
        return

    # Convert all star data to NumPy arrays
    n_stars = len(sky_data.stars)
    alts = np.array([star.alt for star in sky_data.stars], dtype=np.float64)
    azs = np.array([star.az for star in sky_data.stars], dtype=np.float64)
    vmags = np.array([star.vmag for star in sky_data.stars], dtype=np.float64)
    bvs = [star.bv for star in sky_data.stars]

    # Vectorized field-of-view check
    fov_mask = _vectorized_is_in_fov(alts, azs, viewer_data.view_center)

    # Extract only stars within the field of view
    if not np.any(fov_mask):
        return

    visible_alts = alts[fov_mask]
    visible_azs = azs[fov_mask]
    visible_vmags = vmags[fov_mask]
    visible_bvs = [bvs[i] for i in np.where(fov_mask)[0]]

    # Vectorized coordinate transformation and size calculation
    positions = _vectorized_altaz_to_screen(visible_alts, visible_azs, viewer_data, geometry)
    sizes = np.maximum(0.1, star_base_radius * (10.0 ** (-0.2 * visible_vmags))) * geometry.radius / 500.0

    # Calculate colors
    colors = [bv_to_qcolor(bv) for bv in visible_bvs]

    # Vectorized calculation of alpha values for small stars
    small_star_mask = sizes < 4.0
    alphas = np.where(small_star_mask, np.clip(sizes / 4.0, 0.25, 1.0), 1.0)

    # Perform drawing
    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)

    # Separate small and large stars for efficient rendering
    _draw_vectorized_stars(painter, positions, sizes, colors, alphas)

    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)


def _vectorized_is_in_fov(alts: np.ndarray, azs: np.ndarray, view_center: Tuple[float, float]) -> np.ndarray:
    """Vectorized field-of-view check."""
    center_alt, center_az = view_center
    center_alt_rad = np.radians(center_alt)
    center_az_rad = np.radians(center_az)

    # Convert angles to radians in bulk
    alt_rads = np.radians(alts)
    az_rads = np.radians(azs)

    # Vectorized spherical distance calculation
    sin_center_alt = np.sin(center_alt_rad)
    cos_center_alt = np.cos(center_alt_rad)

    cos_theta = sin_center_alt * np.sin(alt_rads) + cos_center_alt * np.cos(alt_rads) * np.cos(az_rads - center_az_rad)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)
    theta_deg = np.degrees(theta)

    # Check if within field of view
    return theta_deg <= (FIELD_OF_VIEW_DEG / 2.0)


def _vectorized_altaz_to_screen(alts: np.ndarray, azs: np.ndarray, viewer_data: ViewerData, geometry: ScreenGeometry) -> List[QPointF]:
    """Vectorized coordinate transformation (altaz → screen coordinates)."""

    center_alt, center_az = viewer_data.view_center
    center_alt_rad = np.radians(center_alt)
    center_az_rad = np.radians(center_az)

    # Convert angles to radians in bulk
    alt_rads = np.radians(alts)
    az_rads = np.radians(azs)

    # Vectorized spherical coordinate calculations
    sin_center_alt = np.sin(center_alt_rad)
    cos_center_alt = np.cos(center_alt_rad)

    # Vectorized spherical distance calculation
    cos_theta = sin_center_alt * np.sin(alt_rads) + cos_center_alt * np.cos(alt_rads) * np.cos(az_rads - center_az_rad)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)
    r = theta / (np.pi / 2.0)

    # Vectorized normalized coordinate calculation
    dx = np.cos(alt_rads) * np.sin(az_rads - center_az_rad)
    dy = cos_center_alt * np.sin(alt_rads) - sin_center_alt * np.cos(alt_rads) * np.cos(az_rads - center_az_rad)

    # Vectorized length normalization
    length = np.sqrt(dx**2 + dy**2)

    # Avoid division by zero using a mask
    non_zero_mask = length > 1e-10
    dx_normalized = np.zeros_like(dx)
    dy_normalized = np.zeros_like(dy)
    dx_normalized[non_zero_mask] = dx[non_zero_mask] / length[non_zero_mask]
    dy_normalized[non_zero_mask] = dy[non_zero_mask] / length[non_zero_mask]

    # Normalized coordinates
    nx = r * dx_normalized
    ny = -r * dy_normalized

    # Vectorized conversion to screen coordinates
    screen_x = geometry.center[0] + nx * geometry.radius
    screen_y = geometry.center[1] + ny * geometry.radius

    # Return as a list of QPointF objects
    return [QPointF(float(x), float(y)) for x, y in zip(screen_x, screen_y)]


def _draw_vectorized_stars(painter: QPainter, positions: List[QPointF], sizes: np.ndarray, colors: List[QColor], alphas: np.ndarray):
    """Efficient drawing using vectorized data."""

    small_star_mask = sizes < 4.0

    # Draw large stars
    painter.setPen(Qt.PenStyle.NoPen)
    for i, (pos, size, color) in enumerate(zip(positions, sizes, colors)):
        if not small_star_mask[i]:
            r = math.sqrt(float(size))
            gradient = QRadialGradient(pos, r)
            gradient.setColorAt(0, color)
            color_transparent = QColor(color)
            color_transparent.setAlpha(0)
            gradient.setColorAt(1, color_transparent)
            painter.setBrush(gradient)
            painter.drawEllipse(pos, r, r)

    # Draw small stars in bulk
    painter.setPen(Qt.PenStyle.NoPen)
    for i, (pos, size, color, alpha) in enumerate(zip(positions, sizes, colors, alphas)):
        if small_star_mask[i]:
            color.setAlphaF(float(alpha))
            painter.fillRect(QRectF(pos.x() - 1, pos.y() - 1, 2, 2), color)


# # Standalone vectorized field-of-view check function
# def is_in_fov_vectorized(alts: np.ndarray, azs: np.ndarray,
#                         view_center: Tuple[float, float]) -> np.ndarray:
#     """
#     Run a vectorized field-of-view check for multiple coordinates.
# 
#     Args:
#         alts: Array of altitudes (degrees)
#         azs: Array of azimuths (degrees)
#         view_center: Center coordinates of the view (alt, az)
#         field_of_view_deg: Field of view in degrees
# 
#     Returns:
#         Boolean array indicating whether each coordinate is within the field of view.
#     """
#     return _vectorized_is_in_fov(alts, azs, view_center)
# 
# # Wrapper for compatibility with the original single-coordinate function
# def is_in_fov_single(alt: float, az: float, view_center: Tuple[float, float]) -> bool:
#     """Field-of-view check for a single coordinate (compatible with the original function)."""
#     alts = np.array([alt])
#     azs = np.array([az])
#     result = _vectorized_is_in_fov(alts, azs, view_center)
#     return bool(result[0])
