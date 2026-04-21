from typing import Any, Dict, List, Optional

from PySide6.QtCore import QPoint, QPointF
from PySide6.QtGui import QColor, QPainter

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import FIELD_OF_VIEW_DEG
from ..satellite_constants import (
    SATELLITE_OVERLAY_MARKER_COLOR_RGB,
    SATELLITE_OVERLAY_MARKER_MAX_ALPHA,
)
from ..satellites.types import SatelliteOverlayPoint
from ..types import ScreenGeometry
from .geometry import normalized_to_screen_xy
from .guides import draw_gauge_cross

_SATELLITE_HOVER_MIN_RADIUS_PX = 12.0
_SATELLITE_HOVER_RADIUS_SCALE = 20.0


def _satellite_hover_radius_px(point: SatelliteOverlayPoint) -> float:
    return max(
        _SATELLITE_HOVER_MIN_RADIUS_PX,
        float(point.marker_scale) * _SATELLITE_HOVER_RADIUS_SCALE,
    )


def find_highlighted_satellite(
    satellite_points: list[SatelliteOverlayPoint] | None,
    mouse_pos: QPoint | QPointF | None,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    *,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> tuple[SatelliteOverlayPoint, QPointF] | None:
    if not satellite_points or mouse_pos is None:
        return None

    mouse_x = float(mouse_pos.x())
    mouse_y = float(mouse_pos.y())
    best_point: tuple[SatelliteOverlayPoint, QPointF] | None = None
    best_dist_sq = float("inf")
    for point in satellite_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        if not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        dist_sq = (mouse_x - float(px)) ** 2 + (mouse_y - float(py)) ** 2
        hover_radius = _satellite_hover_radius_px(point)
        if dist_sq <= hover_radius * hover_radius and dist_sq < best_dist_sq:
            best_dist_sq = dist_sq
            best_point = (point, QPointF(float(px), float(py)))
    return best_point


def draw_satellite_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    satellite_points: list[SatelliteOverlayPoint] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 1.0,
    highlighted_satellite: SatelliteOverlayPoint | None = None,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    preset: str = "night",
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
    marker_scale: float = 1.0,
) -> None:
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if not satellite_points or layer_opacity <= 0.0:
        return

    painter.save()
    width_scale = max(1.0, float(marker_scale))
    marker_color = QColor(
        *SATELLITE_OVERLAY_MARKER_COLOR_RGB,
        max(
            0, min(255, int(round(SATELLITE_OVERLAY_MARKER_MAX_ALPHA * layer_opacity)))
        ),
    )
    for point in satellite_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        if not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pos = QPointF(float(px), float(py))
        draw_gauge_cross(
            painter,
            marker_color,
            pos,
            scale=float(point.marker_scale) * width_scale,
            pen_width=2.0 if point is highlighted_satellite else 1.0,
        )
    painter.restore()
