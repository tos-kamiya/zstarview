from datetime import datetime
from typing import Any, Dict, List, Optional

import astropy.time
from PySide6.QtCore import QPoint, QPointF
from PySide6.QtGui import QColor, QPainter

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import FIELD_OF_VIEW_DEG, ThemeStyle
from ..satellite_constants import (
    SATELLITE_OVERLAY_MARKER_COLOR_RGB,
    SATELLITE_OVERLAY_MARKER_MAX_ALPHA,
)
from ..satellites.types import SatelliteOverlayPoint
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy
from .guides import draw_gauge_cross
from ..satellites import project_satellite_records

_SATELLITE_HOVER_MIN_RADIUS_PX = 12.0
_SATELLITE_HOVER_RADIUS_SCALE = 20.0


def _project_altaz_to_normalized_xy(
    alt_deg: float,
    az_deg: float,
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
) -> tuple[float, float]:
    return altaz_to_normalized_xy(
        alt_deg,
        az_deg,
        view_center,
        edge_fov_deg=edge_fov_deg,
    )


def _satellite_hover_radius_px(point: SatelliteOverlayPoint) -> float:
    return max(
        _SATELLITE_HOVER_MIN_RADIUS_PX,
        float(point.marker_scale) * _SATELLITE_HOVER_RADIUS_SCALE,
    )


def find_highlighted_satellite(
    satellite_records_by_group: object | None = None,
    mouse_pos: QPoint | QPointF | None = None,
    geometry: ScreenGeometry | None = None,
    view_center: tuple[float, float] | None = None,
    *,
    viewer_data: ViewerData | None = None,
    satellite_points: list[SatelliteOverlayPoint] | tuple[SatelliteOverlayPoint, ...] | None = None,
    observer_lat: float | None = None,
    observer_lon: float | None = None,
    observer_height_m: float | None = None,
    time_obj: astropy.time.Time | None = None,
    element_epoch_utc: datetime | None = None,
    edge_fov_deg: float = FIELD_OF_VIEW_DEG,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> tuple[SatelliteOverlayPoint, QPointF] | None:
    if viewer_data is not None:
        view_center = viewer_data.view_center
        observer_lat = viewer_data.lat_deg
        observer_lon = viewer_data.lon_deg
        observer_height_m = viewer_data.observer_height_m
        edge_fov_deg = float(viewer_data.edge_fov_deg)
        content_fov_deg = float(viewer_data.content_fov_deg)
    if view_center is None or geometry is None:
        return None
    if satellite_points is None:
        if isinstance(satellite_records_by_group, (list, tuple)) and all(
            isinstance(point, SatelliteOverlayPoint) for point in satellite_records_by_group
        ):
            satellite_points = satellite_records_by_group
        else:
            if (
                time_obj is None
                or observer_lat is None
                or observer_lon is None
                or observer_height_m is None
            ):
                return None
            satellite_points = project_satellite_records(
                satellite_records_by_group or {},
                observer_lat=float(observer_lat),
                observer_lon=float(observer_lon),
                observer_height_m=float(observer_height_m),
                time_obj=time_obj,
            )
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
        nx, ny = _project_altaz_to_normalized_xy(alt, az, view_center, edge_fov_deg=edge_fov_deg)
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
    satellite_records_by_group: object | None = None,
    view_center: tuple[float, float] | None = None,
    *,
    viewer_data: ViewerData | None = None,
    satellite_points: list[SatelliteOverlayPoint] | tuple[SatelliteOverlayPoint, ...] | None = None,
    observer_lat: float | None = None,
    observer_lon: float | None = None,
    observer_height_m: float | None = None,
    time_obj: astropy.time.Time | None = None,
    element_epoch_utc: datetime | None = None,
    opacity: float = 1.0,
    highlighted_satellite: SatelliteOverlayPoint | None = None,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    theme: ThemeStyle,
    edge_fov_deg: float = FIELD_OF_VIEW_DEG,
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
    marker_scale: float = 1.0,
) -> None:
    if viewer_data is not None:
        view_center = viewer_data.view_center
        observer_lat = viewer_data.lat_deg
        observer_lon = viewer_data.lon_deg
        observer_height_m = viewer_data.observer_height_m
        edge_fov_deg = float(viewer_data.edge_fov_deg)
        content_fov_deg = float(viewer_data.content_fov_deg)
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if layer_opacity <= 0.0:
        return
    if view_center is None:
        return

    if satellite_points is None:
        if isinstance(satellite_records_by_group, (list, tuple)) and all(
            isinstance(point, SatelliteOverlayPoint) for point in satellite_records_by_group
        ):
            satellite_points = satellite_records_by_group
        else:
            if (
                time_obj is None
                or observer_lat is None
                or observer_lon is None
                or observer_height_m is None
            ):
                return
            satellite_points = project_satellite_records(
                satellite_records_by_group or {},
                observer_lat=float(observer_lat),
                observer_lon=float(observer_lon),
                observer_height_m=float(observer_height_m),
                time_obj=time_obj,
            )
    if not satellite_points:
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
        nx, ny = _project_altaz_to_normalized_xy(alt, az, view_center, edge_fov_deg=edge_fov_deg)
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
