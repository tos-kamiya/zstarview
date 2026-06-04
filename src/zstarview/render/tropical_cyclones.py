from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone

from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..location_resolver.place_projection import project_place_targets_to_altaz
from ..paths import ThemeStyle
from ..types import ScreenGeometry, ViewerData
from ..tropical_cyclones.models import TropicalCycloneSnapshot
from ..tropical_cyclones.models import project_tropical_cyclone_snapshot
from .asterisms import (
    ASTERISM_HIGHLIGHT_CORE_WIDTH,
)
from .geometry import normalized_to_screen_xy

TROPICAL_CYCLONE_TARGET_HEIGHT_M = 0.0
TROPICAL_CYCLONE_MARKER_HEIGHT_M = 5000.0
TROPICAL_CYCLONE_TETHER_STEP_M = 1000.0
TROPICAL_CYCLONE_MAX_DISTANCE_KM = 400.0
TROPICAL_CYCLONE_TETHER_PASSES = (
    (ASTERISM_HIGHLIGHT_CORE_WIDTH, 0.80),
)
TROPICAL_CYCLONE_MARKER_ALPHA_SCALE = 0.40
TROPICAL_CYCLONE_FAR_MARKER_HALF_WIDTH_PX = 6.0
TROPICAL_CYCLONE_FAR_MARKER_HALF_HEIGHT_PX = 4.5
TROPICAL_CYCLONE_FAR_HOVER_RADIUS_PX = 16.0
TROPICAL_CYCLONE_FAR_LABEL_OFFSET_X_PX = 8.0
TROPICAL_CYCLONE_FAR_LABEL_OFFSET_Y_PX = 8.0
TROPICAL_CYCLONE_COLOR_RGB = (240, 122, 122)
TROPICAL_CYCLONE_LABEL_RGBA = (240, 122, 122, 255)
TROPICAL_CYCLONE_FAR_LABEL_RGBA = (240, 122, 122, 153)


@dataclass(frozen=True, slots=True)
class _RenderPoint:
    nx: float
    ny: float
    alt_deg: float
    az_deg: float
    distance_km: float


def _project_point_no_cutoff(
    lat_deg: float,
    lon_deg: float,
    *,
    viewer: ViewerData,
    height_m: float,
) -> _RenderPoint | None:
    projections = project_place_targets_to_altaz(
        observer_latitude_deg=float(viewer.lat_deg),
        observer_longitude_deg=float(viewer.lon_deg),
        observer_height_m=float(viewer.ground_elevation_m),
        target_latitude_deg=[float(lat_deg)],
        target_longitude_deg=[float(lon_deg)],
        target_height_m=[float(height_m)],
    )
    if not projections:
        return None
    projection = projections[0]
    view_center = tuple(float(value) for value in viewer.view_center)
    nx, ny = altaz_to_normalized_xy(
        float(projection.alt_deg),
        float(projection.az_deg),
        view_center,
        edge_fov_deg=float(viewer.edge_fov_deg),
    )
    return _RenderPoint(
        nx=float(nx),
        ny=float(ny),
        alt_deg=float(projection.alt_deg),
        az_deg=float(projection.az_deg),
        distance_km=float(projection.distance_km),
    )


def _far_marker_polygon(point: _RenderPoint, *, geometry: ScreenGeometry) -> QPolygonF:
    center_x, center_y = normalized_to_screen_xy(point.nx, point.ny, geometry)
    half_width_px = float(TROPICAL_CYCLONE_FAR_MARKER_HALF_WIDTH_PX)
    half_height_px = float(TROPICAL_CYCLONE_FAR_MARKER_HALF_HEIGHT_PX)
    return QPolygonF(
        [
            QPointF(center_x - half_width_px, center_y - half_height_px),
            QPointF(center_x + half_width_px, center_y - half_height_px),
            QPointF(center_x, center_y + half_height_px),
        ]
    )


def _filled_marker_polygon(point: _RenderPoint, *, geometry: ScreenGeometry) -> QPolygonF:
    center_x, center_y = normalized_to_screen_xy(point.nx, point.ny, geometry)
    half_width_px = float(TROPICAL_CYCLONE_FAR_MARKER_HALF_WIDTH_PX)
    height_px = float(TROPICAL_CYCLONE_FAR_MARKER_HALF_HEIGHT_PX) * 2.0
    return QPolygonF(
        [
            QPointF(center_x - half_width_px, center_y - height_px),
            QPointF(center_x + half_width_px, center_y - height_px),
            QPointF(center_x, center_y),
        ]
    )


def _project_marker_tether_points(
    *,
    ground_point: _RenderPoint,
    lat_deg: float,
    lon_deg: float,
    viewer: ViewerData,
) -> tuple[_RenderPoint, ...] | None:
    max_height_m = int(round(float(TROPICAL_CYCLONE_MARKER_HEIGHT_M)))
    step_m = int(round(float(TROPICAL_CYCLONE_TETHER_STEP_M)))
    if step_m <= 0:
        return None

    tether_points: list[_RenderPoint] = [ground_point]
    height_m = step_m
    while height_m < max_height_m:
        point = _project_point_no_cutoff(
            lat_deg,
            lon_deg,
            viewer=viewer,
            height_m=float(height_m),
        )
        if point is None:
            return None
        tether_points.append(point)
        height_m += step_m

    return tuple(tether_points)


def _draw_marker_tether(
    painter: QPainter,
    tether_points: tuple[_RenderPoint, ...],
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
) -> None:
    if len(tether_points) < 2:
        return
    width_px, pass_alpha_scale = TROPICAL_CYCLONE_TETHER_PASSES[0]
    tether_rgba = (
        int(color_rgba[0]),
        int(color_rgba[1]),
        int(color_rgba[2]),
        max(1, int(round(255.0 * float(pass_alpha_scale)))),
    )
    pen = QPen(QColor(*tether_rgba), float(width_px), Qt.PenStyle.SolidLine)
    pen.setCosmetic(True)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    for start_point, end_point in zip(tether_points, tether_points[1:]):
        start = QPointF(*normalized_to_screen_xy(start_point.nx, start_point.ny, geometry))
        end = QPointF(*normalized_to_screen_xy(end_point.nx, end_point.ny, geometry))
        painter.drawLine(start, end)


def _draw_far_cyclone_marker(
    painter: QPainter,
    point: _RenderPoint,
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
    highlighted: bool = False,
) -> QPointF:
    pen_width = 2.0 if highlighted else 1.0
    pen = QPen(QColor(*color_rgba), float(pen_width), Qt.PenStyle.SolidLine)
    pen.setCosmetic(True)
    painter.setPen(pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    marker_polygon = _far_marker_polygon(point, geometry=geometry)
    painter.drawPolygon(marker_polygon)
    bounds = marker_polygon.boundingRect()
    return QPointF(
        float(bounds.right()) + float(TROPICAL_CYCLONE_FAR_LABEL_OFFSET_X_PX),
        float(bounds.top()) - float(TROPICAL_CYCLONE_FAR_LABEL_OFFSET_Y_PX),
    )


def _draw_filled_cyclone_marker(
    painter: QPainter,
    point: _RenderPoint,
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
) -> None:
    pen = QPen(QColor(*color_rgba), 1.5, Qt.PenStyle.SolidLine)
    pen.setCosmetic(True)
    painter.setPen(pen)
    painter.setBrush(QColor(*color_rgba))
    painter.drawPolygon(_filled_marker_polygon(point, geometry=geometry))


def find_highlighted_tropical_cyclone(
    snapshots: tuple[TropicalCycloneSnapshot, ...] | list[TropicalCycloneSnapshot] | object | None = None,
    mouse_pos: QPoint | QPointF | None = None,
    geometry: ScreenGeometry | None = None,
    *,
    viewer_data: ViewerData | None = None,
    time_obj: datetime | None = None,
) -> tuple[TropicalCycloneSnapshot, QPointF] | None:
    if viewer_data is None or time_obj is None or mouse_pos is None or geometry is None:
        return None
    if not isinstance(snapshots, (list, tuple)) or not snapshots:
        return None

    mouse_x = float(mouse_pos.x())
    mouse_y = float(mouse_pos.y())
    best_match: tuple[TropicalCycloneSnapshot, QPointF] | None = None
    best_dist_sq = float("inf")
    for snapshot in snapshots:
        if not isinstance(snapshot, TropicalCycloneSnapshot):
            continue
        projected_snapshot = project_tropical_cyclone_snapshot(snapshot, time_obj)
        observed = projected_snapshot.observed_position
        center_point = _project_point_no_cutoff(
            observed.lat_deg,
            observed.lon_deg,
            viewer=viewer_data,
            height_m=TROPICAL_CYCLONE_TARGET_HEIGHT_M,
        )
        if center_point is None or float(center_point.distance_km) <= float(TROPICAL_CYCLONE_MAX_DISTANCE_KM):
            continue
        pos_x, pos_y = normalized_to_screen_xy(center_point.nx, center_point.ny, geometry)
        dist_sq = (mouse_x - float(pos_x)) ** 2 + (mouse_y - float(pos_y)) ** 2
        if dist_sq <= float(TROPICAL_CYCLONE_FAR_HOVER_RADIUS_PX) ** 2 and dist_sq < best_dist_sq:
            best_dist_sq = dist_sq
            best_match = (snapshot, QPointF(float(pos_x), float(pos_y)))
    return best_match


def draw_tropical_cyclone_overlay(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    snapshot: TropicalCycloneSnapshot | None,
    when_utc: datetime | None,
    theme: ThemeStyle,
    opacity: float,
    highlighted: bool = False,
    enabled: bool = True,
) -> None:
    if not enabled:
        return
    if snapshot is None:
        return
    del theme

    if when_utc is None:
        when_utc = datetime.now(timezone.utc)
    projected_snapshot = project_tropical_cyclone_snapshot(snapshot, when_utc)
    observed = projected_snapshot.observed_position
    center_point = _project_point_no_cutoff(
        observed.lat_deg,
        observed.lon_deg,
        viewer=viewer,
        height_m=TROPICAL_CYCLONE_TARGET_HEIGHT_M,
    )
    if center_point is None:
        return

    marker_rgba = (
        int(TROPICAL_CYCLONE_COLOR_RGB[0]),
        int(TROPICAL_CYCLONE_COLOR_RGB[1]),
        int(TROPICAL_CYCLONE_COLOR_RGB[2]),
        int(round(255.0 * min(1.0, max(0.0, float(opacity))))),
    )
    if float(center_point.distance_km) > float(TROPICAL_CYCLONE_MAX_DISTANCE_KM):
        marker_point = None
        label_rgba = TROPICAL_CYCLONE_FAR_LABEL_RGBA
    else:
        marker_point = _project_point_no_cutoff(
            observed.lat_deg,
            observed.lon_deg,
            viewer=viewer,
            height_m=TROPICAL_CYCLONE_MARKER_HEIGHT_M,
        )
        if marker_point is None:
            label_rgba = TROPICAL_CYCLONE_FAR_LABEL_RGBA
        else:
            label_rgba = TROPICAL_CYCLONE_LABEL_RGBA
    marker_fill_rgba = (
        int(TROPICAL_CYCLONE_COLOR_RGB[0]),
        int(TROPICAL_CYCLONE_COLOR_RGB[1]),
        int(TROPICAL_CYCLONE_COLOR_RGB[2]),
        max(1, int(round(255.0 * float(TROPICAL_CYCLONE_MARKER_ALPHA_SCALE)))),
    )

    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    if marker_point is None:
        label_pos = _draw_far_cyclone_marker(
            painter,
            center_point,
            geometry=geometry,
            color_rgba=marker_rgba,
            highlighted=highlighted,
        )
        if not highlighted:
            label_pos = None
    else:
        tether_points = _project_marker_tether_points(
            ground_point=center_point,
            lat_deg=observed.lat_deg,
            lon_deg=observed.lon_deg,
            viewer=viewer,
        )
        if tether_points is None:
            label_pos = _draw_far_cyclone_marker(
                painter,
                center_point,
                geometry=geometry,
                color_rgba=marker_rgba,
                highlighted=highlighted,
            )
        else:
            tether_points = (*tether_points, marker_point)
            _draw_marker_tether(
                painter,
                tether_points,
                geometry=geometry,
                color_rgba=marker_rgba,
            )
            _draw_filled_cyclone_marker(
                painter,
                marker_point,
                geometry=geometry,
                color_rgba=marker_fill_rgba,
            )
            screen_x, screen_y = normalized_to_screen_xy(
                marker_point.nx,
                marker_point.ny,
                geometry,
            )
            label_pos = QPointF(float(screen_x + 8.0), float(screen_y - 8.0))
    painter.setPen(QColor(*label_rgba))
    if label_pos is not None:
        painter.drawText(
            label_pos,
            projected_snapshot.storm_name,
        )
    painter.restore()
