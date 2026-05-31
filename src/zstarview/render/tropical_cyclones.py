from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..location_resolver.place_projection import project_place_targets_to_altaz
from ..paths import ThemeStyle
from ..types import ScreenGeometry, ViewerData
from ..tropical_cyclones.models import TropicalCycloneSnapshot
from ..tropical_cyclones.models import project_tropical_cyclone_snapshot
from .geometry import normalized_to_screen_xy
from .terrain import WATER_OVERLAY_POINT_RADIUS_PX

TROPICAL_CYCLONE_TARGET_HEIGHT_M = 0.0
TROPICAL_CYCLONE_MAX_DISTANCE_KM = 400.0
TROPICAL_CYCLONE_COLUMN_HEIGHTS_KM = (0.0, 5.0, 10.0, 15.0)
TROPICAL_CYCLONE_COLUMN_WIDTH_PX = WATER_OVERLAY_POINT_RADIUS_PX
TROPICAL_CYCLONE_FAR_MARKER_HALF_WIDTH_PX = 6.0
TROPICAL_CYCLONE_FAR_MARKER_HALF_HEIGHT_PX = 4.5
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


def _project_cyclone_column_points(
    lat_deg: float,
    lon_deg: float,
    *,
    viewer: ViewerData,
) -> tuple[_RenderPoint, ...] | None:
    points: list[_RenderPoint] = []
    for height_km in TROPICAL_CYCLONE_COLUMN_HEIGHTS_KM:
        point = _project_point_no_cutoff(
            lat_deg,
            lon_deg,
            viewer=viewer,
            height_m=float(height_km) * 1000.0,
        )
        if point is None:
            return None
        points.append(point)
    return tuple(points)


def _draw_column_line(
    painter: QPainter,
    points: tuple[_RenderPoint, ...],
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
    width_px: float,
) -> None:
    pen = QPen(QColor(*color_rgba), float(width_px), Qt.PenStyle.SolidLine)
    pen.setCosmetic(True)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    if len(points) < 2:
        return
    painter.drawPolyline(
        QPolygonF(
            [
                QPointF(*normalized_to_screen_xy(point.nx, point.ny, geometry))
                for point in points
            ]
        )
    )


def _draw_far_cyclone_marker(
    painter: QPainter,
    point: _RenderPoint,
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
) -> QPointF:
    pen = QPen(QColor(*color_rgba), 1.5, Qt.PenStyle.SolidLine)
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


def draw_tropical_cyclone_overlay(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    viewer: ViewerData,
    snapshot: TropicalCycloneSnapshot | None,
    when_utc: datetime | None,
    theme: ThemeStyle,
    opacity: float,
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
        column_points = None
        label_rgba = TROPICAL_CYCLONE_FAR_LABEL_RGBA
    else:
        column_points = _project_cyclone_column_points(
            observed.lat_deg,
            observed.lon_deg,
            viewer=viewer,
        )
        label_rgba = TROPICAL_CYCLONE_LABEL_RGBA

    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    if column_points is None:
        label_pos = _draw_far_cyclone_marker(
            painter,
            center_point,
            geometry=geometry,
            color_rgba=marker_rgba,
        )
    else:
        if len(column_points) >= 2:
            _draw_filled_cyclone_marker(
                painter,
                center_point,
                geometry=geometry,
                color_rgba=marker_rgba,
            )
            _draw_column_line(
                painter,
                column_points,
                geometry=geometry,
                color_rgba=marker_rgba,
                width_px=TROPICAL_CYCLONE_COLUMN_WIDTH_PX,
            )
            top_point = column_points[-1]
            screen_x, screen_y = normalized_to_screen_xy(
                top_point.nx,
                top_point.ny,
                geometry,
            )
            label_pos = QPointF(float(screen_x + 8.0), float(screen_y - 8.0))
        else:
            label_pos = None
    painter.setPen(QColor(*label_rgba))
    if label_pos is not None:
        painter.drawText(
            label_pos,
            projected_snapshot.storm_name,
        )
    painter.restore()
