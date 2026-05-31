from __future__ import annotations

import math
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

TROPICAL_CYCLONE_TARGET_HEIGHT_M = 0.0
TROPICAL_CYCLONE_MAX_DISTANCE_KM = 400.0
TROPICAL_CYCLONE_CYLINDER_HEIGHT_KM = 15.0
TROPICAL_CYCLONE_CYLINDER_LAYER_HEIGHTS_KM = (0.0, 5.0, 10.0, 15.0)
TROPICAL_CYCLONE_CYLINDER_SIDE_RADIUS_KM = 0.5
TROPICAL_CYCLONE_CYLINDER_RING_SAMPLES = 8
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


@dataclass(frozen=True, slots=True)
class _CycloneCylinderGeometry:
    top_center_point: _RenderPoint
    ring_points_by_height: tuple[tuple[_RenderPoint, ...], ...]
    radius_km: float
    height_km: float


def _point_key(point: QPointF) -> tuple[float, float]:
    return (float(point.x()), float(point.y()))


def _cross(
    origin: QPointF,
    a: QPointF,
    b: QPointF,
) -> float:
    return (
        (float(a.x()) - float(origin.x())) * (float(b.y()) - float(origin.y()))
        - (float(a.y()) - float(origin.y())) * (float(b.x()) - float(origin.x()))
    )


def _convex_hull(points: list[QPointF]) -> list[QPointF]:
    unique_points = sorted({(_point_key(point)): point for point in points}.values(), key=_point_key)
    if len(unique_points) <= 2:
        return unique_points

    lower: list[QPointF] = []
    for point in unique_points:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], point) <= 0.0:
            lower.pop()
        lower.append(point)

    upper: list[QPointF] = []
    for point in reversed(unique_points):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], point) <= 0.0:
            upper.pop()
        upper.append(point)

    return lower[:-1] + upper[:-1]


def _destination_point(
    lat_deg: float,
    lon_deg: float,
    *,
    distance_km: float,
    bearing_deg: float,
) -> tuple[float, float]:
    earth_radius_km = 6371.0
    angular_distance = float(distance_km) / earth_radius_km
    bearing = math.radians(float(bearing_deg))
    lat1 = math.radians(float(lat_deg))
    lon1 = math.radians(float(lon_deg))
    sin_lat1 = math.sin(lat1)
    cos_lat1 = math.cos(lat1)
    sin_ang = math.sin(angular_distance)
    cos_ang = math.cos(angular_distance)
    lat2 = math.asin(
        (sin_lat1 * cos_ang) + (cos_lat1 * sin_ang * math.cos(bearing))
    )
    lon2 = lon1 + math.atan2(
        math.sin(bearing) * sin_ang * cos_lat1,
        cos_ang - (sin_lat1 * math.sin(lat2)),
    )
    lon2_deg = ((math.degrees(lon2) + 180.0) % 360.0) - 180.0
    return math.degrees(lat2), lon2_deg


def _project_point(
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
    if float(projection.distance_km) > float(TROPICAL_CYCLONE_MAX_DISTANCE_KM):
        return None
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


def _cyclone_side_radius_km() -> float:
    return float(TROPICAL_CYCLONE_CYLINDER_SIDE_RADIUS_KM)


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


def _project_cyclone_cylinder(
    lat_deg: float,
    lon_deg: float,
    *,
    viewer: ViewerData,
    maxwind_kt: float | None,
) -> _CycloneCylinderGeometry | None:
    base_center_point = _project_point(
        lat_deg,
        lon_deg,
        viewer=viewer,
        height_m=TROPICAL_CYCLONE_TARGET_HEIGHT_M,
    )
    if base_center_point is None:
        return None
    height_km = float(TROPICAL_CYCLONE_CYLINDER_HEIGHT_KM)
    top_center_point: _RenderPoint | None = None
    ring_points_by_height: list[tuple[_RenderPoint, ...]] = []
    side_radius_km = _cyclone_side_radius_km()
    for layer_height_km in TROPICAL_CYCLONE_CYLINDER_LAYER_HEIGHTS_KM:
        layer_center_point = _project_point(
            lat_deg,
            lon_deg,
            viewer=viewer,
            height_m=float(layer_height_km) * 1000.0,
        )
        if layer_center_point is None:
            return None
        if math.isclose(float(layer_height_km), height_km):
            top_center_point = layer_center_point
        ring_radius_km = side_radius_km
        ring_points: list[_RenderPoint] = []
        for index in range(TROPICAL_CYCLONE_CYLINDER_RING_SAMPLES):
            bearing_deg = 360.0 * float(index) / float(TROPICAL_CYCLONE_CYLINDER_RING_SAMPLES)
            ring_lat_deg, ring_lon_deg = _destination_point(
                lat_deg,
                lon_deg,
                distance_km=ring_radius_km,
                bearing_deg=bearing_deg,
            )
            ring_point = _project_point_no_cutoff(
                ring_lat_deg,
                ring_lon_deg,
                viewer=viewer,
                height_m=float(layer_height_km) * 1000.0,
            )
            if ring_point is not None:
                ring_points.append(ring_point)
        if len(ring_points) < 3:
            return None
        ring_points_by_height.append(tuple(ring_points))
    if top_center_point is None:
        return None
    return _CycloneCylinderGeometry(
        top_center_point=top_center_point,
        ring_points_by_height=tuple(ring_points_by_height),
        radius_km=side_radius_km,
        height_km=height_km,
    )


def _draw_line(
    painter: QPainter,
    start: tuple[float, float],
    end: tuple[float, float],
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
    painter.drawLine(
        QPointF(*normalized_to_screen_xy(start[0], start[1], geometry)),
        QPointF(*normalized_to_screen_xy(end[0], end[1], geometry)),
    )


def _draw_cyclone_marker(
    painter: QPainter,
    cylinder: _CycloneCylinderGeometry,
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
) -> None:
    pen = QPen(QColor(*color_rgba), 1.5, Qt.PenStyle.SolidLine)
    pen.setCosmetic(True)
    painter.setPen(pen)
    if len(cylinder.ring_points_by_height) < 4:
        return
    rings = [
        [
            QPointF(*normalized_to_screen_xy(point.nx, point.ny, geometry))
            for point in ring_points
        ]
        for ring_points in cylinder.ring_points_by_height
    ]
    top_points = rings[-1]
    if len(top_points) < 3:
        return
    painter.setBrush(QColor(*color_rgba))
    painter.drawPolygon(QPolygonF(top_points))
    for lower_ring, upper_ring in zip(rings[:-1], rings[1:], strict=False):
        for lower_point, upper_point in zip(lower_ring, upper_ring, strict=False):
            painter.drawLine(lower_point, upper_point)


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
    if not enabled or snapshot is None:
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
        cylinder = None
        line_top_point = None
        label_rgba = TROPICAL_CYCLONE_FAR_LABEL_RGBA
    else:
        cylinder = _project_cyclone_cylinder(
            observed.lat_deg,
            observed.lon_deg,
            viewer=viewer,
            maxwind_kt=observed.maxwind_kt,
        )
        if cylinder is None:
            line_top_point = _project_point_no_cutoff(
                observed.lat_deg,
                observed.lon_deg,
                viewer=viewer,
                height_m=float(TROPICAL_CYCLONE_CYLINDER_HEIGHT_KM) * 1000.0,
            )
            label_rgba = TROPICAL_CYCLONE_LABEL_RGBA
        else:
            line_top_point = None
            label_rgba = TROPICAL_CYCLONE_LABEL_RGBA

    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    if cylinder is None:
        if line_top_point is not None:
            _draw_line(
                painter,
                (center_point.nx, center_point.ny),
                (line_top_point.nx, line_top_point.ny),
                geometry=geometry,
                color_rgba=marker_rgba,
                width_px=1.0,
            )
            screen_x, screen_y = normalized_to_screen_xy(
                line_top_point.nx,
                line_top_point.ny,
                geometry,
            )
            label_pos = QPointF(float(screen_x + 8.0), float(screen_y - 8.0))
        else:
            label_pos = _draw_far_cyclone_marker(
                painter,
                center_point,
                geometry=geometry,
                color_rgba=marker_rgba,
            )
    else:
        _draw_cyclone_marker(
            painter,
            cylinder,
            geometry=geometry,
            color_rgba=marker_rgba,
        )
        _draw_filled_cyclone_marker(
            painter,
            center_point,
            geometry=geometry,
            color_rgba=marker_rgba,
        )
        screen_x, screen_y = normalized_to_screen_xy(
            cylinder.top_center_point.nx,
            cylinder.top_center_point.ny,
            geometry,
        )
        label_pos = QPointF(float(screen_x + 14.0), float(screen_y - 14.0))
    painter.setPen(QColor(*label_rgba))
    if label_pos is not None:
        painter.drawText(
            label_pos,
            projected_snapshot.storm_name,
        )
    painter.restore()
