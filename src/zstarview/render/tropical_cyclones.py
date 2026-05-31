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
TROPICAL_CYCLONE_MARKER_BASE_HEIGHT_M = 15_000.0
TROPICAL_CYCLONE_MAX_DISTANCE_KM = 400.0
TROPICAL_CYCLONE_BASE_RADIUS_KM_MIN = 6.0
TROPICAL_CYCLONE_BASE_RADIUS_KM_MAX = 24.0
TROPICAL_CYCLONE_BASE_RADIUS_KM_PER_KT = 0.12
TROPICAL_CYCLONE_BASE_RING_SAMPLES = 16
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
class _CycloneConeGeometry:
    tip_point: _RenderPoint
    base_center_point: _RenderPoint
    base_ring_points: tuple[_RenderPoint, ...]
    base_radius_km: float


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


def _project_cyclone_cone(
    lat_deg: float,
    lon_deg: float,
    *,
    viewer: ViewerData,
    maxwind_kt: float | None,
) -> _CycloneConeGeometry | None:
    tip_point = _project_point(
        lat_deg,
        lon_deg,
        viewer=viewer,
        height_m=TROPICAL_CYCLONE_TARGET_HEIGHT_M,
    )
    if tip_point is None:
        return None
    base_center_point = _project_point(
        lat_deg,
        lon_deg,
        viewer=viewer,
        height_m=TROPICAL_CYCLONE_MARKER_BASE_HEIGHT_M,
    )
    if base_center_point is None:
        return None
    base_radius_km = max(
        TROPICAL_CYCLONE_BASE_RADIUS_KM_MIN,
        min(
            TROPICAL_CYCLONE_BASE_RADIUS_KM_MAX,
            float(maxwind_kt) * TROPICAL_CYCLONE_BASE_RADIUS_KM_PER_KT
            if maxwind_kt is not None
            else 35.0 * TROPICAL_CYCLONE_BASE_RADIUS_KM_PER_KT,
        ),
    )
    base_ring_points: list[_RenderPoint] = []
    for index in range(TROPICAL_CYCLONE_BASE_RING_SAMPLES):
        bearing_deg = 360.0 * float(index) / float(TROPICAL_CYCLONE_BASE_RING_SAMPLES)
        ring_lat_deg, ring_lon_deg = _destination_point(
            lat_deg,
            lon_deg,
            distance_km=base_radius_km,
            bearing_deg=bearing_deg,
        )
        ring_point = _project_point_no_cutoff(
            ring_lat_deg,
            ring_lon_deg,
            viewer=viewer,
            height_m=TROPICAL_CYCLONE_MARKER_BASE_HEIGHT_M,
        )
        if ring_point is not None:
            base_ring_points.append(ring_point)
    if len(base_ring_points) < 3:
        return None
    return _CycloneConeGeometry(
        tip_point=tip_point,
        base_center_point=base_center_point,
        base_ring_points=tuple(base_ring_points),
        base_radius_km=base_radius_km,
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
    cone: _CycloneConeGeometry,
    *,
    geometry: ScreenGeometry,
    color_rgba: tuple[int, int, int, int],
) -> None:
    pen = QPen(QColor(*color_rgba), 1.5, Qt.PenStyle.SolidLine)
    pen.setCosmetic(True)
    painter.setPen(pen)
    projected_points = [
        QPointF(*normalized_to_screen_xy(point.nx, point.ny, geometry))
        for point in cone.base_ring_points
    ]
    projected_points.append(QPointF(*normalized_to_screen_xy(cone.tip_point.nx, cone.tip_point.ny, geometry)))
    hull_points = _convex_hull(projected_points)
    if len(hull_points) < 3:
        return
    painter.setBrush(QColor(*color_rgba))
    painter.drawPolygon(QPolygonF(hull_points))


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
        cone = None
    else:
        cone = _project_cyclone_cone(
            observed.lat_deg,
            observed.lon_deg,
            viewer=viewer,
            maxwind_kt=observed.maxwind_kt,
        )
        if cone is None:
            return

    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    if cone is None:
        label_pos = _draw_far_cyclone_marker(
            painter,
            center_point,
            geometry=geometry,
            color_rgba=marker_rgba,
        )
    else:
        _draw_cyclone_marker(
            painter,
            cone,
            geometry=geometry,
            color_rgba=marker_rgba,
        )
        screen_x, screen_y = normalized_to_screen_xy(cone.tip_point.nx, cone.tip_point.ny, geometry)
        label_pos = QPointF(
            float(screen_x),
            float(screen_y + 14.0),
    )
    if label_pos is not None:
        painter.setPen(
            QColor(
                *(TROPICAL_CYCLONE_FAR_LABEL_RGBA if cone is None else TROPICAL_CYCLONE_LABEL_RGBA)
            )
        )
        painter.drawText(
            label_pos,
            projected_snapshot.storm_name,
        )
    painter.restore()
