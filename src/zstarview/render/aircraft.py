from dataclasses import dataclass
from typing import Any, Dict, List, Optional
import math

import astropy.time
from astropy import units as u
from astropy.coordinates import EarthLocation
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..aircraft import project_aircraft_snapshots
from ..aircraft.types import AircraftOverlayPoint
from ..aircraft_constants import (
    AIRCRAFT_FADE_START_SECONDS,
    AIRCRAFT_OVERLAY_LINE_COLOR_RGB,
)
from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import ThemeStyle
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy
from .text import resolve_text_style

_AIRCRAFT_CALLSIGN_MAX_DISTANCE_KM = 10.0
_AIRCRAFT_MAX_DRAW_DISTANCE_KM = 50.0
_AIRCRAFT_TRAIL_GAP_PX = 160.0
_AIRCRAFT_RIBBON_FILL_ALPHA_SCALE = 0.72
_AIRCRAFT_RIBBON_MIN_FULL_WIDTH_M = 120.0


@dataclass(frozen=True, slots=True)
class _ObserverProjectionState:
    obs_x: float
    obs_y: float
    obs_z: float
    sin_lat: float
    cos_lat: float
    sin_lon: float
    cos_lon: float


def draw_aircraft_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData | None = None,
    aircraft_snapshots: object | None = None,
    *,
    time_obj: astropy.time.Time | None = None,
    opacity: float = 1.0,
    line_width_scale: float = 1.0,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    theme: ThemeStyle,
) -> None:
    if viewer_data is None or time_obj is None:
        return
    view_center = viewer_data.view_center
    edge_fov_deg = float(viewer_data.edge_fov_deg)
    content_fov_deg = float(viewer_data.content_fov_deg)
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if layer_opacity <= 0.0:
        return

    aircraft_points = _project_aircraft_snapshots(
        aircraft_snapshots,
        viewer_data=viewer_data,
        time_obj=time_obj,
    )
    if not aircraft_points:
        return

    painter.save()
    observer_location = EarthLocation(
        lat=float(viewer_data.lat_deg) * u.deg,
        lon=float(viewer_data.lon_deg) * u.deg,
        height=float(viewer_data.observer_height_m) * u.m,
    )
    observer_xyz = observer_location.to_geocentric()
    observer_state = _ObserverProjectionState(
        obs_x=float(observer_xyz[0].to_value(u.m)),
        obs_y=float(observer_xyz[1].to_value(u.m)),
        obs_z=float(observer_xyz[2].to_value(u.m)),
        sin_lat=math.sin(math.radians(float(viewer_data.lat_deg))),
        cos_lat=math.cos(math.radians(float(viewer_data.lat_deg))),
        sin_lon=math.sin(math.radians(float(viewer_data.lon_deg))),
        cos_lon=math.cos(math.radians(float(viewer_data.lon_deg))),
    )
    width_scale = max(1.0, float(line_width_scale))
    line_color = QColor(*AIRCRAFT_OVERLAY_LINE_COLOR_RGB, 255)
    label_style = resolve_text_style(theme, painter.font(), opacity=layer_opacity)
    label_color = QColor(*AIRCRAFT_OVERLAY_LINE_COLOR_RGB)
    label_color.setAlpha(label_style.text_color.alpha())
    label_style = type(label_style)(
        font=label_style.font,
        text_color=label_color,
        outline_color=label_style.outline_color,
        outline_width=label_style.outline_width,
    )
    line_pen = QPen(line_color, 1.0, Qt.PenStyle.SolidLine)
    line_pen.setCosmetic(True)
    line_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    line_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(line_pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    for point in aircraft_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        distance_km = float(point.distance_km)
        if distance_km > _AIRCRAFT_MAX_DRAW_DISTANCE_KM:
            continue
        if alt <= 0.0 or not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        min_line_alpha = int(round(30.0 * layer_opacity))
        line_alpha = max(
            min_line_alpha,
            min(255, int(round(255.0 * float(point.alpha_scale) * layer_opacity))),
        )
        line_color.setAlpha(line_alpha)
        ribbon_width_px = _aircraft_line_width_px(distance_km, width_scale=width_scale)
        trail_points = tuple(
            (float(sample_alt_deg), float(sample_az_deg))
            for sample_alt_deg, sample_az_deg in point.trail_alt_az_points
        )
        if any(
            is_in_fov(sample_alt_deg, sample_az_deg, view_center, fov_deg=content_fov_deg)
            for sample_alt_deg, sample_az_deg in trail_points
        ):
            trail_screen_points = _project_aircraft_trail_screen_points(
                trail_points,
                geometry=geometry,
                view_center=view_center,
                edge_fov_deg=edge_fov_deg,
            )
            trail_geodetic_points = tuple(
                tuple(float(value) for value in geodetic_point)
                for geodetic_point in getattr(point, "trail_geodetic_points", ())
            )
            ribbon_polygons = _aircraft_ribbon_polygons(
                trail_alt_az_points=trail_points,
                trail_geodetic_points=trail_geodetic_points,
                geometry=geometry,
                view_center=view_center,
                edge_fov_deg=edge_fov_deg,
                full_width_m=_aircraft_full_width_m(point=point),
                observer_state=observer_state,
            )
            if ribbon_polygons:
                fill_alpha = max(
                    1, min(255, int(round(float(line_alpha) * _AIRCRAFT_RIBBON_FILL_ALPHA_SCALE)))
                )
                fill_color = QColor(*AIRCRAFT_OVERLAY_LINE_COLOR_RGB, fill_alpha)
                outline_color = QColor(*AIRCRAFT_OVERLAY_LINE_COLOR_RGB, line_alpha)
                outline_pen = QPen(outline_color, max(1.0, ribbon_width_px * 0.32), Qt.PenStyle.SolidLine)
                outline_pen.setCosmetic(True)
                outline_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                outline_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(fill_color)
                for ribbon_polygon in ribbon_polygons:
                    painter.drawPolygon(ribbon_polygon)
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.setPen(outline_pen)
                for ribbon_polygon in ribbon_polygons:
                    painter.drawPolygon(ribbon_polygon)
            else:
                line_pen.setColor(line_color)
                line_pen.setWidthF(ribbon_width_px)
                painter.setPen(line_pen)
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.drawPolyline(QPolygonF(trail_screen_points))
        nx, ny = altaz_to_normalized_xy(alt, az, view_center, edge_fov_deg=edge_fov_deg)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pos = QPointF(float(px), float(py))
        callsign = (point.callsign or "").strip()
        if (
            callsign
            and float(point.distance_km) <= _AIRCRAFT_CALLSIGN_MAX_DISTANCE_KM
            and float(point.age_seconds) <= AIRCRAFT_FADE_START_SECONDS
            and label_candidates is not None
        ):
            label_candidates.append(
                {
                    "text": callsign,
                    "pos": QPointF(pos.x() + 8.0, pos.y() - 6.0),
                    "style": label_style,
                    "priority": 45,
                    "hide_on_overlap": True,
                }
            )
    painter.restore()


def _project_aircraft_snapshots(
    aircraft_snapshots: object | None,
    *,
    viewer_data: ViewerData | None = None,
    time_obj: astropy.time.Time | None = None,
) -> list[AircraftOverlayPoint]:
    if aircraft_snapshots is None:
        return []
    if viewer_data is None or time_obj is None:
        return []
    return project_aircraft_snapshots(
        aircraft_snapshots,
        observer_lat=float(viewer_data.lat_deg),
        observer_lon=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        time_obj=time_obj,
    )


def _aircraft_line_width_px(distance_km: float, *, width_scale: float = 1.0) -> float:
    d = max(0.0, float(distance_km))
    scale = max(1.0, float(width_scale))
    aircraft_scale = 2.4 * scale
    if d <= 1.0:
        return 3.0 * aircraft_scale
    if d <= 3.0:
        return 2.2 * aircraft_scale
    if d <= 5.0:
        return 1.6 * aircraft_scale
    if d <= 10.0:
        return 1.0 * aircraft_scale
    if d <= 20.0:
        return 0.8 * aircraft_scale
    return 0.6 * aircraft_scale


def _aircraft_ribbon_polygons(
    *,
    trail_alt_az_points: tuple[tuple[float, float], ...],
    trail_geodetic_points: tuple[tuple[float, float, float], ...],
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    full_width_m: float,
    observer_state: _ObserverProjectionState,
) -> tuple[QPolygonF, ...]:
    if len(trail_alt_az_points) < 2:
        return ()
    if len(trail_geodetic_points) != len(trail_alt_az_points) or len(trail_geodetic_points) < 2:
        return ()

    center_screen_points = _project_aircraft_trail_screen_points(
        trail_alt_az_points,
        geometry=geometry,
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
    )
    run_indices_list = _split_trail_run_indices(center_screen_points, gap_px=_AIRCRAFT_TRAIL_GAP_PX)
    if not run_indices_list:
        return ()

    polygons: list[QPolygonF] = []
    for run_indices in run_indices_list:
        if len(run_indices) < 2:
            continue
        left_points: list[QPointF] = []
        right_points: list[QPointF] = []
        for sample_index in run_indices:
            sample_lat_deg, sample_lon_deg, sample_alt_m = trail_geodetic_points[sample_index]
            direction_east_m, direction_north_m = _aircraft_trail_direction_m(
                trail_geodetic_points,
                sample_index,
            )
            if direction_east_m == 0.0 and direction_north_m == 0.0:
                continue
            perp_east_m = -direction_north_m
            perp_north_m = direction_east_m
            norm = (perp_east_m * perp_east_m + perp_north_m * perp_north_m) ** 0.5
            if norm <= 0.0:
                continue
            perp_east_m /= norm
            perp_north_m /= norm
            half_width_m = max(_AIRCRAFT_RIBBON_MIN_FULL_WIDTH_M * 0.5, float(full_width_m) * 0.5)
            left_lat_deg, left_lon_deg = _offset_latlon_by_local_m(
                sample_lat_deg,
                sample_lon_deg,
                east_m=perp_east_m * half_width_m,
                north_m=perp_north_m * half_width_m,
            )
            right_lat_deg, right_lon_deg = _offset_latlon_by_local_m(
                sample_lat_deg,
                sample_lon_deg,
                east_m=-perp_east_m * half_width_m,
                north_m=-perp_north_m * half_width_m,
            )
            left_alt_deg, left_az_deg, _ = _project_geodetic_to_altaz(
                left_lat_deg,
                left_lon_deg,
                sample_alt_m,
                observer_state=observer_state,
            )
            right_alt_deg, right_az_deg, _ = _project_geodetic_to_altaz(
                right_lat_deg,
                right_lon_deg,
                sample_alt_m,
                observer_state=observer_state,
            )
            left_nx, left_ny = altaz_to_normalized_xy(
                left_alt_deg,
                left_az_deg,
                view_center,
                edge_fov_deg=edge_fov_deg,
            )
            right_nx, right_ny = altaz_to_normalized_xy(
                right_alt_deg,
                right_az_deg,
                view_center,
                edge_fov_deg=edge_fov_deg,
            )
            left_px, left_py = normalized_to_screen_xy(left_nx, left_ny, geometry)
            right_px, right_py = normalized_to_screen_xy(right_nx, right_ny, geometry)
            left_points.append(QPointF(float(left_px), float(left_py)))
            right_points.append(QPointF(float(right_px), float(right_py)))
        if len(left_points) >= 2 and len(right_points) >= 2:
            ribbon_points = [*left_points, *reversed(right_points)]
            ribbon_polygon = QPolygonF(ribbon_points)
            if not ribbon_polygon.isEmpty():
                polygons.append(ribbon_polygon)

    if polygons:
        return tuple(polygons)

    fallback_screen_points = center_screen_points
    # Fall back to a compact marker when the trail collapses to a point.
    center = fallback_screen_points[len(fallback_screen_points) // 2]
    radius = max(1.0, float(full_width_m) * 0.01)
    fallback_polygon = QPolygonF(
        [
            QPointF(center.x() - radius, center.y() - radius),
            QPointF(center.x() + radius, center.y() - radius),
            QPointF(center.x() + radius, center.y() + radius),
            QPointF(center.x() - radius, center.y() + radius),
        ]
    )
    if fallback_polygon.isEmpty():
        return ()
    return (fallback_polygon,)


def _project_aircraft_trail_screen_points(
    trail_alt_az_points: tuple[tuple[float, float], ...],
    *,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    edge_fov_deg: float,
) -> list[QPointF]:
    screen_points: list[QPointF] = []
    for sample_alt_deg, sample_az_deg in trail_alt_az_points:
        sample_nx, sample_ny = altaz_to_normalized_xy(
            float(sample_alt_deg),
            float(sample_az_deg),
            view_center,
            edge_fov_deg=edge_fov_deg,
        )
        sample_px, sample_py = normalized_to_screen_xy(sample_nx, sample_ny, geometry)
        screen_points.append(QPointF(float(sample_px), float(sample_py)))
    return screen_points


def _split_trail_run_indices(
    screen_points: list[QPointF],
    *,
    gap_px: float,
) -> tuple[tuple[int, ...], ...]:
    if len(screen_points) < 2:
        return ()
    runs: list[list[int]] = []
    current_run: list[int] = [0]
    prev_point = screen_points[0]
    for index, point in enumerate(screen_points[1:], start=1):
        delta_x = float(point.x()) - float(prev_point.x())
        delta_y = float(point.y()) - float(prev_point.y())
        if (delta_x * delta_x + delta_y * delta_y) ** 0.5 > float(gap_px):
            if len(current_run) >= 2:
                runs.append(current_run)
            current_run = [index]
        else:
            current_run.append(index)
        prev_point = point
    if len(current_run) >= 2:
        runs.append(current_run)
    return tuple(tuple(run) for run in runs)


def _aircraft_trail_direction_m(
    trail_geodetic_points: tuple[tuple[float, float, float], ...],
    index: int,
) -> tuple[float, float]:
    if index < 0 or index >= len(trail_geodetic_points):
        return (0.0, 0.0)
    prev_index = max(0, index - 1)
    next_index = min(len(trail_geodetic_points) - 1, index + 1)
    if prev_index == next_index:
        return (0.0, 0.0)
    prev_lat_deg, prev_lon_deg, _ = trail_geodetic_points[prev_index]
    next_lat_deg, next_lon_deg, _ = trail_geodetic_points[next_index]
    base_lat_deg, base_lon_deg, _ = trail_geodetic_points[index]
    prev_east_m, prev_north_m = _geodetic_delta_to_local_m(
        base_lat_deg,
        base_lon_deg,
        prev_lat_deg,
        prev_lon_deg,
    )
    next_east_m, next_north_m = _geodetic_delta_to_local_m(
        base_lat_deg,
        base_lon_deg,
        next_lat_deg,
        next_lon_deg,
    )
    east_m = next_east_m - prev_east_m
    north_m = next_north_m - prev_north_m
    norm = (east_m * east_m + north_m * north_m) ** 0.5
    if norm <= 0.0:
        east_m, north_m = next_east_m, next_north_m
        norm = (east_m * east_m + north_m * north_m) ** 0.5
    if norm <= 0.0:
        east_m, north_m = -prev_east_m, -prev_north_m
        norm = (east_m * east_m + north_m * north_m) ** 0.5
    if norm <= 0.0:
        return (0.0, 0.0)
    return (east_m / norm, north_m / norm)


def _aircraft_full_width_m(point: AircraftOverlayPoint) -> float:
    trail_points = tuple(point.trail_geodetic_points)
    if len(trail_points) < 2:
        return _AIRCRAFT_RIBBON_MIN_FULL_WIDTH_M
    center_index = len(trail_points) // 2
    current = trail_points[center_index]
    if center_index + 1 < len(trail_points):
        neighbour = trail_points[center_index + 1]
    else:
        neighbour = trail_points[center_index - 1]
    east_m, north_m = _geodetic_delta_to_local_m(
        current[0],
        current[1],
        neighbour[0],
        neighbour[1],
    )
    full_width_m = (east_m * east_m + north_m * north_m) ** 0.5
    return max(_AIRCRAFT_RIBBON_MIN_FULL_WIDTH_M, full_width_m)


def _geodetic_delta_to_local_m(
    base_lat_deg: float,
    base_lon_deg: float,
    target_lat_deg: float,
    target_lon_deg: float,
) -> tuple[float, float]:
    lat_scale_m = 111_320.0
    lon_scale_m = max(1.0, lat_scale_m * abs(math.cos(math.radians(float(base_lat_deg)))))
    north_m = (float(target_lat_deg) - float(base_lat_deg)) * lat_scale_m
    east_m = (float(target_lon_deg) - float(base_lon_deg)) * lon_scale_m
    return east_m, north_m


def _offset_latlon_by_local_m(
    lat_deg: float,
    lon_deg: float,
    *,
    east_m: float,
    north_m: float,
) -> tuple[float, float]:
    lat_scale_m = 111_320.0
    lon_scale_m = max(1.0, lat_scale_m * abs(math.cos(math.radians(float(lat_deg)))))
    offset_lat = float(lat_deg) + (float(north_m) / lat_scale_m)
    offset_lon = float(lon_deg) + (float(east_m) / lon_scale_m)
    return offset_lat, offset_lon


def _project_geodetic_to_altaz(
    target_lat_deg: float,
    target_lon_deg: float,
    target_alt_m: float,
    *,
    observer_state: _ObserverProjectionState,
) -> tuple[float, float, float]:
    from astropy import units as u
    from astropy.coordinates import EarthLocation

    target_location = EarthLocation(
        lat=float(target_lat_deg) * u.deg,
        lon=float(target_lon_deg) * u.deg,
        height=float(target_alt_m) * u.m,
    )
    target_xyz = target_location.to_geocentric()
    dx = float(target_xyz[0].to_value(u.m)) - float(observer_state.obs_x)
    dy = float(target_xyz[1].to_value(u.m)) - float(observer_state.obs_y)
    dz = float(target_xyz[2].to_value(u.m)) - float(observer_state.obs_z)

    east_m = (-float(observer_state.sin_lon) * dx) + (float(observer_state.cos_lon) * dy)
    north_m = (
        (-float(observer_state.sin_lat) * float(observer_state.cos_lon) * dx)
        - (float(observer_state.sin_lat) * float(observer_state.sin_lon) * dy)
        + (float(observer_state.cos_lat) * dz)
    )
    up_m = (
        (float(observer_state.cos_lat) * float(observer_state.cos_lon) * dx)
        + (float(observer_state.cos_lat) * float(observer_state.sin_lon) * dy)
        + (float(observer_state.sin_lat) * dz)
    )

    horizontal_m = math.hypot(east_m, north_m)
    distance_km = math.sqrt((horizontal_m * horizontal_m) + (up_m * up_m)) / 1000.0
    alt_deg = math.degrees(math.atan2(up_m, horizontal_m))
    az_deg = math.degrees(math.atan2(east_m, north_m)) % 360.0
    return alt_deg, az_deg, distance_km
