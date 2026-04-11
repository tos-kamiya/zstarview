from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPainterPathStroker, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import EARTH_GUIDE_LAND_FILE, TERRAIN_HORIZON_LINE_COLOR
from ..types import ScreenGeometry
from .geometry import normalized_to_screen_xy
from .terrain import terrain_horizon_line_alpha


EARTH_MEAN_RADIUS_M = 6_371_008.8

EARTH_GUIDE_DEAD_ZONE_MIN_KM = 20.0
EARTH_GUIDE_DEAD_ZONE_MAX_KM = 80.0
EARTH_GUIDE_DEAD_ZONE_SCALE = 0.25
EARTH_GUIDE_HORIZON_MARGIN_DEG = 1.0
EARTH_GUIDE_FILL_WIDTH = 3.6
EARTH_GUIDE_FILL_ALPHA_SCALE = 0.22
EARTH_GUIDE_FOREGROUND_WIDTH = 1.5


@dataclass(frozen=True)
class EarthGuideRing:
    source_name: str
    label_name: str | None
    points_lonlat_deg: np.ndarray
    points_xyz: np.ndarray


def _lonlat_to_unit_xyz(lon_deg: float, lat_deg: float) -> tuple[float, float, float]:
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)
    cos_lat = math.cos(lat_rad)
    return (
        cos_lat * math.cos(lon_rad),
        cos_lat * math.sin(lon_rad),
        math.sin(lat_rad),
    )


def _observer_basis(
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lat_rad = math.radians(observer_lat_deg)
    lon_rad = math.radians(observer_lon_deg)
    up = np.array(
        (
            math.cos(lat_rad) * math.cos(lon_rad),
            math.cos(lat_rad) * math.sin(lon_rad),
            math.sin(lat_rad),
        ),
        dtype=np.float64,
    )
    east = np.array((-math.sin(lon_rad), math.cos(lon_rad), 0.0), dtype=np.float64)
    north = np.array(
        (
            -math.sin(lat_rad) * math.cos(lon_rad),
            -math.sin(lat_rad) * math.sin(lon_rad),
            math.cos(lat_rad),
        ),
        dtype=np.float64,
    )
    radius_scale = 1.0 + (max(0.0, float(observer_height_m)) / EARTH_MEAN_RADIUS_M)
    origin = up * radius_scale
    return origin, east, north, up


def _observer_dead_zone_km(observer_height_m: float) -> float:
    height_km = max(0.0, float(observer_height_m)) / 1000.0
    horizon_km = math.sqrt(max(0.0, 2.0 * (EARTH_MEAN_RADIUS_M / 1000.0) * height_km))
    return max(
        EARTH_GUIDE_DEAD_ZONE_MIN_KM,
        min(EARTH_GUIDE_DEAD_ZONE_SCALE * horizon_km, EARTH_GUIDE_DEAD_ZONE_MAX_KM),
    )


def _observer_horizon_dip_deg(observer_height_m: float) -> float:
    height_m = max(0.0, float(observer_height_m))
    if height_m <= 0.0:
        return 0.0
    ratio = EARTH_MEAN_RADIUS_M / (EARTH_MEAN_RADIUS_M + height_m)
    return math.degrees(math.acos(max(-1.0, min(1.0, ratio))))


def _observer_visible_altitude_limit_deg(
    az_deg: float,
    observer_height_m: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> float:
    terrain_limit = _effective_visible_altitude_limit_deg(az_deg, terrain_profile_altaz)
    horizon_limit = -(_observer_horizon_dip_deg(observer_height_m) + EARTH_GUIDE_HORIZON_MARGIN_DEG)
    return min(terrain_limit, horizon_limit)


def _surface_distance_km(point_xyz: np.ndarray, observer_up: np.ndarray) -> float:
    dot = float(np.dot(point_xyz, observer_up))
    dot = max(-1.0, min(1.0, dot))
    return (EARTH_MEAN_RADIUS_M / 1000.0) * math.acos(dot)


def _fragment_fill_path(fragment: list[tuple[float, float]]) -> QPainterPath:
    path = QPainterPath()
    if not fragment:
        return path
    first_x, first_y = fragment[0]
    path.moveTo(QPointF(first_x, first_y))
    for x, y in fragment[1:]:
        path.lineTo(QPointF(x, y))
    stroker = QPainterPathStroker()
    stroker.setCapStyle(Qt.PenCapStyle.RoundCap)
    stroker.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    stroker.setWidth(EARTH_GUIDE_FILL_WIDTH)
    return stroker.createStroke(path)


def _project_point_altaz(
    point_xyz: np.ndarray,
    *,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
) -> tuple[float, float] | None:
    ray = point_xyz - origin
    norm = float(np.linalg.norm(ray))
    if norm <= 1.0e-12:
        return None
    ray /= norm
    east_c = float(np.dot(ray, east))
    north_c = float(np.dot(ray, north))
    up_c = float(np.dot(ray, up))
    alt_deg = math.degrees(math.asin(max(-1.0, min(1.0, up_c))))
    az_deg = math.degrees(math.atan2(east_c, north_c)) % 360.0
    return alt_deg, az_deg


def _interpolate_horizon_altitude_deg(
    az_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> float:
    if not terrain_profile_altaz:
        return 0.0
    if len(terrain_profile_altaz) == 1:
        return float(terrain_profile_altaz[0][0])
    pts = sorted((float(az) % 360.0, float(alt)) for alt, az in terrain_profile_altaz)
    target = float(az_deg) % 360.0
    wrapped = pts + [(pts[0][0] + 360.0, pts[0][1])]
    prev_az, prev_alt = wrapped[0]
    if target < prev_az:
        target += 360.0
    for next_az, next_alt in wrapped[1:]:
        if target <= next_az:
            span = next_az - prev_az
            if span <= 1.0e-12:
                return prev_alt
            t = (target - prev_az) / span
            return prev_alt + ((next_alt - prev_alt) * t)
        prev_az, prev_alt = next_az, next_alt
    return wrapped[-1][1]


def _effective_visible_altitude_limit_deg(
    az_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> float:
    if not terrain_profile_altaz:
        return 0.0
    terrain_limit = _interpolate_horizon_altitude_deg(az_deg, terrain_profile_altaz)
    return terrain_limit


def _interpolate_xyz_on_sphere(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    point = ((1.0 - t) * a) + (t * b)
    norm = float(np.linalg.norm(point))
    if norm <= 1.0e-12:
        return a
    return point / norm


def _project_xyz_to_screen_point(
    xyz: np.ndarray,
    *,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    observer_height_m: float,
    dead_zone_km: float,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    content_fov_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> tuple[tuple[float, float], bool] | None:
    altaz = _project_point_altaz(xyz, origin=origin, east=east, north=north, up=up)
    if altaz is None:
        return None
    alt_deg, az_deg = altaz
    if not is_in_fov(alt_deg, az_deg, view_center, fov_deg=content_fov_deg):
        return None
    nx, ny = altaz_to_normalized_xy(alt_deg, az_deg, view_center)
    screen_xy = normalized_to_screen_xy(nx, ny, geometry)
    limit_deg = _observer_visible_altitude_limit_deg(
        az_deg,
        observer_height_m,
        terrain_profile_altaz,
    )
    visible = _surface_distance_km(xyz, up) >= dead_zone_km and alt_deg <= limit_deg
    return (screen_xy, visible)


def _segment_screen_fragments(
    xyz0: np.ndarray,
    xyz1: np.ndarray,
    *,
    depth: int,
    max_depth: int,
    threshold_px: float,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    observer_height_m: float,
    dead_zone_km: float,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    content_fov_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> list[list[tuple[float, float]]]:
    start = _project_xyz_to_screen_point(
        xyz0,
        origin=origin,
        east=east,
        north=north,
        up=up,
        observer_height_m=observer_height_m,
        dead_zone_km=dead_zone_km,
        geometry=geometry,
        view_center=view_center,
        content_fov_deg=content_fov_deg,
        terrain_profile_altaz=terrain_profile_altaz,
    )
    end = _project_xyz_to_screen_point(
        xyz1,
        origin=origin,
        east=east,
        north=north,
        up=up,
        observer_height_m=observer_height_m,
        dead_zone_km=dead_zone_km,
        geometry=geometry,
        view_center=view_center,
        content_fov_deg=content_fov_deg,
        terrain_profile_altaz=terrain_profile_altaz,
    )
    if start is None or end is None:
        return []

    start_xy, start_visible = start
    end_xy, end_visible = end
    midpoint_xyz = _interpolate_xyz_on_sphere(xyz0, xyz1, 0.5)
    midpoint = _project_xyz_to_screen_point(
        midpoint_xyz,
        origin=origin,
        east=east,
        north=north,
        up=up,
        observer_height_m=observer_height_m,
        dead_zone_km=dead_zone_km,
        geometry=geometry,
        view_center=view_center,
        content_fov_deg=content_fov_deg,
        terrain_profile_altaz=terrain_profile_altaz,
    )
    if midpoint is None:
        return []

    midpoint_xy, midpoint_visible = midpoint
    screen_span = math.hypot(end_xy[0] - start_xy[0], end_xy[1] - start_xy[1])
    midpoint_distance = float(np.linalg.norm(midpoint_xyz - origin))
    effective_threshold_px = threshold_px * min(3.0, 1.0 + max(0.0, midpoint_distance - 0.35) * 1.5)
    should_split = (
        screen_span > effective_threshold_px
        or start_visible != end_visible
        or midpoint_visible != start_visible
        or midpoint_visible != end_visible
    )
    if should_split and depth < max_depth:
        left = _segment_screen_fragments(
            xyz0,
            midpoint_xyz,
            depth=depth + 1,
            max_depth=max_depth,
            threshold_px=threshold_px,
            origin=origin,
            east=east,
            north=north,
            up=up,
            observer_height_m=observer_height_m,
            dead_zone_km=dead_zone_km,
            geometry=geometry,
            view_center=view_center,
            content_fov_deg=content_fov_deg,
            terrain_profile_altaz=terrain_profile_altaz,
        )
        right = _segment_screen_fragments(
            midpoint_xyz,
            xyz1,
            depth=depth + 1,
            max_depth=max_depth,
            threshold_px=threshold_px,
            origin=origin,
            east=east,
            north=north,
            up=up,
            observer_height_m=observer_height_m,
            dead_zone_km=dead_zone_km,
            geometry=geometry,
            view_center=view_center,
            content_fov_deg=content_fov_deg,
            terrain_profile_altaz=terrain_profile_altaz,
        )
        if len(left) == 1 and len(right) == 1:
            left_fragment = left[0]
            right_fragment = right[0]
            if left_fragment and right_fragment and left_fragment[-1] == right_fragment[0]:
                return [left_fragment[:-1] + right_fragment]
        if left and right and left[-1][-1] == right[0][0]:
            return left[:-1] + right
        return left + right

    if start_visible and end_visible:
        return [[start_xy, end_xy]]
    if start_visible and midpoint_visible:
        return [[start_xy, midpoint_xy]]
    if midpoint_visible and end_visible:
        return [[midpoint_xy, end_xy]]
    return []


def _ring_fragments_altaz(
    ring: EarthGuideRing,
    *,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    observer_height_m: float,
    dead_zone_km: float,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    content_fov_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> list[list[tuple[float, float]]]:
    xyz_points = ring.points_xyz
    fragments: list[list[tuple[float, float]]] = []
    count = len(xyz_points)
    for index in range(count):
        next_index = (index + 1) % count
        sampled = _segment_screen_fragments(
            xyz_points[index],
            xyz_points[next_index],
            depth=0,
            max_depth=12,
            threshold_px=24.0,
            origin=origin,
            east=east,
            north=north,
            up=up,
            observer_height_m=observer_height_m,
            dead_zone_km=dead_zone_km,
            geometry=geometry,
            view_center=view_center,
            content_fov_deg=content_fov_deg,
            terrain_profile_altaz=terrain_profile_altaz,
        )
        fragments.extend(sampled)

    stitched: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []
    close_threshold_px = 0.75
    for fragment in fragments:
        if len(fragment) < 2:
            continue
        if not current:
            current = fragment[:]
            continue
        last_x, last_y = current[-1]
        first_x, first_y = fragment[0]
        if math.hypot(first_x - last_x, first_y - last_y) <= close_threshold_px:
            current.extend(fragment[1:])
            continue
        if len(current) >= 2:
            stitched.append(current)
        current = fragment[:]
    if len(current) >= 2:
        stitched.append(current)
    return stitched


@lru_cache(maxsize=1)
def load_earth_guide_rings(path_str: str = EARTH_GUIDE_LAND_FILE) -> tuple[EarthGuideRing, ...]:
    path = Path(path_str)
    payload = json.loads(path.read_text(encoding="utf-8"))
    raw_rings = payload.get("rings", [])
    rings: list[EarthGuideRing] = []
    for item in raw_rings:
        if not isinstance(item, dict):
            continue
        raw_points = item.get("points_lonlat_deg")
        if not isinstance(raw_points, list):
            continue
        lonlat: list[tuple[float, float]] = []
        xyz: list[tuple[float, float, float]] = []
        for pair in raw_points:
            if not isinstance(pair, list) or len(pair) < 2:
                continue
            lon_deg = float(pair[0])
            lat_deg = float(pair[1])
            lonlat.append((lon_deg, lat_deg))
            xyz.append(_lonlat_to_unit_xyz(lon_deg, lat_deg))
        if len(lonlat) < 3:
            continue
        rings.append(
            EarthGuideRing(
                source_name=str(item.get("source_name", "")),
                label_name=(str(item["label_name"]) if item.get("label_name") else None),
                points_lonlat_deg=np.asarray(lonlat, dtype=np.float64),
                points_xyz=np.asarray(xyz, dtype=np.float64),
            )
        )
    return tuple(rings)


def draw_earth_guide(
    painter: QPainter,
    *,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
    terrain_horizon_opacity: float = 0.028,
    content_fov_deg: float = 90.0,
) -> None:
    rings = load_earth_guide_rings()
    if not rings:
        return
    if float(terrain_horizon_opacity) <= 0.0:
        return
    origin, east, north, up = _observer_basis(observer_lat_deg, observer_lon_deg, observer_height_m)
    dead_zone_km = _observer_dead_zone_km(observer_height_m)
    painter.save()
    try:
        alpha = terrain_horizon_line_alpha(terrain_horizon_opacity)
        fill_color = QColor(*TERRAIN_HORIZON_LINE_COLOR)
        fill_color.setAlphaF(alpha * EARTH_GUIDE_FILL_ALPHA_SCALE)
        fill_brush = fill_color

        line_color = QColor(*TERRAIN_HORIZON_LINE_COLOR)
        line_color.setAlphaF(alpha)
        line_pen = QPen(line_color, EARTH_GUIDE_FOREGROUND_WIDTH, Qt.PenStyle.SolidLine)
        line_pen.setCosmetic(True)
        line_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        line_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

        for ring in rings:
            for fragment in _ring_fragments_altaz(
                ring,
                origin=origin,
                east=east,
                north=north,
                up=up,
                observer_height_m=observer_height_m,
                dead_zone_km=dead_zone_km,
                geometry=geometry,
                view_center=view_center,
                content_fov_deg=content_fov_deg,
                terrain_profile_altaz=terrain_profile_altaz,
            ):
                screen_points = [QPointF(x, y) for x, y in fragment]
                if len(screen_points) < 2:
                    continue
                poly = QPolygonF(screen_points)
                painter.fillPath(_fragment_fill_path(fragment), fill_brush)
                painter.setPen(line_pen)
                painter.drawPolyline(poly)
    finally:
        painter.restore()
