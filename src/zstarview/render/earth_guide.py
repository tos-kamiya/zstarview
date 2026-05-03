from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import EARTH_GUIDE_LAND_FILE, EARTH_GUIDE_LINE_COLOR
from ..types import ScreenGeometry
from .geometry import normalized_to_screen_xy

EARTH_MEAN_RADIUS_M = 6_371_008.8

EARTH_GUIDE_DEAD_ZONE_MIN_KM = 20.0
EARTH_GUIDE_DEAD_ZONE_MAX_KM = 80.0
EARTH_GUIDE_DEAD_ZONE_SCALE = 0.25
EARTH_GUIDE_HORIZON_MARGIN_DEG = 1.0
EARTH_GUIDE_FILL_DEAD_ZONE_KM = 600.0
EARTH_GUIDE_FILL_CELL_AREA_DEG2 = 10.0
EARTH_GUIDE_FILL_LINE_WIDTH_PX = 1.35
EARTH_GUIDE_FILL_ALPHA = 0.12
EARTH_GUIDE_FILL_FAST_MODE_THINNING = 3
EARTH_GUIDE_FILL_SAMPLER_MODE = "equal_area"
EARTH_GUIDE_FILL_SOUTH_POLE_CUTOFF_LAT_DEG = -60.0
EARTH_GUIDE_FILL_LAT_BAND_DEG = 0.5
EARTH_GUIDE_FILL_MAX_LON_GAP_DEG = 8.0
EARTH_GUIDE_UNDERLAY_WIDTH = 12.0
EARTH_GUIDE_FOREGROUND_WIDTH = 1.5


@dataclass(frozen=True)
class EarthGuideRing:
    source_name: str
    label_name: str | None
    points_lonlat_deg: np.ndarray
    points_xyz: np.ndarray
    approx_area_deg2: float | None = None
    fill_points_lonlat_deg: np.ndarray | None = None
    fill_points_xyz: np.ndarray | None = None


def _lonlat_to_unit_xyz(lon_deg: float, lat_deg: float) -> tuple[float, float, float]:
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)
    cos_lat = math.cos(lat_rad)
    return (
        cos_lat * math.cos(lon_rad),
        cos_lat * math.sin(lon_rad),
        math.sin(lat_rad),
    )


def _wrap_lon_deg(lon_deg: float) -> float:
    return ((float(lon_deg) + 180.0) % 360.0) - 180.0


def _unwrap_ring_lonlat_deg(points_lonlat_deg: np.ndarray) -> np.ndarray:
    lon = np.asarray(points_lonlat_deg[:, 0], dtype=np.float64)
    lat = np.asarray(points_lonlat_deg[:, 1], dtype=np.float64)
    lon_unwrapped = np.degrees(np.unwrap(np.radians(lon)))
    return np.column_stack((lon_unwrapped, lat))


def _point_in_polygon_2d(x: float, y: float, polygon_xy: np.ndarray) -> bool:
    count = len(polygon_xy)
    if count < 3:
        return False
    inside = False
    prev_x = float(polygon_xy[-1][0])
    prev_y = float(polygon_xy[-1][1])
    for index in range(count):
        curr_x = float(polygon_xy[index][0])
        curr_y = float(polygon_xy[index][1])
        if ((curr_y > y) != (prev_y > y)) and (
            x < ((prev_x - curr_x) * (y - curr_y) / ((prev_y - curr_y) or 1.0e-12)) + curr_x
        ):
            inside = not inside
        prev_x = curr_x
        prev_y = curr_y
    return inside


def _fill_lat_band_key(lat_deg: float, band_deg: float = EARTH_GUIDE_FILL_LAT_BAND_DEG) -> int:
    band = max(1.0e-3, float(band_deg))
    return int(round(float(lat_deg) / band))


def _build_ring_fill_points(
    ring: EarthGuideRing,
    *,
    sampler_mode: str = EARTH_GUIDE_FILL_SAMPLER_MODE,
    target_cell_area_deg2: float = EARTH_GUIDE_FILL_CELL_AREA_DEG2,
) -> tuple[np.ndarray, np.ndarray]:
    points_lonlat = np.asarray(ring.points_lonlat_deg, dtype=np.float64)
    if len(points_lonlat) < 3:
        return (
            np.empty((0, 2), dtype=np.float64),
            np.empty((0, 3), dtype=np.float64),
        )

    lonlat_unwrapped = _unwrap_ring_lonlat_deg(points_lonlat)
    lon_min = float(np.min(lonlat_unwrapped[:, 0]))
    lon_max = float(np.max(lonlat_unwrapped[:, 0]))
    lat_min = float(np.min(lonlat_unwrapped[:, 1]))
    lat_max = float(np.max(lonlat_unwrapped[:, 1]))

    if lat_max < EARTH_GUIDE_FILL_SOUTH_POLE_CUTOFF_LAT_DEG:
        return (
            np.empty((0, 2), dtype=np.float64),
            np.empty((0, 3), dtype=np.float64),
        )

    lat_span = max(1.0e-6, lat_max - lat_min)
    approx_area = float(ring.approx_area_deg2 or max(1.0, lat_span * max(1.0e-6, lon_max - lon_min)))
    cell_area = max(6.0, float(target_cell_area_deg2))
    target_points = max(8, int(round(approx_area / cell_area)))
    sampler_mode = str(sampler_mode)

    points: list[tuple[float, float]] = []
    if sampler_mode == "latlon":
        lat_step = max(1.5, math.sqrt(cell_area))
        lon_step = lat_step
        lat_positions = np.arange(lat_min, lat_max + (lat_step * 0.5), lat_step, dtype=np.float64)
        for row_index, lat in enumerate(lat_positions):
            if not (lat_min <= float(lat) <= lat_max):
                continue
            lon_offset = 0.0 if (row_index % 2 == 0) else lon_step * 0.5
            lon_positions = np.arange(
                lon_min + lon_offset,
                lon_max + lon_step,
                lon_step,
                dtype=np.float64,
            )
            for lon in lon_positions:
                if _point_in_polygon_2d(float(lon), float(lat), lonlat_unwrapped):
                    points.append((_wrap_lon_deg(float(lon)), float(lat)))
    else:
        band_count = max(2, int(round(math.sqrt(target_points))))
        sin_min = math.sin(math.radians(lat_min))
        sin_max = math.sin(math.radians(lat_max))
        sin_span = max(1.0e-6, sin_max - sin_min)
        row_height_deg = lat_span / float(band_count)
        for row_index in range(band_count):
            row_sin_min = sin_min + (sin_span * (row_index / band_count))
            row_sin_max = sin_min + (sin_span * ((row_index + 1) / band_count))
            row_lat = math.degrees(math.asin(max(-1.0, min(1.0, (row_sin_min + row_sin_max) * 0.5))))
            cos_lat = max(0.2, math.cos(math.radians(row_lat)))
            lon_step = max(1.5, cell_area / max(1.0e-6, row_height_deg * cos_lat))
            lon_offset = 0.0 if (row_index % 2 == 0) else lon_step * 0.5
            lon_positions = np.arange(
                lon_min + lon_offset,
                lon_max + lon_step,
                lon_step,
                dtype=np.float64,
            )
            if sampler_mode == "jitter":
                jitter = min(0.35 * lon_step, 0.65)
                lon_positions = lon_positions + ((row_index % 3) - 1) * jitter * 0.5
            for lon in lon_positions:
                lat = math.degrees(math.asin(max(-1.0, min(1.0, (row_sin_min + row_sin_max) * 0.5))))
                if _point_in_polygon_2d(float(lon), float(lat), lonlat_unwrapped):
                    points.append((_wrap_lon_deg(float(lon)), float(lat)))

    if not points:
        return (
            np.empty((0, 2), dtype=np.float64),
            np.empty((0, 3), dtype=np.float64),
        )

    lonlat = np.asarray(points, dtype=np.float64)
    xyz = np.asarray([_lonlat_to_unit_xyz(float(lon), float(lat)) for lon, lat in lonlat], dtype=np.float64)
    return lonlat, xyz


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
    horizon_limit = -(
        _observer_horizon_dip_deg(observer_height_m) + EARTH_GUIDE_HORIZON_MARGIN_DEG
    )
    return min(terrain_limit, horizon_limit)


def earth_guide_line_alpha(opacity: float) -> float:
    """Return the alpha curve used by the Earth guide foreground line."""
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(opacity, min(1.0, 0.32 + (opacity * 0.78)))


def earth_guide_underlay_line_alpha(opacity: float) -> float:
    """Return the softer alpha curve used by the Earth guide underlay line."""
    opacity = max(0.0, min(1.0, float(opacity)))
    return max(0.0, min(1.0, 0.03 + (opacity * 0.12)))


def _earth_guide_underlay_pass_specs(opacity: float) -> tuple[tuple[float, float], ...]:
    base_alpha = earth_guide_underlay_line_alpha(opacity)
    return (
        (EARTH_GUIDE_UNDERLAY_WIDTH, base_alpha * 0.55),
        (EARTH_GUIDE_UNDERLAY_WIDTH * 0.68, base_alpha * 0.78),
        (EARTH_GUIDE_UNDERLAY_WIDTH * 0.42, base_alpha),
    )


def _surface_distance_km(point_xyz: np.ndarray, observer_up: np.ndarray) -> float:
    dot = float(np.dot(point_xyz, observer_up))
    dot = max(-1.0, min(1.0, dot))
    return (EARTH_MEAN_RADIUS_M / 1000.0) * math.acos(dot)


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
    edge_fov_deg: float,
    content_fov_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> tuple[tuple[float, float], bool] | None:
    altaz = _project_point_altaz(xyz, origin=origin, east=east, north=north, up=up)
    if altaz is None:
        return None
    alt_deg, az_deg = altaz
    if not is_in_fov(alt_deg, az_deg, view_center, fov_deg=content_fov_deg):
        return None
    nx, ny = altaz_to_normalized_xy(alt_deg, az_deg, view_center, edge_fov_deg=edge_fov_deg)
    screen_xy = normalized_to_screen_xy(nx, ny, geometry)
    limit_deg = _observer_visible_altitude_limit_deg(
        az_deg,
        observer_height_m,
        terrain_profile_altaz,
    )
    visible = _surface_distance_km(xyz, up) >= dead_zone_km and alt_deg <= limit_deg
    return (screen_xy, visible)


def _project_fill_xyz_to_screen_point(
    xyz: np.ndarray,
    *,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    observer_height_m: float,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> tuple[tuple[float, float], bool] | None:
    projected = _project_xyz_to_screen_point(
        xyz,
        origin=origin,
        east=east,
        north=north,
        up=up,
        observer_height_m=observer_height_m,
        dead_zone_km=EARTH_GUIDE_FILL_DEAD_ZONE_KM,
        geometry=geometry,
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
        terrain_profile_altaz=terrain_profile_altaz,
    )
    return projected


def _draw_fill_segments_for_ring(
    painter: QPainter,
    ring: EarthGuideRing,
    *,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    observer_height_m: float,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    fill_pen: QPen,
    fill_thinning: int,
    fast_mode: bool,
) -> None:
    fill_points_xyz = ring.fill_points_xyz
    fill_points_lonlat = ring.fill_points_lonlat_deg
    if fill_points_xyz is None or fill_points_lonlat is None:
        return
    if len(fill_points_xyz) == 0 or len(fill_points_lonlat) == 0:
        return

    grouped: dict[int, list[tuple[float, tuple[float, float]]]] = {}
    for point_index, ((lon_deg, lat_deg), xyz) in enumerate(zip(fill_points_lonlat, fill_points_xyz)):
        if fill_thinning > 1 and ((point_index + _fill_lat_band_key(float(lat_deg))) % fill_thinning != 0):
            continue
        projected = _project_fill_xyz_to_screen_point(
            xyz,
            origin=origin,
            east=east,
            north=north,
            up=up,
            observer_height_m=observer_height_m,
            geometry=geometry,
            view_center=view_center,
            edge_fov_deg=edge_fov_deg,
            content_fov_deg=content_fov_deg,
            terrain_profile_altaz=terrain_profile_altaz,
        )
        if projected is None:
            continue
        screen_xy, visible = projected
        if not visible:
            continue
        band_key = _fill_lat_band_key(float(lat_deg))
        grouped.setdefault(band_key, []).append((float(lon_deg), screen_xy))

    if not grouped:
        return

    painter.setPen(fill_pen)
    lon_gap_limit = EARTH_GUIDE_FILL_MAX_LON_GAP_DEG * (1.15 if fast_mode else 1.0)
    for band_key in sorted(grouped):
        band_points = grouped[band_key]
        if len(band_points) < 2:
            continue
        band_points.sort(key=lambda item: item[0])
        prev_lon, prev_xy = band_points[0]
        for lon_deg, xy in band_points[1:]:
            if abs(float(lon_deg) - prev_lon) <= lon_gap_limit:
                painter.drawLine(
                    QPointF(prev_xy[0], prev_xy[1]),
                    QPointF(xy[0], xy[1]),
                )
            prev_lon = lon_deg
            prev_xy = xy


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
    edge_fov_deg: float,
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
        edge_fov_deg=edge_fov_deg,
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
        edge_fov_deg=edge_fov_deg,
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
        edge_fov_deg=edge_fov_deg,
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
            edge_fov_deg=edge_fov_deg,
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
            edge_fov_deg=edge_fov_deg,
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
    edge_fov_deg: float,
    content_fov_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    max_depth: int,
    threshold_px: float,
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
            edge_fov_deg=edge_fov_deg,
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
        approx_area_deg2: float | None = None
        raw_area = item.get("approx_area_deg2")
        if isinstance(raw_area, (int, float)):
            approx_area_deg2 = float(raw_area)
        ring = EarthGuideRing(
            source_name=str(item.get("source_name", "")),
            label_name=(str(item["label_name"]) if item.get("label_name") else None),
            points_lonlat_deg=np.asarray(lonlat, dtype=np.float64),
            points_xyz=np.asarray(xyz, dtype=np.float64),
            approx_area_deg2=approx_area_deg2,
        )
        fill_lonlat, fill_xyz = _build_ring_fill_points(ring)
        rings.append(
            EarthGuideRing(
                source_name=ring.source_name,
                label_name=ring.label_name,
                points_lonlat_deg=ring.points_lonlat_deg,
                points_xyz=ring.points_xyz,
                approx_area_deg2=ring.approx_area_deg2,
                fill_points_lonlat_deg=fill_lonlat,
                fill_points_xyz=fill_xyz,
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
    earth_guide_opacity: float = 0.028,
    visibility_boost: float = 1.0,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = 90.0,
    fast_mode: bool = False,
) -> None:
    rings = load_earth_guide_rings()
    if not rings:
        return
    if float(earth_guide_opacity) <= 0.0:
        return
    origin, east, north, up = _observer_basis(
        observer_lat_deg,
        observer_lon_deg,
        observer_height_m,
    )
    dead_zone_km = _observer_dead_zone_km(observer_height_m)
    boost_scale = max(1.0, float(visibility_boost))
    painter.save()
    try:
        if fast_mode:
            rings = rings[::2] or rings[:1]
            max_depth = 8
            threshold_px = 40.0
        else:
            max_depth = 12
            threshold_px = 24.0
            fill_alpha = max(
                0.0,
                min(1.0, EARTH_GUIDE_FILL_ALPHA + (float(earth_guide_opacity) * 0.35)),
            )
            fill_color = QColor(*EARTH_GUIDE_LINE_COLOR)
            fill_color.setAlphaF(fill_alpha)
            fill_line_width = 1.35
            fill_pen = QPen(fill_color, fill_line_width, Qt.PenStyle.SolidLine)
            fill_pen.setCosmetic(True)
            fill_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            fill_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            fill_thinning = 1
            painter.save()
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setPen(fill_pen)
            for ring in rings:
                _draw_fill_segments_for_ring(
                    painter,
                    ring,
                    origin=origin,
                    east=east,
                    north=north,
                    up=up,
                    observer_height_m=observer_height_m,
                    geometry=geometry,
                    view_center=view_center,
                    edge_fov_deg=edge_fov_deg,
                    content_fov_deg=content_fov_deg,
                    terrain_profile_altaz=terrain_profile_altaz,
                    fill_pen=fill_pen,
                    fill_thinning=fill_thinning,
                    fast_mode=fast_mode,
                )
            painter.restore()
        if fast_mode:
            line_alpha = max(0.0, min(1.0, 0.18 + (earth_guide_opacity * 0.25)))
            line_width = max(0.7, EARTH_GUIDE_FOREGROUND_WIDTH * 0.75)
            underlay_pens: list[QPen] = []
        else:
            underlay_pens = []
            for pass_index, (width, alpha) in enumerate(_earth_guide_underlay_pass_specs(earth_guide_opacity)):
                pass_boost = boost_scale if pass_index == 2 else 1.0
                underlay_color = QColor(*EARTH_GUIDE_LINE_COLOR)
                underlay_color.setAlphaF(max(0.0, min(1.0, alpha * pass_boost)))
                pen = QPen(underlay_color, width * pass_boost, Qt.PenStyle.SolidLine)
                pen.setCosmetic(True)
                pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
                underlay_pens.append(pen)
            line_alpha = earth_guide_line_alpha(earth_guide_opacity)
            line_width = EARTH_GUIDE_FOREGROUND_WIDTH

        line_color = QColor(*EARTH_GUIDE_LINE_COLOR)
        line_color.setAlphaF(line_alpha)
        line_pen = QPen(
            line_color,
            line_width,
            Qt.PenStyle.SolidLine,
        )
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
                edge_fov_deg=edge_fov_deg,
                content_fov_deg=content_fov_deg,
                terrain_profile_altaz=terrain_profile_altaz,
                max_depth=max_depth,
                threshold_px=threshold_px,
            ):
                screen_points = [QPointF(x, y) for x, y in fragment]
                if len(screen_points) < 2:
                    continue
                poly = QPolygonF(screen_points)
                for underlay_pen in underlay_pens:
                    painter.setPen(underlay_pen)
                    painter.drawPolyline(poly)
                painter.setPen(line_pen)
                painter.drawPolyline(poly)
    finally:
        painter.restore()
