from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import EARTH_GUIDE_LAND_FILE
from ..types import ScreenGeometry
from .geometry import normalized_to_screen_xy


EARTH_MEAN_RADIUS_M = 6_371_008.8
EARTH_GUIDE_LINE_RGBA = (245, 248, 252, 150)
EARTH_GUIDE_OUTLINE_RGBA = (70, 54, 26, 120)
EARTH_GUIDE_MAX_VISIBLE_ALT_DEG = -5.0


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
        return float(EARTH_GUIDE_MAX_VISIBLE_ALT_DEG)
    terrain_limit = _interpolate_horizon_altitude_deg(az_deg, terrain_profile_altaz)
    return max(float(EARTH_GUIDE_MAX_VISIBLE_ALT_DEG), terrain_limit)


def _signed_visible_distance_deg(
    alt_deg: float,
    az_deg: float,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> float:
    return float(alt_deg) - _effective_visible_altitude_limit_deg(az_deg, terrain_profile_altaz)


def _interpolate_xyz_on_sphere(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    point = ((1.0 - t) * a) + (t * b)
    norm = float(np.linalg.norm(point))
    if norm <= 1.0e-12:
        return a
    return point / norm


def _edge_crossing_altaz(
    xyz0: np.ndarray,
    xyz1: np.ndarray,
    *,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    distance0: float,
    distance1: float,
) -> tuple[float, float] | None:
    if distance0 == distance1:
        return None
    lo = 0.0
    hi = 1.0
    dlo = float(distance0)
    dhi = float(distance1)
    for _ in range(20):
        mid = (lo + hi) * 0.5
        altaz = _project_point_altaz(
            _interpolate_xyz_on_sphere(xyz0, xyz1, mid),
            origin=origin,
            east=east,
            north=north,
            up=up,
        )
        if altaz is None:
            return None
        mid_dist = _signed_visible_distance_deg(altaz[0], altaz[1], terrain_profile_altaz)
        if (dlo < 0.0) == (mid_dist < 0.0):
            lo = mid
            dlo = mid_dist
        else:
            hi = mid
            dhi = mid_dist
    _ = dhi
    return _project_point_altaz(
        _interpolate_xyz_on_sphere(xyz0, xyz1, (lo + hi) * 0.5),
        origin=origin,
        east=east,
        north=north,
        up=up,
    )


def _ring_fragments_altaz(
    ring: EarthGuideRing,
    *,
    origin: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> list[list[tuple[float, float]]]:
    xyz_points = ring.points_xyz
    projected: list[tuple[float, float] | None] = []
    distances: list[float | None] = []
    for point_xyz in xyz_points:
        altaz = _project_point_altaz(point_xyz, origin=origin, east=east, north=north, up=up)
        projected.append(altaz)
        if altaz is None:
            distances.append(None)
        else:
            distances.append(_signed_visible_distance_deg(altaz[0], altaz[1], terrain_profile_altaz))

    fragments: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []
    count = len(xyz_points)
    for index in range(count):
        next_index = (index + 1) % count
        altaz0 = projected[index]
        altaz1 = projected[next_index]
        if altaz0 is None or altaz1 is None:
            if len(current) >= 2:
                fragments.append(current)
            current = []
            continue
        dist0 = float(distances[index])
        dist1 = float(distances[next_index])
        inside0 = dist0 < 0.0
        inside1 = dist1 < 0.0
        if inside0 and not current:
            current = [altaz0]
        if inside0 and inside1:
            current.append(altaz1)
            continue
        crossing = None
        if inside0 != inside1:
            crossing = _edge_crossing_altaz(
                xyz_points[index],
                xyz_points[next_index],
                origin=origin,
                east=east,
                north=north,
                up=up,
                terrain_profile_altaz=terrain_profile_altaz,
                distance0=dist0,
                distance1=dist1,
            )
        if inside0 and not inside1:
            if crossing is not None:
                current.append(crossing)
            if len(current) >= 2:
                fragments.append(current)
            current = []
        elif not inside0 and inside1:
            current = [crossing] if crossing is not None else []
            current.append(altaz1)
        elif len(current) >= 2:
            fragments.append(current)
            current = []
    if len(current) >= 2:
        fragments.append(current)
    return fragments


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
    content_fov_deg: float = 90.0,
) -> None:
    rings = load_earth_guide_rings()
    if not rings:
        return
    origin, east, north, up = _observer_basis(observer_lat_deg, observer_lon_deg, observer_height_m)
    painter.save()
    try:
        outline_pen = QPen(QColor(*EARTH_GUIDE_OUTLINE_RGBA), 3.0, Qt.PenStyle.SolidLine)
        outline_pen.setCosmetic(True)
        outline_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        outline_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        line_pen = QPen(QColor(*EARTH_GUIDE_LINE_RGBA), 1.3, Qt.PenStyle.SolidLine)
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
                terrain_profile_altaz=terrain_profile_altaz,
            ):
                screen_points: list[QPointF] = []
                for alt_deg, az_deg in fragment:
                    if not is_in_fov(alt_deg, az_deg, view_center, fov_deg=content_fov_deg):
                        continue
                    nx, ny = altaz_to_normalized_xy(alt_deg, az_deg, view_center)
                    screen_points.append(QPointF(*normalized_to_screen_xy(nx, ny, geometry)))
                if len(screen_points) < 2:
                    continue
                poly = QPolygonF(screen_points)
                painter.setPen(outline_pen)
                painter.drawPolyline(poly)
                painter.setPen(line_pen)
                painter.drawPolyline(poly)
    finally:
        painter.restore()
