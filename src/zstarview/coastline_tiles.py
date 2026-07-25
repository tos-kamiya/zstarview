"""Small local reader for the experimental coastline vector tiles.

This is intentionally a read-only preview path.  It does not download tiles,
validate manifests, or populate the application cache.
"""

from __future__ import annotations

import gzip
import json
import logging
import math
import os
from pathlib import Path
from typing import Iterable

from pyproj.enums import TransformDirection

from .astro import is_in_fov
from .location_resolver.place_projection import project_place_targets_to_altaz
from .paths import CACHE_PATH
from .water_overlay import WaterOverlayPoint, WaterOverlayPolyline
from .water_surface_mesh import make_local_transformer, project_ring_xy

logger = logging.getLogger(__name__)

GRID_COLS = 32
GRID_ROWS = 16
PREVIEW_RADIUS_KM = 10.0
PREVIEW_ROOT_ENV = "ZSTARVIEW_COASTLINE_TILE_DIR"


def _repository_preview_roots() -> tuple[Path, ...]:
    repository_root = Path(__file__).resolve().parents[2]
    return (
        repository_root / "raw-data" / "coastline_vector_tiles_32x16",
        repository_root / "raw-data" / "coastline_vector_tiles_32x16_top8_split4x4",
    )


def _tile_roots() -> tuple[Path, ...]:
    configured = os.environ.get(PREVIEW_ROOT_ENV, "").strip()
    if configured:
        configured_root = Path(configured)
        return (configured_root,) if configured_root.exists() else ()
    cache_root = (
        Path(CACHE_PATH)
        / "coastline"
        / "osm-water-polygons"
        / "20260725"
        / "schema-1"
        / "grid-32x16"
    )
    preview_roots = _repository_preview_roots()
    # The second experimental directory contains the expensive parent-tile
    # subdivisions.  Prefer it before opening one of the 200+ MB parents.
    roots = (cache_root, preview_roots[1], preview_roots[0])
    return tuple(root for root in roots if root.exists())


def _parent_bounds(row: int, col: int) -> tuple[float, float, float, float]:
    width = 360.0 / GRID_COLS
    height = 180.0 / GRID_ROWS
    min_lon = -180.0 + col * width
    max_lat = 90.0 - row * height
    return min_lon, max_lat - height, min_lon + width, max_lat


def _bbox_overlaps(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    return a[0] <= b[2] and a[2] >= b[0] and a[1] <= b[3] and a[3] >= b[1]


def _read_geojson(path: Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict) or payload.get("type") != "FeatureCollection":
        raise ValueError(f"invalid coastline GeoJSON: {path}")
    return payload


def _line_coordinates(geometry: object) -> Iterable[tuple[tuple[float, float], ...]]:
    if not isinstance(geometry, dict):
        return ()
    kind = geometry.get("type")
    coordinates = geometry.get("coordinates")
    if kind == "LineString" and isinstance(coordinates, list):
        return (_normalize_line(coordinates),)
    if kind == "MultiLineString" and isinstance(coordinates, list):
        return tuple(_normalize_line(line) for line in coordinates)
    if kind == "GeometryCollection":
        geometries = geometry.get("geometries")
        if isinstance(geometries, list):
            return tuple(line for item in geometries for line in _line_coordinates(item))
    return ()


def _normalize_line(values: object) -> tuple[tuple[float, float], ...]:
    if not isinstance(values, list):
        return ()
    result = []
    for value in values:
        if isinstance(value, (list, tuple)) and len(value) >= 2:
            result.append((float(value[0]), float(value[1])))
    return tuple(result)


def _clip_line_to_radius(
    points: Iterable[tuple[float, float]], radius_m: float
) -> tuple[tuple[tuple[float, float], ...], ...]:
    points = tuple(points)
    if len(points) < 2:
        return ()
    radius_sq = float(radius_m) ** 2
    fragments: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []

    def finish() -> None:
        nonlocal current
        if len(current) >= 2:
            fragments.append(current)
        current = []

    for start, end in zip(points, points[1:]):
        x0, y0 = start
        dx, dy = end[0] - x0, end[1] - y0
        cuts = [0.0, 1.0]
        a = dx * dx + dy * dy
        if a > 0.0:
            b = 2.0 * (x0 * dx + y0 * dy)
            c = x0 * x0 + y0 * y0 - radius_sq
            discriminant = b * b - 4.0 * a * c
            if discriminant > 0.0:
                root = discriminant**0.5
                cuts.extend(
                    t for t in ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a))
                    if 0.0 < t < 1.0
                )
        for left, right in zip(sorted(set(cuts)), sorted(set(cuts))[1:]):
            p0 = (x0 + dx * left, y0 + dy * left)
            p1 = (x0 + dx * right, y0 + dy * right)
            midpoint = ((p0[0] + p1[0]) * 0.5, (p0[1] + p1[1]) * 0.5)
            if midpoint[0] * midpoint[0] + midpoint[1] * midpoint[1] <= radius_sq:
                if not current:
                    current.append(p0)
                elif current[-1] != p0:
                    current.append(p0)
                current.append(p1)
            else:
                finish()
    finish()
    return tuple(tuple(fragment) for fragment in fragments)


def _candidate_paths(root: Path, row: int, col: int, bbox: tuple[float, float, float, float]) -> list[Path]:
    parent_bbox = _parent_bounds(row, col)
    split_pattern = root / f"tile_y{row}_x{col}_q*.geojson"
    split_paths = sorted(root.glob(split_pattern.name))
    if split_paths:
        result = []
        parent_width = (parent_bbox[2] - parent_bbox[0]) / 4.0
        parent_height = (parent_bbox[3] - parent_bbox[1]) / 4.0
        for path in split_paths:
            stem = path.stem
            try:
                qrow, qcol = int(stem[-2]), int(stem[-1])
            except ValueError:
                continue
            child_bbox = (
                parent_bbox[0] + qcol * parent_width,
                parent_bbox[1] + qrow * parent_height,
                parent_bbox[0] + (qcol + 1) * parent_width,
                parent_bbox[1] + (qrow + 1) * parent_height,
            )
            if _bbox_overlaps(child_bbox, bbox):
                result.append(path)
        return result
    for suffix in (".geojson", ".geojson.gz"):
        path = root / f"tile_y{row}_x{col}{suffix}"
        if path.exists():
            return [path]
    return []


def load_coastline_overlay_polylines(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    max_distance_km: float = PREVIEW_RADIUS_KM,
    view_center: tuple[float, float] | None = None,
    fov_deg: float | None = None,
) -> tuple[WaterOverlayPolyline, ...]:
    """Load nearby experimental vector tiles and project their linework."""
    request_bbox = _bbox_from_point(observer_lat_deg, observer_lon_deg, max_distance_km)
    transformer = make_local_transformer(observer_lat_deg, observer_lon_deg)
    radius_m = max_distance_km * 1000.0
    result: list[WaterOverlayPolyline] = []
    seen: set[tuple[str, int, tuple[tuple[float, float], ...]]] = set()
    for root in _tile_roots():
        for row in range(GRID_ROWS):
            for col in range(GRID_COLS):
                if not _bbox_overlaps(_parent_bounds(row, col), request_bbox):
                    continue
                for path in _candidate_paths(root, row, col, request_bbox):
                    try:
                        payload = _read_geojson(path)
                    except (OSError, ValueError, json.JSONDecodeError) as exc:
                        logger.warning("Cannot read coastline preview tile %s: %s", path, exc)
                        continue
                    for feature_index, feature in enumerate(payload.get("features", ())):
                        if not isinstance(feature, dict):
                            continue
                        for line_index, line in enumerate(_line_coordinates(feature.get("geometry"))):
                            if len(line) < 2:
                                continue
                            line_xy = project_ring_xy(line, transformer)
                            for fragment_index, fragment in enumerate(_clip_line_to_radius(line_xy, radius_m)):
                                if len(fragment) < 2:
                                    continue
                                xs = [point[0] for point in fragment]
                                ys = [point[1] for point in fragment]
                                lon_values, lat_values = transformer.transform(
                                    xs, ys, direction=TransformDirection.INVERSE
                                )
                                projections = project_place_targets_to_altaz(
                                    observer_latitude_deg=observer_lat_deg,
                                    observer_longitude_deg=observer_lon_deg,
                                    observer_height_m=observer_height_m,
                                    target_latitude_deg=[float(value) for value in lat_values],
                                    target_longitude_deg=[float(value) for value in lon_values],
                                    target_height_m=[0.0] * len(fragment),
                                )
                                points = tuple(
                                    WaterOverlayPoint(
                                        water_id=f"coastline/{path.name}/{feature_index}/{line_index}/{fragment_index}",
                                        alt_deg=float(item.alt_deg),
                                        az_deg=float(item.az_deg),
                                        distance_km=float(item.distance_km),
                                        water_category="coastline",
                                    )
                                    for item in projections
                                    if view_center is None
                                    or fov_deg is None
                                    or is_in_fov(float(item.alt_deg), float(item.az_deg), view_center, fov_deg=fov_deg)
                                )
                                if len(points) >= 2:
                                    key = (path.name, feature_index, tuple((point.alt_deg, point.az_deg) for point in points))
                                    if key in seen:
                                        continue
                                    seen.add(key)
                                    result.append(
                                        WaterOverlayPolyline(
                                            water_id=f"coastline/{path.name}/{feature_index}/{line_index}/{fragment_index}",
                                            water_category="coastline",
                                            points=points,
                                        )
                                    )
        if result:
            break
    if result:
        logger.info("Coastline preview: %d projected line fragments", len(result))
    return tuple(result)


def _bbox_from_point(lat_deg: float, lon_deg: float, radius_km: float) -> tuple[float, float, float, float]:
    lat_delta = radius_km / 111.2
    cos_lat = abs(math.cos(math.radians(lat_deg)))
    lon_delta = 180.0 if cos_lat < 1.0e-9 else radius_km / (111.2 * cos_lat)
    return (
        max(-180.0, lon_deg - lon_delta),
        max(-90.0, lat_deg - lat_delta),
        min(180.0, lon_deg + lon_delta),
        min(90.0, lat_deg + lat_delta),
    )
