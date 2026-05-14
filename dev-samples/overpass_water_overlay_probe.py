#!/usr/bin/env python3
"""Fetch and normalize nearby OSM water polygons for overlay prototyping.

This sample keeps the pipeline intentionally simple:

1. Query Overpass around a latitude/longitude center point.
2. Reconstruct closed polygon rings from OSM nodes and ways.
3. Keep only water polygons:
   - outer rings
   - inner rings
4. Write a normalized JSON payload and, optionally, a quick SVG preview.
5. Optionally triangulate the polygon rings into SVG triangles for shape checks through the shared mesh helper.
6. Optionally split polygons into a meter grid before triangulation to probe thin-shape behavior.

The script is meant for data-structure exploration, not production caching.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from zstarview.water_overlay import WaterPolygonFootprint
from zstarview.water_surface_mesh import build_water_surface_mesh


EARTH_RADIUS_KM = 6371.0088
EARTH_RADIUS_M = EARTH_RADIUS_KM * 1000.0
DEFAULT_ENDPOINT = "https://overpass-api.de/api/interpreter"
DEFAULT_USER_AGENT = "zstarview-water-overlay-probe/0.1"
DEFAULT_TIMEOUT_S = 60.0

POLYGON_WATER_KEYS = {
    ("natural", "water"),
    ("waterway", "riverbank"),
}


@dataclass(frozen=True)
class WaterPolygon:
    osm_id: str
    kind: str
    outer_rings: tuple[tuple[tuple[float, float], ...], ...]
    inner_rings: tuple[tuple[tuple[float, float], ...], ...]
    source: str
    tags: dict[str, str]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch and normalize OSM water polygons around a point. "
            "Coordinates are EPSG:4326 lon/lat degrees."
        )
    )
    parser.add_argument("--lat", type=float, required=True, help="Center latitude in degrees.")
    parser.add_argument("--lon", type=float, required=True, help="Center longitude in degrees.")
    parser.add_argument(
        "--radius-km",
        type=float,
        default=5.0,
        help="Search radius in kilometers around the center point (default: 5.0).",
    )
    parser.add_argument(
        "--endpoint",
        default=DEFAULT_ENDPOINT,
        help=f"Overpass endpoint URL (default: {DEFAULT_ENDPOINT}).",
    )
    parser.add_argument(
        "--user-agent",
        default=DEFAULT_USER_AGENT,
        help=f"HTTP User-Agent string (default: {DEFAULT_USER_AGENT}).",
    )
    parser.add_argument(
        "--timeout-s",
        type=float,
        default=DEFAULT_TIMEOUT_S,
        help=f"HTTP timeout in seconds (default: {DEFAULT_TIMEOUT_S:.0f}).",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        help="Optional output path for the normalized JSON payload. If omitted, stdout is used.",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        help="Optional output path for a quick SVG preview.",
    )
    parser.add_argument(
        "--output-triangulated-svg",
        type=Path,
        help="Optional output path for a triangulated SVG preview.",
    )
    parser.add_argument(
        "--svg-width",
        type=int,
        default=1200,
        help="SVG width in pixels (default: 1200).",
    )
    parser.add_argument(
        "--svg-height",
        type=int,
        default=1200,
        help="SVG height in pixels (default: 1200).",
    )
    parser.add_argument(
        "--svg-padding",
        type=float,
        default=0.05,
        help="Extra padding around the query bbox in the SVG preview (default: 0.05).",
    )
    parser.add_argument(
        "--triangulation-grid-m",
        type=float,
        default=1.0,
        help="Precision grid in meters before triangulation (default: 1.0).",
    )
    parser.add_argument(
        "--triangulation-simplify-m",
        type=float,
        default=20.0,
        help="Topology-preserving simplify tolerance in meters before triangulation (default: 20.0).",
    )
    parser.add_argument(
        "--triangulation-cell-m",
        type=float,
        default=300.0,
        help=(
            "Optional grid cell size in meters for splitting polygons before triangulation "
            "(default: 300.0)."
        ),
    )
    return parser


def bbox_from_point(lat_deg: float, lon_deg: float, radius_km: float) -> tuple[float, float, float, float]:
    if radius_km <= 0.0:
        raise ValueError("radius_km must be positive")
    angular_distance = radius_km / EARTH_RADIUS_KM
    lat_delta_deg = math.degrees(angular_distance)
    cos_lat = math.cos(math.radians(lat_deg))
    if abs(cos_lat) < 1.0e-9:
        lon_delta_deg = 180.0
    else:
        lon_delta_deg = math.degrees(angular_distance / abs(cos_lat))
    min_lon = max(-180.0, lon_deg - lon_delta_deg)
    max_lon = min(180.0, lon_deg + lon_delta_deg)
    min_lat = max(-90.0, lat_deg - lat_delta_deg)
    max_lat = min(90.0, lat_deg + lat_delta_deg)
    return min_lon, min_lat, max_lon, max_lat


def build_overpass_query(bbox: tuple[float, float, float, float]) -> str:
    west, south, east, north = bbox
    return "\n".join(
        [
            "[out:json][timeout:60];",
            "(",
            f'  way["natural"="water"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});',
            f'  relation["natural"="water"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});',
            f'  way["waterway"="riverbank"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});',
            f'  relation["waterway"="riverbank"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});',
            f'  way["waterway"~"^(river|stream|canal|drain)$"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});',
            ");",
            "out body;",
            ">;",
            "out skel qt;",
        ]
    )


def fetch_overpass_json(*, endpoint: str, query: str, user_agent: str, timeout_s: float) -> dict[str, Any]:
    request_body = urllib.parse.urlencode({"data": query}).encode("utf-8")
    request = urllib.request.Request(
        endpoint,
        data=request_body,
        headers={
            "Content-Type": "application/x-www-form-urlencoded; charset=utf-8",
            "User-Agent": user_agent,
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_s) as response:
            payload = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Overpass HTTP error: {exc.code} {exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Overpass network error: {exc.reason}") from exc
    try:
        loaded = json.loads(payload.decode("utf-8"))
    except Exception as exc:
        raise RuntimeError("Overpass returned invalid JSON") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError("Overpass returned a non-object JSON payload")
    return loaded


def tags_to_dict(tags: object) -> dict[str, str]:
    if not isinstance(tags, dict):
        return {}
    result: dict[str, str] = {}
    for key, value in tags.items():
        if isinstance(key, str) and isinstance(value, str):
            result[key] = value
    return result


def is_closed_ring(points: tuple[tuple[float, float], ...]) -> bool:
    return len(points) >= 4 and points[0] == points[-1]


def almost_same_point(a: tuple[float, float], b: tuple[float, float]) -> bool:
    return a == b


def normalize_ring(points: Iterable[tuple[float, float]]) -> tuple[tuple[float, float], ...]:
    ring = tuple(points)
    if not ring:
        return ()
    if ring[0] != ring[-1]:
        ring = ring + (ring[0],)
    return ring


def assemble_rings_from_segments(
    segments: list[tuple[tuple[float, float], ...]],
) -> list[tuple[tuple[float, float], ...]]:
    remaining = [segment for segment in segments if len(segment) >= 2]
    rings: list[tuple[tuple[float, float], ...]] = []
    while remaining:
        ring = list(remaining.pop())
        progress = True
        while progress and remaining:
            progress = False
            start = ring[0]
            end = ring[-1]
            for index, segment in enumerate(remaining):
                seg_start = segment[0]
                seg_end = segment[-1]
                if almost_same_point(end, seg_start):
                    ring.extend(segment[1:])
                    remaining.pop(index)
                    progress = True
                    break
                if almost_same_point(end, seg_end):
                    ring.extend(reversed(segment[:-1]))
                    remaining.pop(index)
                    progress = True
                    break
                if almost_same_point(start, seg_end):
                    ring[:0] = list(segment[:-1])
                    remaining.pop(index)
                    progress = True
                    break
                if almost_same_point(start, seg_start):
                    ring[:0] = list(reversed(segment[1:]))
                    remaining.pop(index)
                    progress = True
                    break
        normalized = normalize_ring(ring)
        if is_closed_ring(normalized):
            rings.append(normalized)
    return rings


def node_lookup(elements: list[dict[str, Any]]) -> dict[int, tuple[float, float]]:
    nodes: dict[int, tuple[float, float]] = {}
    for element in elements:
        if element.get("type") != "node":
            continue
        try:
            node_id = int(element["id"])
            lon = float(element["lon"])
            lat = float(element["lat"])
        except (KeyError, TypeError, ValueError):
            continue
        nodes[node_id] = (lon, lat)
    return nodes


def way_lookup(elements: list[dict[str, Any]], nodes: dict[int, tuple[float, float]]) -> dict[int, tuple[tuple[float, float], ...]]:
    ways: dict[int, tuple[tuple[float, float], ...]] = {}
    for element in elements:
        if element.get("type") != "way":
            continue
        try:
            way_id = int(element["id"])
        except (KeyError, TypeError, ValueError):
            continue
        coords: list[tuple[float, float]] = []
        for node_id in element.get("nodes", []):
            try:
                coords.append(nodes[int(node_id)])
            except (KeyError, TypeError, ValueError):
                continue
        ways[way_id] = tuple(coords)
    return ways


def relation_lookup(elements: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    relations: dict[int, dict[str, Any]] = {}
    for element in elements:
        if element.get("type") != "relation":
            continue
        try:
            relation_id = int(element["id"])
        except (KeyError, TypeError, ValueError):
            continue
        relations[relation_id] = element
    return relations


def is_polygon_water_way(tags: dict[str, str], points: tuple[tuple[float, float], ...]) -> bool:
    if not is_closed_ring(points):
        return False
    for key, value in POLYGON_WATER_KEYS:
        if tags.get(key) == value:
            return True
    return False


def polygon_kind_from_tags(tags: dict[str, str]) -> str:
    if tags.get("natural") == "water":
        return "natural_water"
    if tags.get("waterway") == "riverbank":
        return "riverbank"
    return "water_polygon"


def relation_is_water_relation(tags: dict[str, str]) -> bool:
    for key, value in POLYGON_WATER_KEYS:
        if tags.get(key) == value:
            return True
    return False


def extract_water_features(elements: list[dict[str, Any]]) -> list[WaterPolygon]:
    nodes = node_lookup(elements)
    ways = way_lookup(elements, nodes)
    relations = relation_lookup(elements)

    polygons: list[WaterPolygon] = []

    way_tags: dict[int, dict[str, str]] = {}
    for element in elements:
        if element.get("type") != "way":
            continue
        try:
            way_id = int(element["id"])
        except (KeyError, TypeError, ValueError):
            continue
        way_tags[way_id] = tags_to_dict(element.get("tags"))

    for way_id, points in ways.items():
        tags = way_tags.get(way_id, {})
        if is_polygon_water_way(tags, points):
            polygons.append(
                WaterPolygon(
                    osm_id=f"way/{way_id}",
                    kind=polygon_kind_from_tags(tags),
                    outer_rings=(points,),
                    inner_rings=(),
                    source="way",
                    tags=tags,
                )
            )

    for relation_id, relation in relations.items():
        tags = tags_to_dict(relation.get("tags"))
        if not relation_is_water_relation(tags):
            continue
        if tags.get("type", "multipolygon") != "multipolygon":
            continue

        outer_segments: list[tuple[tuple[float, float], ...]] = []
        inner_segments: list[tuple[tuple[float, float], ...]] = []
        for member in relation.get("members", []):
            if not isinstance(member, dict):
                continue
            if member.get("type") != "way":
                continue
            role = member.get("role")
            try:
                member_way_id = int(member["ref"])
            except (KeyError, TypeError, ValueError):
                continue
            member_points = ways.get(member_way_id)
            if not member_points:
                continue
            if role == "inner":
                inner_segments.append(member_points)
            else:
                outer_segments.append(member_points)

        outer_rings = assemble_rings_from_segments(outer_segments)
        inner_rings = assemble_rings_from_segments(inner_segments)
        if not outer_rings:
            continue
        polygons.append(
            WaterPolygon(
                osm_id=f"relation/{relation_id}",
                kind=polygon_kind_from_tags(tags),
                outer_rings=tuple(outer_rings),
                inner_rings=tuple(inner_rings),
                source="relation",
                tags=tags,
            )
        )

    return polygons


def geometry_bounds(
    polygons: list[WaterPolygon],
    fallback_bbox: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    for polygon in polygons:
        for ring in list(polygon.outer_rings) + list(polygon.inner_rings):
            for lon, lat in ring:
                xs.append(lon)
                ys.append(lat)
    if not xs or not ys:
        return fallback_bbox
    return min(xs), min(ys), max(xs), max(ys)


def project_lonlat(
    lon_deg: float,
    lat_deg: float,
    *,
    bbox: tuple[float, float, float, float],
    width: int,
    height: int,
    padding: float,
) -> tuple[float, float]:
    west, south, east, north = bbox
    lon_span = max(1.0e-9, east - west)
    lat_span = max(1.0e-9, north - south)
    pad_x = lon_span * float(padding)
    pad_y = lat_span * float(padding)
    west -= pad_x
    east += pad_x
    south -= pad_y
    north += pad_y
    lon_span = max(1.0e-9, east - west)
    lat_span = max(1.0e-9, north - south)
    x = (lon_deg - west) / lon_span * float(width)
    y = (north - lat_deg) / lat_span * float(height)
    return x, y


def project_lonlat_to_local_m(
    lon_deg: float,
    lat_deg: float,
    *,
    center_lon_deg: float,
    center_lat_deg: float,
) -> tuple[float, float]:
    center_lat_rad = math.radians(center_lat_deg)
    x = math.radians(lon_deg - center_lon_deg) * EARTH_RADIUS_M * math.cos(center_lat_rad)
    y = math.radians(lat_deg - center_lat_deg) * EARTH_RADIUS_M
    return x, y


def project_local_m_to_lonlat(
    x_m: float,
    y_m: float,
    *,
    center_lon_deg: float,
    center_lat_deg: float,
) -> tuple[float, float]:
    center_lat_rad = math.radians(center_lat_deg)
    lon_scale = EARTH_RADIUS_M * math.cos(center_lat_rad)
    if abs(lon_scale) < 1.0e-9:
        lon_deg = center_lon_deg
    else:
        lon_deg = center_lon_deg + math.degrees(x_m / lon_scale)
    lat_deg = center_lat_deg + math.degrees(y_m / EARTH_RADIUS_M)
    return lon_deg, lat_deg


def ring_to_svg_points(
    ring: tuple[tuple[float, float], ...],
    *,
    bbox: tuple[float, float, float, float],
    width: int,
    height: int,
    padding: float,
) -> list[tuple[float, float]]:
    return [
        project_lonlat(lon, lat, bbox=bbox, width=width, height=height, padding=padding)
        for lon, lat in ring
    ]


def ring_to_local_points(
    ring: tuple[tuple[float, float], ...],
    *,
    center_lon_deg: float,
    center_lat_deg: float,
) -> list[tuple[float, float]]:
    return [
        project_lonlat_to_local_m(lon, lat, center_lon_deg=center_lon_deg, center_lat_deg=center_lat_deg)
        for lon, lat in ring
    ]


def local_points_to_ring(
    points: Iterable[tuple[float, float]],
    *,
    center_lon_deg: float,
    center_lat_deg: float,
) -> tuple[tuple[float, float], ...]:
    return tuple(
        project_local_m_to_lonlat(x, y, center_lon_deg=center_lon_deg, center_lat_deg=center_lat_deg)
        for x, y in points
    )


def path_for_points(points: list[tuple[float, float]]) -> str:
    if not points:
        return ""
    return " ".join(
        [
            f"M {points[0][0]:.2f} {points[0][1]:.2f}",
            *[f"L {x:.2f} {y:.2f}" for x, y in points[1:]],
            "Z",
        ]
    )


def local_points_bounds(points_groups: list[list[tuple[float, float]]]) -> tuple[float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    for points in points_groups:
        for x, y in points:
            xs.append(x)
            ys.append(y)
    if not xs or not ys:
        return 0.0, 0.0, 1.0, 1.0
    return min(xs), min(ys), max(xs), max(ys)


def _import_shapely_geometry() -> tuple[Any, Any]:
    try:
        from shapely.geometry import Polygon, box
    except ImportError as exc:  # pragma: no cover - optional dev dependency
        raise RuntimeError(
            "shapely is required for grid-split triangulation; install the dev extra and rerun uv sync"
        ) from exc
    return Polygon, box


def _geometry_polygons(geometry: Any) -> list[Any]:
    geom_type = getattr(geometry, "geom_type", "")
    if geom_type == "Polygon":
        return [geometry]
    if geom_type in {"MultiPolygon", "GeometryCollection"}:
        parts: list[Any] = []
        for part in getattr(geometry, "geoms", []):
            parts.extend(_geometry_polygons(part))
        return parts
    return []


def _build_local_polygon_geometry(
    outer_ring: list[tuple[float, float]],
    hole_rings: list[list[tuple[float, float]]],
) -> Any | None:
    polygon_cls, _ = _import_shapely_geometry()
    if len(outer_ring) < 4:
        return None
    polygon = polygon_cls(outer_ring, holes=hole_rings or None)
    if polygon.is_empty or float(polygon.area) <= 0.0:
        return None
    if not polygon.is_valid:
        repaired = polygon.buffer(0.0)
        if not repaired.is_empty:
            polygon = repaired
    if polygon.is_empty or float(polygon.area) <= 0.0:
        return None
    return polygon


def _ring_body(points_xy: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
    ring = list(points_xy)
    if len(ring) >= 2 and ring[0] == ring[-1]:
        ring = ring[:-1]
    if len(ring) < 3:
        return []
    return ring


def _ring_signed_area(points_xy: Iterable[tuple[float, float]]) -> float:
    body = _ring_body(points_xy)
    if len(body) < 3:
        return 0.0
    area = 0.0
    for index in range(len(body)):
        x0, y0 = body[index]
        x1, y1 = body[(index + 1) % len(body)]
        area += x0 * y1 - x1 * y0
    return 0.5 * area


def _point_in_ring(point: tuple[float, float], ring: Iterable[tuple[float, float]]) -> bool:
    body = _ring_body(ring)
    if len(body) < 3:
        return False
    x, y = point
    inside = False
    for index in range(len(body)):
        x0, y0 = body[index]
        x1, y1 = body[(index + 1) % len(body)]
        if ((y0 > y) != (y1 > y)) and (y1 != y0):
            x_cross = ((x1 - x0) * (y - y0) / (y1 - y0)) + x0
            if x < x_cross:
                inside = not inside
    return inside


def _split_polygon_by_grid(polygon: Any, cell_m: float) -> list[Any]:
    if cell_m <= 0.0:
        return [polygon]
    _, box = _import_shapely_geometry()
    min_x, min_y, max_x, max_y = polygon.bounds
    if min_x >= max_x or min_y >= max_y:
        return []
    start_x = math.floor(min_x / cell_m)
    end_x = math.floor(max_x / cell_m)
    start_y = math.floor(min_y / cell_m)
    end_y = math.floor(max_y / cell_m)
    pieces: list[Any] = []
    seen_keys: set[tuple[float, float, float, float]] = set()
    for cell_x in range(int(start_x), int(end_x) + 1):
        for cell_y in range(int(start_y), int(end_y) + 1):
            min_cell_x = float(cell_x) * cell_m
            min_cell_y = float(cell_y) * cell_m
            max_cell_x = min_cell_x + cell_m
            max_cell_y = min_cell_y + cell_m
            cell_key = (min_cell_x, min_cell_y, max_cell_x, max_cell_y)
            if cell_key in seen_keys:
                continue
            seen_keys.add(cell_key)
            cell = box(min_cell_x, min_cell_y, max_cell_x, max_cell_y)
            clipped = polygon.intersection(cell)
            for part in _geometry_polygons(clipped):
                if part.is_empty or float(part.area) <= 0.0:
                    continue
                pieces.append(part)
    return pieces


def _assign_holes_to_shells(
    outer_rings: list[list[tuple[float, float]]],
    inner_rings: list[list[tuple[float, float]]],
) -> list[list[list[tuple[float, float]]]]:
    shell_holes: list[list[list[tuple[float, float]]]] = [[] for _ in outer_rings]
    shell_areas = []
    for shell in outer_rings:
        shell_areas.append(abs(_ring_signed_area(shell)))
    for hole in inner_rings:
        if len(hole) < 4:
            continue
        anchor = hole[0]
        matches: list[int] = []
        for shell_index, shell in enumerate(outer_rings):
            if _point_in_ring(anchor, shell):
                matches.append(shell_index)
        if not matches:
            continue
        target_index = min(matches, key=lambda index: shell_areas[index])
        shell_holes[target_index].append(hole)
    return shell_holes


def build_local_triangulation_pieces(
    polygon: WaterPolygon,
    *,
    center_lon_deg: float,
    center_lat_deg: float,
    cell_m: float,
) -> list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]]:
    outer_rings = [
        ring_to_local_points(
            ring,
            center_lon_deg=center_lon_deg,
            center_lat_deg=center_lat_deg,
        )
        for ring in polygon.outer_rings
    ]
    inner_rings = [
        ring_to_local_points(
            ring,
            center_lon_deg=center_lon_deg,
            center_lat_deg=center_lat_deg,
        )
        for ring in polygon.inner_rings
    ]
    shell_holes = _assign_holes_to_shells(outer_rings, inner_rings)
    pieces: list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]] = []
    for shell, holes in zip(outer_rings, shell_holes):
        local_polygon = _build_local_polygon_geometry(shell, holes)
        if local_polygon is None:
            continue
        for piece in _split_polygon_by_grid(local_polygon, cell_m):
            if piece.is_empty or float(piece.area) <= 0.0:
                continue
            piece_shell = list(piece.exterior.coords)
            piece_holes = [list(interior.coords) for interior in piece.interiors if len(interior.coords) >= 4]
            if len(piece_shell) >= 4:
                pieces.append((piece_shell, piece_holes))
    return pieces


def project_local_points_to_svg(
    points: list[tuple[float, float]],
    *,
    bounds: tuple[float, float, float, float],
    width: int,
    height: int,
    padding: float,
) -> list[tuple[float, float]]:
    west, south, east, north = bounds
    x_span = max(1.0e-9, east - west)
    y_span = max(1.0e-9, north - south)
    pad_x = x_span * float(padding)
    pad_y = y_span * float(padding)
    west -= pad_x
    east += pad_x
    south -= pad_y
    north += pad_y
    x_span = max(1.0e-9, east - west)
    y_span = max(1.0e-9, north - south)
    projected: list[tuple[float, float]] = []
    for x, y in points:
        svg_x = (x - west) / x_span * float(width)
        svg_y = (north - y) / y_span * float(height)
        projected.append((svg_x, svg_y))
    return projected


def svg_path_for_ring(
    ring: tuple[tuple[float, float], ...],
    *,
    bbox: tuple[float, float, float, float],
    width: int,
    height: int,
    padding: float,
) -> str:
    if not ring:
        return ""
    points = ring_to_svg_points(ring, bbox=bbox, width=width, height=height, padding=padding)
    return path_for_points(points)


def build_triangulated_svg_preview(
    polygons: list[WaterPolygon],
    *,
    bbox: tuple[float, float, float, float],
    width: int,
    height: int,
    padding: float,
    center_lon_deg: float,
    center_lat_deg: float,
    grid_m: float,
    simplify_m: float,
    cell_m: float,
) -> str:
    background_fill = "#f7fbff"
    water_fill = "#8ecae6"
    water_stroke = "#4a90c2"

    local_rings: list[list[tuple[float, float]]] = []
    for polygon in polygons:
        for ring in polygon.outer_rings:
            local_rings.append(
                ring_to_local_points(
                    ring,
                    center_lon_deg=center_lon_deg,
                    center_lat_deg=center_lat_deg,
                )
            )
        for ring in polygon.inner_rings:
            local_rings.append(
                ring_to_local_points(
                    ring,
                    center_lon_deg=center_lon_deg,
                    center_lat_deg=center_lat_deg,
                )
            )

    bounds = local_points_bounds(local_rings)

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">'
        ),
        f'<rect x="0" y="0" width="100%" height="100%" fill="{background_fill}"/>',
    ]
    for polygon in polygons:
        if cell_m > 0.0:
            pieces = build_local_triangulation_pieces(
                polygon,
                center_lon_deg=center_lon_deg,
                center_lat_deg=center_lat_deg,
                cell_m=cell_m,
            )
            for shell_points, hole_points_list in pieces:
                footprint = WaterPolygonFootprint(
                    water_id=f"{polygon.osm_id}@cell",
                    kind=polygon.kind,
                    outer_rings_lonlat=(
                        local_points_to_ring(
                            shell_points,
                            center_lon_deg=center_lon_deg,
                            center_lat_deg=center_lat_deg,
                        ),
                    ),
                    inner_rings_lonlat=tuple(
                        local_points_to_ring(
                            hole_points,
                            center_lon_deg=center_lon_deg,
                            center_lat_deg=center_lat_deg,
                        )
                        for hole_points in hole_points_list
                    ),
                    source=polygon.source,
                    tags=polygon.tags,
                )
                mesh = build_water_surface_mesh(
                    footprint,
                    center_lat_deg=center_lat_deg,
                    center_lon_deg=center_lon_deg,
                    grid_m=grid_m,
                    simplify_tolerance_m=simplify_m,
                )
                if mesh is None:
                    continue
                for triangle in mesh.triangles_xy_m:
                    points = project_local_points_to_svg(
                        list(triangle),
                        bounds=bounds,
                        width=width,
                        height=height,
                        padding=padding,
                    )
                    triangle_path_d = path_for_points(points)
                    if not triangle_path_d:
                        continue
                    parts.append(
                        (
                            '<path d="{d}" fill="{fill}" fill-opacity="0.50" '
                            'stroke="{stroke}" stroke-width="0.8" vector-effect="non-scaling-stroke"/>'
                        ).format(
                            d=triangle_path_d,
                            fill=water_fill,
                            stroke=water_stroke,
                        )
                    )
        else:
            footprint = WaterPolygonFootprint(
                water_id=polygon.osm_id,
                kind=polygon.kind,
                outer_rings_lonlat=polygon.outer_rings,
                inner_rings_lonlat=polygon.inner_rings,
                source=polygon.source,
                tags=polygon.tags,
            )
            mesh = build_water_surface_mesh(
                footprint,
                center_lat_deg=center_lat_deg,
                center_lon_deg=center_lon_deg,
                grid_m=grid_m,
                simplify_tolerance_m=simplify_m,
            )
            if mesh is None:
                continue
            for triangle in mesh.triangles_xy_m:
                points = project_local_points_to_svg(
                    list(triangle),
                    bounds=bounds,
                    width=width,
                    height=height,
                    padding=padding,
                )
                triangle_path_d = path_for_points(points)
                if not triangle_path_d:
                    continue
                parts.append(
                    (
                        '<path d="{d}" fill="{fill}" fill-opacity="0.50" '
                        'stroke="{stroke}" stroke-width="0.8" vector-effect="non-scaling-stroke"/>'
                    ).format(
                        d=triangle_path_d,
                        fill=water_fill,
                        stroke=water_stroke,
                    )
                )
    parts.append("</svg>")
    return "\n".join(parts)


def build_svg_preview(
    polygons: list[WaterPolygon],
    *,
    bbox: tuple[float, float, float, float],
    width: int,
    height: int,
    padding: float,
) -> str:
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">'
        ),
        '<rect x="0" y="0" width="100%" height="100%" fill="#f7fbff"/>',
    ]
    for polygon in polygons:
        if not polygon.outer_rings:
            continue
        path_d = " ".join(
            svg_path_for_ring(ring, bbox=bbox, width=width, height=height, padding=padding)
            for ring in polygon.outer_rings
        )
        path_d += "".join(
            f" {svg_path_for_ring(ring, bbox=bbox, width=width, height=height, padding=padding)}"
            for ring in polygon.inner_rings
        )
        parts.append(
            (
                '<path d="{d}" fill="#8ecae6" fill-opacity="0.45" '
                'stroke="#4a90c2" stroke-width="1.2" vector-effect="non-scaling-stroke" '
                'fill-rule="evenodd"/>'
            ).format(d=path_d)
        )
    parts.append("</svg>")
    return "\n".join(parts)


def serialize_water_layer(
    polygons: list[WaterPolygon],
    *,
    bbox: tuple[float, float, float, float],
    radius_km: float,
    endpoint: str,
    center_lat_deg: float,
    center_lon_deg: float,
) -> dict[str, Any]:
    return {
        "query": {
            "center": {"lat_deg": float(center_lat_deg), "lon_deg": float(center_lon_deg)},
            "radius_km": float(radius_km),
            "bbox": {
                "west": bbox[0],
                "south": bbox[1],
                "east": bbox[2],
                "north": bbox[3],
            },
            "source": endpoint,
        },
        "water_polygons": [
            {
                "osm_id": polygon.osm_id,
                "kind": polygon.kind,
                "source": polygon.source,
                "tags": polygon.tags,
                "outer_rings": [[list(point) for point in ring] for ring in polygon.outer_rings],
                "inner_rings": [[list(point) for point in ring] for ring in polygon.inner_rings],
            }
            for polygon in polygons
        ],
        "summary": {
            "polygon_count": len(polygons),
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    bbox = bbox_from_point(args.lat, args.lon, float(args.radius_km))
    query = build_overpass_query(bbox)
    payload = fetch_overpass_json(
        endpoint=args.endpoint,
        query=query,
        user_agent=args.user_agent,
        timeout_s=float(args.timeout_s),
    )
    elements = payload.get("elements")
    if not isinstance(elements, list):
        raise RuntimeError("Overpass payload is missing the elements array")

    polygons = extract_water_features(elements)
    normalized = serialize_water_layer(
        polygons,
        bbox=bbox,
        radius_km=float(args.radius_km),
        endpoint=args.endpoint,
        center_lat_deg=float(args.lat),
        center_lon_deg=float(args.lon),
    )

    json_text = json.dumps(normalized, ensure_ascii=True, indent=2)
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json_text + "\n", encoding="utf-8")
        print(f"wrote_json={args.output_json}")
    else:
        print(json_text)

    if args.output_svg is not None:
        args.output_svg.parent.mkdir(parents=True, exist_ok=True)
        svg_text = build_svg_preview(
            polygons,
            bbox=geometry_bounds(polygons, bbox),
            width=int(args.svg_width),
            height=int(args.svg_height),
            padding=float(args.svg_padding),
        )
        args.output_svg.write_text(svg_text + "\n", encoding="utf-8")
        print(f"wrote_svg={args.output_svg}")

    if args.output_triangulated_svg is not None:
        args.output_triangulated_svg.parent.mkdir(parents=True, exist_ok=True)
        svg_text = build_triangulated_svg_preview(
            polygons,
            bbox=geometry_bounds(polygons, bbox),
            width=int(args.svg_width),
            height=int(args.svg_height),
            padding=float(args.svg_padding),
            center_lon_deg=float(args.lon),
            center_lat_deg=float(args.lat),
            grid_m=float(args.triangulation_grid_m),
            simplify_m=float(args.triangulation_simplify_m),
            cell_m=float(args.triangulation_cell_m),
        )
        args.output_triangulated_svg.write_text(svg_text + "\n", encoding="utf-8")
        print(f"wrote_triangulated_svg={args.output_triangulated_svg}")

    print(
        f"summary polygons={len(polygons)} bbox={bbox[0]:.6f},{bbox[1]:.6f},{bbox[2]:.6f},{bbox[3]:.6f}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
