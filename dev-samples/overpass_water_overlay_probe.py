#!/usr/bin/env python3
"""Fetch and normalize nearby OSM water polygons for overlay prototyping.

This sample keeps the pipeline intentionally simple and dependency-free:

1. Query Overpass around a latitude/longitude center point.
2. Reconstruct closed polygon rings from OSM nodes and ways.
3. Keep only water polygons:
   - outer rings
   - inner rings
4. Write a normalized JSON payload and, optionally, a quick SVG preview.
5. Optionally triangulate the polygon rings into SVG triangles for shape checks.

The script is meant for data-structure exploration, not production caching.
It uses only the Python standard library so it can be run directly.
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

try:
    from shapely import Point, Polygon, constrained_delaunay_triangles, set_precision, simplify
except ImportError:  # pragma: no cover - optional dev dependency
    Point = None  # type: ignore[assignment]
    Polygon = None  # type: ignore[assignment]
    constrained_delaunay_triangles = None  # type: ignore[assignment]
    set_precision = None  # type: ignore[assignment]
    simplify = None  # type: ignore[assignment]


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


def triangle_geometry_to_svg_path(
    triangle: Any,
    *,
    bounds: tuple[float, float, float, float],
    width: int,
    height: int,
    padding: float,
) -> str:
    coords = list(triangle.exterior.coords)
    if len(coords) < 4:
        return ""
    if coords[0] == coords[-1]:
        coords = coords[:-1]
    points = project_local_points_to_svg(coords, bounds=bounds, width=width, height=height, padding=padding)
    return path_for_points(points)


def build_local_polygon(
    shell_points: list[tuple[float, float]],
    hole_points_list: list[list[tuple[float, float]]],
    *,
    grid_m: float,
    simplify_m: float,
) -> Any:
    if Polygon is None:
        raise RuntimeError("shapely is required for triangulated SVG output")
    polygon = Polygon(shell_points, holes=hole_points_list or None)
    if polygon.is_empty or not polygon.is_valid or polygon.area <= 0.0:
        return None
    if set_precision is not None and grid_m > 0.0:
        polygon = set_precision(polygon, grid_m)
    if simplify is not None and simplify_m > 0.0:
        polygon = simplify(polygon, simplify_m, preserve_topology=True)
    if polygon.is_empty or not polygon.is_valid or polygon.area <= 0.0:
        return None
    return polygon


def ring_to_triangle_paths(
    polygon_points: list[tuple[float, float]],
    hole_points_list: list[list[tuple[float, float]]],
    *,
    bounds: tuple[float, float, float, float],
    width: int,
    height: int,
    padding: float,
    grid_m: float,
    simplify_m: float,
) -> list[str]:
    if Polygon is None or constrained_delaunay_triangles is None:
        raise RuntimeError("shapely is required for triangulated SVG output")
    if len(polygon_points) < 3:
        return []
    polygon = build_local_polygon(
        polygon_points,
        hole_points_list,
        grid_m=grid_m,
        simplify_m=simplify_m,
    )
    if polygon is None:
        return []
    triangles = constrained_delaunay_triangles(polygon)
    paths: list[str] = []
    for triangle in triangles.geoms:
        if triangle.is_empty or triangle.geom_type != "Polygon":
            continue
        path_d = triangle_geometry_to_svg_path(
            triangle,
            bounds=bounds,
            width=width,
            height=height,
            padding=padding,
        )
        if path_d:
            paths.append(path_d)
    return paths


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
) -> str:
    if Polygon is None or constrained_delaunay_triangles is None:
        raise RuntimeError(
            "shapely is required for triangulated SVG output; install the dev extra and rerun uv sync"
        )
    background_fill = "#f7fbff"
    water_fill = "#8ecae6"
    water_stroke = "#4a90c2"

    shell_entries: list[tuple[list[tuple[float, float]], Any]] = []
    hole_entries: list[list[tuple[float, float]]] = []
    for polygon in polygons:
        for ring in polygon.outer_rings:
            points = ring_to_local_points(
                ring,
                center_lon_deg=center_lon_deg,
                center_lat_deg=center_lat_deg,
            )
            shell_geom = build_local_polygon(points, [], grid_m=grid_m, simplify_m=simplify_m)
            if shell_geom is None:
                continue
            shell_entries.append((points, shell_geom))
        for ring in polygon.inner_rings:
            points = ring_to_local_points(
                ring,
                center_lon_deg=center_lon_deg,
                center_lat_deg=center_lat_deg,
            )
            if len(points) >= 3:
                hole_entries.append(points)

    bounds = local_points_bounds([points for points, _ in shell_entries] + hole_entries)

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">'
        ),
        f'<rect x="0" y="0" width="100%" height="100%" fill="{background_fill}"/>',
    ]
    shell_holes: list[list[list[tuple[float, float]]]] = [[] for _ in shell_entries]
    for points in hole_entries:
        if len(points) < 3:
            continue
        assigned = False
        hole_anchor = Point(points[0])
        best_index = -1
        best_area = float("inf")
        for index, (_, shell_geom) in enumerate(shell_entries):
            if not shell_geom.covers(hole_anchor):
                continue
            shell_area = float(shell_geom.area)
            if shell_area < best_area:
                best_index = index
                best_area = shell_area
        if best_index >= 0:
            shell_holes[best_index].append(points)
            assigned = True
        if not assigned:
            fallback_polygon = build_local_polygon(points, [], grid_m=grid_m, simplify_m=simplify_m)
            if fallback_polygon is None:
                continue
            for triangle in constrained_delaunay_triangles(fallback_polygon).geoms:
                if triangle.is_empty or triangle.geom_type != "Polygon":
                    continue
                path_d = triangle_geometry_to_svg_path(
                    triangle,
                    bounds=bounds,
                    width=width,
                    height=height,
                    padding=padding,
                )
                if path_d:
                    parts.append(
                        (
                            '<path d="{d}" fill="{fill}" fill-opacity="1.0" '
                            'stroke="{stroke}" stroke-width="0.8" vector-effect="non-scaling-stroke"/>'
                        ).format(
                            d=path_d,
                            fill=background_fill,
                            stroke="#dbe9f6",
                        )
                    )
    for index, (points, _) in enumerate(shell_entries):
        hole_list = shell_holes[index]
        for triangle_path_d in ring_to_triangle_paths(
            points,
            hole_list,
            bounds=bounds,
            width=width,
            height=height,
            padding=padding,
            grid_m=grid_m,
            simplify_m=simplify_m,
        ):
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
