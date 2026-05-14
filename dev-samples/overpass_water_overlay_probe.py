#!/usr/bin/env python3
"""Fetch and normalize nearby OSM water polygons for overlay prototyping.

This sample keeps the pipeline intentionally simple:

1. Query Overpass around a latitude/longitude center point.
2. Reconstruct closed polygon rings from OSM nodes and ways.
3. Keep only water polygons:
   - outer rings
   - inner rings
4. Write a normalized JSON payload and, optionally, a quick SVG preview.
   The SVG preview can also show the observer marker and distance rings.

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

from zstarview.water_overlay import (
    classify_water_surface_category,
    extract_water_polygons as extract_core_water_polygons,
)

EARTH_RADIUS_KM = 6371.0088
EARTH_RADIUS_M = EARTH_RADIUS_KM * 1000.0
DEFAULT_ENDPOINT = "https://overpass-api.de/api/interpreter"
DEFAULT_USER_AGENT = "zstarview-water-overlay-probe/0.1"
DEFAULT_TIMEOUT_S = 60.0
QUERY_MARGIN_KM = 2.0
OVERPASS_SAMPLE_ATTEMPTS = 8

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
        help="Optional output path for the simplified SVG preview.",
    )
    parser.add_argument(
        "--output-svg-raw",
        type=Path,
        help="Optional output path for the unsimplified SVG preview.",
    )
    parser.add_argument(
        "--input-cache",
        type=Path,
        help=(
            "Optional input cache JSON written by zstarview. "
            "When set, the script skips Overpass and renders polygons from the cache."
        ),
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
        "--simplify-min-distance-km",
        type=float,
        default=1.0,
        help="Minimum distance before nearby polygons can be merged in the simplified SVG (default: 1.0).",
    )
    parser.add_argument(
        "--simplify-max-distance-gap-km",
        type=float,
        default=0.5,
        help="Maximum distance gap between adjacent far polygons for simplification (default: 0.5).",
    )
    parser.add_argument(
        "--simplify-max-offset-m",
        type=float,
        default=400.0,
        help="Maximum local offset between adjacent far polygons for simplification (default: 400.0).",
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


def extract_water_features(
    elements: list[dict[str, Any]],
    *,
    bbox: tuple[float, float, float, float] | None,
) -> list[WaterPolygon]:
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

    if bbox is not None:
        core_polygons = extract_core_water_polygons(elements, bbox=bbox)
        for footprint in core_polygons:
            if footprint.kind != "coastline":
                continue
            polygons.append(
                WaterPolygon(
                    osm_id=str(footprint.water_id),
                    kind=str(footprint.kind),
                    outer_rings=tuple(
                        tuple(tuple(point) for point in ring)
                        for ring in footprint.outer_rings_lonlat
                    ),
                    inner_rings=tuple(
                        tuple(tuple(point) for point in ring)
                        for ring in footprint.inner_rings_lonlat
                    ),
                    source=str(footprint.source),
                    tags=dict(footprint.tags),
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


def _circle_ring_lonlat(
    *,
    center_lon_deg: float,
    center_lat_deg: float,
    radius_km: float,
    sample_count: int = 96,
) -> tuple[tuple[float, float], ...]:
    radius_m = max(0.0, float(radius_km)) * 1000.0
    if radius_m <= 0.0:
        return ((float(center_lon_deg), float(center_lat_deg)),)
    sample_count = max(8, int(sample_count))
    points: list[tuple[float, float]] = []
    for index in range(sample_count):
        angle = (2.0 * math.pi * index) / float(sample_count)
        x_m = math.cos(angle) * radius_m
        y_m = math.sin(angle) * radius_m
        lon_deg, lat_deg = project_local_m_to_lonlat(
            x_m,
            y_m,
            center_lon_deg=float(center_lon_deg),
            center_lat_deg=float(center_lat_deg),
        )
        points.append((lon_deg, lat_deg))
    points.append(points[0])
    return tuple(points)


def _build_distance_ring_values(radius_km: float) -> list[float]:
    if radius_km <= 0.0:
        return []
    max_ring = max(1, int(math.floor(float(radius_km))))
    rings = [float(index) for index in range(1, max_ring + 1)]
    if not rings or abs(rings[-1] - float(radius_km)) > 1.0e-6:
        rings.append(float(radius_km))
    return rings


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


def _ring_bounds(ring: tuple[tuple[float, float], ...]) -> tuple[float, float, float, float]:
    min_x = float("inf")
    min_y = float("inf")
    max_x = float("-inf")
    max_y = float("-inf")
    for x, y in ring:
        min_x = min(min_x, float(x))
        min_y = min(min_y, float(y))
        max_x = max(max_x, float(x))
        max_y = max(max_y, float(y))
    if not math.isfinite(min_x):
        return 0.0, 0.0, 0.0, 0.0
    return min_x, min_y, max_x, max_y


def _polygon_overlaps_bbox(polygon: WaterPolygon, bbox: tuple[float, float, float, float]) -> bool:
    west, south, east, north = bbox
    for ring in polygon.outer_rings + polygon.inner_rings:
        min_x, min_y, max_x, max_y = _ring_bounds(ring)
        if max_x < west or max_y < south or min_x > east or min_y > north:
            continue
        return True
    return False


def _water_style(polygon: WaterPolygon) -> tuple[str, str, float]:
    category = classify_water_surface_category(polygon.tags, kind=polygon.kind)
    if category == "sea":
        return "#f59e0b", "#c2410c", 0.42
    if category == "river":
        return "#4db6ff", "#2b7fbf", 0.55
    return "#8ecae6", "#4a90c2", 0.50


def _bbox_union(
    a: tuple[float, float, float, float],
    b: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    return (
        min(float(a[0]), float(b[0])),
        min(float(a[1]), float(b[1])),
        max(float(a[2]), float(b[2])),
        max(float(a[3]), float(b[3])),
    )


def _polygon_preview_anchor(
    polygon: WaterPolygon,
    *,
    center_lon_deg: float,
    center_lat_deg: float,
) -> tuple[float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    for ring in polygon.outer_rings:
        for lon, lat in ring:
            x_m, y_m = project_lonlat_to_local_m(
                float(lon),
                float(lat),
                center_lon_deg=center_lon_deg,
                center_lat_deg=center_lat_deg,
            )
            xs.append(x_m)
            ys.append(y_m)
    if not xs or not ys:
        return float("inf"), 0.0, 0.0
    anchor_x = sum(xs) / len(xs)
    anchor_y = sum(ys) / len(ys)
    return math.hypot(anchor_x, anchor_y) / 1000.0, anchor_x, anchor_y


def _polygon_vertex_count(polygon: WaterPolygon) -> int:
    return sum(len(ring) for ring in polygon.outer_rings) + sum(len(ring) for ring in polygon.inner_rings)


def _simplification_thresholds_for_distance_km(distance_km: float) -> tuple[float, float, float]:
    if distance_km < 1.0:
        return float("inf"), 0.0, 0.0
    if distance_km < 2.0:
        return 1.0, 0.20, 120.0
    if distance_km < 3.5:
        return 1.2, 0.35, 220.0
    return 1.5, 0.55, 420.0


def _dedupe_adjacent_far_polygons(
    polygons: list[WaterPolygon],
    *,
    center_lon_deg: float,
    center_lat_deg: float,
    min_distance_km: float = 1.0,
    max_neighbor_distance_gap_km: float = 0.5,
    max_neighbor_offset_m: float = 400.0,
) -> tuple[list[WaterPolygon], dict[str, int]]:
    if len(polygons) < 2:
        kept = list(polygons)
        vertex_count = sum(_polygon_vertex_count(polygon) for polygon in kept)
        return kept, {
            "raw_polygons": len(kept),
            "kept_polygons": len(kept),
            "removed_polygons": 0,
            "raw_vertices": vertex_count,
            "kept_vertices": vertex_count,
            "removed_vertices": 0,
        }

    scored: list[tuple[float, float, float, int, WaterPolygon]] = []
    for index, polygon in enumerate(polygons):
        distance_km, anchor_x, anchor_y = _polygon_preview_anchor(
            polygon,
            center_lon_deg=center_lon_deg,
            center_lat_deg=center_lat_deg,
        )
        scored.append((distance_km, anchor_x, anchor_y, index, polygon))
    scored.sort(key=lambda item: (item[0], item[1], item[2], item[3]))

    deduped: list[WaterPolygon] = []
    previous: tuple[float, float, float, int, WaterPolygon] | None = None
    dropped = 0
    removed_vertices = 0
    for item in scored:
        if previous is not None:
            previous_distance_km, previous_x_m, previous_y_m, _, _ = previous
            distance_km, x_m, y_m, _, _ = item
            effective_distance_km = max(previous_distance_km, distance_km)
            tier_min_distance_km, tier_gap_km, tier_offset_m = _simplification_thresholds_for_distance_km(
                effective_distance_km
            )
            tier_min_distance_km = max(tier_min_distance_km, float(min_distance_km))
            tier_gap_km = max(tier_gap_km, float(max_neighbor_distance_gap_km))
            tier_offset_m = max(tier_offset_m, float(max_neighbor_offset_m))
            is_far = previous_distance_km >= tier_min_distance_km and distance_km >= tier_min_distance_km
            close_in_range = abs(distance_km - previous_distance_km) <= tier_gap_km
            close_in_space = math.hypot(x_m - previous_x_m, y_m - previous_y_m) <= tier_offset_m
            if is_far and close_in_range and close_in_space:
                dropped += 1
                removed_vertices += _polygon_vertex_count(item[4])
                continue
        deduped.append(item[4])
        previous = item

    raw_vertices = sum(_polygon_vertex_count(polygon) for polygon in polygons)
    kept_vertices = sum(_polygon_vertex_count(polygon) for polygon in deduped)
    if dropped:
        print(
            "deduped_far_polygons={dropped} raw_polygons={raw_polygons} kept_polygons={kept_polygons} "
            "raw_vertices={raw_vertices} kept_vertices={kept_vertices} removed_vertices={removed_vertices}".format(
                dropped=dropped,
                raw_polygons=len(polygons),
                kept_polygons=len(deduped),
                raw_vertices=raw_vertices,
                kept_vertices=kept_vertices,
                removed_vertices=removed_vertices,
            ),
            file=sys.stderr,
        )
    return deduped, {
        "raw_polygons": len(polygons),
        "kept_polygons": len(deduped),
        "removed_polygons": dropped,
        "raw_vertices": raw_vertices,
        "kept_vertices": kept_vertices,
        "removed_vertices": removed_vertices,
    }


def filter_water_polygons_to_bbox(
    polygons: list[WaterPolygon],
    *,
    bbox: tuple[float, float, float, float],
) -> list[WaterPolygon]:
    return [polygon for polygon in polygons if _polygon_overlaps_bbox(polygon, bbox)]


def load_water_polygons_from_cache(path: Path) -> list[WaterPolygon]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError("Cache payload is not a JSON object")
    if int(payload.get("cache_format_version", 0)) != 2:
        raise RuntimeError("Unsupported cache format version")
    footprints = payload.get("footprints", [])
    if not isinstance(footprints, list):
        raise RuntimeError("Cache payload missing footprints")

    polygons: list[WaterPolygon] = []
    for item in footprints:
        if not isinstance(item, dict):
            continue
        outer_rings_raw = item.get("outer_rings_lonlat", [])
        inner_rings_raw = item.get("inner_rings_lonlat", [])
        outer_rings: list[tuple[tuple[float, float], ...]] = []
        inner_rings: list[tuple[tuple[float, float], ...]] = []
        if isinstance(outer_rings_raw, list):
            for ring in outer_rings_raw:
                if not isinstance(ring, list):
                    continue
                points = tuple((float(lon), float(lat)) for lon, lat in ring)
                if len(points) >= 4:
                    outer_rings.append(points)
        if isinstance(inner_rings_raw, list):
            for ring in inner_rings_raw:
                if not isinstance(ring, list):
                    continue
                points = tuple((float(lon), float(lat)) for lon, lat in ring)
                if len(points) >= 4:
                    inner_rings.append(points)
        if not outer_rings:
            continue
        polygons.append(
            WaterPolygon(
                osm_id=str(item.get("water_id", "water")),
                kind=str(item.get("kind", "water_polygon")),
                outer_rings=tuple(outer_rings),
                inner_rings=tuple(inner_rings),
                source=str(item.get("source", "")),
                tags={
                    str(key): str(value)
                    for key, value in dict(item.get("tags", {})).items()
                    if isinstance(key, str) and isinstance(value, str)
                },
            )
        )
    return polygons


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


def build_svg_preview(
    polygons: list[WaterPolygon],
    *,
    bbox: tuple[float, float, float, float],
    center_lat_deg: float,
    center_lon_deg: float,
    radius_km: float,
    simplify: bool,
    simplify_min_distance_km: float,
    simplify_max_distance_gap_km: float,
    simplify_max_offset_m: float,
    width: int,
    height: int,
    padding: float,
) -> str:
    preview_bbox = _bbox_union(
        bbox,
        bbox_from_point(float(center_lat_deg), float(center_lon_deg), float(radius_km)),
    )
    raw_vertex_count = sum(_polygon_vertex_count(polygon) for polygon in polygons)
    simplification_stats: dict[str, int] | None = None
    if simplify:
        polygons, simplification_stats = _dedupe_adjacent_far_polygons(
            list(polygons),
            center_lon_deg=float(center_lon_deg),
            center_lat_deg=float(center_lat_deg),
            min_distance_km=float(simplify_min_distance_km),
            max_neighbor_distance_gap_km=float(simplify_max_distance_gap_km),
            max_neighbor_offset_m=float(simplify_max_offset_m),
        )
    else:
        polygons = list(polygons)
    kept_vertex_count = sum(_polygon_vertex_count(polygon) for polygon in polygons)
    if simplification_stats is None:
        simplification_stats = {
            "raw_polygons": len(polygons),
            "kept_polygons": len(polygons),
            "removed_polygons": 0,
            "raw_vertices": raw_vertex_count,
            "kept_vertices": kept_vertex_count,
            "removed_vertices": 0,
        }
    print(
        "svg_preview simplify={simplify} raw_polygons={raw_polygons} kept_polygons={kept_polygons} "
        "raw_vertices={raw_vertices} kept_vertices={kept_vertices} removed_vertices={removed_vertices}".format(
            simplify=str(bool(simplify)).lower(),
            raw_polygons=int(simplification_stats["raw_polygons"]),
            kept_polygons=int(simplification_stats["kept_polygons"]),
            raw_vertices=int(simplification_stats["raw_vertices"]),
            kept_vertices=int(simplification_stats["kept_vertices"]),
            removed_vertices=int(simplification_stats["removed_vertices"]),
        ),
        file=sys.stderr,
    )
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">'
        ),
        '<rect x="0" y="0" width="100%" height="100%" fill="#f7fbff"/>',
    ]
    ring_values = _build_distance_ring_values(radius_km)
    if ring_values:
        center_x, center_y = project_lonlat(
            float(center_lon_deg),
            float(center_lat_deg),
            bbox=preview_bbox,
            width=width,
            height=height,
            padding=padding,
        )
        for ring_km in ring_values:
            ring_points = ring_to_svg_points(
                _circle_ring_lonlat(
                    center_lon_deg=float(center_lon_deg),
                    center_lat_deg=float(center_lat_deg),
                    radius_km=float(ring_km),
                ),
                bbox=preview_bbox,
                width=width,
                height=height,
                padding=padding,
            )
            parts.append(
                (
                    '<path d="{d}" fill="none" stroke="#6f8597" stroke-opacity="0.42" '
                    'stroke-width="1.0" stroke-dasharray="5 4" vector-effect="non-scaling-stroke"/>'
                ).format(d=path_for_points(ring_points))
            )
            label_x = min(float(width) - 24.0, center_x + 8.0)
            label_y = max(16.0, center_y - 6.0)
            parts.append(
                (
                    '<text x="{x:.1f}" y="{y:.1f}" font-family="sans-serif" font-size="12" '
                    'fill="#5d7284" fill-opacity="0.92">{label:.1f} km</text>'
                ).format(x=label_x, y=label_y, label=float(ring_km))
            )
        parts.append(
            (
                '<circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="#111827" fill-opacity="0.9" '
                'stroke="#f8fafc" stroke-width="1.4" vector-effect="non-scaling-stroke"/>'
            ).format(x=center_x, y=center_y)
        )
        parts.append(
            (
                '<path d="M {x1:.1f} {y1:.1f} L {x2:.1f} {y2:.1f} M {x3:.1f} {y3:.1f} L {x4:.1f} {y4:.1f}" '
                'fill="none" stroke="#f8fafc" stroke-width="1.2" vector-effect="non-scaling-stroke"/>'
            ).format(
                x1=center_x - 9.0,
                y1=center_y,
                x2=center_x + 9.0,
                y2=center_y,
                x3=center_x,
                y3=center_y - 9.0,
                x4=center_x,
                y4=center_y + 9.0,
            )
        )
    for polygon in polygons:
        if not polygon.outer_rings:
            continue
        fill, stroke, fill_opacity = _water_style(polygon)
        path_d = " ".join(
            svg_path_for_ring(ring, bbox=preview_bbox, width=width, height=height, padding=padding)
            for ring in polygon.outer_rings
        )
        path_d += "".join(
            f" {svg_path_for_ring(ring, bbox=preview_bbox, width=width, height=height, padding=padding)}"
            for ring in polygon.inner_rings
        )
        parts.append(
            (
                '<path d="{d}" fill="{fill}" fill-opacity="{fill_opacity:.2f}" '
                'stroke="{stroke}" stroke-width="1.2" vector-effect="non-scaling-stroke" '
                'fill-rule="evenodd"/>'
            ).format(d=path_d, fill=fill, stroke=stroke, fill_opacity=fill_opacity)
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

    view_bbox = bbox_from_point(args.lat, args.lon, float(args.radius_km))
    query_bbox = bbox_from_point(
        args.lat,
        args.lon,
        float(args.radius_km) + QUERY_MARGIN_KM,
    )
    polygons: list[WaterPolygon]
    if args.input_cache is not None:
        polygons = load_water_polygons_from_cache(args.input_cache)
        polygons = filter_water_polygons_to_bbox(polygons, bbox=view_bbox)
        query = build_overpass_query(query_bbox)
    else:
        query = build_overpass_query(query_bbox)
        polygons = []
        best_coastline_count = -1
        best_polygon_count = -1
        for _ in range(OVERPASS_SAMPLE_ATTEMPTS):
            payload = fetch_overpass_json(
                endpoint=args.endpoint,
                query=query,
                user_agent=args.user_agent,
                timeout_s=float(args.timeout_s),
            )
            elements = payload.get("elements")
            if not isinstance(elements, list):
                raise RuntimeError("Overpass payload is missing the elements array")

            candidate_polygons = extract_water_features(elements, bbox=query_bbox)
            candidate_polygons = filter_water_polygons_to_bbox(candidate_polygons, bbox=view_bbox)
            coastline_count = sum(1 for polygon in candidate_polygons if polygon.kind == "coastline")
            polygon_count = len(candidate_polygons)
            if (
                coastline_count > best_coastline_count
                or (
                    coastline_count == best_coastline_count
                    and polygon_count > best_polygon_count
                )
            ):
                polygons = candidate_polygons
                best_coastline_count = coastline_count
                best_polygon_count = polygon_count
            if coastline_count > 0:
                break
    normalized = serialize_water_layer(
        polygons,
        bbox=view_bbox,
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
            bbox=geometry_bounds(polygons, view_bbox),
            center_lat_deg=float(args.lat),
            center_lon_deg=float(args.lon),
            radius_km=float(args.radius_km),
            simplify=True,
            simplify_min_distance_km=float(args.simplify_min_distance_km),
            simplify_max_distance_gap_km=float(args.simplify_max_distance_gap_km),
            simplify_max_offset_m=float(args.simplify_max_offset_m),
            width=int(args.svg_width),
            height=int(args.svg_height),
            padding=float(args.svg_padding),
        )
        args.output_svg.write_text(svg_text + "\n", encoding="utf-8")
        print(f"wrote_svg={args.output_svg}")

    if args.output_svg_raw is not None:
        args.output_svg_raw.parent.mkdir(parents=True, exist_ok=True)
        svg_text = build_svg_preview(
            polygons,
            bbox=geometry_bounds(polygons, view_bbox),
            center_lat_deg=float(args.lat),
            center_lon_deg=float(args.lon),
            radius_km=float(args.radius_km),
            simplify=False,
            simplify_min_distance_km=float(args.simplify_min_distance_km),
            simplify_max_distance_gap_km=float(args.simplify_max_distance_gap_km),
            simplify_max_offset_m=float(args.simplify_max_offset_m),
            width=int(args.svg_width),
            height=int(args.svg_height),
            padding=float(args.svg_padding),
        )
        args.output_svg_raw.write_text(svg_text + "\n", encoding="utf-8")
        print(f"wrote_svg_raw={args.output_svg_raw}")

    print(
        f"summary polygons={len(polygons)} bbox={view_bbox[0]:.6f},{view_bbox[1]:.6f},{view_bbox[2]:.6f},{view_bbox[3]:.6f}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
