from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, Sequence


EARTH_RADIUS_KM = 6371.0088

POLYGON_WATER_KEYS = {
    ("natural", "water"),
    ("waterway", "riverbank"),
}


@dataclass(frozen=True, slots=True)
class WaterPolygonFootprint:
    water_id: str
    kind: str
    outer_rings_lonlat: tuple[tuple[tuple[float, float], ...], ...]
    inner_rings_lonlat: tuple[tuple[tuple[float, float], ...], ...]
    source: str
    tags: dict[str, str]


@dataclass(frozen=True, slots=True)
class WaterSurfacePatch:
    patch_id: str
    polygon_id: str
    anchor_points_lonlat: tuple[
        tuple[float, float],
        tuple[float, float],
        tuple[float, float],
    ]
    anchor_elevations_m: tuple[float, float, float]
    flat_threshold_m: float = 1.0

    @property
    def surface_mode(self) -> str:
        return classify_water_surface_mode(
            self.anchor_elevations_m,
            flat_threshold_m=self.flat_threshold_m,
        )

    @property
    def height_span_m(self) -> float:
        values = tuple(float(value) for value in self.anchor_elevations_m)
        return max(values) - min(values)


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
    body = ring[:-1]
    if not body:
        return ring
    start_index = min(range(len(body)), key=lambda index: body[index])
    rotated = body[start_index:] + body[:start_index]
    return rotated + (rotated[0],)


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


def way_lookup(
    elements: list[dict[str, Any]],
    nodes: dict[int, tuple[float, float]],
) -> dict[int, tuple[tuple[float, float], ...]]:
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


def extract_water_polygons(elements: list[dict[str, Any]]) -> tuple[WaterPolygonFootprint, ...]:
    nodes = node_lookup(elements)
    ways = way_lookup(elements, nodes)
    relations = relation_lookup(elements)

    polygons: list[WaterPolygonFootprint] = []

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
                WaterPolygonFootprint(
                    water_id=f"way/{way_id}",
                    kind=polygon_kind_from_tags(tags),
                    outer_rings_lonlat=(points,),
                    inner_rings_lonlat=(),
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
            WaterPolygonFootprint(
                water_id=f"relation/{relation_id}",
                kind=polygon_kind_from_tags(tags),
                outer_rings_lonlat=tuple(outer_rings),
                inner_rings_lonlat=tuple(inner_rings),
                source="relation",
                tags=tags,
            )
        )

    return tuple(polygons)


def classify_water_surface_mode(
    anchor_elevations_m: Sequence[float],
    *,
    flat_threshold_m: float = 1.0,
) -> str:
    if len(anchor_elevations_m) != 3:
        raise ValueError("anchor_elevations_m must contain exactly three values")
    threshold = float(flat_threshold_m)
    if threshold < 0.0:
        raise ValueError("flat_threshold_m must be non-negative")
    values = tuple(float(value) for value in anchor_elevations_m)
    if not all(math.isfinite(value) for value in values):
        return "sloped"
    span = max(values) - min(values)
    return "flat" if span <= threshold else "sloped"
