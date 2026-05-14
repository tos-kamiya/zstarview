from __future__ import annotations

import json
import math
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Any, Callable, Iterable, Sequence

import numpy as np

from .location_resolver.place_projection import project_place_target_to_altaz
from .terrain import WGS84_GEOD, build_ray_scan_grid


EARTH_RADIUS_KM = 6371.0088
DEFAULT_WATER_OVERPASS_ENDPOINT = "https://overpass-api.de/api/interpreter"
DEFAULT_WATER_USER_AGENT = "zstarview-water-overlay/0.1"
DEFAULT_WATER_TIMEOUT_S = 60.0
DEFAULT_WATER_RADIUS_KM = 2.0
DEFAULT_WATER_HORIZON_MARGIN_KM = 1.0
DEFAULT_WATER_SAMPLE_STEP_M = 1.25**5
DEFAULT_WATER_AZIMUTH_STEP_DEG = 2.0
DEFAULT_WATER_SAMPLE_GROWTH_FACTOR = 1.15
DEFAULT_WATER_ALPHA_MIN = 0.04

POLYGON_WATER_KEYS = {
    ("natural", "water"),
    ("waterway", "riverbank"),
}
COASTLINE_WATER_TAG = {"natural": "coastline"}


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


@dataclass(frozen=True, slots=True)
class WaterOverlayPoint:
    water_id: str
    alt_deg: float
    az_deg: float
    distance_km: float
    alpha_scale: float = 1.0
    scan_azimuth_index: int | None = None
    scan_distance_index: int | None = None


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


def horizon_distance_km_from_height(height_m: float) -> float:
    height_m = max(0.0, float(height_m))
    if height_m <= 0.0:
        return 0.0
    earth_radius_m = EARTH_RADIUS_KM * 1000.0
    return math.sqrt(((earth_radius_m + height_m) ** 2) - (earth_radius_m**2)) / 1000.0


def resolve_water_scan_radius_km(
    observer_height_m: float,
    *,
    minimum_distance_km: float = DEFAULT_WATER_RADIUS_KM,
    horizon_margin_km: float = DEFAULT_WATER_HORIZON_MARGIN_KM,
) -> float:
    if minimum_distance_km <= 0.0:
        raise ValueError("minimum_distance_km must be positive")
    if horizon_margin_km < 0.0:
        raise ValueError("horizon_margin_km must be non-negative")
    return max(
        float(minimum_distance_km),
        horizon_distance_km_from_height(observer_height_m) + float(horizon_margin_km),
    )


def build_geometric_distance_samples(
    max_distance_km: float,
    sample_start_m: float,
    *,
    growth_factor: float = DEFAULT_WATER_SAMPLE_GROWTH_FACTOR,
) -> np.ndarray:
    max_distance_m = max_distance_km * 1000.0
    if max_distance_m <= 0.0:
        raise ValueError("max_distance_km must be positive")
    if sample_start_m <= 0.0:
        raise ValueError("sample_start_m must be positive")
    if growth_factor <= 1.0:
        raise ValueError("growth_factor must be greater than 1.0")

    samples: list[float] = []
    distance_m = float(sample_start_m)
    while distance_m <= max_distance_m:
        samples.append(distance_m)
        distance_m *= float(growth_factor)
    return np.asarray(samples, dtype=np.float64)


def water_overlay_alpha_scale(distance_m: float, max_distance_m: float) -> float:
    if max_distance_m <= 0.0:
        raise ValueError("max_distance_m must be positive")
    distance_ratio = max(0.0, min(1.0, float(distance_m) / float(max_distance_m)))
    alpha = DEFAULT_WATER_ALPHA_MIN + (1.0 - DEFAULT_WATER_ALPHA_MIN) * math.exp(
        -3.6 * distance_ratio
    )
    return max(DEFAULT_WATER_ALPHA_MIN, min(1.0, alpha))


def build_overpass_query(bbox: tuple[float, float, float, float]) -> str:
    west, south, east, north = bbox
    return "\n".join(
        [
            "[out:json][timeout:60];",
            "(",
            f'  way["natural"="water"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});',
            f'  relation["natural"="water"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});',
            f'  way["natural"="coastline"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});',
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


def _distance_to_segment(
    point: tuple[float, float],
    segment: tuple[tuple[float, float], tuple[float, float]],
) -> float:
    px, py = point
    (ax, ay), (bx, by) = segment
    dx = bx - ax
    dy = by - ay
    if dx == 0.0 and dy == 0.0:
        return math.hypot(px - ax, py - ay)
    t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
    t = max(0.0, min(1.0, t))
    cx = ax + t * dx
    cy = ay + t * dy
    return math.hypot(px - cx, py - cy)


def _point_is_right_of_segment(
    point: tuple[float, float],
    segment: tuple[tuple[float, float], tuple[float, float]],
) -> bool:
    px, py = point
    (ax, ay), (bx, by) = segment
    cross = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
    return cross < 0.0


def _collect_coastline_segments(
    elements: list[dict[str, Any]],
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    nodes = node_lookup(elements)
    ways = way_lookup(elements, nodes)
    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for element in elements:
        if element.get("type") != "way":
            continue
        tags = tags_to_dict(element.get("tags"))
        if tags.get("natural") != "coastline":
            continue
        try:
            way_id = int(element["id"])
        except (KeyError, TypeError, ValueError):
            continue
        points = ways.get(way_id)
        if not points or len(points) < 2:
            continue
        segments.extend(tuple(pair) for pair in zip(points, points[1:]))
    return segments


def _polygon_is_water_face(
    polygon: "Any",
    coastline_segments: Sequence[tuple[tuple[float, float], tuple[float, float]]],
) -> bool:
    if not coastline_segments:
        return False
    rep = polygon.representative_point()
    point = (float(rep.x), float(rep.y))
    nearest_segments = sorted(
        coastline_segments,
        key=lambda segment: _distance_to_segment(point, segment),
    )[:3]
    if not nearest_segments:
        return False
    return _point_is_right_of_segment(point, nearest_segments[0])


def _build_coastline_water_polygons(
    elements: list[dict[str, Any]],
    *,
    bbox: tuple[float, float, float, float],
) -> tuple[WaterPolygonFootprint, ...]:
    try:
        from shapely.geometry import LineString, box
        from shapely.ops import polygonize, unary_union
    except Exception:
        return _build_coastline_water_polygons_fallback(elements, bbox=bbox)

    coastline_segments = _collect_coastline_segments(elements)
    if not coastline_segments:
        return ()

    west, south, east, north = bbox
    coastline_lines = [
        LineString([segment[0], segment[1]])
        for segment in coastline_segments
    ]
    coastline_lines.append(LineString(list(box(west, south, east, north).exterior.coords)))
    merged = unary_union(coastline_lines)
    polygons = list(polygonize(merged))

    footprints: list[WaterPolygonFootprint] = []
    for polygon_index, polygon in enumerate(polygons):
        if polygon.is_empty or polygon.area <= 0.0:
            continue
        if not _polygon_is_water_face(polygon, coastline_segments):
            continue
        outer_ring = normalize_ring(
            tuple((float(lon), float(lat)) for lon, lat in polygon.exterior.coords)
        )
        inner_rings = tuple(
            normalize_ring(tuple((float(lon), float(lat)) for lon, lat in interior.coords))
            for interior in polygon.interiors
            if len(interior.coords) >= 4
        )
        footprints.append(
            WaterPolygonFootprint(
                water_id=f"coastline/{polygon_index}",
                kind="coastline",
                outer_rings_lonlat=(outer_ring,),
                inner_rings_lonlat=inner_rings,
                source="coastline",
                tags=dict(COASTLINE_WATER_TAG),
            )
        )
    return tuple(footprints)


def _build_coastline_water_polygons_fallback(
    elements: list[dict[str, Any]],
    *,
    bbox: tuple[float, float, float, float],
) -> tuple[WaterPolygonFootprint, ...]:
    coastline_segments = _collect_coastline_segments(elements)
    if not coastline_segments:
        return ()

    coastline_rings = assemble_rings_from_segments(
        [tuple(segment) for segment in coastline_segments]
    )
    if not coastline_rings:
        return ()

    west, south, east, north = bbox
    outer_ring = normalize_ring(
        (
            (west, south),
            (east, south),
            (east, north),
            (west, north),
            (west, south),
        )
    )
    inner_rings = tuple(
        normalize_ring(ring)
        for ring in coastline_rings
        if is_closed_ring(normalize_ring(ring))
    )
    if not inner_rings:
        return ()

    return (
        WaterPolygonFootprint(
            water_id="coastline/0",
            kind="coastline",
            outer_rings_lonlat=(outer_ring,),
            inner_rings_lonlat=inner_rings,
            source="coastline",
            tags=dict(COASTLINE_WATER_TAG),
        ),
    )


def fetch_overpass_json(
    *,
    bbox: tuple[float, float, float, float],
    endpoint: str = DEFAULT_WATER_OVERPASS_ENDPOINT,
    user_agent: str = DEFAULT_WATER_USER_AGENT,
    timeout_s: float = DEFAULT_WATER_TIMEOUT_S,
) -> dict[str, Any]:
    request_body = urllib.parse.urlencode({"data": build_overpass_query(bbox)}).encode("utf-8")
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
        with urllib.request.urlopen(request, timeout=float(timeout_s)) as response:
            payload = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code}") from exc
    except urllib.error.URLError as exc:
        reason = str(exc.reason).strip().casefold()
        if "timed out" in reason:
            raise RuntimeError("timeout") from exc
        raise RuntimeError("network error") from exc

    try:
        loaded = json.loads(payload.decode("utf-8"))
    except Exception as exc:
        raise RuntimeError("invalid JSON") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError("invalid payload")
    return loaded


def extract_water_polygons(
    elements: list[dict[str, Any]],
    *,
    bbox: tuple[float, float, float, float] | None = None,
) -> tuple[WaterPolygonFootprint, ...]:
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

    if bbox is not None:
        polygons.extend(_build_coastline_water_polygons(elements, bbox=bbox))

    return tuple(polygons)


def _point_in_ring(lon_deg: float, lat_deg: float, ring: Sequence[tuple[float, float]]) -> bool:
    if len(ring) < 4:
        return False
    inside = False
    x = float(lon_deg)
    y = float(lat_deg)
    for index in range(len(ring) - 1):
        x0, y0 = ring[index]
        x1, y1 = ring[index + 1]
        if ((y0 > y) != (y1 > y)) and (y1 != y0):
            x_cross = ((x1 - x0) * (y - y0) / (y1 - y0)) + x0
            if x < x_cross:
                inside = not inside
    return inside


def _footprint_bounds(footprint: WaterPolygonFootprint) -> tuple[float, float, float, float]:
    min_lon = float("inf")
    min_lat = float("inf")
    max_lon = float("-inf")
    max_lat = float("-inf")
    for ring in footprint.outer_rings_lonlat:
        for lon_deg, lat_deg in ring:
            lon = float(lon_deg)
            lat = float(lat_deg)
            min_lon = min(min_lon, lon)
            min_lat = min(min_lat, lat)
            max_lon = max(max_lon, lon)
            max_lat = max(max_lat, lat)
    if not math.isfinite(min_lon):
        return (0.0, 0.0, 0.0, 0.0)
    return (min_lon, min_lat, max_lon, max_lat)


def _point_in_footprint(lon_deg: float, lat_deg: float, footprint: WaterPolygonFootprint) -> bool:
    for outer_ring in footprint.outer_rings_lonlat:
        if not _point_in_ring(lon_deg, lat_deg, outer_ring):
            continue
        if any(_point_in_ring(lon_deg, lat_deg, inner_ring) for inner_ring in footprint.inner_rings_lonlat):
            return False
        return True
    return False


def _footprint_explicit_surface_height_m(footprint: WaterPolygonFootprint) -> float | None:
    tags = footprint.tags
    raw_ele = tags.get("ele")
    if raw_ele is not None:
        try:
            return float(raw_ele)
        except (TypeError, ValueError):
            pass
    raw_water_level = tags.get("water_level")
    if raw_water_level is not None:
        try:
            return float(raw_water_level)
        except (TypeError, ValueError):
            pass
    return None


def _footprint_is_sea_like(footprint: WaterPolygonFootprint) -> bool:
    if footprint.kind == "coastline":
        return True
    tags = footprint.tags
    if tags.get("natural") == "coastline":
        return True
    water_tag = tags.get("water")
    return water_tag in {"sea", "ocean"}


def sample_water_overlay_points(
    footprints: Sequence[WaterPolygonFootprint],
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    fallback_surface_height_m: float = 0.0,
    target_ground_elevation_m_sampler: Callable[[float, float], float] | None = None,
    max_distance_km: float = DEFAULT_WATER_RADIUS_KM,
    sample_step_m: float = DEFAULT_WATER_SAMPLE_STEP_M,
    azimuth_step_deg: float = DEFAULT_WATER_AZIMUTH_STEP_DEG,
) -> tuple[WaterOverlayPoint, ...]:
    if max_distance_km <= 0.0:
        raise ValueError("max_distance_km must be positive")
    if sample_step_m <= 0.0:
        raise ValueError("sample_step_m must be positive")
    if azimuth_step_deg <= 0.0:
        raise ValueError("azimuth_step_deg must be positive")

    distances_m = build_geometric_distance_samples(
        float(max_distance_km),
        float(sample_step_m),
    )
    observer_lon = float(observer_lon_deg)
    observer_lat = float(observer_lat_deg)
    ray_scan = build_ray_scan_grid(
        geod=WGS84_GEOD,
        observer_latitude_deg=observer_lat,
        observer_longitude_deg=observer_lon,
        azimuth_step_deg=float(azimuth_step_deg),
        distance_samples_m=distances_m,
    )
    points: list[WaterOverlayPoint] = []
    observer_height = float(observer_height_m)
    footprint_bounds = tuple((_footprint_bounds(footprint), footprint) for footprint in footprints)

    for row_index, azimuth_deg in enumerate(ray_scan.azimuths_deg):
        for col_index, distance_m in enumerate(ray_scan.distance_grid_m[row_index]):
            lon = float(ray_scan.ray_lon_deg[row_index, col_index])
            lat = float(ray_scan.ray_lat_deg[row_index, col_index])
            matched_footprint = None
            for bounds, footprint in footprint_bounds:
                if not ((bounds[0] <= lon <= bounds[2]) and (bounds[1] <= lat <= bounds[3])):
                    continue
                if _point_in_footprint(lon, lat, footprint):
                    matched_footprint = footprint
                    break
            if matched_footprint is None:
                continue
            explicit_target_height_m = _footprint_explicit_surface_height_m(matched_footprint)
            if explicit_target_height_m is not None:
                target_height_m = explicit_target_height_m
            elif _footprint_is_sea_like(matched_footprint):
                target_height_m = 0.0
            elif target_ground_elevation_m_sampler is not None:
                try:
                    target_height_m = float(target_ground_elevation_m_sampler(lat, lon))
                except Exception:
                    target_height_m = float(fallback_surface_height_m)
            else:
                target_height_m = float(fallback_surface_height_m)
            projection = project_place_target_to_altaz(
                observer_latitude_deg=observer_lat,
                observer_longitude_deg=observer_lon,
                observer_height_m=observer_height,
                target_latitude_deg=lat,
                target_longitude_deg=lon,
                target_height_m=target_height_m,
            )
            points.append(
                WaterOverlayPoint(
                    water_id="water",
                    alt_deg=float(projection.alt_deg),
                    az_deg=float(projection.az_deg),
                    distance_km=float(projection.distance_km),
                    alpha_scale=water_overlay_alpha_scale(float(distance_m), float(max_distance_km) * 1000.0),
                    scan_azimuth_index=int(row_index),
                    scan_distance_index=int(col_index),
                )
            )
    return tuple(points)


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
