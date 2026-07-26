from __future__ import annotations

import json
import math
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, replace
from typing import Any, Callable, Iterable, Sequence

import numpy as np
from pyproj import Transformer
from pyproj.enums import TransformDirection

from .clouddisc.types import DownloadCancelledError
from .astro import is_in_fov
from .location_resolver.place_projection import project_place_targets_to_altaz
from .terrain import WGS84_GEOD, build_ray_scan_grid
from .user_agent import build_user_agent
from .water_surface_mesh import make_local_transformer, project_ring_xy

EARTH_RADIUS_KM = 6371.0088
DEFAULT_WATER_OVERPASS_ENDPOINT = "https://overpass-api.de/api/interpreter"
DEFAULT_WATER_USER_AGENT = build_user_agent("water-overlay")
DEFAULT_WATER_TIMEOUT_S = 60.0
DEFAULT_WATER_RADIUS_KM = 2.0
DEFAULT_WATER_HORIZON_MARGIN_KM = 1.0
DEFAULT_WATER_SAMPLE_STEP_M = 1.25**5
DEFAULT_WATER_INLAND_SAMPLE_MIN_DISTANCE_M = 5.0
DEFAULT_WATER_AZIMUTH_STEP_DEG = 2.0
DEFAULT_WATER_SAMPLE_GROWTH_FACTOR = 1.15
DEFAULT_WATER_SIMPLIFICATION_APPARENT_ANGLE_DEG = 0.5
DEFAULT_WATER_SIMPLIFICATION_MIN_GRID_M = 1.0
DEFAULT_WATER_QUERY_BBOX_SCALE = 1.2
DEFAULT_WATER_SCAN_RADIUS_MAX_KM = 128.0
DEFAULT_WATER_BOUNDARY_RADIUS_KM = 10.0

POLYGON_WATER_KEYS = {
    ("natural", "water"),
    ("waterway", "riverbank"),
}
COASTLINE_WATER_TAG = {"natural": "coastline"}
RIVER_LIKE_WATER_TAG_VALUES = {"river", "stream", "canal", "drain"}


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
    scan_distance_m: float | None = None
    scan_azimuth_index: int | None = None
    scan_distance_index: int | None = None
    water_category: str = "lake"


@dataclass(frozen=True, slots=True)
class WaterOverlayPolyline:
    water_id: str
    water_category: str
    points: tuple[WaterOverlayPoint, ...]


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


def expanded_query_bbox_from_point(
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    *,
    scale: float = DEFAULT_WATER_QUERY_BBOX_SCALE,
) -> tuple[float, float, float, float]:
    if scale <= 0.0:
        raise ValueError("scale must be positive")
    return bbox_from_point(lat_deg, lon_deg, float(radius_km) * float(scale))


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
    scan_radius_km = max(
        float(minimum_distance_km),
        horizon_distance_km_from_height(observer_height_m) + float(horizon_margin_km),
    )
    return min(float(DEFAULT_WATER_SCAN_RADIUS_MAX_KM), float(scan_radius_km))


def resolve_water_surface_azimuth_step_deg(
    surface_width_px: int,
    surface_height_px: int,
) -> float:
    """Return the fixed water-surface azimuth step used for all render sizes."""
    return DEFAULT_WATER_AZIMUTH_STEP_DEG


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


def classify_water_surface_category(
    tags: dict[str, str],
    *,
    kind: str | None = None,
) -> str:
    if kind == "coastline" or tags.get("natural") == "coastline":
        return "sea"
    water_tag = tags.get("water")
    waterway_tag = tags.get("waterway")
    if water_tag in RIVER_LIKE_WATER_TAG_VALUES or waterway_tag == "riverbank":
        return "river"
    return "lake"


def relation_is_water_relation(tags: dict[str, str]) -> bool:
    for key, value in POLYGON_WATER_KEYS:
        if tags.get(key) == value:
            return True
    return False


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


def _build_coastline_water_polygons(
    elements: list[dict[str, Any]],
    *,
    bbox: tuple[float, float, float, float],
) -> tuple[WaterPolygonFootprint, ...]:
    return _build_coastline_water_polygons_fallback(elements, bbox=bbox)


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
    abort_event: threading.Event | None = None,
) -> dict[str, Any]:
    _raise_if_abort_requested(abort_event)
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
            _raise_if_abort_requested(abort_event)
            payload = response.read()
            _raise_if_abort_requested(abort_event)
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


def fetch_water_overlay_footprints(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    max_distance_km: float,
    endpoint: str = DEFAULT_WATER_OVERPASS_ENDPOINT,
    user_agent: str = DEFAULT_WATER_USER_AGENT,
    timeout_s: float = DEFAULT_WATER_TIMEOUT_S,
    abort_event: threading.Event | None = None,
) -> tuple[WaterPolygonFootprint, ...]:
    bbox = expanded_query_bbox_from_point(
        float(observer_lat_deg),
        float(observer_lon_deg),
        float(max_distance_km),
    )
    payload = fetch_overpass_json(
        bbox=bbox,
        endpoint=endpoint,
        user_agent=user_agent,
        timeout_s=timeout_s,
        abort_event=abort_event,
    )
    elements = payload.get("elements")
    if not isinstance(elements, list):
        raise RuntimeError("invalid payload")
    footprints = extract_water_polygons(elements, bbox=bbox, abort_event=abort_event)
    return tuple(
        footprint
        for footprint in footprints
        if not _footprint_is_sea_like(footprint)
    )


def sample_water_overlay_points_for_observer(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    max_distance_km: float = DEFAULT_WATER_RADIUS_KM,
    fallback_surface_height_m: float = 0.0,
    target_ground_elevation_m_sampler: Callable[[float, float], float] | None = None,
    endpoint: str = DEFAULT_WATER_OVERPASS_ENDPOINT,
    user_agent: str = DEFAULT_WATER_USER_AGENT,
    timeout_s: float = DEFAULT_WATER_TIMEOUT_S,
    abort_event: threading.Event | None = None,
    front_hemisphere_view_center: tuple[float, float] | None = None,
    front_hemisphere_fov_deg: float = 90.0,
) -> tuple[WaterOverlayPoint, ...]:
    footprints = fetch_water_overlay_footprints(
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        max_distance_km=max_distance_km,
        endpoint=endpoint,
        user_agent=user_agent,
        timeout_s=timeout_s,
        abort_event=abort_event,
    )
    simplified_footprints = simplify_water_footprints_for_observer(
        footprints,
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
    )
    return sample_water_overlay_points(
        simplified_footprints,
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        observer_height_m=observer_height_m,
        fallback_surface_height_m=fallback_surface_height_m,
        target_ground_elevation_m_sampler=target_ground_elevation_m_sampler,
        max_distance_km=max_distance_km,
        abort_event=abort_event,
        front_hemisphere_view_center=front_hemisphere_view_center,
        front_hemisphere_fov_deg=front_hemisphere_fov_deg,
    )


def _clip_ring_xy_to_radius(
    ring_xy: Sequence[tuple[float, float]],
    radius_m: float,
) -> tuple[tuple[tuple[float, float], ...], ...]:
    """Clip a closed local-XY ring to an observer-centered distance circle."""
    body = list(ring_xy[:-1] if ring_xy and ring_xy[0] == ring_xy[-1] else ring_xy)
    if len(body) < 2:
        return ()
    radius_sq = float(radius_m) * float(radius_m)
    fragments: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []

    def inside(point: tuple[float, float]) -> bool:
        return (point[0] * point[0]) + (point[1] * point[1]) <= radius_sq

    def finish() -> None:
        nonlocal current
        if len(current) >= 2:
            deduped: list[tuple[float, float]] = []
            for point in current:
                if not deduped or point != deduped[-1]:
                    deduped.append(point)
            if len(deduped) >= 2:
                fragments.append(deduped)
        current = []

    for index, start in enumerate(body):
        end = body[(index + 1) % len(body)]
        x0, y0 = float(start[0]), float(start[1])
        dx = float(end[0]) - x0
        dy = float(end[1]) - y0
        a = (dx * dx) + (dy * dy)
        cuts = [0.0, 1.0]
        if a > 0.0:
            b = 2.0 * ((x0 * dx) + (y0 * dy))
            c = (x0 * x0) + (y0 * y0) - radius_sq
            discriminant = (b * b) - (4.0 * a * c)
            if discriminant > 0.0:
                root = math.sqrt(discriminant)
                for t in ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a)):
                    if 0.0 < t < 1.0:
                        cuts.append(float(t))
        cuts = sorted(set(cuts))
        for left, right in zip(cuts, cuts[1:]):
            p0 = (x0 + (dx * left), y0 + (dy * left))
            p1 = (x0 + (dx * right), y0 + (dy * right))
            midpoint = (
                0.5 * (p0[0] + p1[0]),
                0.5 * (p0[1] + p1[1]),
            )
            if not inside(midpoint):
                finish()
                continue
            if not current:
                current.append(p0)
            elif current[-1] != p0:
                current.append(p0)
            current.append(p1)
    finish()
    if len(fragments) >= 2 and fragments[0][0] == fragments[-1][-1]:
        fragments[0] = fragments[-1] + fragments[0][1:]
        fragments.pop()
    return tuple(tuple(fragment) for fragment in fragments)


def build_water_overlay_polylines(
    footprints: Sequence[WaterPolygonFootprint],
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    fallback_surface_height_m: float = 0.0,
    target_ground_elevation_m_sampler: Callable[[float, float], float] | None = None,
    max_distance_km: float = DEFAULT_WATER_RADIUS_KM,
) -> tuple[WaterOverlayPolyline, ...]:
    """Project simplified water-ring boundaries into observer-relative Alt/Az."""
    if max_distance_km <= 0.0:
        raise ValueError("max_distance_km must be positive")

    transformer = make_local_transformer(
        float(observer_lat_deg), float(observer_lon_deg)
    )
    result: list[WaterOverlayPolyline] = []
    max_distance_m = float(max_distance_km) * 1000.0
    for footprint in footprints:
        for ring_index, ring in enumerate(
            (*footprint.outer_rings_lonlat, *footprint.inner_rings_lonlat)
        ):
            if len(ring) < 4:
                continue
            ring_xy = project_ring_xy(ring, transformer)
            for fragment_index, clipped_ring in enumerate(
                _clip_ring_xy_to_radius(ring_xy, max_distance_m)
            ):
                xs = [point[0] for point in clipped_ring]
                ys = [point[1] for point in clipped_ring]
                lon_values, lat_values = transformer.transform(
                    xs,
                    ys,
                    direction=TransformDirection.INVERSE,
                )
                latitudes = [float(value) for value in lat_values]
                longitudes = [float(value) for value in lon_values]
                target_heights = [
                    _water_surface_height_m(
                        footprint,
                        fallback_surface_height_m=fallback_surface_height_m,
                        target_ground_elevation_m_sampler=target_ground_elevation_m_sampler,
                        latitude_deg=latitude_deg,
                        longitude_deg=longitude_deg,
                    )
                    for latitude_deg, longitude_deg in zip(latitudes, longitudes)
                ]
                projections = project_place_targets_to_altaz(
                    observer_latitude_deg=float(observer_lat_deg),
                    observer_longitude_deg=float(observer_lon_deg),
                    observer_height_m=float(observer_height_m),
                    target_latitude_deg=latitudes,
                    target_longitude_deg=longitudes,
                    target_height_m=target_heights,
                )
                points = tuple(
                    WaterOverlayPoint(
                        water_id=f"{footprint.water_id}/ring-{ring_index}-{fragment_index}",
                        alt_deg=float(projection.alt_deg),
                        az_deg=float(projection.az_deg),
                        distance_km=float(projection.distance_km),
                        water_category=classify_water_surface_category(
                            footprint.tags, kind=footprint.kind
                        ),
                    )
                    for projection in projections
                )
                if len(points) >= 2:
                    result.append(
                        WaterOverlayPolyline(
                            water_id=f"{footprint.water_id}/ring-{ring_index}-{fragment_index}",
                            water_category=classify_water_surface_category(
                                footprint.tags, kind=footprint.kind
                            ),
                            points=points,
                        )
                    )
    return tuple(result)


def extract_water_polygons(
    elements: list[dict[str, Any]],
    *,
    bbox: tuple[float, float, float, float] | None = None,
    abort_event: threading.Event | None = None,
) -> tuple[WaterPolygonFootprint, ...]:
    _raise_if_abort_requested(abort_event)
    nodes = node_lookup(elements)
    ways = way_lookup(elements, nodes)
    relations = relation_lookup(elements)

    polygons: list[WaterPolygonFootprint] = []

    way_tags: dict[int, dict[str, str]] = {}
    for element_index, element in enumerate(elements):
        _cooperative_yield(abort_event, interval=256, iteration_index=element_index)
        if element.get("type") != "way":
            continue
        try:
            way_id = int(element["id"])
        except (KeyError, TypeError, ValueError):
            continue
        way_tags[way_id] = tags_to_dict(element.get("tags"))

    for way_index, (way_id, points) in enumerate(ways.items()):
        _cooperative_yield(abort_event, interval=256, iteration_index=way_index)
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

    for relation_index, (relation_id, relation) in enumerate(relations.items()):
        _cooperative_yield(abort_event, interval=128, iteration_index=relation_index)
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


def _water_surface_height_m(
    footprint: WaterPolygonFootprint,
    *,
    fallback_surface_height_m: float,
    target_ground_elevation_m_sampler: Callable[[float, float], float] | None = None,
    latitude_deg: float,
    longitude_deg: float,
) -> float:
    explicit_target_height_m = _footprint_explicit_surface_height_m(footprint)
    if explicit_target_height_m is not None:
        return explicit_target_height_m
    if _footprint_is_sea_like(footprint):
        return 0.0
    if target_ground_elevation_m_sampler is not None:
        try:
            return float(target_ground_elevation_m_sampler(latitude_deg, longitude_deg))
        except Exception:
            return float(fallback_surface_height_m)
    return float(fallback_surface_height_m)


def _footprint_is_sea_like(footprint: WaterPolygonFootprint) -> bool:
    tags = footprint.tags
    return classify_water_surface_category(tags, kind=footprint.kind) == "sea"


def _raise_if_abort_requested(abort_event: threading.Event | None) -> None:
    if abort_event is not None and abort_event.is_set():
        raise DownloadCancelledError("water overlay update cancelled")


def _cooperative_yield(
    abort_event: threading.Event | None,
    *,
    interval: int,
    iteration_index: int,
) -> None:
    if interval <= 0 or iteration_index % interval != 0:
        return
    _raise_if_abort_requested(abort_event)
    time.sleep(0)


def _ring_body_xy(points_xy: Sequence[tuple[float, float]]) -> list[tuple[float, float]]:
    ring = list(points_xy)
    if len(ring) >= 2 and ring[0] == ring[-1]:
        ring = ring[:-1]
    if len(ring) < 3:
        return []
    return ring


def _ring_area_xy(ring_xy: Sequence[tuple[float, float]]) -> float:
    points = list(ring_xy)
    if len(points) < 3:
        return 0.0
    if points[0] != points[-1]:
        points.append(points[0])
    if len(points) < 4:
        return 0.0
    area2 = 0.0
    for (x0, y0), (x1, y1) in zip(points, points[1:], strict=False):
        area2 += float(x0) * float(y1) - float(x1) * float(y0)
    return abs(area2) * 0.5


def _snap_value_to_grid_m(value_m: float, grid_size_m: float) -> float:
    if grid_size_m <= 0.0:
        raise ValueError("grid_size_m must be positive")
    cell = float(grid_size_m)
    return math.floor(float(value_m) / cell + 0.5) * cell


def _grid_size_for_distance_m(distance_m: float) -> int:
    distance_m = max(0.0, float(distance_m))
    apparent_cell_m = max(
        DEFAULT_WATER_SIMPLIFICATION_MIN_GRID_M,
        distance_m * math.tan(math.radians(DEFAULT_WATER_SIMPLIFICATION_APPARENT_ANGLE_DEG)),
    )
    exponent = max(0, int(math.floor(math.log2(apparent_cell_m))))
    return max(int(DEFAULT_WATER_SIMPLIFICATION_MIN_GRID_M), 1 << exponent)


def _simplify_closed_ring_xy_by_distance_grid(
    ring_xy: Sequence[tuple[float, float]],
) -> tuple[list[tuple[float, float]] | None, bool]:
    body = _ring_body_xy(ring_xy)
    if len(body) < 3:
        return list(ring_xy), False

    snapped_body: list[tuple[float, float]] = []
    for x, y in body:
        distance_m = math.hypot(float(x), float(y))
        grid_size_m = float(_grid_size_for_distance_m(distance_m))
        snapped = (
            _snap_value_to_grid_m(float(x), grid_size_m),
            _snap_value_to_grid_m(float(y), grid_size_m),
        )
        if snapped_body and snapped_body[-1] == snapped:
            continue
        snapped_body.append(snapped)

    if len(snapped_body) < 3:
        return None, True

    closed_snapped = snapped_body + [snapped_body[0]]
    if _ring_area_xy(closed_snapped) <= 0.0:
        return None, True
    return closed_snapped, True


def _simplify_ring_lonlat_by_distance_grid(
    ring_lonlat: Sequence[tuple[float, float]],
    *,
    transformer: Transformer,
) -> tuple[tuple[tuple[float, float], ...] | None, bool]:
    ring_xy = project_ring_xy(tuple(ring_lonlat), transformer)
    simplified_xy, changed = _simplify_closed_ring_xy_by_distance_grid(ring_xy)
    if simplified_xy is None:
        return None, True
    if not changed:
        return tuple(tuple(point) for point in ring_lonlat), False
    xs = [point[0] for point in simplified_xy]
    ys = [point[1] for point in simplified_xy]
    lon_values, lat_values = transformer.transform(
        xs,
        ys,
        direction=TransformDirection.INVERSE,
    )
    simplified_ring = tuple((float(lon), float(lat)) for lon, lat in zip(lon_values, lat_values))
    if len(simplified_ring) < 4 or simplified_ring[0] != simplified_ring[-1]:
        simplified_ring = simplified_ring + (simplified_ring[0],)
    return simplified_ring, simplified_ring != tuple(tuple(point) for point in ring_lonlat)


def _footprint_distance_km(
    footprint: WaterPolygonFootprint,
    *,
    transformer: Transformer,
) -> float:
    xs: list[float] = []
    ys: list[float] = []
    for ring in footprint.outer_rings_lonlat:
        ring_xy = project_ring_xy(ring, transformer)  # type: ignore[arg-type]
        for x, y in ring_xy:
            xs.append(float(x))
            ys.append(float(y))
    if not xs or not ys:
        return float("inf")
    anchor_x = sum(xs) / len(xs)
    anchor_y = sum(ys) / len(ys)
    return math.hypot(anchor_x, anchor_y) / 1000.0


def simplify_water_footprints_for_observer(
    footprints: Sequence[WaterPolygonFootprint],
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> tuple[WaterPolygonFootprint, ...]:
    transformer = make_local_transformer(float(observer_lat_deg), float(observer_lon_deg))
    simplified: list[WaterPolygonFootprint] = []
    for footprint in footprints:
        footprint_distance_km = _footprint_distance_km(footprint, transformer=transformer)
        if footprint_distance_km <= 0.0:
            simplified.append(footprint)
            continue

        outer_rings: list[tuple[tuple[float, float], ...]] = []
        inner_rings: list[tuple[tuple[float, float], ...]] = []
        changed = False

        for ring in footprint.outer_rings_lonlat:
            simplified_ring, ring_changed = _simplify_ring_lonlat_by_distance_grid(
                ring,
                transformer=transformer,
            )
            if simplified_ring is None:
                changed = True
                continue
            outer_rings.append(simplified_ring)
            changed = changed or ring_changed

        for ring in footprint.inner_rings_lonlat:
            simplified_ring, ring_changed = _simplify_ring_lonlat_by_distance_grid(
                ring,
                transformer=transformer,
            )
            if simplified_ring is None:
                changed = True
                continue
            inner_rings.append(simplified_ring)
            changed = changed or ring_changed

        if not outer_rings:
            continue
        if not changed:
            simplified.append(footprint)
            continue
        simplified.append(
            replace(
                footprint,
                outer_rings_lonlat=tuple(outer_rings),
                inner_rings_lonlat=tuple(inner_rings),
            )
        )
    return tuple(simplified)


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
    abort_event: threading.Event | None = None,
    front_hemisphere_view_center: tuple[float, float] | None = None,
    front_hemisphere_fov_deg: float = 90.0,
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
    distances_m = distances_m[
        distances_m >= float(DEFAULT_WATER_INLAND_SAMPLE_MIN_DISTANCE_M)
    ]
    _raise_if_abort_requested(abort_event)
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
    front_hemisphere_view_center = (
        tuple(float(value) for value in front_hemisphere_view_center)
        if front_hemisphere_view_center is not None
        else None
    )
    front_hemisphere_fov_deg = float(front_hemisphere_fov_deg)
    front_hemisphere_enabled = front_hemisphere_view_center is not None

    for row_index, azimuth_deg in enumerate(ray_scan.azimuths_deg):
        _cooperative_yield(abort_event, interval=8, iteration_index=row_index)
        if front_hemisphere_enabled and not is_in_fov(
            0.0,
            float(azimuth_deg),
            front_hemisphere_view_center,
            fov_deg=front_hemisphere_fov_deg,
        ):
            continue
        matched_records: list[tuple[int, float, float, float, float, WaterPolygonFootprint]] = []
        for col_index, distance_m in enumerate(ray_scan.distance_grid_m[row_index]):
            if col_index % 128 == 0:
                _raise_if_abort_requested(abort_event)
                time.sleep(0)
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
            target_height_m = _water_surface_height_m(
                matched_footprint,
                fallback_surface_height_m=fallback_surface_height_m,
                target_ground_elevation_m_sampler=target_ground_elevation_m_sampler,
                latitude_deg=lat,
                longitude_deg=lon,
            )
            matched_records.append(
                (
                    int(col_index),
                    float(distance_m),
                    float(lat),
                    float(lon),
                    float(target_height_m),
                    matched_footprint,
                )
            )
        if not matched_records:
            continue
        projections = project_place_targets_to_altaz(
            observer_latitude_deg=observer_lat,
            observer_longitude_deg=observer_lon,
            observer_height_m=observer_height,
            target_latitude_deg=[record[2] for record in matched_records],
            target_longitude_deg=[record[3] for record in matched_records],
            target_height_m=[record[4] for record in matched_records],
        )
        for (col_index, distance_m, lat, lon, _target_height_m, matched_footprint), projection in zip(
            matched_records,
            projections,
            strict=False,
        ):
            if front_hemisphere_enabled and not is_in_fov(
                float(projection.alt_deg),
                float(projection.az_deg),
                front_hemisphere_view_center,
                fov_deg=front_hemisphere_fov_deg,
            ):
                continue
            points.append(
                WaterOverlayPoint(
                    water_id="water",
                    alt_deg=float(projection.alt_deg),
                    az_deg=float(projection.az_deg),
                    distance_km=float(projection.distance_km),
                    alpha_scale=1.0,
                    scan_distance_m=float(distance_m),
                    scan_azimuth_index=int(row_index),
                    scan_distance_index=int(col_index),
                    water_category=classify_water_surface_category(
                        matched_footprint.tags,
                        kind=matched_footprint.kind,
                    ),
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
