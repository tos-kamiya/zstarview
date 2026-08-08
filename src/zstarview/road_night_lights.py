"""OSM road geometry used by the Road Night Lights overlay."""

from __future__ import annotations

import json
import math
import threading
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import pairwise
from pathlib import Path

from pyproj import Transformer

from .location_resolver.place_projection import project_place_targets_to_altaz
from .paths import CACHE_PATH
from .user_agent import build_user_agent
from .water_surface_mesh import make_local_transformer

ROAD_NIGHT_LIGHT_HIGHWAY_TYPES = (
    "motorway",
    "trunk",
    "primary",
    "secondary",
    "tertiary",
)
ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM = 10.0
ROAD_NIGHT_LIGHT_MIN_DISTANCE_KM = 0.5
ROAD_NIGHT_LIGHT_ENDPOINT = "https://overpass-api.de/api/interpreter"
ROAD_NIGHT_LIGHT_USER_AGENT = build_user_agent("road-night-lights")
ROAD_NIGHT_LIGHT_TIMEOUT_S = 60.0
ROAD_NIGHT_LIGHT_CACHE_ROOT = Path(CACHE_PATH) / "road_night_lights"
ROAD_NIGHT_LIGHT_CACHE_FORMAT_VERSION = 1
ROAD_NIGHT_LIGHT_SIMPLIFICATION_APPARENT_ANGLE_DEG = 0.5
ROAD_NIGHT_LIGHT_SIMPLIFICATION_MIN_GRID_M = 1.0


@dataclass(frozen=True, slots=True)
class RoadNightLightWay:
    """An OSM way with a supported highway classification."""

    way_id: int
    highway: str
    coordinates_lonlat: tuple[tuple[float, float], ...]


@dataclass(frozen=True, slots=True)
class RoadNightLightCacheSnapshot:
    ways: tuple[RoadNightLightWay, ...]
    fetched_at_utc: datetime | None = None


@dataclass(frozen=True, slots=True)
class RoadNightLightPoint:
    way_id: int
    highway: str
    alt_deg: float
    az_deg: float
    distance_km: float


@dataclass(frozen=True, slots=True)
class RoadNightLightPolyline:
    way_id: int
    highway: str
    points: tuple[RoadNightLightPoint, ...]


def road_night_lights_scope_key(
    *, observer_lat_deg: float, observer_lon_deg: float, radius_km: float
) -> str:
    """Return a stable cache key with a modestly rounded observer location."""
    return (
        f"earth_{float(observer_lat_deg):+.4f}_{float(observer_lon_deg):+.4f}"
        f"_r{float(radius_km):.2f}"
    )


def road_night_lights_cache_path(
    scope_key: str, *, cache_root: str | Path = ROAD_NIGHT_LIGHT_CACHE_ROOT
) -> Path:
    return Path(cache_root) / f"{scope_key}.json"


def build_road_night_lights_query(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float = ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
) -> str:
    """Build one Overpass query for all supported road classes."""
    if radius_km <= 0.0:
        raise ValueError("radius_km must be positive")
    highway_pattern = "|".join(ROAD_NIGHT_LIGHT_HIGHWAY_TYPES)
    return "\n".join(
        (
            "[out:json][timeout:60];",
            "way(around:{radius:.1f},{lat:.8f},{lon:.8f})".format(
                radius=float(radius_km) * 1000.0,
                lat=float(observer_lat_deg),
                lon=float(observer_lon_deg),
            ),
            f'["highway"~"^({highway_pattern})$"];',
            "out tags geom;",
        )
    )


def parse_road_night_lights_payload(payload: object) -> tuple[RoadNightLightWay, ...]:
    """Parse an Overpass ``out tags geom`` response without requiring a network."""
    if not isinstance(payload, dict) or not isinstance(payload.get("elements"), list):
        raise TypeError("invalid Overpass road payload")
    result: list[RoadNightLightWay] = []
    for element in payload["elements"]:
        if not isinstance(element, dict) or element.get("type") != "way":
            continue
        try:
            way_id = int(element["id"])
        except (KeyError, TypeError, ValueError):
            continue
        tags = element.get("tags")
        highway = tags.get("highway") if isinstance(tags, dict) else None
        if highway not in ROAD_NIGHT_LIGHT_HIGHWAY_TYPES:
            continue
        geometry = element.get("geometry")
        if not isinstance(geometry, list):
            continue
        coordinates: list[tuple[float, float]] = []
        for point in geometry:
            if not isinstance(point, dict):
                continue
            try:
                coordinates.append((float(point["lon"]), float(point["lat"])))
            except (KeyError, TypeError, ValueError):
                continue
        if len(coordinates) >= 2:
            result.append(RoadNightLightWay(way_id, highway, tuple(coordinates)))
    return tuple(result)


def _parse_optional_utc(value: object) -> datetime | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (TypeError, ValueError):
        return None
    return (
        parsed.astimezone(timezone.utc)
        if parsed.tzinfo is not None
        else parsed.replace(tzinfo=timezone.utc)
    )


def save_road_night_lights_cache(
    scope_key: str,
    snapshot: RoadNightLightCacheSnapshot,
    *,
    cache_root: str | Path = ROAD_NIGHT_LIGHT_CACHE_ROOT,
) -> Path:
    path = road_night_lights_cache_path(scope_key, cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    fetched_at = snapshot.fetched_at_utc
    if fetched_at is not None:
        fetched_at = fetched_at.astimezone(timezone.utc)
    payload = {
        "cache_format_version": ROAD_NIGHT_LIGHT_CACHE_FORMAT_VERSION,
        "scope_key": scope_key,
        "fetched_at_utc": fetched_at.isoformat().replace("+00:00", "Z")
        if fetched_at
        else None,
        "ways": [
            {
                "id": way.way_id,
                "tags": {"highway": way.highway},
                "geometry": [
                    {"lon": lon, "lat": lat} for lon, lat in way.coordinates_lonlat
                ],
            }
            for way in snapshot.ways
        ],
    }
    path.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True),
        encoding="utf-8",
    )
    return path


def load_road_night_lights_cache(
    scope_key: str, *, cache_root: str | Path = ROAD_NIGHT_LIGHT_CACHE_ROOT
) -> RoadNightLightCacheSnapshot | None:
    path = road_night_lights_cache_path(scope_key, cache_root=cache_root)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if (
        not isinstance(payload, dict)
        or payload.get("cache_format_version") != ROAD_NIGHT_LIGHT_CACHE_FORMAT_VERSION
    ):
        return None
    if payload.get("scope_key") != scope_key:
        return None
    try:
        cached_ways = payload.get("ways", [])
        if not isinstance(cached_ways, list):
            return None
        ways = parse_road_night_lights_payload(
            {
                "elements": [
                    dict(item, type="way")
                    for item in cached_ways
                    if isinstance(item, dict)
                ]
            }
        )
    except (TypeError, ValueError):
        return None
    return RoadNightLightCacheSnapshot(
        ways=ways, fetched_at_utc=_parse_optional_utc(payload.get("fetched_at_utc"))
    )


def fetch_road_night_lights(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float = ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
    endpoint: str = ROAD_NIGHT_LIGHT_ENDPOINT,
    timeout_s: float = ROAD_NIGHT_LIGHT_TIMEOUT_S,
    abort_event: threading.Event | None = None,
) -> RoadNightLightCacheSnapshot:
    if abort_event is not None and abort_event.is_set():
        raise RuntimeError("cancelled")
    query = build_road_night_lights_query(
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        radius_km=radius_km,
    )
    request = urllib.request.Request(
        endpoint,
        data=urllib.parse.urlencode({"data": query}).encode("utf-8"),
        headers={
            "Content-Type": "application/x-www-form-urlencoded; charset=utf-8",
            "User-Agent": ROAD_NIGHT_LIGHT_USER_AGENT,
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=float(timeout_s)) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code}") from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise RuntimeError("road data request failed") from exc
    return RoadNightLightCacheSnapshot(
        ways=parse_road_night_lights_payload(payload),
        fetched_at_utc=datetime.now(timezone.utc),
    )


def load_or_fetch_road_night_lights(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float = ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
    endpoint: str = ROAD_NIGHT_LIGHT_ENDPOINT,
    timeout_s: float = ROAD_NIGHT_LIGHT_TIMEOUT_S,
    cache_root: str | Path = ROAD_NIGHT_LIGHT_CACHE_ROOT,
    abort_event: threading.Event | None = None,
) -> RoadNightLightCacheSnapshot:
    """Reuse a scope cache and fetch only when that scope is absent."""
    scope_key = road_night_lights_scope_key(
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        radius_km=radius_km,
    )
    cached = load_road_night_lights_cache(scope_key, cache_root=cache_root)
    if cached is not None:
        return cached
    fresh = fetch_road_night_lights(
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        radius_km=radius_km,
        endpoint=endpoint,
        timeout_s=timeout_s,
        abort_event=abort_event,
    )
    save_road_night_lights_cache(scope_key, fresh, cache_root=cache_root)
    return fresh


def _circle_intersections(
    start: tuple[float, float], end: tuple[float, float], radius_m: float
) -> tuple[float, ...]:
    x0, y0 = start
    dx, dy = end[0] - x0, end[1] - y0
    a = dx * dx + dy * dy
    if a <= 0.0:
        return ()
    b = 2.0 * (x0 * dx + y0 * dy)
    c = x0 * x0 + y0 * y0 - radius_m * radius_m
    discriminant = b * b - 4.0 * a * c
    if discriminant <= 0.0:
        return ()
    root = discriminant**0.5
    return tuple(
        t for t in ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a)) if 0.0 < t < 1.0
    )


def _road_simplification_grid_size_m(distance_m: float) -> float:
    apparent_grid_m = max(
        ROAD_NIGHT_LIGHT_SIMPLIFICATION_MIN_GRID_M,
        float(distance_m)
        * math.tan(math.radians(ROAD_NIGHT_LIGHT_SIMPLIFICATION_APPARENT_ANGLE_DEG)),
    )
    exponent = max(0, int(math.floor(math.log2(apparent_grid_m))))
    return float(max(ROAD_NIGHT_LIGHT_SIMPLIFICATION_MIN_GRID_M, 1 << exponent))


def simplify_road_night_light_way_for_observer(
    way: RoadNightLightWay,
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> RoadNightLightWay:
    """Reduce an open way using a distance-dependent local metric grid."""
    if len(way.coordinates_lonlat) < 3:
        return way
    forward = make_local_transformer(observer_lat_deg, observer_lon_deg)
    xs, ys = forward.transform(
        [point[0] for point in way.coordinates_lonlat],
        [point[1] for point in way.coordinates_lonlat],
    )
    snapped: list[tuple[float, float]] = []
    for x, y in zip(xs, ys):
        x_float, y_float = float(x), float(y)
        grid_size = _road_simplification_grid_size_m(math.hypot(x_float, y_float))
        point = (
            math.floor(x_float / grid_size + 0.5) * grid_size,
            math.floor(y_float / grid_size + 0.5) * grid_size,
        )
        if not snapped or snapped[-1] != point:
            snapped.append(point)
    if len(snapped) == len(way.coordinates_lonlat):
        return way
    inverse = Transformer.from_crs(forward.target_crs, "EPSG:4326", always_xy=True)
    lons, lats = inverse.transform(
        [point[0] for point in snapped], [point[1] for point in snapped]
    )
    return RoadNightLightWay(
        way.way_id,
        way.highway,
        tuple((float(lon), float(lat)) for lon, lat in zip(lons, lats)),
    )


def _clip_segment_to_annulus(
    start: tuple[float, float], end: tuple[float, float], inner_m: float, outer_m: float
) -> tuple[tuple[tuple[float, float], tuple[float, float]], ...]:
    cuts = {0.0, 1.0}
    cuts.update(_circle_intersections(start, end, inner_m))
    cuts.update(_circle_intersections(start, end, outer_m))
    ordered = sorted(cuts)
    fragments = []
    for left, right in pairwise(ordered):
        midpoint_t = (left + right) * 0.5
        midpoint = (
            start[0] + (end[0] - start[0]) * midpoint_t,
            start[1] + (end[1] - start[1]) * midpoint_t,
        )
        distance = (midpoint[0] ** 2 + midpoint[1] ** 2) ** 0.5
        if inner_m <= distance <= outer_m:
            p0 = (
                start[0] + (end[0] - start[0]) * left,
                start[1] + (end[1] - start[1]) * left,
            )
            p1 = (
                start[0] + (end[0] - start[0]) * right,
                start[1] + (end[1] - start[1]) * right,
            )
            fragments.append((p0, p1))
    return tuple(fragments)


def clip_road_night_light_way_to_annulus(
    way: RoadNightLightWay,
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    inner_distance_km: float = ROAD_NIGHT_LIGHT_MIN_DISTANCE_KM,
    outer_distance_km: float = ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
) -> tuple[RoadNightLightWay, ...]:
    if inner_distance_km < 0.0 or outer_distance_km <= inner_distance_km:
        raise ValueError("invalid road night light annulus")
    transformer = make_local_transformer(observer_lat_deg, observer_lon_deg)
    xs, ys = transformer.transform(
        [point[0] for point in way.coordinates_lonlat],
        [point[1] for point in way.coordinates_lonlat],
    )
    local = [(float(x), float(y)) for x, y in zip(xs, ys)]
    inner_m, outer_m = inner_distance_km * 1000.0, outer_distance_km * 1000.0
    result: list[RoadNightLightWay] = []
    current: list[tuple[float, float]] = []
    for start, end in pairwise(local):
        segments = _clip_segment_to_annulus(start, end, inner_m, outer_m)
        if not segments:
            if len(current) >= 2:
                result.append(
                    RoadNightLightWay(way.way_id, way.highway, tuple(current))
                )
            current = []
            continue
        for segment in segments:
            if not current or current[-1] != segment[0]:
                if len(current) >= 2:
                    result.append(
                        RoadNightLightWay(way.way_id, way.highway, tuple(current))
                    )
                current = [segment[0]]
            current.append(segment[1])
    if len(current) >= 2:
        result.append(RoadNightLightWay(way.way_id, way.highway, tuple(current)))
    forward = make_local_transformer(observer_lat_deg, observer_lon_deg)
    inverse = Transformer.from_crs(forward.target_crs, "EPSG:4326", always_xy=True)
    converted: list[RoadNightLightWay] = []
    for fragment in result:
        x_values = [point[0] for point in fragment.coordinates_lonlat]
        y_values = [point[1] for point in fragment.coordinates_lonlat]
        lons, lats = inverse.transform(x_values, y_values)
        converted.append(
            RoadNightLightWay(
                fragment.way_id,
                fragment.highway,
                tuple((float(lon), float(lat)) for lon, lat in zip(lons, lats)),
            )
        )
    return tuple(converted)


def clip_road_night_lights_to_annulus(
    ways: tuple[RoadNightLightWay, ...],
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    inner_distance_km: float = ROAD_NIGHT_LIGHT_MIN_DISTANCE_KM,
    outer_distance_km: float = ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
) -> tuple[RoadNightLightWay, ...]:
    return tuple(
        fragment
        for way in ways
        for fragment in clip_road_night_light_way_to_annulus(
            way,
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
            inner_distance_km=inner_distance_km,
            outer_distance_km=outer_distance_km,
        )
    )


def project_road_night_lights(
    ways: tuple[RoadNightLightWay, ...],
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float = 0.0,
) -> tuple[RoadNightLightPolyline, ...]:
    """Project clipped road centerlines into the viewer's Alt/Az coordinates."""
    result: list[RoadNightLightPolyline] = []
    for way in ways:
        latitudes = [point[1] for point in way.coordinates_lonlat]
        longitudes = [point[0] for point in way.coordinates_lonlat]
        projections = project_place_targets_to_altaz(
            observer_latitude_deg=float(observer_lat_deg),
            observer_longitude_deg=float(observer_lon_deg),
            observer_height_m=float(observer_height_m),
            target_latitude_deg=latitudes,
            target_longitude_deg=longitudes,
            target_height_m=[0.0] * len(latitudes),
        )
        points = tuple(
            RoadNightLightPoint(
                way_id=way.way_id,
                highway=way.highway,
                alt_deg=float(projection.alt_deg),
                az_deg=float(projection.az_deg),
                distance_km=float(projection.distance_km),
            )
            for projection in projections
        )
        if len(points) >= 2:
            result.append(RoadNightLightPolyline(way.way_id, way.highway, points))
    return tuple(result)
