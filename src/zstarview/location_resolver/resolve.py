from __future__ import annotations

import logging
import math
import re
import urllib.error
import urllib.parse
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, List

from ..config import load_last_city, save_last_city
from ..paths import CITY_ADMIN1_CODES_FILE, CITY_COORD_FILE
from ..paths import OVERTURE_DERIVED_ROOT_DIR
from ..utils.resolve_city import (
    CityRec,
    load_admin1_names,
    resolve_city,
    resolve_city_by_geonameid,
    resolve_city_by_name,
)
from .mountains import resolve_mountain_viewpoint
from .nominatim import search_nominatim
from .towers import resolve_tower_viewpoint
from .viewpoints import Viewpoint, prefixed_viewpoint_name, split_prefixed_viewpoint

if TYPE_CHECKING:
    from ..data.urban_outline_common import BuildingFootprint

logger = logging.getLogger(__name__)
DEFAULT_OBSERVER_HEIGHT_M = 1.7
BUILDING_TOP_FETCH_RADIUS_KM = 0.15
BUILDING_TOP_MATCH_RADIUS_M = 5.0


class LocationResolveError(Exception):
    """Abort launch because location resolution failed."""


@dataclass(frozen=True)
class ResolvedLocation:
    display_name: str
    lat: float
    lon: float
    tz: str
    persistence_key: str
    observer_height_m: float
    kind: str
    persistence_value: str | dict[str, Any] | None = None
    location_height_label: str | None = None
    location_height_m: float | None = None
    cc: str = ""


def _point_in_ring(point_lonlat: tuple[float, float], ring_lonlat: tuple[tuple[float, float], ...]) -> bool:
    px, py = point_lonlat
    if len(ring_lonlat) < 3:
        return False
    inside = False
    points = ring_lonlat if ring_lonlat[0] == ring_lonlat[-1] else ring_lonlat + (ring_lonlat[0],)
    for (x0, y0), (x1, y1) in zip(points[:-1], points[1:]):
        intersects = ((y0 > py) != (y1 > py)) and (
            px < (x1 - x0) * (py - y0) / ((y1 - y0) or 1e-12) + x0
        )
        if intersects:
            inside = not inside
    return inside


def _building_contains_lonlat(
    building: "BuildingFootprint",
    *,
    lon_deg: float,
    lat_deg: float,
) -> bool:
    if not building.rings_lonlat:
        return False
    point = (lon_deg, lat_deg)
    outer_ring = building.rings_lonlat[0]
    if not _point_in_ring(point, outer_ring):
        return False
    return not any(_point_in_ring(point, hole_ring) for hole_ring in building.rings_lonlat[1:])


def _lonlat_to_local_xy_m(
    lon_deg: float,
    lat_deg: float,
    *,
    origin_lon_deg: float,
    origin_lat_deg: float,
) -> tuple[float, float]:
    mean_lat_rad = math.radians((lat_deg + origin_lat_deg) * 0.5)
    x_m = (lon_deg - origin_lon_deg) * 111_320.0 * math.cos(mean_lat_rad)
    y_m = (lat_deg - origin_lat_deg) * 110_540.0
    return (x_m, y_m)


def _point_to_segment_distance_m(
    px: float,
    py: float,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
) -> float:
    dx = x1 - x0
    dy = y1 - y0
    denom = dx * dx + dy * dy
    if denom <= 1e-12:
        return math.hypot(px - x0, py - y0)
    t = max(0.0, min(1.0, ((px - x0) * dx + (py - y0) * dy) / denom))
    nearest_x = x0 + t * dx
    nearest_y = y0 + t * dy
    return math.hypot(px - nearest_x, py - nearest_y)


def _ring_min_distance_to_lonlat_m(
    ring_lonlat: tuple[tuple[float, float], ...],
    *,
    lon_deg: float,
    lat_deg: float,
) -> float:
    if len(ring_lonlat) < 2:
        return float("inf")
    px, py = 0.0, 0.0
    points = ring_lonlat if ring_lonlat[0] == ring_lonlat[-1] else ring_lonlat + (ring_lonlat[0],)
    best = float("inf")
    for (lon0, lat0), (lon1, lat1) in zip(points[:-1], points[1:]):
        x0, y0 = _lonlat_to_local_xy_m(lon0, lat0, origin_lon_deg=lon_deg, origin_lat_deg=lat_deg)
        x1, y1 = _lonlat_to_local_xy_m(lon1, lat1, origin_lon_deg=lon_deg, origin_lat_deg=lat_deg)
        best = min(best, _point_to_segment_distance_m(px, py, x0, y0, x1, y1))
    return best


def _building_is_near_lonlat(
    building: "BuildingFootprint",
    *,
    lon_deg: float,
    lat_deg: float,
    max_distance_m: float,
) -> bool:
    if _building_contains_lonlat(building, lon_deg=lon_deg, lat_deg=lat_deg):
        return True
    if not building.rings_lonlat:
        return False
    outer_distance_m = _ring_min_distance_to_lonlat_m(
        building.rings_lonlat[0],
        lon_deg=lon_deg,
        lat_deg=lat_deg,
    )
    return outer_distance_m <= max_distance_m


def _find_building_top_height_m(
    buildings: tuple["BuildingFootprint", ...],
    *,
    lon_deg: float,
    lat_deg: float,
    max_distance_m: float = BUILDING_TOP_MATCH_RADIUS_M,
) -> float | None:
    nearby = tuple(
        building
        for building in buildings
        if _building_is_near_lonlat(
            building,
            lon_deg=lon_deg,
            lat_deg=lat_deg,
            max_distance_m=max_distance_m,
        )
    )
    if not nearby:
        return None
    root_ids = {building.parent_building_id or building.building_id for building in nearby}
    related = tuple(
        building
        for building in buildings
        if building.building_id in root_ids or building.parent_building_id in root_ids
    )
    candidates = related or nearby
    return max(float(building.height_m) for building in candidates)


def _resolve_building_top_height_m(
    *,
    lat_deg: float,
    lon_deg: float,
    derived_root_dir: Path | str = OVERTURE_DERIVED_ROOT_DIR,
    fetch_radius_km: float = BUILDING_TOP_FETCH_RADIUS_KM,
    overturemaps_bin: str = "overturemaps",
) -> float | None:
    derived_root = Path(derived_root_dir)
    all_buildings: list["BuildingFootprint"] = []
    from ..data.derived_tile_cache import parse_derived_tile_buildings, select_derived_tile_envelopes
    from ..data.import_overture_buildings import import_overture_buildings

    for feature_type in ("building", "building_part"):
        try:
            derived_dir = import_overture_buildings(
                lat_deg=lat_deg,
                lon_deg=lon_deg,
                radius_km=fetch_radius_km,
                derived_root_dir=derived_root,
                min_building_height_m=0.0,
                feature_type=feature_type,
                fmt="geojsonseq",
                overturemaps_bin=overturemaps_bin,
                dataset_name=None,
                keep_download=None,
                no_stac=False,
            )
        except Exception:
            logger.warning(
                "Building-top viewpoint fetch failed for feature_type=%s at lat=%.6f lon=%.6f",
                feature_type,
                lat_deg,
                lon_deg,
                exc_info=True,
            )
            continue
        try:
            envelopes = select_derived_tile_envelopes(
                derived_dir,
                observer_lat_deg=lat_deg,
                observer_lon_deg=lon_deg,
                radius_km=fetch_radius_km,
            )
        except ValueError:
            continue
        for envelope in envelopes:
            all_buildings.extend(parse_derived_tile_buildings(envelope.path))
    if not all_buildings:
        return None
    return _find_building_top_height_m(
        tuple(all_buildings),
        lon_deg=lon_deg,
        lat_deg=lat_deg,
    )


def _maybe_apply_building_top_viewpoint(
    location: ResolvedLocation,
    *,
    enabled: bool,
) -> ResolvedLocation:
    if not enabled or location.kind in {"tower", "mountain"}:
        return location
    building_top_height_m = _resolve_building_top_height_m(
        lat_deg=float(location.lat),
        lon_deg=float(location.lon),
    )
    if building_top_height_m is None or building_top_height_m <= 0.0:
        return location
    logger.info(
        "Using building-top viewpoint height %.1f m for %s",
        building_top_height_m,
        location.display_name,
    )
    return ResolvedLocation(
        display_name=location.display_name,
        lat=location.lat,
        lon=location.lon,
        tz=location.tz,
        persistence_key=location.persistence_key,
        observer_height_m=float(building_top_height_m) + DEFAULT_OBSERVER_HEIGHT_M,
        kind=location.kind,
        persistence_value=location.persistence_value,
        location_height_label="Building height",
        location_height_m=float(building_top_height_m),
        cc=location.cc,
    )


def format_splash_location(city: ResolvedLocation) -> str:
    return f"Location: {city.display_name}"


def _great_circle_distance_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    lat1_rad = math.radians(lat1)
    lat2_rad = math.radians(lat2)
    dlat = lat2_rad - lat1_rad
    dlon = math.radians(lon2 - lon1)
    a = (
        math.sin(dlat / 2.0) ** 2
        + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2.0) ** 2
    )
    return 6371.0088 * 2.0 * math.asin(min(1.0, math.sqrt(a)))


def _resolve_nearest_city(lat: float, lon: float, admin1_map: dict[tuple[str, str], str]) -> CityRec | None:
    best_city: CityRec | None = None
    best_distance_km = float("inf")
    with open(CITY_COORD_FILE, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 19:
                continue
            rec = CityRec.from_cols(cols, admin1_name=admin1_map.get((cols[8], cols[10])))
            distance_km = _great_circle_distance_km(lat, lon, rec.lat, rec.lon)
            if distance_km < best_distance_km:
                best_distance_km = distance_km
                best_city = rec
    return best_city


def _viewpoint_to_location(
    viewpoint: Viewpoint | None,
    admin1_map: dict[tuple[str, str], str],
) -> ResolvedLocation | None:
    if viewpoint is None:
        return None
    nearest_city = _resolve_nearest_city(viewpoint.latitude_deg, viewpoint.longitude_deg, admin1_map)
    timezone_name = nearest_city.tz if nearest_city is not None else "UTC"
    viewpoint_height_m = 0.0 if viewpoint.viewpoint_height_m is None else float(viewpoint.viewpoint_height_m)
    location_height_label: str | None = None
    location_height_m: float | None = None
    if viewpoint.kind == "tower" and viewpoint.height_m > 0.0:
        location_height_label = "Tower height"
        location_height_m = float(viewpoint.height_m)
    elif viewpoint.kind == "mountain":
        raw_elevation_m = viewpoint.meta.get("elevation_m")
        if isinstance(raw_elevation_m, (int, float)) and float(raw_elevation_m) > 0.0:
            location_height_label = "Elevation"
            location_height_m = float(raw_elevation_m)
    return ResolvedLocation(
        display_name=prefixed_viewpoint_name(viewpoint.kind, viewpoint.name),
        lat=viewpoint.latitude_deg,
        lon=viewpoint.longitude_deg,
        tz=timezone_name,
        persistence_key=prefixed_viewpoint_name(viewpoint.kind, viewpoint.name),
        observer_height_m=viewpoint_height_m + DEFAULT_OBSERVER_HEIGHT_M,
        kind=viewpoint.kind,
        location_height_label=location_height_label,
        location_height_m=location_height_m,
        cc=nearest_city.cc if nearest_city is not None else "",
    )


def _nominatim_result_to_location(
    query: str,
    result: dict[str, Any],
    admin1_map: dict[tuple[str, str], str],
) -> ResolvedLocation:
    try:
        lat = float(result["lat"])
        lon = float(result["lon"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("invalid Nominatim result coordinates") from exc

    name = result.get("name")
    if not isinstance(name, str) or not name:
        raise ValueError("invalid Nominatim result name")

    nearest_city = _resolve_nearest_city(lat, lon, admin1_map)
    timezone_name = nearest_city.tz if nearest_city is not None else "UTC"
    persistence_value = {
        "resolver": "nominatim",
        "query": query,
        "result": dict(result),
    }
    return ResolvedLocation(
        display_name=name,
        lat=lat,
        lon=lon,
        tz=timezone_name,
        persistence_key=name,
        observer_height_m=DEFAULT_OBSERVER_HEIGHT_M,
        kind="place",
        persistence_value=persistence_value,
        location_height_label=None,
        location_height_m=None,
        cc=nearest_city.cc if nearest_city is not None else "",
    )


def _restore_persisted_location(
    stored_location: dict[str, Any],
    admin1_map: dict[tuple[str, str], str],
) -> ResolvedLocation | None:
    if stored_location.get("resolver") != "nominatim":
        return None
    query = stored_location.get("query")
    result = stored_location.get("result")
    if not isinstance(query, str) or not isinstance(result, dict):
        logger.warning("Ignoring malformed persisted Nominatim location payload")
        return None
    try:
        location = _nominatim_result_to_location(query, result, admin1_map)
    except ValueError:
        logger.warning("Ignoring invalid persisted Nominatim location payload")
        return None
    logger.info("Restored saved place: %s", location.display_name)
    return location


def _parse_direct_coordinate_location(raw_value: str) -> tuple[float, float] | None:
    text = raw_value.strip()
    if not text:
        return None

    lat_token: str | None = None
    lon_token: str | None = None

    if ";" in text:
        parts = [part.strip() for part in text.split(";")]
        if len(parts) != 2:
            raise ValueError("Expected 'lat;lon'")
        lat_token, lon_token = parts
    elif text.startswith("@"):
        match = re.fullmatch(
            r"@\s*([+-]?\d+(?:\.\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?)(?:\s*[,/?#&].*)?",
            text,
        )
        if match is None:
            raise ValueError("Expected '@lat,lon'")
        lat_token, lon_token = match.group(1), match.group(2)
    else:
        coord_match = re.fullmatch(
            r"([+-]?\d+(?:\.\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?)",
            text,
        )
        if coord_match is not None:
            lat_token, lon_token = coord_match.group(1), coord_match.group(2)
        else:
            candidate = text
            if candidate.startswith(("maps.google.com/", "www.google.com/maps/")):
                candidate = "https://" + candidate
            parsed = urllib.parse.urlparse(candidate)
            host = parsed.netloc.lower()
            if host in {"maps.google.com", "www.google.com"}:
                full_path = parsed.path or ""
                if host == "maps.google.com":
                    if not full_path.startswith("/maps/"):
                        return None
                elif not full_path.startswith("/maps/"):
                    return None
                pin_match = re.search(r"!3d([+-]?\d+(?:\.\d+)?)!4d([+-]?\d+(?:\.\d+)?)", candidate)
                if pin_match is not None:
                    lat_token, lon_token = pin_match.group(1), pin_match.group(2)
                else:
                    center_match = re.search(
                        r"@\s*([+-]?\d+(?:\.\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?)(?:\s*[,/?#&].*)?",
                        candidate,
                    )
                    if center_match is None:
                        raise ValueError(
                            "Google Maps URL does not contain '!3d...!4d...' or '@lat,lon'"
                        )
                    lat_token, lon_token = center_match.group(1), center_match.group(2)
            else:
                return None

    def _parse_coord(token: str, dirs: str) -> float:
        token_upper = token.strip().upper()
        found = {ch for ch in token_upper if ch in "NSEW"}
        allowed = set(dirs)
        if found and not found.issubset(allowed):
            raise ValueError(f"Invalid direction in '{token}' (expected one of {sorted(allowed)}).")
        sign = -1.0 if (("S" in found) or ("W" in found)) else 1.0
        value_text = re.sub(r"[^0-9.-]", "", token)
        if not value_text:
            raise ValueError(f"No numeric value found in '{token}'")
        value = float(value_text)
        return value if value < 0 else value * sign

    assert lat_token is not None
    assert lon_token is not None
    lat = _parse_coord(lat_token, "NS")
    lon = _parse_coord(lon_token, "EW")
    if not (-90 <= lat <= 90):
        raise ValueError(f"Latitude out of range (-90 to 90): {lat}")
    if not (-180 <= lon <= 180):
        raise ValueError(f"Longitude out of range (-180 to 180): {lon}")
    return lat, lon


def _resolve_place_query(
    query: str,
    countrycode: str | None,
    language: str,
    admin1_map: dict[tuple[str, str], str],
) -> ResolvedLocation:
    logger.info("Searching Nominatim for '%s'...", query)
    try:
        results = search_nominatim(query, limit=5, countrycode=countrycode, language=language)
    except urllib.error.HTTPError as exc:
        logger.error("Nominatim HTTP error for '%s': %s %s", query, exc.code, exc.reason)
        raise LocationResolveError() from exc
    except urllib.error.URLError as exc:
        logger.error("Nominatim network error for '%s': %s", query, exc.reason)
        raise LocationResolveError() from exc
    except ValueError as exc:
        logger.error("Nominatim response error for '%s': %s", query, exc)
        raise LocationResolveError() from exc
    except Exception as exc:
        logger.error("Nominatim search failed for '%s': %s", query, exc)
        raise LocationResolveError() from exc

    if not results:
        logger.error("No Nominatim result for '%s'", query)
        raise LocationResolveError()

    logger.info("Nominatim found %d match(es) for '%s':", len(results), query)
    for index, result in enumerate(results, start=1):
        logger.info(
            "[%d] %s (lat=%.6f lon=%.6f category=%s type=%s importance=%.6f)",
            index,
            result["name"],
            result["lat"],
            result["lon"],
            result.get("category") or "-",
            result.get("type") or "-",
            float(result.get("importance") or 0.0),
        )

    try:
        location = _nominatim_result_to_location(query, results[0], admin1_map)
    except ValueError as exc:
        logger.error("Invalid top Nominatim result for '%s': %s", query, exc)
        raise LocationResolveError() from exc

    logger.info("Using top Nominatim result: %s", location.display_name)
    return location


def _tower_to_location(args_city: str, admin1_map: dict[tuple[str, str], str]) -> ResolvedLocation | None:
    return _viewpoint_to_location(resolve_tower_viewpoint(args_city), admin1_map)


def _mountain_to_location(args_city: str, admin1_map: dict[tuple[str, str], str]) -> ResolvedLocation | None:
    return _viewpoint_to_location(resolve_mountain_viewpoint(args_city), admin1_map)


def resolve_launch_location(
    args_city: str | None,
    place_query: str | None = None,
    place_countrycode: str | None = None,
    place_lang: str = "en",
    use_building_top: bool = False,
) -> ResolvedLocation:
    last_city = load_last_city()
    stored_location: dict[str, Any] | None = None
    if not args_city and place_query is None:
        if isinstance(last_city, str):
            args_city = last_city or "Tokyo"
        elif isinstance(last_city, dict):
            stored_location = last_city
        else:
            args_city = "Tokyo"

    resolved_location: ResolvedLocation | None = None
    persist_location = False

    parsed_coords: tuple[float, float] | None = None
    if args_city:
        try:
            parsed_coords = _parse_direct_coordinate_location(args_city)
        except ValueError as exc:
            logger.error("Invalid latitude/longitude format: '%s'. %s", args_city, exc)
            raise LocationResolveError() from exc

    logger.info("Loading city data...")
    try:
        admin1_map = load_admin1_names(CITY_ADMIN1_CODES_FILE)
    except FileNotFoundError as exc:
        logger.error("Fail to load admin1CodesASCII.txt.")
        raise LocationResolveError() from exc

    recs: List[CityRec] = []
    try:
        if parsed_coords is not None:
            lat, lon = parsed_coords
            nearest_city = _resolve_nearest_city(lat, lon, admin1_map)
            timezone_name = nearest_city.tz if nearest_city is not None else "UTC"
            logger.info("Parsed location: Lat=%.6f, Lon=%.6f, Timezone=%s", lat, lon, timezone_name)
            return _maybe_apply_building_top_viewpoint(
                ResolvedLocation(
                display_name=f"Lat: {lat:.2f}, Lon: {lon:.2f}",
                lat=lat,
                lon=lon,
                tz=timezone_name,
                persistence_key=f"{lat:.6f};{lon:.6f}",
                observer_height_m=DEFAULT_OBSERVER_HEIGHT_M,
                kind="coords",
                location_height_label=None,
                location_height_m=None,
                cc=nearest_city.cc if nearest_city is not None else "",
                ),
                enabled=use_building_top,
            )
        if stored_location is not None:
            resolved_location = _restore_persisted_location(stored_location, admin1_map)
            if resolved_location is None:
                args_city = "Tokyo"

        if place_query is not None:
            resolved_location = _resolve_place_query(place_query, place_countrycode, place_lang, admin1_map)
            persist_location = True
        elif resolved_location is None:
            assert args_city is not None
            explicit_viewpoint = split_prefixed_viewpoint(args_city)
            tower_query = args_city.startswith("wikidata:") or re.match(r"^Q\d+$", args_city) is not None

            if explicit_viewpoint is not None:
                explicit_kind, explicit_name = explicit_viewpoint
                if explicit_kind == "tower":
                    tower_location = _tower_to_location(explicit_name, admin1_map)
                    if tower_location is None:
                        logger.error("No tower found for '%s'", explicit_name)
                        raise LocationResolveError()
                    resolved_location = tower_location
                else:
                    mountain_location = _mountain_to_location(explicit_name, admin1_map)
                    if mountain_location is None:
                        logger.error("No mountain found for '%s'", explicit_name)
                        raise LocationResolveError()
                    resolved_location = mountain_location
                persist_location = True
            elif tower_query:
                tower_location = _tower_to_location(args_city, admin1_map)
                mountain_location = _mountain_to_location(args_city, admin1_map)
                if tower_location is not None:
                    resolved_location = tower_location
                elif mountain_location is not None:
                    resolved_location = mountain_location
                else:
                    logger.error("No tower or mountain found for '%s'", args_city)
                    raise LocationResolveError()
                persist_location = True
            elif re.match(r"^\d+$", args_city):
                rec = resolve_city_by_geonameid(int(args_city), CITY_COORD_FILE)
                if rec:
                    recs.append(rec)
                else:
                    logger.error("No city found for geonameid %s", args_city)
                    raise LocationResolveError()
            else:
                if "/" not in args_city:
                    recs = resolve_city_by_name(args_city, CITY_COORD_FILE, admin1_map)
                else:
                    recs = resolve_city(args_city, CITY_COORD_FILE, admin1_map)
                if recs:
                    logger.info("Found %d match(es) for '%s':", len(recs), args_city)
                    for rec in recs:
                        logger.info(
                            "- %s/%s, lat: %.6f, lon: %.6f, tz: %s  (geonameid=%s)",
                            rec.cc,
                            rec.name,
                            rec.lat,
                            rec.lon,
                            rec.tz,
                            rec.geonameid,
                        )
                    if len(recs) > 1:
                        logger.warning("Multiple matches found for '%s'", args_city)
                else:
                    tower_location = _tower_to_location(args_city, admin1_map)
                    mountain_location = _mountain_to_location(args_city, admin1_map)
                    if tower_location is not None:
                        resolved_location = tower_location
                        persist_location = True
                    elif mountain_location is not None:
                        resolved_location = mountain_location
                        persist_location = True
                    else:
                        logger.error("No match for '%s'", args_city)
                        raise LocationResolveError()
    except FileNotFoundError as exc:
        logger.error("Fail to load cities1000.txt.")
        raise LocationResolveError() from exc

    if resolved_location is None:
        city = recs[0]
        city_str = f"{city.cc}/{city.name}"
        resolved_location = ResolvedLocation(
            display_name=city_str,
            lat=city.lat,
            lon=city.lon,
            tz=city.tz,
            persistence_key=city_str,
            observer_height_m=DEFAULT_OBSERVER_HEIGHT_M,
            kind="city",
            location_height_label=None,
            location_height_m=None,
            cc=city.cc,
        )
        persist_location = True

    resolved_location = _maybe_apply_building_top_viewpoint(
        resolved_location,
        enabled=use_building_top,
    )

    if persist_location:
        save_last_city(
            resolved_location.persistence_value
            if resolved_location.persistence_value is not None
            else resolved_location.persistence_key
        )
        if resolved_location.kind == "tower":
            logger.info("Tower: %s", resolved_location.persistence_key)
        elif resolved_location.kind == "mountain":
            logger.info("Mountain: %s", resolved_location.persistence_key)
        elif resolved_location.kind == "place":
            logger.info("Place: %s", resolved_location.persistence_key)
        else:
            logger.info("City: %s", resolved_location.persistence_key)
    return resolved_location
