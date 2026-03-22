import logging
import math
import re
import sys
import urllib.error
import urllib.parse
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, List, Optional, Tuple

import polars as pl

from .astro import _starfield_load
from .catalog import load_dso_catalog, load_star_catalog
from .config import load_last_city, save_last_city
from .nominatim_search import search as search_nominatim
from .viewpoints import Viewpoint
from .viewpoints import prefixed_viewpoint_name, split_prefixed_viewpoint
from .mountain_viewpoints import resolve_mountain_viewpoint
from .paths import (
    CITY_ADMIN1_CODES_FILE,
    CITY_COORD_FILE,
    DSO_CSV_FILE,
    EPHEMERIS_FILENAME,
    LOG_PATH,
    STARS_CSV_FILE,
)
from .tower_viewpoints import resolve_tower_viewpoint
from .utils.resolve_city import (
    CityRec,
    load_admin1_names,
    resolve_city,
    resolve_city_by_geonameid,
    resolve_city_by_name,
)
from .utils.timezone_parser import parse_tz_string

logger = logging.getLogger(__name__)
DEFAULT_OBSERVER_HEIGHT_M = 1.7


class StartupAbortError(Exception):
    """Abort the startup sequence (handled by main to show splash for 3s)."""


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


def setup_root_logger() -> logging.Logger:
    """Configure and return the root logger for the application."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stderr,
    )
    root_logger = logging.getLogger()

    log_dir = Path(LOG_PATH)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "app.log"

    try:
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
    except OSError as exc:
        logger.warning("Failed to open log file %s: %s", log_path, exc)
        return root_logger

    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
    )
    root_logger.addHandler(file_handler)

    logger.info("Logging to file: %s", log_path)
    return root_logger


def _format_splash_location(city: ResolvedLocation) -> str:
    """Create a concise location label for splash screen context."""
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
    nearest_city = _resolve_nearest_city(
        viewpoint.latitude_deg,
        viewpoint.longitude_deg,
        admin1_map,
    )
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
            r"@([+-]?\d+(?:\.\d+)?),([+-]?\d+(?:\.\d+)?)(?:[,/?#&].*)?",
            text,
        )
        if match is None:
            raise ValueError("Expected '@lat,lon'")
        lat_token, lon_token = match.group(1), match.group(2)
    else:
        candidate = text
        if candidate.startswith(("maps.google.com/", "www.google.com/maps/", "google.com/maps/")):
            candidate = "https://" + candidate
        parsed = urllib.parse.urlparse(candidate)
        host = parsed.netloc.lower()
        if host in {"maps.google.com", "www.google.com", "google.com"}:
            full_path = parsed.path or ""
            if host == "maps.google.com":
                if not full_path.startswith("/maps/"):
                    return None
            elif not full_path.startswith("/maps/"):
                return None
            match = re.search(
                r"/maps/@([+-]?\d+(?:\.\d+)?),([+-]?\d+(?:\.\d+)?)(?:[,/?#&].*)?$",
                full_path,
            )
            if match is None:
                raise ValueError("Google Maps URL does not contain '/maps/@lat,lon'")
            lat_token, lon_token = match.group(1), match.group(2)
        else:
            return None

    def _parse_coord(token: str, dirs: str) -> float:
        token_upper = token.strip().upper()
        found = {ch for ch in token_upper if ch in "NSEW"}
        allowed = set(dirs)
        if found and not found.issubset(allowed):
            raise ValueError(
                f"Invalid direction in '{token}' (expected one of {sorted(allowed)})."
            )
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
        results = search_nominatim(
            query,
            limit=5,
            countrycode=countrycode,
            language=language,
        )
    except urllib.error.HTTPError as exc:
        logger.error("Nominatim HTTP error for '%s': %s %s", query, exc.code, exc.reason)
        raise StartupAbortError() from exc
    except urllib.error.URLError as exc:
        logger.error("Nominatim network error for '%s': %s", query, exc.reason)
        raise StartupAbortError() from exc
    except ValueError as exc:
        logger.error("Nominatim response error for '%s': %s", query, exc)
        raise StartupAbortError() from exc
    except Exception as exc:
        logger.error("Nominatim search failed for '%s': %s", query, exc)
        raise StartupAbortError() from exc

    if not results:
        logger.error("No Nominatim result for '%s'", query)
        raise StartupAbortError()

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
        raise StartupAbortError() from exc

    logger.info("Using top Nominatim result: %s", location.display_name)
    return location


def _tower_to_location(args_city: str, admin1_map: dict[tuple[str, str], str]) -> ResolvedLocation | None:
    return _viewpoint_to_location(resolve_tower_viewpoint(args_city), admin1_map)


def _mountain_to_location(args_city: str, admin1_map: dict[tuple[str, str], str]) -> ResolvedLocation | None:
    return _viewpoint_to_location(resolve_mountain_viewpoint(args_city), admin1_map)


def _startup_resolve_city(
    args_city: Optional[str],
    place_query: str | None = None,
    place_countrycode: str | None = None,
    place_lang: str = "en",
) -> ResolvedLocation:
    """
    Resolve target city from CLI or last used city.

    Also handles direct latitude/longitude input like "35.68;139.76" or "N35.68;E139.76".

    Kept importable despite the leading underscore because focused startup tests
    exercise this resolver directly.
    """
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
            raise StartupAbortError() from exc

    logger.info("Loading city data...")
    try:
        admin1_map = load_admin1_names(CITY_ADMIN1_CODES_FILE)
    except FileNotFoundError as exc:
        logger.error("Fail to load admin1CodesASCII.txt.")
        raise StartupAbortError() from exc

    recs: List[CityRec] = []
    try:
        if parsed_coords is not None:
            lat, lon = parsed_coords
            nearest_city = _resolve_nearest_city(lat, lon, admin1_map)
            timezone_name = nearest_city.tz if nearest_city is not None else "UTC"
            logger.info(
                "Parsed location: Lat=%.6f, Lon=%.6f, Timezone=%s",
                lat,
                lon,
                timezone_name,
            )
            return ResolvedLocation(
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
            )
        if stored_location is not None:
            resolved_location = _restore_persisted_location(stored_location, admin1_map)
            if resolved_location is None:
                args_city = "Tokyo"

        if place_query is not None:
            resolved_location = _resolve_place_query(
                place_query,
                place_countrycode,
                place_lang,
                admin1_map,
            )
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
                        raise StartupAbortError()
                    resolved_location = tower_location
                else:
                    mountain_location = _mountain_to_location(explicit_name, admin1_map)
                    if mountain_location is None:
                        logger.error("No mountain found for '%s'", explicit_name)
                        raise StartupAbortError()
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
                    raise StartupAbortError()
                persist_location = True
            elif re.match(r"^\d+$", args_city):
                rec = resolve_city_by_geonameid(int(args_city), CITY_COORD_FILE)
                if rec:
                    recs.append(rec)
                else:
                    logger.error("No city found for geonameid %s", args_city)
                    raise StartupAbortError()
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
                        raise StartupAbortError()
    except FileNotFoundError as exc:
        logger.error("Fail to load cities1000.txt.")
        raise StartupAbortError() from exc

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


def _parse_flexible_time(time_str: str) -> Tuple[int, int, int]:
    """Parse time string that may omit minutes/seconds."""
    match = re.fullmatch(r"\s*(\d{1,2})(?::(\d{1,2}))?(?::(\d{1,2}))?\s*", time_str)
    if not match:
        raise ValueError(f"Invalid time: {time_str!r}. Use HH, HH:MM, or HH:MM:SS.")
    hour = int(match.group(1))
    minute = int(match.group(2)) if match.group(2) is not None else 0
    second = int(match.group(3)) if match.group(3) is not None else 0

    if not (0 <= hour <= 23):
        raise ValueError(f"Hour out of range (0-23): {hour}")
    if not (0 <= minute <= 59):
        raise ValueError(f"Minute out of range (0-59): {minute}")
    if not (0 <= second <= 59):
        raise ValueError(f"Second out of range (0-59): {second}")
    return hour, minute, second


def _startup_parse_time_arguments(
    args_datetime: Optional[str],
    args_days: int,
    args_hours: int,
    *,
    timezone_name: str = "UTC",
    timezone_override: str | None = None,
) -> timedelta:
    """
    Parse time-related arguments and return a timedelta from now.

    Supports --datetime in 'YYYY-MM-DD HH[:MM[:SS]] [TZ]' (TZ optional).
    """
    if not args_datetime:
        return timedelta(days=args_days, hours=args_hours)

    if args_hours != 0 or args_days != 0:
        logger.error("Invalid option: --datetime cannot be used with --hours or --days.")
        raise StartupAbortError()

    try:
        parts = args_datetime.split()
        if len(parts) < 2:
            raise ValueError("Missing time part. Use 'YYYY-MM-DD HH[:MM[:SS]] [TZ]'.")

        date_str = parts[0]
        time_str = parts[1]
        tz_str = parts[2] if len(parts) >= 3 else None

        try:
            date_only = datetime.strptime(date_str, "%Y-%m-%d")
        except ValueError as exc:
            raise ValueError(f"Invalid date: {date_str!r}. Use YYYY-MM-DD.") from exc

        hour, minute, second = _parse_flexible_time(time_str)
        dt_naive = date_only.replace(hour=hour, minute=minute, second=second)

        effective_tz_str = timezone_override if timezone_override is not None else tz_str

        if effective_tz_str:
            try:
                tz = parse_tz_string(effective_tz_str)
                dt_local = dt_naive.replace(tzinfo=tz)
                target_time_utc = dt_local.astimezone(timezone.utc)
            except Exception as exc:
                logger.error("Invalid timezone '%s'. %s", effective_tz_str, exc)
                raise StartupAbortError() from exc
        else:
            try:
                tz = parse_tz_string(timezone_name)
            except Exception as exc:
                logger.error("Invalid timezone '%s'. %s", timezone_name, exc)
                raise StartupAbortError() from exc
            target_time_utc = dt_naive.replace(tzinfo=tz).astimezone(timezone.utc)

        now_utc = datetime.now(timezone.utc)
        return target_time_utc - now_utc
    except ValueError as exc:
        logger.error("%s Input was: %r", exc, args_datetime)
        raise StartupAbortError() from exc


def _startup_load_stars(args_vmag_limit: Optional[float]) -> pl.DataFrame:
    """Load the star catalog from the source file."""
    logger.info("Loading city and star data...")

    try:
        star_catalog = load_star_catalog(STARS_CSV_FILE, vmag_threshold=args_vmag_limit)
    except FileNotFoundError as exc:
        logger.error("Fail to load file: %s", STARS_CSV_FILE)
        raise StartupAbortError() from exc

    limit_str = args_vmag_limit if args_vmag_limit is not None else "no limit"
    logger.info("Loaded %d stars (Vmag <= %s)", len(star_catalog), limit_str)
    return star_catalog


def _startup_verify_ephemeris() -> None:
    """Ensure the Skyfield ephemeris is available before background sky updates begin."""
    logger.info("Checking ephemeris cache...")
    try:
        _starfield_load(EPHEMERIS_FILENAME)
    except OSError as exc:
        logger.error(
            "Failed to load ephemeris %s. The app cannot continue until this file "
            "is available in the cache. %s",
            EPHEMERIS_FILENAME,
            exc,
        )
        raise StartupAbortError() from exc
    except Exception as exc:
        logger.error("Unexpected ephemeris load failure for %s: %s", EPHEMERIS_FILENAME, exc)
        raise StartupAbortError() from exc
    logger.info("Ephemeris ready: %s", EPHEMERIS_FILENAME)


def _startup_load_dso() -> Optional[pl.DataFrame]:
    """Load optional DSO catalog used for shape overlays."""
    try:
        dso_catalog = load_dso_catalog(DSO_CSV_FILE)
    except FileNotFoundError:
        logger.info("DSO catalog not found: %s (shape overlay disabled)", DSO_CSV_FILE)
        return None
    except Exception:
        logger.warning("Failed to load DSO catalog: %s", DSO_CSV_FILE, exc_info=True)
        return None
    logger.info("Loaded %d DSO rows", len(dso_catalog))
    return dso_catalog
