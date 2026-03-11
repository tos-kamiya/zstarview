import logging
import math
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import List, Optional, Tuple

import polars as pl

from .catalog import load_dso_catalog, load_star_catalog
from .config import load_last_city, save_last_city
from .viewpoints import Viewpoint
from .viewpoints import prefixed_viewpoint_name, split_prefixed_viewpoint
from .mountain_viewpoints import resolve_mountain_viewpoint
from .paths import (
    CITY_ADMIN1_CODES_FILE,
    CITY_COORD_FILE,
    DSO_CSV_FILE,
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

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
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
    return ResolvedLocation(
        display_name=prefixed_viewpoint_name(viewpoint.kind, viewpoint.name),
        lat=viewpoint.latitude_deg,
        lon=viewpoint.longitude_deg,
        tz=timezone_name,
        persistence_key=prefixed_viewpoint_name(viewpoint.kind, viewpoint.name),
        observer_height_m=viewpoint.observer_height_m + DEFAULT_OBSERVER_HEIGHT_M,
        kind=viewpoint.kind,
        cc=nearest_city.cc if nearest_city is not None else "",
    )


def _tower_to_location(args_city: str, admin1_map: dict[tuple[str, str], str]) -> ResolvedLocation | None:
    return _viewpoint_to_location(resolve_tower_viewpoint(args_city), admin1_map)


def _mountain_to_location(args_city: str, admin1_map: dict[tuple[str, str], str]) -> ResolvedLocation | None:
    return _viewpoint_to_location(resolve_mountain_viewpoint(args_city), admin1_map)


def _startup_resolve_city(args_city: Optional[str]) -> ResolvedLocation:
    """
    Resolve target city from CLI or last used city.

    Also handles direct latitude/longitude input like "35.68;139.76" or "N35.68;E139.76".
    """
    last_city = load_last_city()
    if not args_city:
        args_city = last_city or "Tokyo"

    resolved_location: ResolvedLocation | None = None
    persist_location = False

    if ";" in args_city:
        try:
            lat_str, lon_str = [s.strip() for s in args_city.split(";")]

            def _parse_coord(s: str, dirs: str) -> float:
                s_upper = s.strip().upper()
                found = {ch for ch in s_upper if ch in "NSEW"}
                allowed = set(dirs)
                if found and not found.issubset(allowed):
                    raise ValueError(
                        f"Invalid direction in '{s}' (expected one of {sorted(allowed)})."
                    )
                sign = -1.0 if (("S" in found) or ("W" in found)) else 1.0
                val_str = re.sub(r"[^0-9.-]", "", s)
                if not val_str:
                    raise ValueError(f"No numeric value found in '{s}'")
                val = float(val_str)
                return val if val < 0 else val * sign

            lat = _parse_coord(lat_str, "NS")
            lon = _parse_coord(lon_str, "EW")

            if not (-90 <= lat <= 90):
                raise ValueError(f"Latitude out of range (-90 to 90): {lat}")
            if not (-180 <= lon <= 180):
                raise ValueError(f"Longitude out of range (-180 to 180): {lon}")

            logger.info("Parsed location: Lat=%.6f, Lon=%.6f, Timezone=UTC", lat, lon)

            name = f"Lat: {lat:.2f}, Lon: {lon:.2f}"
            resolved_location = ResolvedLocation(
                display_name=name,
                lat=lat,
                lon=lon,
                tz="UTC",
                persistence_key=f"{lat:.6f};{lon:.6f}",
                observer_height_m=DEFAULT_OBSERVER_HEIGHT_M,
                kind="coords",
            )
            return resolved_location
        except (ValueError, IndexError) as exc:
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
            cc=city.cc,
        )
        persist_location = True

    if persist_location:
        save_last_city(resolved_location.persistence_key)
        if resolved_location.kind == "tower":
            logger.info("Tower: %s", resolved_location.persistence_key)
        elif resolved_location.kind == "mountain":
            logger.info("Mountain: %s", resolved_location.persistence_key)
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
    args_datetime: Optional[str], args_days: int, args_hours: int
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

        if tz_str:
            try:
                tz = parse_tz_string(tz_str)
                dt_local = dt_naive.replace(tzinfo=tz)
                target_time_utc = dt_local.astimezone(timezone.utc)
            except Exception as exc:
                logger.error("Invalid timezone '%s'. %s", tz_str, exc)
                raise StartupAbortError() from exc
        else:
            target_time_utc = dt_naive.replace(tzinfo=timezone.utc)

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
    logger.info("Loaded %d stars (Vmag ≤ %s)", len(star_catalog), limit_str)
    return star_catalog


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
