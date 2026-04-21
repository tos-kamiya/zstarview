from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from collections.abc import Iterable
import logging
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen
from urllib.error import URLError

from skyfield.api import EarthSatellite
import skyfield.api

from ..satellite_constants import (
    SATELLITE_FETCH_TIMEOUT_SECONDS,
    SATELLITE_HORIZONS_CACHE_KEY,
    SATELLITE_ISS_CACHE_KEY,
)
from .types import SatelliteOmmRecord

logger = logging.getLogger(__name__)

CELESTRAK_GP_JSON_URL = "https://celestrak.org/NORAD/elements/gp.php"
WHERETHEISS_API_URL = "https://api.wheretheiss.at/v1"
HORIZONS_LOOKUP_API_URL = "https://ssd.jpl.nasa.gov/api/horizons_lookup.api"
HORIZONS_API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
CELESTRAK_GROUP_BY_KEY = {
    SATELLITE_ISS_CACHE_KEY: "stations",
}
HORIZONS_TARGETS_BY_KEY = {
    SATELLITE_HORIZONS_CACHE_KEY: (
        ("JWST", ("JWST", "James Webb Space Telescope", "James Webb")),
        ("Voyager 1", ("Voyager 1", "Voyager-1")),
        ("Voyager 2", ("Voyager 2", "Voyager-2")),
        ("Parker", ("Parker Solar Probe", "Parker", "Solar Probe Plus")),
    ),
}
_ISS_NORAD_CAT_ID = "25544"
_SOURCE_KEY = "_SOURCE"
_TLE_LINE1_KEY = "_TLE_LINE1"
_TLE_LINE2_KEY = "_TLE_LINE2"


def build_celestrak_group_url(
    group_name: str,
    *,
    base_url: str = CELESTRAK_GP_JSON_URL,
    data_format: str = "json",
) -> str:
    return f"{base_url}?{urlencode({'GROUP': group_name, 'FORMAT': data_format})}"


def build_horizons_lookup_url(
    search_text: str,
    *,
    base_url: str = HORIZONS_LOOKUP_API_URL,
    group: str = "sct",
) -> str:
    return f"{base_url}?{urlencode({'sstr': search_text, 'group': group, 'format': 'json'}, quote_via=quote)}"


def build_horizons_observer_url(
    command: str,
    *,
    target_time_utc: datetime,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    base_url: str = HORIZONS_API_URL,
) -> str:
    observer_height_km = float(observer_height_m) / 1000.0
    site_coord = f"{float(observer_lon):.6f},{float(observer_lat):.6f},{observer_height_km:.6f}"
    tlist = target_time_utc.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    params = {
        "format": "json",
        "COMMAND": f"'{command}'",
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "OBSERVER",
        "CENTER": "'coord'",
        "COORD_TYPE": "GEODETIC",
        "SITE_COORD": f"'{site_coord}'",
        "TIME_TYPE": "UT",
        "TIME_DIGITS": "SECONDS",
        "EXTRA_PREC": "YES",
        "CSV_FORMAT": "YES",
        "QUANTITIES": "'4'",
        "TLIST": f"'{tlist}'",
    }
    return f"{base_url}?{urlencode(params, quote_via=quote)}"


def fetch_celestrak_group_omm(
    group_name: str,
    *,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    base_url: str = CELESTRAK_GP_JSON_URL,
) -> list[SatelliteOmmRecord]:
    request = Request(
        build_celestrak_group_url(group_name, base_url=base_url),
        headers={"Accept": "application/json"},
    )
    with urlopen(request, timeout=float(timeout_s)) as response:
        payload = json.load(response)
    return normalize_celestrak_omm_payload(payload)


def fetch_horizons_lookup(
    search_text: str,
    *,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    base_url: str = HORIZONS_LOOKUP_API_URL,
    group: str = "sct",
) -> dict[str, object]:
    request = Request(
        build_horizons_lookup_url(search_text, base_url=base_url, group=group),
        headers={"Accept": "application/json"},
    )
    with urlopen(request, timeout=float(timeout_s)) as response:
        payload = json.load(response)
    if not isinstance(payload, dict):
        raise RuntimeError("Horizons lookup returned an invalid payload")
    return payload


def fetch_horizons_observer_csv(
    command: str,
    *,
    target_time_utc: datetime,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    base_url: str = HORIZONS_API_URL,
) -> list[list[str]]:
    request = Request(
        build_horizons_observer_url(
            command,
            target_time_utc=target_time_utc,
            observer_lat=observer_lat,
            observer_lon=observer_lon,
            observer_height_m=observer_height_m,
            base_url=base_url,
        ),
        headers={"Accept": "application/json"},
    )
    with urlopen(request, timeout=float(timeout_s)) as response:
        payload = json.load(response)
    if not isinstance(payload, dict):
        return []
    error = payload.get("error")
    if error:
        raise RuntimeError(str(error).strip())
    result_text = str(payload.get("result", ""))
    return _parse_horizons_observer_csv(result_text)


def build_wheretheiss_tle_url(
    norad_cat_id: str = _ISS_NORAD_CAT_ID,
    *,
    base_url: str = WHERETHEISS_API_URL,
) -> str:
    return f"{base_url}/satellites/{norad_cat_id}/tles"


def fetch_wheretheiss_iss_tle(
    *,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    base_url: str = WHERETHEISS_API_URL,
) -> list[SatelliteOmmRecord]:
    request = Request(
        build_wheretheiss_tle_url(base_url=base_url),
        headers={"Accept": "application/json"},
    )
    with urlopen(request, timeout=float(timeout_s)) as response:
        payload = json.load(response)
    return normalize_wheretheiss_tle_payload(payload)


def fetch_celestrak_group_by_key(
    group_key: str,
    *,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    base_url: str = CELESTRAK_GP_JSON_URL,
) -> list[SatelliteOmmRecord]:
    group_name = CELESTRAK_GROUP_BY_KEY[group_key]
    records = fetch_celestrak_group_omm(group_name, timeout_s=timeout_s, base_url=base_url)
    return filter_records_for_group(group_key, records)


def fetch_iss_records(
    group_key: str,
    *,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    target_time_utc: datetime | None = None,
    observer_lat: float | None = None,
    observer_lon: float | None = None,
    observer_height_m: float | None = None,
    wheretheiss_base_url: str = WHERETHEISS_API_URL,
    celestrak_base_url: str = CELESTRAK_GP_JSON_URL,
) -> list[SatelliteOmmRecord]:
    del target_time_utc, observer_lat, observer_lon, observer_height_m
    if group_key != SATELLITE_ISS_CACHE_KEY:
        raise KeyError(group_key)
    try:
        return fetch_wheretheiss_iss_tle(timeout_s=timeout_s, base_url=wheretheiss_base_url)
    except Exception as primary_exc:
        _log_fetch_attempt_failure("wheretheiss.at", primary_exc)
    try:
        return fetch_celestrak_group_by_key(
            group_key,
            timeout_s=timeout_s,
            base_url=celestrak_base_url,
        )
    except Exception as fallback_exc:
        _log_fetch_attempt_failure("CelesTrak fallback", fallback_exc)
        raise


def fetch_horizons_records(
    group_key: str,
    *,
    target_time_utc: datetime,
    observer_lat: float | None,
    observer_lon: float | None,
    observer_height_m: float | None,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    lookup_base_url: str = HORIZONS_LOOKUP_API_URL,
    horizons_base_url: str = HORIZONS_API_URL,
) -> list[SatelliteOmmRecord]:
    if group_key != SATELLITE_HORIZONS_CACHE_KEY:
        raise KeyError(group_key)
    if observer_lat is None or observer_lon is None or observer_height_m is None:
        raise ValueError("Horizons spacecraft fetch requires observer coordinates")
    records: list[SatelliteOmmRecord] = []
    for label, aliases in HORIZONS_TARGETS_BY_KEY[SATELLITE_HORIZONS_CACHE_KEY]:
        resolved = _resolve_horizons_target(
            label,
            aliases,
            timeout_s=timeout_s,
            lookup_base_url=lookup_base_url,
        )
        if resolved is None:
            continue
        command = str(resolved.get("spkid", "")).strip()
        if not command:
            continue
        try:
            rows = fetch_horizons_observer_csv(
                command,
                target_time_utc=target_time_utc,
                observer_lat=observer_lat,
                observer_lon=observer_lon,
                observer_height_m=observer_height_m,
                timeout_s=timeout_s,
                base_url=horizons_base_url,
            )
        except Exception as exc:
            _log_fetch_attempt_failure(f"Horizons {label}", exc)
            continue
        parsed_altaz: tuple[float, float] | None = None
        for row in rows:
            parsed_altaz = _extract_altaz_from_csv_row(row)
            if parsed_altaz[0] is not None and parsed_altaz[1] is not None:
                break
        if parsed_altaz is None or parsed_altaz[0] is None or parsed_altaz[1] is None:
            continue
        records.append(
            {
                "OBJECT_NAME": label,
                "HORIZONS_TARGET_NAME": str(resolved.get("name", label)).strip() or label,
                "HORIZONS_SPKID": command,
                "EPOCH": target_time_utc.astimezone(timezone.utc).isoformat(),
                "ALT_DEG": float(parsed_altaz[0]),
                "AZ_DEG": float(parsed_altaz[1]),
                _SOURCE_KEY: "horizons",
            }
        )
    if not records:
        raise RuntimeError("Horizons fetch returned no spacecraft records")
    return records


def normalize_celestrak_omm_payload(payload: object) -> list[SatelliteOmmRecord]:
    if not isinstance(payload, list):
        return []
    records: list[SatelliteOmmRecord] = []
    for row in payload:
        if isinstance(row, dict):
            normalized = dict(row)
            normalized.setdefault(_SOURCE_KEY, "celestrak")
            records.append(normalized)
    return records


def normalize_wheretheiss_tle_payload(payload: object) -> list[SatelliteOmmRecord]:
    if not isinstance(payload, dict):
        return []
    line1 = str(payload.get("line1", "")).strip()
    line2 = str(payload.get("line2", "")).strip()
    if not line1 or not line2:
        return []
    header = str(payload.get("header", "")).strip() or "ISS (ZARYA)"
    record: SatelliteOmmRecord = {
        "OBJECT_NAME": header,
        "NORAD_CAT_ID": str(payload.get("id", _ISS_NORAD_CAT_ID)).strip() or _ISS_NORAD_CAT_ID,
        _TLE_LINE1_KEY: line1,
        _TLE_LINE2_KEY: line2,
        _SOURCE_KEY: "wheretheiss",
    }
    tle_timestamp = payload.get("tle_timestamp")
    if isinstance(tle_timestamp, (int, float)):
        record["TLE_TIMESTAMP"] = float(tle_timestamp)
    return [record]


def _resolve_horizons_target(
    label: str,
    aliases: tuple[str, ...],
    *,
    timeout_s: float,
    lookup_base_url: str,
) -> dict[str, object] | None:
    for alias in aliases:
        try:
            payload = fetch_horizons_lookup(alias, timeout_s=timeout_s, base_url=lookup_base_url)
        except Exception as exc:
            _log_fetch_attempt_failure(f"Horizons lookup {label}", exc)
            continue
        if not isinstance(payload, dict):
            continue
        result = payload.get("result")
        if not isinstance(result, list) or not result:
            continue
        exact = _pick_exact_horizons_match(alias, result)
        if exact is not None:
            return exact
        first = result[0]
        if isinstance(first, dict):
            return first
    return None


def _pick_exact_horizons_match(search_text: str, results: list[object]) -> dict[str, object] | None:
    query = _normalize_horizons_name(search_text)
    if not query:
        return None
    for item in results:
        if not isinstance(item, dict):
            continue
        name = _normalize_horizons_name(str(item.get("name", "")))
        raw_aliases = item.get("alias", []) or []
        aliases = {
            _normalize_horizons_name(alias)
            for alias in raw_aliases
            if isinstance(alias, str)
        }
        if query == name or query in aliases:
            return item
    return None


def _parse_horizons_observer_csv(result_text: str) -> list[list[str]]:
    lines: list[str] = []
    inside = False
    for raw_line in result_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line == "$$SOE":
            inside = True
            continue
        if line == "$$EOE":
            break
        if inside:
            lines.append(line)
    if not lines:
        return []
    return list(csv.reader(lines))


def _extract_altaz_from_csv_row(row: list[str]) -> tuple[float | None, float | None]:
    numeric_values: list[float] = []
    for value in row:
        parsed = _parse_float(value)
        if parsed is not None:
            numeric_values.append(parsed)
    if len(numeric_values) < 2:
        return None, None
    # Horizons observer CSV reports azimuth first and elevation second.
    return numeric_values[1], numeric_values[0]


def _parse_float(value: object) -> float | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def _normalize_horizons_name(value: str) -> str:
    return " ".join(part for part in value.casefold().replace("(", " ").replace(")", " ").split() if part)


def extract_element_epoch_utc(records: Iterable[SatelliteOmmRecord]) -> datetime | None:
    epochs: list[datetime] = []
    for record in records:
        raw = str(record.get("EPOCH", "")).strip()
        if raw:
            try:
                parsed = datetime.fromisoformat(raw.replace("Z", "+00:00"))
            except ValueError:
                parsed = None
            if parsed is not None:
                if parsed.tzinfo is None:
                    parsed = parsed.replace(tzinfo=timezone.utc)
                else:
                    parsed = parsed.astimezone(timezone.utc)
                epochs.append(parsed)
                continue
        tle_timestamp = record.get("TLE_TIMESTAMP")
        if isinstance(tle_timestamp, (int, float)):
            epochs.append(datetime.fromtimestamp(float(tle_timestamp), tz=timezone.utc))
    if not epochs:
        return None
    return max(epochs)


def extract_record_source(records: Iterable[SatelliteOmmRecord], *, default: str = "celestrak") -> str:
    for record in records:
        source = str(record.get(_SOURCE_KEY, "")).strip()
        if source:
            return source
    return default


def filter_records_for_group(
    group_key: str,
    records: Iterable[SatelliteOmmRecord],
) -> list[SatelliteOmmRecord]:
    normalized = [dict(record) for record in records]
    if group_key != SATELLITE_ISS_CACHE_KEY:
        return normalized
    return [
        record
        for record in normalized
        if str(record.get("NORAD_CAT_ID", "")).strip() == _ISS_NORAD_CAT_ID
    ]


def build_earth_satellites(
    records: Iterable[SatelliteOmmRecord],
    *,
    ts: object | None = None,
) -> list[EarthSatellite]:
    timescale = ts or skyfield.api.load.timescale()
    satellites: list[EarthSatellite] = []
    for record in records:
        line1 = str(record.get(_TLE_LINE1_KEY, "")).strip()
        line2 = str(record.get(_TLE_LINE2_KEY, "")).strip()
        if line1 and line2:
            satellites.append(
                EarthSatellite(
                    line1,
                    line2,
                    str(record.get("OBJECT_NAME", "")).strip() or "ISS",
                    timescale,
                )
            )
            continue
        satellites.append(EarthSatellite.from_omm(timescale, record))
    return satellites


def _log_fetch_attempt_failure(source_name: str, exc: Exception) -> None:
    detail = f"{exc.__class__.__name__}: {exc!s}" if str(exc).strip() else exc.__class__.__name__
    if _is_timeout_error(exc):
        logger.warning("Satellite fetch failed via %s: %s", source_name, detail)
        return
    logger.info("Satellite fetch failed via %s: %s", source_name, detail)


def _is_timeout_error(exc: Exception) -> bool:
    if isinstance(exc, TimeoutError):
        return True
    if isinstance(exc, URLError):
        reason = getattr(exc, "reason", None)
        if isinstance(reason, TimeoutError):
            return True
    return "timed out" in str(exc).lower()
