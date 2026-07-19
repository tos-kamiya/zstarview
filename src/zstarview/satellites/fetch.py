from __future__ import annotations

import csv
import json
import logging
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from datetime import datetime, timezone
from time import monotonic, sleep
from urllib.error import URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

import skyfield.api
from skyfield.api import EarthSatellite

from ..user_agent import build_user_agent
from ..satellite_constants import (
    SATELLITE_FETCH_TIMEOUT_SECONDS,
    SATELLITE_HORIZONS_CACHE_KEY,
    SATELLITE_ISS_CACHE_KEY,
)
from .types import SatelliteOmmRecord

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class HorizonsTargetSpec:
    label: str
    aliases: tuple[str, ...]
    spkid: str
    pdes: str = ""
    target_name: str = ""


CELESTRAK_GP_JSON_URL = "https://celestrak.org/NORAD/elements/gp.php"
WHERETHEISS_API_URL = "https://api.wheretheiss.at/v1"
HORIZONS_LOOKUP_API_URL = "https://ssd.jpl.nasa.gov/api/horizons_lookup.api"
HORIZONS_API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
CELESTRAK_GROUP_BY_KEY = {
    SATELLITE_ISS_CACHE_KEY: "stations",
}
HORIZONS_TARGETS_BY_KEY = {
    SATELLITE_HORIZONS_CACHE_KEY: (
        HorizonsTargetSpec(
            label="JWST",
            aliases=("JWST", "James Webb Space Telescope", "James Webb"),
            spkid="-170",
            pdes="2021-130A",
            target_name="James Webb Space Telescope (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="Voyager 1",
            aliases=("Voyager 1", "Voyager-1"),
            spkid="-31",
            pdes="1977-084A",
            target_name="Voyager 1 (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="Voyager 2",
            aliases=("Voyager 2", "Voyager-2"),
            spkid="-32",
            pdes="1977-076A",
            target_name="Voyager 2 (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="Parker",
            aliases=("Parker Solar Probe", "Parker", "Solar Probe Plus"),
            spkid="-96",
            pdes="2018-065A",
            target_name="Parker Solar Probe (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="Europa Clipper",
            aliases=("Europa Clipper",),
            spkid="-159",
            pdes="2024-182A",
            target_name="Europa Clipper (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="Lucy",
            aliases=("Lucy",),
            spkid="-49",
            pdes="2021-093A",
            target_name="Lucy (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="Psyche",
            aliases=("Psyche", "Psyche spacecraft"),
            spkid="-255",
            pdes="2023-157A",
            target_name="Psyche (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="JUICE",
            aliases=("JUICE",),
            spkid="-28",
            pdes="2023-053A",
            target_name="JUICE (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="Solar Orbiter",
            aliases=("Solar Orbiter", "Solo"),
            spkid="-144",
            pdes="2020-010A",
            target_name="Solar Orbiter (spacecraft)",
        ),
        HorizonsTargetSpec(
            label="BepiColombo",
            aliases=("BepiColombo", "MPO MMO"),
            spkid="-121",
            pdes="2018-080A",
            target_name="BepiColombo (spacecraft)",
        ),
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


def build_horizons_vector_url(
    command: str,
    *,
    target_time_utc: datetime,
    base_url: str = HORIZONS_API_URL,
) -> str:
    tlist = target_time_utc.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    params = {
        "format": "json",
        "COMMAND": f"'{command}'",
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "VECTORS",
        "CENTER": "'500@399'",
        "REF_SYSTEM": "ICRF",
        "REF_PLANE": "FRAME",
        "OUT_UNITS": "KM-S",
        "VEC_TABLE": "'2'",
        "TIME_TYPE": "UT",
        "TIME_DIGITS": "SECONDS",
        "EXTRA_PREC": "YES",
        "CSV_FORMAT": "YES",
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
        headers={
            "Accept": "application/json",
            "User-Agent": build_user_agent("satellites-celestrak"),
        },
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
        headers={
            "Accept": "application/json",
            "User-Agent": build_user_agent("satellites-horizons"),
        },
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
        headers={
            "Accept": "application/json",
            "User-Agent": build_user_agent("satellites-horizons"),
        },
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


def fetch_horizons_vector_csv(
    command: str,
    *,
    target_time_utc: datetime,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    base_url: str = HORIZONS_API_URL,
) -> list[list[str]]:
    request = Request(
        build_horizons_vector_url(
            command,
            target_time_utc=target_time_utc,
            base_url=base_url,
        ),
        headers={
            "Accept": "application/json",
            "User-Agent": build_user_agent("satellites-horizons"),
        },
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
        headers={
            "Accept": "application/json",
            "User-Agent": build_user_agent("satellites-wheretheiss"),
        },
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
    request_interval_s: float = 0.0,
    record_callback: Callable[[SatelliteOmmRecord], None] | None = None,
) -> list[SatelliteOmmRecord]:
    del observer_lat, observer_lon, observer_height_m
    if group_key != SATELLITE_HORIZONS_CACHE_KEY:
        raise KeyError(group_key)
    records: list[SatelliteOmmRecord] = []
    last_request_completed_at: float | None = None

    def wait_for_request_slot() -> None:
        nonlocal last_request_completed_at
        interval_s = max(0.0, float(request_interval_s))
        if last_request_completed_at is not None and interval_s > 0.0:
            remaining_s = interval_s - (monotonic() - last_request_completed_at)
            if remaining_s > 0.0:
                sleep(remaining_s)

    def mark_request_completed() -> None:
        nonlocal last_request_completed_at
        last_request_completed_at = monotonic()

    for target_spec in HORIZONS_TARGETS_BY_KEY[SATELLITE_HORIZONS_CACHE_KEY]:
        resolved = _resolve_horizons_target(
            target_spec,
            timeout_s=timeout_s,
            lookup_base_url=lookup_base_url,
            wait_for_request_slot=wait_for_request_slot,
            mark_request_completed=mark_request_completed,
        )
        if resolved is None:
            continue
        command = str(resolved.get("spkid", "")).strip()
        if not command:
            command = target_spec.spkid
        if not command:
            continue
        try:
            wait_for_request_slot()
            rows = fetch_horizons_vector_csv(
                command,
                target_time_utc=target_time_utc,
                timeout_s=timeout_s,
                base_url=horizons_base_url,
            )
            mark_request_completed()
        except Exception as exc:
            mark_request_completed()
            _log_fetch_attempt_failure(f"Horizons {target_spec.label}", exc)
            continue
        parsed_state_vector: tuple[float, float, float, float, float, float] | None = None
        for row in rows:
            parsed_state_vector = _extract_state_vector_from_csv_row(row)
            if parsed_state_vector is not None:
                break
        if parsed_state_vector is None:
            continue
        record: SatelliteOmmRecord = {
            "OBJECT_NAME": target_spec.label,
            "HORIZONS_TARGET_NAME": str(resolved.get("name", target_spec.label)).strip() or target_spec.label,
            "HORIZONS_SPKID": command,
            "HORIZONS_PDES": target_spec.pdes,
            "EPOCH": target_time_utc.astimezone(timezone.utc).isoformat(),
            "HORIZONS_X_KM": float(parsed_state_vector[0]),
            "HORIZONS_Y_KM": float(parsed_state_vector[1]),
            "HORIZONS_Z_KM": float(parsed_state_vector[2]),
            "HORIZONS_VX_KM_S": float(parsed_state_vector[3]),
            "HORIZONS_VY_KM_S": float(parsed_state_vector[4]),
            "HORIZONS_VZ_KM_S": float(parsed_state_vector[5]),
            "HORIZONS_CENTER": "500@399",
            "HORIZONS_REF_SYSTEM": "ICRF",
            _SOURCE_KEY: "horizons",
        }
        records.append(record)
        if record_callback is not None:
            record_callback(dict(record))
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
    target_spec: HorizonsTargetSpec,
    *,
    timeout_s: float,
    lookup_base_url: str,
    wait_for_request_slot: Callable[[], None] | None = None,
    mark_request_completed: Callable[[], None] | None = None,
) -> dict[str, object] | None:
    for alias in target_spec.aliases:
        try:
            if wait_for_request_slot is not None:
                wait_for_request_slot()
            payload = fetch_horizons_lookup(alias, timeout_s=timeout_s, base_url=lookup_base_url)
            if mark_request_completed is not None:
                mark_request_completed()
        except Exception as exc:
            if mark_request_completed is not None:
                mark_request_completed()
            _log_fetch_attempt_failure(f"Horizons lookup {target_spec.label}", exc)
            continue
        if not isinstance(payload, dict):
            continue
        result = payload.get("result")
        if not isinstance(result, list) or not result:
            continue
        exact = _pick_expected_horizons_spacecraft(target_spec, result)
        if exact is not None:
            return exact
        exact = _pick_exact_horizons_match(alias, result)
        if exact is not None:
            return exact
        for item in result:
            if isinstance(item, dict) and _is_horizons_spacecraft(item):
                return item
    return None


def _pick_expected_horizons_spacecraft(
    target_spec: HorizonsTargetSpec,
    results: list[object],
) -> dict[str, object] | None:
    expected_name = _normalize_horizons_name(target_spec.target_name)
    for item in results:
        if not isinstance(item, dict):
            continue
        if not _is_horizons_spacecraft(item):
            continue
        spkid = str(item.get("spkid", "")).strip()
        pdes = str(item.get("pdes", "")).strip()
        name = _normalize_horizons_name(str(item.get("name", "")))
        if target_spec.spkid and spkid == target_spec.spkid:
            return item
        if target_spec.pdes and pdes == target_spec.pdes:
            return item
        if expected_name and name == expected_name:
            return item
    return None


def _is_horizons_spacecraft(item: dict[str, object]) -> bool:
    return str(item.get("type", "")).strip().casefold() == "spacecraft"


def _pick_exact_horizons_match(search_text: str, results: list[object]) -> dict[str, object] | None:
    query = _normalize_horizons_name(search_text)
    if not query:
        return None
    for item in results:
        if not isinstance(item, dict):
            continue
        if not _is_horizons_spacecraft(item):
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


def _extract_state_vector_from_csv_row(
    row: list[str],
) -> tuple[float, float, float, float, float, float] | None:
    numeric_values: list[float] = []
    for value in row:
        parsed = _parse_float(value)
        if parsed is not None:
            numeric_values.append(parsed)
    if len(numeric_values) >= 7 and _looks_like_julian_date(numeric_values[0]):
        numeric_values = numeric_values[1:]
    if len(numeric_values) < 6:
        return None
    return (
        numeric_values[0],
        numeric_values[1],
        numeric_values[2],
        numeric_values[3],
        numeric_values[4],
        numeric_values[5],
    )


def _parse_float(value: object) -> float | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def _looks_like_julian_date(value: float) -> bool:
    return 2_000_000.0 <= float(value) <= 3_000_000.0


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
