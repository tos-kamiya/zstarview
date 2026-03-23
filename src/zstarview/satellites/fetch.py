from __future__ import annotations

from datetime import datetime, timezone
import json
from collections.abc import Iterable
import logging
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import URLError

from skyfield.api import EarthSatellite
import skyfield.api

from ..satellite_constants import SATELLITE_FETCH_TIMEOUT_SECONDS, SATELLITE_ISS_CACHE_KEY
from .types import SatelliteOmmRecord

logger = logging.getLogger(__name__)

CELESTRAK_GP_JSON_URL = "https://celestrak.org/NORAD/elements/gp.php"
WHERETHEISS_API_URL = "https://api.wheretheiss.at/v1"
CELESTRAK_GROUP_BY_KEY = {
    SATELLITE_ISS_CACHE_KEY: "stations",
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
    wheretheiss_base_url: str = WHERETHEISS_API_URL,
    celestrak_base_url: str = CELESTRAK_GP_JSON_URL,
) -> list[SatelliteOmmRecord]:
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
    if _is_timeout_error(exc):
        logger.warning("ISS fetch failed via %s: %s", source_name, exc)
        return
    logger.info("ISS fetch failed via %s: %s", source_name, exc)


def _is_timeout_error(exc: Exception) -> bool:
    if isinstance(exc, TimeoutError):
        return True
    if isinstance(exc, URLError):
        reason = getattr(exc, "reason", None)
        if isinstance(reason, TimeoutError):
            return True
    return "timed out" in str(exc).lower()
