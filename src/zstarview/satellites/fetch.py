from __future__ import annotations

import json
from collections.abc import Iterable
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from skyfield.api import EarthSatellite
import skyfield.api

from .types import SatelliteOmmRecord

CELESTRAK_GP_JSON_URL = "https://celestrak.org/NORAD/elements/gp.php"
CELESTRAK_GROUP_BY_KEY = {
    "station": "stations",
    "starlink": "starlink",
}
_ISS_NORAD_CAT_ID = "25544"
_CSS_TIANHE_NORAD_CAT_ID = "48274"


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
    timeout_s: float = 20.0,
    base_url: str = CELESTRAK_GP_JSON_URL,
) -> list[SatelliteOmmRecord]:
    request = Request(
        build_celestrak_group_url(group_name, base_url=base_url),
        headers={"Accept": "application/json"},
    )
    with urlopen(request, timeout=float(timeout_s)) as response:
        payload = json.load(response)
    return normalize_celestrak_omm_payload(payload)


def fetch_celestrak_group_by_key(
    group_key: str,
    *,
    timeout_s: float = 20.0,
    base_url: str = CELESTRAK_GP_JSON_URL,
) -> list[SatelliteOmmRecord]:
    group_name = CELESTRAK_GROUP_BY_KEY[group_key]
    records = fetch_celestrak_group_omm(group_name, timeout_s=timeout_s, base_url=base_url)
    return filter_records_for_group(group_key, records)


def normalize_celestrak_omm_payload(payload: object) -> list[SatelliteOmmRecord]:
    if not isinstance(payload, list):
        return []
    records: list[SatelliteOmmRecord] = []
    for row in payload:
        if isinstance(row, dict):
            records.append(dict(row))
    return records


def filter_records_for_group(
    group_key: str,
    records: Iterable[SatelliteOmmRecord],
) -> list[SatelliteOmmRecord]:
    normalized = [dict(record) for record in records]
    if group_key != "station":
        return normalized
    allowed_cat_ids = {_ISS_NORAD_CAT_ID, _CSS_TIANHE_NORAD_CAT_ID}
    return [
        record
        for record in normalized
        if str(record.get("NORAD_CAT_ID", "")).strip() in allowed_cat_ids
    ]


def build_earth_satellites(
    records: Iterable[SatelliteOmmRecord],
    *,
    ts: object | None = None,
) -> list[EarthSatellite]:
    timescale = ts or skyfield.api.load.timescale()
    satellites: list[EarthSatellite] = []
    for record in records:
        satellites.append(EarthSatellite.from_omm(timescale, record))
    return satellites
