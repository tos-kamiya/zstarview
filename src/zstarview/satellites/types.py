from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Any

SatelliteOmmRecord = dict[str, Any]


@dataclass(frozen=True)
class CachedSatelliteElementSet:
    group_key: str
    element_epoch_utc: datetime
    fetched_at_utc: datetime
    source: str
    records: list[SatelliteOmmRecord]
    last_fetch_attempt_utc: datetime | None = None
    last_fetch_failed: bool = False
    last_fetch_error: str | None = None
    last_fetch_failure_utc: datetime | None = None
    failure_backoff_until_utc: datetime | None = None


@dataclass(frozen=True)
class SatelliteOverlayPoint:
    group_key: str
    satellite_name: str
    alt_deg: float
    az_deg: float
    marker_scale: float
    show_label: bool


def satellite_altaz_from_record(record: SatelliteOmmRecord) -> tuple[float | None, float | None]:
    alt_keys = ("ALT_DEG", "ALTITUDE_DEG", "ELEV_DEG", "ELEVATION_DEG")
    az_keys = ("AZ_DEG", "AZIMUTH_DEG")

    alt_deg: float | None = None
    az_deg: float | None = None

    for key in alt_keys:
        raw = record.get(key)
        if raw is None:
            continue
        try:
            alt_deg = float(raw)
            break
        except (TypeError, ValueError):
            continue

    for key in az_keys:
        raw = record.get(key)
        if raw is None:
            continue
        try:
            az_deg = float(raw)
            break
        except (TypeError, ValueError):
            continue

    return alt_deg, az_deg
