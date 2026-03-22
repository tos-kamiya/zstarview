from __future__ import annotations

import json
import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable

from ..paths import SATELLITE_CACHE_ROOT_DIR
from ..satellite_constants import (
    SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS,
)
from .fetch import fetch_celestrak_group_by_key, filter_records_for_group
from .types import CachedSatelliteElementSet, SatelliteOmmRecord

logger = logging.getLogger(__name__)

SatelliteFetcher = Callable[..., list[SatelliteOmmRecord]]


def satellite_group_cache_path(
    group_key: str,
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / f"{group_key}.json"


def load_satellite_cache(
    group_key: str,
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
) -> CachedSatelliteElementSet | None:
    path = satellite_group_cache_path(group_key, cache_root=cache_root)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return _cached_set_from_payload(payload)
    except Exception:
        logger.warning("Failed to read satellite cache: %s", path, exc_info=True)
        return None


def save_satellite_cache(
    group_key: str,
    records: list[SatelliteOmmRecord],
    *,
    fetched_at_utc: datetime,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    source: str = "celestrak",
) -> Path:
    path = satellite_group_cache_path(group_key, cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "group_key": str(group_key),
        "fetched_at_utc": _normalize_utc(fetched_at_utc).isoformat(),
        "source": str(source),
        "records": filter_records_for_group(group_key, records),
    }
    path.write_text(json.dumps(payload, separators=(",", ":"), sort_keys=True), encoding="utf-8")
    return path


def cleanup_satellite_cache(
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    now_utc: datetime | None = None,
    max_age_seconds: int = SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS,
) -> None:
    root = Path(cache_root)
    if not root.is_dir():
        return
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    cutoff = now - timedelta(seconds=max(0, int(max_age_seconds)))
    for path in root.glob("*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            cached = _cached_set_from_payload(payload)
        except Exception:
            logger.warning("Removing unreadable satellite cache file: %s", path, exc_info=True)
            path.unlink(missing_ok=True)
            continue
        if cached.fetched_at_utc < cutoff:
            path.unlink(missing_ok=True)


def fetch_cached_satellite_elements(
    group_key: str,
    *,
    fetcher: SatelliteFetcher = fetch_celestrak_group_by_key,
    timeout_s: float = 20.0,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    fresh_ttl_seconds: int = SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS,
    now_utc: datetime | None = None,
) -> CachedSatelliteElementSet:
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    cached = load_satellite_cache(group_key, cache_root=cache_root)
    if cached is not None:
        age_seconds = (now - cached.fetched_at_utc).total_seconds()
        if age_seconds <= max(0, int(fresh_ttl_seconds)):
            return CachedSatelliteElementSet(
                group_key=group_key,
                fetched_at_utc=cached.fetched_at_utc,
                source="cache-fresh",
                records=filter_records_for_group(group_key, cached.records),
            )
    try:
        records = fetcher(group_key, timeout_s=timeout_s)
    except Exception:
        raise
    fetched_at_utc = now
    save_satellite_cache(
        group_key,
        records,
        fetched_at_utc=fetched_at_utc,
        cache_root=cache_root,
        source="celestrak",
    )
    cleanup_satellite_cache(cache_root=cache_root, now_utc=now, max_age_seconds=fresh_ttl_seconds)
    return CachedSatelliteElementSet(
        group_key=group_key,
        fetched_at_utc=fetched_at_utc,
        source="celestrak",
        records=filter_records_for_group(group_key, records),
    )


def _cached_set_from_payload(payload: object) -> CachedSatelliteElementSet:
    if not isinstance(payload, dict):
        raise ValueError("cache payload must be a dict")
    group_key = str(payload["group_key"])
    fetched_at_raw = payload.get("fetched_at_utc")
    records_payload = payload.get("records")
    if not isinstance(fetched_at_raw, str):
        raise ValueError("cache payload missing fetched_at_utc")
    if not isinstance(records_payload, list):
        raise ValueError("cache payload missing records")
    records = [dict(item) for item in records_payload if isinstance(item, dict)]
    return CachedSatelliteElementSet(
        group_key=group_key,
        fetched_at_utc=_normalize_utc(datetime.fromisoformat(fetched_at_raw)),
        source=str(payload.get("source", "cache")),
        records=filter_records_for_group(group_key, records),
    )


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)
