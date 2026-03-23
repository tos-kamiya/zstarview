from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable

from ..overlay_time import TimeMode, current_utc_time
from ..paths import SATELLITE_CACHE_ROOT_DIR
from ..satellite_constants import (
    SATELLITE_ELEMENT_VALID_SECONDS,
    SATELLITE_FAILURE_RETRY_SECONDS,
    SATELLITE_FETCH_TIMEOUT_SECONDS,
    SATELLITE_GROUP_VALIDITY_SECONDS,
    SATELLITE_ISS_CACHE_KEY,
)
from .fetch import (
    extract_record_source,
    extract_element_epoch_utc,
    fetch_iss_records,
    filter_records_for_group,
)
from .types import CachedSatelliteElementSet, SatelliteOmmRecord

logger = logging.getLogger(__name__)

SatelliteFetcher = Callable[..., list[SatelliteOmmRecord]]


@dataclass(frozen=True)
class SatelliteFetchMetadata:
    last_fetch_attempt_utc: datetime | None = None
    last_fetch_failed: bool = False
    last_fetch_error: str | None = None
    last_fetch_failure_utc: datetime | None = None
    failure_backoff_until_utc: datetime | None = None


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
    return _load_cached_set_from_path(path, group_key=group_key)


def save_satellite_cache(
    group_key: str,
    records: list[SatelliteOmmRecord],
    *,
    element_epoch_utc: datetime,
    fetched_at_utc: datetime,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    source: str = "celestrak",
    last_fetch_attempt_utc: datetime | None = None,
    last_fetch_failed: bool = False,
    last_fetch_error: str | None = None,
    last_fetch_failure_utc: datetime | None = None,
    failure_backoff_until_utc: datetime | None = None,
) -> Path:
    path = satellite_group_cache_path(group_key, cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "group_key": str(group_key),
        "element_epoch_utc": _normalize_utc(element_epoch_utc).isoformat(),
        "fetched_at_utc": _normalize_utc(fetched_at_utc).isoformat(),
        "source": str(source),
        "records": filter_records_for_group(group_key, records),
        "last_fetch_attempt_utc": _serialize_optional_utc(last_fetch_attempt_utc or fetched_at_utc),
        "last_fetch_failed": bool(last_fetch_failed),
        "last_fetch_error": str(last_fetch_error).strip() if last_fetch_error else None,
        "last_fetch_failure_utc": _serialize_optional_utc(last_fetch_failure_utc),
        "failure_backoff_until_utc": _serialize_optional_utc(failure_backoff_until_utc),
    }
    path.write_text(json.dumps(payload, separators=(",", ":"), sort_keys=True), encoding="utf-8")
    return path


def save_satellite_fetch_failure(
    group_key: str,
    *,
    attempted_at_utc: datetime,
    error_text: str,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    backoff_seconds: int = SATELLITE_FAILURE_RETRY_SECONDS,
) -> Path:
    path = satellite_group_cache_path(group_key, cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _load_cache_payload(path) or {}
    attempted_at = _normalize_utc(attempted_at_utc)
    payload["group_key"] = str(group_key)
    payload["last_fetch_attempt_utc"] = attempted_at.isoformat()
    payload["last_fetch_failed"] = True
    payload["last_fetch_error"] = str(error_text).strip()
    payload["last_fetch_failure_utc"] = attempted_at.isoformat()
    payload["failure_backoff_until_utc"] = (
        attempted_at.replace()
        + timedelta(seconds=max(0, int(backoff_seconds)))
    ).isoformat()
    path.write_text(json.dumps(payload, separators=(",", ":"), sort_keys=True), encoding="utf-8")
    return path


def fetch_cached_satellite_elements(
    group_key: str,
    *,
    fetcher: SatelliteFetcher = fetch_iss_records,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    fresh_ttl_seconds: int | None = None,
    now_utc: datetime | None = None,
) -> CachedSatelliteElementSet:
    now = _normalize_utc(now_utc or current_utc_time())
    ttl_seconds = _group_validity_seconds(group_key, fresh_ttl_seconds)
    path = satellite_group_cache_path(group_key, cache_root=cache_root)
    payload = _load_cache_payload(path)
    metadata = _fetch_metadata_from_payload(payload)
    cached = _load_cached_set_from_payload(payload, group_key=group_key)
    if cached is not None and _within_validity(now, cached.element_epoch_utc, ttl_seconds):
        return CachedSatelliteElementSet(
            group_key=group_key,
            element_epoch_utc=cached.element_epoch_utc,
            fetched_at_utc=cached.fetched_at_utc,
            source="cache-fresh",
            records=filter_records_for_group(group_key, cached.records),
            last_fetch_attempt_utc=metadata.last_fetch_attempt_utc,
            last_fetch_failed=metadata.last_fetch_failed,
            last_fetch_error=metadata.last_fetch_error,
            last_fetch_failure_utc=metadata.last_fetch_failure_utc,
            failure_backoff_until_utc=metadata.failure_backoff_until_utc,
        )
    if (
        cached is not None
        and metadata.failure_backoff_until_utc is not None
        and now < metadata.failure_backoff_until_utc
    ):
        return CachedSatelliteElementSet(
            group_key=group_key,
            element_epoch_utc=cached.element_epoch_utc,
            fetched_at_utc=cached.fetched_at_utc,
            source="cache-backoff",
            records=filter_records_for_group(group_key, cached.records),
            last_fetch_attempt_utc=metadata.last_fetch_attempt_utc,
            last_fetch_failed=metadata.last_fetch_failed,
            last_fetch_error=metadata.last_fetch_error,
            last_fetch_failure_utc=metadata.last_fetch_failure_utc,
            failure_backoff_until_utc=metadata.failure_backoff_until_utc,
        )
    try:
        records = fetcher(group_key, timeout_s=timeout_s)
    except Exception as exc:
        save_satellite_fetch_failure(
            group_key,
            attempted_at_utc=now,
            error_text=str(exc),
            cache_root=cache_root,
        )
        raise
    element_epoch_utc = extract_element_epoch_utc(records) or now
    fetched_at_utc = now
    source = extract_record_source(records)
    save_satellite_cache(
        group_key,
        records,
        element_epoch_utc=element_epoch_utc,
        fetched_at_utc=fetched_at_utc,
        cache_root=cache_root,
        source=source,
        last_fetch_attempt_utc=fetched_at_utc,
    )
    return CachedSatelliteElementSet(
        group_key=group_key,
        element_epoch_utc=element_epoch_utc,
        fetched_at_utc=fetched_at_utc,
        source=source,
        records=filter_records_for_group(group_key, records),
        last_fetch_attempt_utc=fetched_at_utc,
    )


def resolve_satellite_elements_for_time(
    group_key: str,
    *,
    target_time_utc: datetime,
    time_mode: TimeMode,
    fetcher: SatelliteFetcher = fetch_iss_records,
    timeout_s: float = SATELLITE_FETCH_TIMEOUT_SECONDS,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    validity_seconds: int | None = None,
    now_utc: datetime | None = None,
) -> CachedSatelliteElementSet:
    del target_time_utc
    if time_mode != "present":
        raise RuntimeError("Satellites: time-shifted view is not supported")
    return fetch_cached_satellite_elements(
        group_key,
        fetcher=fetcher,
        timeout_s=timeout_s,
        cache_root=cache_root,
        fresh_ttl_seconds=_group_validity_seconds(group_key, validity_seconds),
        now_utc=now_utc,
    )


def _load_cached_set_from_path(path: Path, *, group_key: str | None = None) -> CachedSatelliteElementSet | None:
    payload = _load_cache_payload(path)
    return _load_cached_set_from_payload(payload, group_key=group_key)


def _load_cached_set_from_payload(
    payload: object,
    *,
    group_key: str | None = None,
) -> CachedSatelliteElementSet | None:
    if payload is None:
        return None
    try:
        return _cached_set_from_payload(payload, group_key=group_key)
    except Exception:
        return None


def _load_cache_payload(path: Path) -> object | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        logger.warning("Failed to read satellite cache: %s", path, exc_info=True)
        return None


def _cached_set_from_payload(
    payload: object,
    *,
    group_key: str | None = None,
) -> CachedSatelliteElementSet:
    if not isinstance(payload, dict):
        raise ValueError("cache payload must be a dict")
    stored_group_key = str(payload.get("group_key", group_key or ""))
    effective_group_key = str(group_key or stored_group_key)
    records_payload = payload.get("records")
    if not isinstance(records_payload, list):
        raise ValueError("cache payload missing records")
    records = [dict(item) for item in records_payload if isinstance(item, dict)]
    element_epoch_raw = payload.get("element_epoch_utc")
    element_epoch_utc: datetime | None = None
    if isinstance(element_epoch_raw, str) and element_epoch_raw.strip():
        element_epoch_utc = _normalize_utc(datetime.fromisoformat(element_epoch_raw.replace("Z", "+00:00")))
    if element_epoch_utc is None:
        legacy_element_time_raw = payload.get("element_time_utc")
        if isinstance(legacy_element_time_raw, str) and legacy_element_time_raw.strip():
            element_epoch_utc = _normalize_utc(datetime.fromisoformat(legacy_element_time_raw.replace("Z", "+00:00")))
    if element_epoch_utc is None:
        element_epoch_utc = extract_element_epoch_utc(records)
    if element_epoch_utc is None:
        raise ValueError("cache payload missing element_epoch_utc")
    fetched_at_raw = payload.get("fetched_at_utc")
    fetched_at_utc: datetime | None = None
    if isinstance(fetched_at_raw, str) and fetched_at_raw.strip():
        fetched_at_utc = _normalize_utc(datetime.fromisoformat(fetched_at_raw.replace("Z", "+00:00")))
    if fetched_at_utc is None:
        fetched_at_utc = element_epoch_utc
    metadata = _fetch_metadata_from_payload(payload)
    return CachedSatelliteElementSet(
        group_key=effective_group_key,
        element_epoch_utc=element_epoch_utc,
        fetched_at_utc=fetched_at_utc,
        source=str(payload.get("source", "cache")),
        records=filter_records_for_group(effective_group_key, records),
        last_fetch_attempt_utc=metadata.last_fetch_attempt_utc,
        last_fetch_failed=metadata.last_fetch_failed,
        last_fetch_error=metadata.last_fetch_error,
        last_fetch_failure_utc=metadata.last_fetch_failure_utc,
        failure_backoff_until_utc=metadata.failure_backoff_until_utc,
    )


def _fetch_metadata_from_payload(payload: object) -> SatelliteFetchMetadata:
    if not isinstance(payload, dict):
        return SatelliteFetchMetadata()
    return SatelliteFetchMetadata(
        last_fetch_attempt_utc=_parse_optional_utc(payload.get("last_fetch_attempt_utc")),
        last_fetch_failed=bool(payload.get("last_fetch_failed", False)),
        last_fetch_error=_parse_optional_text(payload.get("last_fetch_error")),
        last_fetch_failure_utc=_parse_optional_utc(payload.get("last_fetch_failure_utc")),
        failure_backoff_until_utc=_parse_optional_utc(payload.get("failure_backoff_until_utc")),
    )


def _within_validity(
    reference_utc: datetime,
    element_epoch_utc: datetime,
    validity_seconds: int,
) -> bool:
    delta = abs((_normalize_utc(reference_utc) - _normalize_utc(element_epoch_utc)).total_seconds())
    return delta <= max(0, int(validity_seconds))


def _group_validity_seconds(group_key: str, override: int | None) -> int:
    if override is not None:
        return int(override)
    return int(SATELLITE_GROUP_VALIDITY_SECONDS.get(group_key, SATELLITE_ELEMENT_VALID_SECONDS))


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _parse_optional_utc(value: object) -> datetime | None:
    if not isinstance(value, str) or not value.strip():
        return None
    return _normalize_utc(datetime.fromisoformat(value.replace("Z", "+00:00")))


def _parse_optional_text(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    text = value.strip()
    return text or None


def _serialize_optional_utc(value: datetime | None) -> str | None:
    if value is None:
        return None
    return _normalize_utc(value).isoformat()
