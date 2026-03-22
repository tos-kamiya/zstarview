from __future__ import annotations

import json
import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable, Iterable

from ..overlay_time import TimeMode, current_utc_time
from ..paths import SATELLITE_CACHE_ROOT_DIR
from ..satellite_constants import (
    SATELLITE_ARCHIVE_RETENTION_SECONDS,
    SATELLITE_ELEMENT_VALID_SECONDS,
)
from .fetch import (
    extract_element_epoch_utc,
    fetch_celestrak_group_by_key,
    filter_records_for_group,
)
from .types import CachedSatelliteElementSet, SatelliteOmmRecord

logger = logging.getLogger(__name__)

SatelliteFetcher = Callable[..., list[SatelliteOmmRecord]]


def satellite_group_cache_path(
    group_key: str,
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / f"{group_key}.json"


def satellite_archive_dir(
    group_key: str,
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / "archive" / str(group_key)


def satellite_archive_cache_path(
    group_key: str,
    element_epoch_utc: datetime,
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
) -> Path:
    normalized = _normalize_utc(element_epoch_utc)
    filename = normalized.strftime("%Y-%m-%dT%H-%M-%SZ.json")
    return satellite_archive_dir(group_key, cache_root=cache_root) / filename


def load_satellite_cache(
    group_key: str,
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
) -> CachedSatelliteElementSet | None:
    path = satellite_group_cache_path(group_key, cache_root=cache_root)
    return _load_cached_set_from_path(path, group_key=group_key)


def iter_satellite_archive(
    group_key: str,
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
) -> list[CachedSatelliteElementSet]:
    archive_dir = satellite_archive_dir(group_key, cache_root=cache_root)
    if not archive_dir.is_dir():
        return []
    loaded: list[CachedSatelliteElementSet] = []
    for path in sorted(archive_dir.glob("*.json")):
        cached = _load_cached_set_from_path(path, group_key=group_key)
        if cached is not None:
            loaded.append(cached)
    return loaded


def save_satellite_cache(
    group_key: str,
    records: list[SatelliteOmmRecord],
    *,
    element_epoch_utc: datetime,
    fetched_at_utc: datetime,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    source: str = "celestrak",
) -> Path:
    path = satellite_group_cache_path(group_key, cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "group_key": str(group_key),
        "element_epoch_utc": _normalize_utc(element_epoch_utc).isoformat(),
        "fetched_at_utc": _normalize_utc(fetched_at_utc).isoformat(),
        "source": str(source),
        "records": filter_records_for_group(group_key, records),
    }
    path.write_text(json.dumps(payload, separators=(",", ":"), sort_keys=True), encoding="utf-8")
    return path


def archive_current_satellite_cache(
    group_key: str,
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
) -> Path | None:
    current_path = satellite_group_cache_path(group_key, cache_root=cache_root)
    if not current_path.is_file():
        return None
    try:
        payload = json.loads(current_path.read_text(encoding="utf-8"))
        cached = _cached_set_from_payload(payload, group_key=group_key)
    except Exception:
        logger.warning("Removing unreadable satellite current cache: %s", current_path, exc_info=True)
        current_path.unlink(missing_ok=True)
        return None

    archive_path = satellite_archive_cache_path(
        group_key,
        cached.element_epoch_utc,
        cache_root=cache_root,
    )
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    if archive_path.exists():
        current_path.unlink(missing_ok=True)
        return archive_path
    current_path.replace(archive_path)
    return archive_path


def cleanup_satellite_cache(
    *,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    now_utc: datetime | None = None,
    max_age_seconds: int = SATELLITE_ARCHIVE_RETENTION_SECONDS,
) -> None:
    root = Path(cache_root) / "archive"
    if not root.is_dir():
        return
    now = _normalize_utc(now_utc or current_utc_time())
    cutoff = now - timedelta(seconds=max(0, int(max_age_seconds)))
    for path in root.glob("*/*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            cached = _cached_set_from_payload(payload)
        except Exception:
            logger.warning("Removing unreadable satellite archive file: %s", path, exc_info=True)
            path.unlink(missing_ok=True)
            continue
        if cached.fetched_at_utc < cutoff:
            path.unlink(missing_ok=True)
    for path in root.glob("*"):
        if path.is_dir():
            try:
                next(path.iterdir())
            except StopIteration:
                path.rmdir()


def fetch_cached_satellite_elements(
    group_key: str,
    *,
    fetcher: SatelliteFetcher = fetch_celestrak_group_by_key,
    timeout_s: float = 20.0,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    fresh_ttl_seconds: int = SATELLITE_ELEMENT_VALID_SECONDS,
    now_utc: datetime | None = None,
) -> CachedSatelliteElementSet:
    now = _normalize_utc(now_utc or current_utc_time())
    cached = load_satellite_cache(group_key, cache_root=cache_root)
    if cached is not None and _within_validity(now, cached.element_epoch_utc, fresh_ttl_seconds):
        return CachedSatelliteElementSet(
            group_key=group_key,
            element_epoch_utc=cached.element_epoch_utc,
            fetched_at_utc=cached.fetched_at_utc,
            source="cache-fresh",
            records=filter_records_for_group(group_key, cached.records),
        )
    records = fetcher(group_key, timeout_s=timeout_s)
    element_epoch_utc = extract_element_epoch_utc(records) or now
    fetched_at_utc = now
    archive_current_satellite_cache(group_key, cache_root=cache_root)
    save_satellite_cache(
        group_key,
        records,
        element_epoch_utc=element_epoch_utc,
        fetched_at_utc=fetched_at_utc,
        cache_root=cache_root,
        source="celestrak",
    )
    cleanup_satellite_cache(cache_root=cache_root, now_utc=now)
    return CachedSatelliteElementSet(
        group_key=group_key,
        element_epoch_utc=element_epoch_utc,
        fetched_at_utc=fetched_at_utc,
        source="celestrak",
        records=filter_records_for_group(group_key, records),
    )


def resolve_satellite_elements_for_time(
    group_key: str,
    *,
    target_time_utc: datetime,
    time_mode: TimeMode,
    fetcher: SatelliteFetcher = fetch_celestrak_group_by_key,
    timeout_s: float = 20.0,
    cache_root: str | Path = SATELLITE_CACHE_ROOT_DIR,
    validity_seconds: int = SATELLITE_ELEMENT_VALID_SECONDS,
    now_utc: datetime | None = None,
) -> CachedSatelliteElementSet:
    if time_mode == "future":
        raise RuntimeError("Satellites: future view is not supported")
    if time_mode == "present":
        return fetch_cached_satellite_elements(
            group_key,
            fetcher=fetcher,
            timeout_s=timeout_s,
            cache_root=cache_root,
            fresh_ttl_seconds=validity_seconds,
            now_utc=now_utc,
        )

    target = _normalize_utc(target_time_utc)
    candidates: list[CachedSatelliteElementSet] = []
    current = load_satellite_cache(group_key, cache_root=cache_root)
    if current is not None:
        candidates.append(current)
    candidates.extend(iter_satellite_archive(group_key, cache_root=cache_root))
    selected = _select_best_cached_set(
        candidates,
        target_time_utc=target,
        validity_seconds=validity_seconds,
    )
    if selected is None:
        raise RuntimeError("Satellites: no cached orbital elements within validity window")
    return CachedSatelliteElementSet(
        group_key=selected.group_key,
        element_epoch_utc=selected.element_epoch_utc,
        fetched_at_utc=selected.fetched_at_utc,
        source="cache-archive" if selected is not current else "cache-current",
        records=filter_records_for_group(group_key, selected.records),
    )


def _load_cached_set_from_path(path: Path, *, group_key: str | None = None) -> CachedSatelliteElementSet | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return _cached_set_from_payload(payload, group_key=group_key)
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
    return CachedSatelliteElementSet(
        group_key=effective_group_key,
        element_epoch_utc=element_epoch_utc,
        fetched_at_utc=fetched_at_utc,
        source=str(payload.get("source", "cache")),
        records=filter_records_for_group(effective_group_key, records),
    )


def _select_best_cached_set(
    candidates: Iterable[CachedSatelliteElementSet],
    *,
    target_time_utc: datetime,
    validity_seconds: int,
) -> CachedSatelliteElementSet | None:
    target = _normalize_utc(target_time_utc)
    best: tuple[float, float, CachedSatelliteElementSet] | None = None
    for candidate in candidates:
        age_seconds = abs((candidate.element_epoch_utc - target).total_seconds())
        if age_seconds > max(0, int(validity_seconds)):
            continue
        score = (age_seconds, -candidate.element_epoch_utc.timestamp())
        if best is None or score < best[:2]:
            best = (score[0], score[1], candidate)
    return None if best is None else best[2]


def _within_validity(now_utc: datetime, element_epoch_utc: datetime, validity_seconds: int) -> bool:
    return abs((now_utc - element_epoch_utc).total_seconds()) <= max(0, int(validity_seconds))


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)
