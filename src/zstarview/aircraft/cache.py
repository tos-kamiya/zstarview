from __future__ import annotations

import json
import logging
from collections.abc import Callable
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from uuid import uuid4

from filelock import FileLock, Timeout

from ..aircraft_constants import (
    AIRCRAFT_CACHE_RETENTION_SECONDS,
    AIRCRAFT_CACHE_STALE_FALLBACK_SECONDS,
    AIRCRAFT_REFRESH_INTERVAL_SECONDS,
)
from ..paths import AIRCRAFT_CACHE_ROOT_DIR
from .opensky import AircraftBoundingBox, fetch_opensky_states
from .types import AircraftSnapshot

logger = logging.getLogger(__name__)

AircraftFetcher = Callable[..., list[AircraftSnapshot]]
RATE_LIMIT_FILE_NAME = "rate_limit.json"
FETCH_LOCK_FILE_NAME = "fetch.lock"
AIRCRAFT_HISTORY_DIR_NAME = "history"


@dataclass(frozen=True)
class CachedAircraftSnapshotSet:
    snapshots: list[AircraftSnapshot]
    bbox: AircraftBoundingBox
    fetched_at_utc: datetime
    source: str
    is_stale: bool = False


def bbox_cache_key(bbox: AircraftBoundingBox) -> str:
    return (
        f"{bbox.min_lat:+08.4f}_{bbox.max_lat:+08.4f}_"
        f"{bbox.min_lon:+09.4f}_{bbox.max_lon:+09.4f}"
    ).replace("+", "p").replace("-", "m")


def aircraft_cache_path(
    bbox: AircraftBoundingBox,
    *,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / f"{bbox_cache_key(bbox)}.json"


def _aircraft_history_path(
    bbox: AircraftBoundingBox,
    fetched_at_utc: datetime,
    *,
    cache_root: str | Path,
) -> Path:
    fetched_at = _normalize_utc(fetched_at_utc)
    timestamp_dir = fetched_at.strftime("%Y/%m/%d/%H/%M")
    filename = f"{bbox_cache_key(bbox)}-{fetched_at.strftime('%S%f')}-{uuid4().hex}.json"
    return Path(cache_root) / AIRCRAFT_HISTORY_DIR_NAME / timestamp_dir / filename


def aircraft_rate_limit_path(
    *,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / RATE_LIMIT_FILE_NAME


def aircraft_fetch_lock_path(
    *,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / FETCH_LOCK_FILE_NAME


def load_aircraft_cache(
    bbox: AircraftBoundingBox,
    *,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
) -> CachedAircraftSnapshotSet | None:
    root = Path(cache_root)
    key = bbox_cache_key(bbox)
    candidates = [aircraft_cache_path(bbox, cache_root=root)]
    history_root = root / AIRCRAFT_HISTORY_DIR_NAME
    if history_root.is_dir():
        candidates.extend(history_root.glob(f"*/*/*/*/*/{key}-*.json"))

    cached_sets: list[CachedAircraftSnapshotSet] = []
    for path in candidates:
        if not path.is_file():
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            cached = _cached_set_from_payload(payload)
        except Exception:
            logger.warning("Failed to read aircraft cache: %s", path, exc_info=True)
            continue
        if cached.bbox == bbox:
            cached_sets.append(cached)
    if not cached_sets:
        return None
    return max(cached_sets, key=lambda cached: cached.fetched_at_utc)


def save_aircraft_cache(
    bbox: AircraftBoundingBox,
    snapshots: list[AircraftSnapshot],
    *,
    fetched_at_utc: datetime,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
    source: str = "opensky",
) -> Path:
    path = aircraft_cache_path(bbox, cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    _archive_aircraft_latest(path, bbox=bbox, cache_root=cache_root)
    payload = {
        "bbox": asdict(bbox),
        "fetched_at_utc": _normalize_utc(fetched_at_utc).isoformat(),
        "source": str(source),
        "snapshots": [asdict(snapshot) for snapshot in snapshots],
    }
    payload_text = json.dumps(payload, separators=(",", ":"), sort_keys=True)
    _write_aircraft_history(
        payload_text,
        bbox=bbox,
        fetched_at_utc=fetched_at_utc,
        cache_root=cache_root,
    )
    path.write_text(payload_text, encoding="utf-8")
    return path


def load_aircraft_rate_limit(
    *,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
) -> datetime | None:
    path = aircraft_rate_limit_path(cache_root=cache_root)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("rate limit payload must be a dict")
        raw = payload.get("last_successful_fetch_at_utc")
        if not isinstance(raw, str):
            return None
        return _normalize_utc(datetime.fromisoformat(raw))
    except Exception:
        logger.warning("Failed to read aircraft rate-limit metadata: %s", path, exc_info=True)
        return None


def save_aircraft_rate_limit(
    fetched_at_utc: datetime,
    *,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
) -> Path:
    path = aircraft_rate_limit_path(cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "last_successful_fetch_at_utc": _normalize_utc(fetched_at_utc).isoformat(),
    }
    path.write_text(json.dumps(payload, separators=(",", ":"), sort_keys=True), encoding="utf-8")
    return path


def cleanup_aircraft_cache(
    *,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
    now_utc: datetime | None = None,
    max_age_seconds: int = AIRCRAFT_CACHE_RETENTION_SECONDS,
) -> None:
    root = Path(cache_root)
    if not root.is_dir():
        return
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    cutoff = now - timedelta(seconds=max(0, int(max_age_seconds)))
    candidates = list(root.glob("*.json"))
    history_root = root / AIRCRAFT_HISTORY_DIR_NAME
    if history_root.is_dir():
        candidates.extend(history_root.rglob("*.json"))
    for path in candidates:
        if path.name == RATE_LIMIT_FILE_NAME:
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            cached = _cached_set_from_payload(payload)
        except Exception:
            logger.warning("Removing unreadable aircraft cache file: %s", path, exc_info=True)
            _remove_aircraft_cache_file(path)
            continue
        if cached.fetched_at_utc < cutoff:
            _remove_aircraft_cache_file(path)
    if history_root.is_dir():
        for directory in sorted(history_root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
            if directory.is_dir():
                try:
                    directory.rmdir()
                except OSError:
                    pass
        try:
            history_root.rmdir()
        except OSError:
            pass


def fetch_cached_opensky_states(
    bbox: AircraftBoundingBox,
    *,
    fetcher: AircraftFetcher = fetch_opensky_states,
    timeout_s: float = 20.0,
    cache_root: str | Path = AIRCRAFT_CACHE_ROOT_DIR,
    fresh_ttl_seconds: int = AIRCRAFT_REFRESH_INTERVAL_SECONDS,
    stale_fallback_seconds: int = AIRCRAFT_CACHE_STALE_FALLBACK_SECONDS,
    enforce_global_rate_limit: bool = True,
    lock_timeout_s: float = 0.0,
    now_utc: datetime | None = None,
) -> CachedAircraftSnapshotSet:
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    cached = load_aircraft_cache(bbox, cache_root=cache_root)
    if cached is not None:
        age_seconds = (now - cached.fetched_at_utc).total_seconds()
        if age_seconds <= max(0, int(fresh_ttl_seconds)):
            return CachedAircraftSnapshotSet(
                snapshots=cached.snapshots,
                bbox=cached.bbox,
                fetched_at_utc=cached.fetched_at_utc,
                source="cache-fresh",
                is_stale=False,
            )
    lock_path = aircraft_fetch_lock_path(cache_root=cache_root)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        lock = FileLock(lock_path, timeout=float(lock_timeout_s))
        lock.acquire()
    except Timeout:
        if enforce_global_rate_limit:
            return CachedAircraftSnapshotSet(
                snapshots=[],
                bbox=bbox,
                fetched_at_utc=now,
                source="rate-limited-skip",
                is_stale=False,
            )
        raise
    try:
        cached = load_aircraft_cache(bbox, cache_root=cache_root)
        if cached is not None:
            age_seconds = (now - cached.fetched_at_utc).total_seconds()
            if age_seconds <= max(0, int(fresh_ttl_seconds)):
                return CachedAircraftSnapshotSet(
                    snapshots=cached.snapshots,
                    bbox=cached.bbox,
                    fetched_at_utc=cached.fetched_at_utc,
                    source="cache-fresh",
                    is_stale=False,
                )
        if enforce_global_rate_limit:
            last_success_utc = load_aircraft_rate_limit(cache_root=cache_root)
            if last_success_utc is not None:
                age_seconds = (now - last_success_utc).total_seconds()
                if age_seconds <= max(0, int(fresh_ttl_seconds)):
                    return CachedAircraftSnapshotSet(
                        snapshots=[],
                        bbox=bbox,
                        fetched_at_utc=last_success_utc,
                        source="rate-limited-skip",
                        is_stale=False,
                    )
        try:
            snapshots = fetcher(bbox, timeout_s=timeout_s)
        except Exception:
            if cached is not None:
                age_seconds = (now - cached.fetched_at_utc).total_seconds()
                if age_seconds <= max(0, int(stale_fallback_seconds)):
                    return CachedAircraftSnapshotSet(
                        snapshots=cached.snapshots,
                        bbox=cached.bbox,
                        fetched_at_utc=cached.fetched_at_utc,
                        source="cache-stale",
                        is_stale=True,
                    )
            raise
        fetched_at_utc = now
        save_aircraft_cache(
            bbox,
            snapshots,
            fetched_at_utc=fetched_at_utc,
            cache_root=cache_root,
            source="opensky",
        )
        save_aircraft_rate_limit(fetched_at_utc, cache_root=cache_root)
        cleanup_aircraft_cache(
            cache_root=cache_root,
            now_utc=now,
            max_age_seconds=AIRCRAFT_CACHE_RETENTION_SECONDS,
        )
        return CachedAircraftSnapshotSet(
            snapshots=snapshots,
            bbox=bbox,
            fetched_at_utc=fetched_at_utc,
            source="opensky",
            is_stale=False,
        )
    finally:
        lock.release()


def _cached_set_from_payload(payload: object) -> CachedAircraftSnapshotSet:
    if not isinstance(payload, dict):
        raise ValueError("cache payload must be a dict")
    bbox_payload = payload.get("bbox")
    fetched_at_raw = payload.get("fetched_at_utc")
    snapshots_payload = payload.get("snapshots")
    if not isinstance(bbox_payload, dict):
        raise ValueError("cache payload missing bbox")
    if not isinstance(fetched_at_raw, str):
        raise ValueError("cache payload missing fetched_at_utc")
    if not isinstance(snapshots_payload, list):
        raise ValueError("cache payload missing snapshots")
    bbox = AircraftBoundingBox(
        min_lat=float(bbox_payload["min_lat"]),
        max_lat=float(bbox_payload["max_lat"]),
        min_lon=float(bbox_payload["min_lon"]),
        max_lon=float(bbox_payload["max_lon"]),
    )
    snapshots = [
        AircraftSnapshot(
            icao24=str(item["icao24"]),
            callsign=None if item.get("callsign") is None else str(item.get("callsign")),
            latitude=float(item["latitude"]),
            longitude=float(item["longitude"]),
            baro_altitude_m=_optional_float(item.get("baro_altitude_m")),
            velocity_mps=_optional_float(item.get("velocity_mps")),
            heading_deg=_optional_float(item.get("heading_deg")),
            vertical_rate_mps=_optional_float(item.get("vertical_rate_mps")),
            on_ground=bool(item.get("on_ground", False)),
            last_contact_unix=_optional_int(item.get("last_contact_unix")),
        )
        for item in snapshots_payload
        if isinstance(item, dict)
    ]
    fetched_at_utc = _normalize_utc(datetime.fromisoformat(fetched_at_raw))
    source = str(payload.get("source", "cache"))
    return CachedAircraftSnapshotSet(
        snapshots=snapshots,
        bbox=bbox,
        fetched_at_utc=fetched_at_utc,
        source=source,
        is_stale=False,
    )


def _archive_aircraft_latest(
    latest_path: Path,
    *,
    bbox: AircraftBoundingBox,
    cache_root: str | Path,
) -> Path | None:
    if not latest_path.is_file():
        return None
    try:
        payload = json.loads(latest_path.read_text(encoding="utf-8"))
        cached = _cached_set_from_payload(payload)
    except Exception:
        logger.warning("Skipping invalid aircraft cache history source: %s", latest_path, exc_info=True)
        return None
    if cached.bbox != bbox:
        return None
    return _write_aircraft_history(
        latest_path.read_text(encoding="utf-8"),
        bbox=bbox,
        fetched_at_utc=cached.fetched_at_utc,
        cache_root=cache_root,
    )


def _write_aircraft_history(
    payload_text: str,
    *,
    bbox: AircraftBoundingBox,
    fetched_at_utc: datetime,
    cache_root: str | Path,
) -> Path | None:
    fetched_at = _normalize_utc(fetched_at_utc)
    history_root = Path(cache_root) / AIRCRAFT_HISTORY_DIR_NAME
    key = bbox_cache_key(bbox)
    for existing_path in history_root.glob(f"*/*/*/*/*/{key}-*.json"):
        try:
            existing_payload = json.loads(existing_path.read_text(encoding="utf-8"))
            existing = _cached_set_from_payload(existing_payload)
        except Exception:
            continue
        if existing.bbox == bbox and existing.fetched_at_utc == fetched_at:
            if existing_path.read_text(encoding="utf-8") == payload_text:
                return existing_path

    history_path = _aircraft_history_path(bbox, fetched_at, cache_root=cache_root)
    temporary_path = history_path.with_name(f".{history_path.name}.{uuid4().hex}.tmp")
    try:
        history_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path.write_text(payload_text, encoding="utf-8")
        temporary_path.replace(history_path)
    except OSError:
        logger.warning("Failed to write aircraft cache history: %s", history_path, exc_info=True)
        _remove_aircraft_cache_file(temporary_path)
        return None
    return history_path


def _remove_aircraft_cache_file(path: Path) -> None:
    try:
        path.unlink(missing_ok=True)
    except OSError:
        logger.warning("Failed to remove aircraft cache file: %s", path, exc_info=True)


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _optional_float(value: object) -> float | None:
    if value is None:
        return None
    return float(value)


def _optional_int(value: object) -> int | None:
    if value is None:
        return None
    return int(value)
