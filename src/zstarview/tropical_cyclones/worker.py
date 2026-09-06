"""Process entry point for tropical-cyclone fetch and normalization jobs."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import traceback
from pathlib import Path
from typing import Any

from .cache import (
    TROPICAL_CYCLONE_CACHE_TTL_SECONDS,
    TROPICAL_CYCLONE_CACHE_VERSION,
    TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS,
    TropicalCycloneCacheEntry,
    is_tropical_cyclone_cache_current,
    is_tropical_cyclone_cache_stale,
    load_tropical_cyclone_cache,
    save_tropical_cyclone_cache,
)
from .client import (
    DEFAULT_SERVICE_URL,
    DEFAULT_TIMEOUT_S,
    DEFAULT_USER_AGENT,
    TropicalCycloneFetchError,
    fetch_active_hurricanes_snapshot,
    fetch_latest_observed_feature,
)
from .models import TropicalCycloneSnapshotCollection

EMPTY_OBSERVED_POSITION_MESSAGE = "No observed position features returned"


def _iso(value: dt.datetime) -> str:
    return value.astimezone(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def _payload(
    collection: TropicalCycloneSnapshotCollection,
    *,
    cached_at_utc: dt.datetime,
    last_checked_utc: dt.datetime,
    next_check_utc: dt.datetime,
    next_refresh_utc: dt.datetime,
    service_url: str,
    banner: str = "",
) -> dict[str, Any]:
    return {
        "snapshot_collection": collection.to_dict(),
        "cached_at_utc": _iso(cached_at_utc),
        "last_checked_utc": _iso(last_checked_utc),
        "next_check_utc": _iso(next_check_utc),
        "next_refresh_utc": _iso(next_refresh_utc),
        "banner": banner,
        "service_url": service_url,
    }


def _save_empty_overlay(*, cache_root: Path, service_url: str, now: dt.datetime) -> dict[str, Any]:
    collection = TropicalCycloneSnapshotCollection(
        snapshots=(),
        source_url=service_url,
        service_name="",
        refreshed_at_utc=now,
    )
    cached_at = dt.datetime.now(dt.timezone.utc)
    entry = TropicalCycloneCacheEntry(
        snapshot_collection=collection,
        cached_at_utc=cached_at,
        cache_version=TROPICAL_CYCLONE_CACHE_VERSION,
    )
    save_tropical_cyclone_cache(entry, cache_root=cache_root)
    return _payload(
        collection,
        cached_at_utc=cached_at,
        last_checked_utc=cached_at,
        next_check_utc=cached_at + dt.timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
        next_refresh_utc=cached_at + dt.timedelta(seconds=TROPICAL_CYCLONE_CACHE_TTL_SECONDS),
        service_url=service_url,
        banner="Typhoon: none",
    )


def build_payload(request: dict[str, Any]) -> dict[str, Any]:
    payload = request.get("payload")
    if not isinstance(payload, dict):
        raise ValueError("tropical cyclone request payload is missing")
    service_url = str(payload.get("service_url", DEFAULT_SERVICE_URL))
    cache_root = Path(str(payload["cache_root"]))
    timeout_s = float(payload.get("timeout_s", DEFAULT_TIMEOUT_S))
    user_agent = str(payload.get("user_agent", DEFAULT_USER_AGENT))
    now = dt.datetime.now(dt.timezone.utc)
    cached_entry = load_tropical_cyclone_cache(cache_root)
    cached_is_stale = cached_entry is not None and is_tropical_cyclone_cache_stale(
        cached_entry, now_utc=now
    )
    if (
        cached_entry is not None
        and not cached_is_stale
        and is_tropical_cyclone_cache_current(cached_entry)
        and now
        < cached_entry.cached_at_utc
        + dt.timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS)
    ):
        cached_at = cached_entry.cached_at_utc
        return _payload(
            cached_entry.snapshot_collection,
            cached_at_utc=cached_at,
            last_checked_utc=now,
            next_check_utc=cached_at + dt.timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
            next_refresh_utc=cached_at + dt.timedelta(seconds=TROPICAL_CYCLONE_CACHE_TTL_SECONDS),
            service_url=service_url,
        )

    latest_feature = fetch_latest_observed_feature(
        service_url=service_url, timeout_s=timeout_s, user_agent=user_agent
    )
    if latest_feature is None:
        return _save_empty_overlay(cache_root=cache_root, service_url=service_url, now=now)
    latest_attrs = latest_feature.get("attributes")
    if not isinstance(latest_attrs, dict):
        raise TropicalCycloneFetchError("Observed position payload missing attributes")
    latest_storm_name = latest_attrs.get("STORMNAME")
    latest_basin = latest_attrs.get("BASIN")
    latest_advdate = latest_attrs.get("ADVDATE")
    latest_advdate_int = (
        int(latest_advdate) if isinstance(latest_advdate, (int, float)) else None
    )
    cached_snapshot = (
        cached_entry.snapshot_collection.snapshots[0]
        if cached_entry is not None and cached_entry.snapshot_collection.snapshots
        else None
    )
    if (
        cached_snapshot is not None
        and not cached_is_stale
        and cached_entry is not None
        and is_tropical_cyclone_cache_current(cached_entry)
        and cached_snapshot.has_projectable_timeline()
        and isinstance(latest_storm_name, str)
        and latest_storm_name == cached_snapshot.storm_name
        and ((latest_basin is None and cached_snapshot.basin is None) or (
            isinstance(latest_basin, str) and latest_basin == cached_snapshot.basin
        ))
        and latest_advdate_int is not None
        and cached_snapshot.advdate_utc is not None
        and int(cached_snapshot.advdate_utc.timestamp() * 1000.0) == latest_advdate_int
    ):
        cached_at = cached_entry.cached_at_utc
        return _payload(
            cached_entry.snapshot_collection,
            cached_at_utc=cached_at,
            last_checked_utc=now,
            next_check_utc=now + dt.timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
            next_refresh_utc=cached_at + dt.timedelta(seconds=TROPICAL_CYCLONE_CACHE_TTL_SECONDS),
            service_url=service_url,
        )

    snapshot_collection = fetch_active_hurricanes_snapshot(
        service_url=service_url, timeout_s=timeout_s, user_agent=user_agent
    )
    cached_at = dt.datetime.now(dt.timezone.utc)
    save_tropical_cyclone_cache(
        TropicalCycloneCacheEntry(
            snapshot_collection=snapshot_collection,
            cached_at_utc=cached_at,
            cache_version=TROPICAL_CYCLONE_CACHE_VERSION,
        ),
        cache_root=cache_root,
    )
    return _payload(
        snapshot_collection,
        cached_at_utc=cached_at,
        last_checked_utc=cached_at,
        next_check_utc=cached_at + dt.timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
        next_refresh_utc=cached_at + dt.timedelta(seconds=TROPICAL_CYCLONE_CACHE_TTL_SECONDS),
        service_url=service_url,
    )


def _write_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=True, sort_keys=True), encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    args = parser.parse_args()
    request = json.loads(args.request.read_text(encoding="utf-8"))
    result_path = args.result
    payload_path = result_path.with_name("payload.json")
    try:
        payload = build_payload(request)
        _write_json(payload_path, payload)
        result = {
            **{key: request[key] for key in (
                "protocol_version", "session_id", "worker_epoch", "request_id",
                "job_kind", "layer_generation", "view_generation",
            )},
            "status": "ok",
            "artifacts": [{
                "relative_path": payload_path.name,
                "schema": "zstarview.tropical-cyclone-payload.v1",
                "size_bytes": payload_path.stat().st_size,
            }],
        }
    except Exception as exc:
        result = {
            **{key: request.get(key) for key in (
                "protocol_version", "session_id", "worker_epoch", "request_id",
                "job_kind", "layer_generation", "view_generation",
            )},
            "status": "failed",
            "artifacts": [],
            "error_kind": type(exc).__name__,
            "error_message": str(exc),
            "error_traceback": traceback.format_exc(),
        }
    _write_json(result_path, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
