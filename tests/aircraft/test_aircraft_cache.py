from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest
from filelock import Timeout

from zstarview.aircraft.cache import (
    aircraft_cache_path,
    aircraft_fetch_lock_path,
    aircraft_rate_limit_path,
    cleanup_aircraft_cache,
    fetch_cached_opensky_states,
    load_aircraft_cache,
    load_aircraft_rate_limit,
    save_aircraft_cache,
    save_aircraft_rate_limit,
)
from zstarview.aircraft.opensky import AircraftBoundingBox
from zstarview.aircraft.types import AircraftSnapshot


def _sample_bbox() -> AircraftBoundingBox:
    return AircraftBoundingBox(min_lat=34.0, max_lat=36.0, min_lon=132.0, max_lon=134.0)


def _sample_snapshot() -> AircraftSnapshot:
    return AircraftSnapshot(
        icao24="abcd01",
        callsign="TEST123",
        latitude=35.0,
        longitude=133.0,
        baro_altitude_m=10000.0,
        velocity_mps=250.0,
        heading_deg=90.0,
        vertical_rate_mps=0.0,
        on_ground=False,
        last_contact_unix=1710000000,
    )


def test_save_and_load_aircraft_cache_round_trip(tmp_path) -> None:
    bbox = _sample_bbox()
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_aircraft_cache(bbox, [_sample_snapshot()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    cached = load_aircraft_cache(bbox, cache_root=tmp_path)

    assert cached is not None
    assert cached.bbox == bbox
    assert cached.fetched_at_utc == fetched_at
    assert len(cached.snapshots) == 1
    assert cached.snapshots[0].icao24 == "abcd01"


def test_fetch_cached_opensky_states_uses_fresh_cache_without_network(tmp_path) -> None:
    bbox = _sample_bbox()
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_aircraft_cache(bbox, [_sample_snapshot()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    called = False

    def fail_fetcher(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("network fetch should not run")

    result = fetch_cached_opensky_states(
        bbox,
        fetcher=fail_fetcher,
        cache_root=tmp_path,
        now_utc=fetched_at + timedelta(seconds=30),
        fresh_ttl_seconds=180,
        stale_fallback_seconds=600,
    )

    assert not called
    assert result.source == "cache-fresh"
    assert not result.is_stale
    assert len(result.snapshots) == 1


def test_fetch_cached_opensky_states_falls_back_to_stale_cache_on_failure(tmp_path) -> None:
    bbox = _sample_bbox()
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_aircraft_cache(bbox, [_sample_snapshot()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    def failing_fetcher(*args, **kwargs):
        raise RuntimeError("opensky unavailable")

    result = fetch_cached_opensky_states(
        bbox,
        fetcher=failing_fetcher,
        cache_root=tmp_path,
        now_utc=fetched_at + timedelta(seconds=240),
        fresh_ttl_seconds=180,
        stale_fallback_seconds=600,
    )

    assert result.source == "cache-stale"
    assert result.is_stale
    assert result.fetched_at_utc == fetched_at


def test_fetch_cached_opensky_states_refreshes_and_overwrites_cache(tmp_path) -> None:
    bbox = _sample_bbox()
    original_fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_aircraft_cache(bbox, [_sample_snapshot()], fetched_at_utc=original_fetched_at, cache_root=tmp_path)
    new_snapshot = AircraftSnapshot(
        icao24="efgh02",
        callsign="NEW456",
        latitude=35.1,
        longitude=133.1,
        baro_altitude_m=11000.0,
        velocity_mps=240.0,
        heading_deg=100.0,
        vertical_rate_mps=1.0,
        on_ground=False,
        last_contact_unix=1710000300,
    )

    def fetcher(_bbox, *, timeout_s):
        assert timeout_s == 12.0
        return [new_snapshot]

    now_utc = original_fetched_at + timedelta(seconds=240)
    result = fetch_cached_opensky_states(
        bbox,
        fetcher=fetcher,
        timeout_s=12.0,
        cache_root=tmp_path,
        now_utc=now_utc,
        fresh_ttl_seconds=180,
        stale_fallback_seconds=600,
    )

    assert result.source == "opensky"
    assert not result.is_stale
    assert result.snapshots[0].icao24 == "efgh02"
    cached = load_aircraft_cache(bbox, cache_root=tmp_path)
    assert cached is not None
    assert cached.snapshots[0].icao24 == "efgh02"
    assert cached.fetched_at_utc == now_utc
    assert load_aircraft_rate_limit(cache_root=tmp_path) == now_utc


def test_fetch_cached_opensky_states_skips_when_global_rate_limit_is_fresh(tmp_path) -> None:
    bbox = _sample_bbox()
    now_utc = datetime(2026, 3, 22, 6, 5, tzinfo=timezone.utc)
    last_success = now_utc - timedelta(seconds=60)
    save_aircraft_rate_limit(last_success, cache_root=tmp_path)

    called = False

    def fail_fetcher(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("network fetch should not run")

    result = fetch_cached_opensky_states(
        bbox,
        fetcher=fail_fetcher,
        cache_root=tmp_path,
        now_utc=now_utc,
        fresh_ttl_seconds=300,
        stale_fallback_seconds=600,
    )

    assert not called
    assert result.source == "rate-limited-skip"
    assert result.snapshots == []
    assert result.fetched_at_utc == last_success
    assert not result.is_stale


def test_fetch_cached_opensky_states_does_not_show_stale_cache_on_global_skip(
    tmp_path,
) -> None:
    bbox = _sample_bbox()
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    now_utc = fetched_at + timedelta(seconds=240)
    save_aircraft_cache(bbox, [_sample_snapshot()], fetched_at_utc=fetched_at, cache_root=tmp_path)
    save_aircraft_rate_limit(now_utc - timedelta(seconds=30), cache_root=tmp_path)

    def fail_fetcher(*args, **kwargs):
        raise AssertionError("network fetch should not run")

    result = fetch_cached_opensky_states(
        bbox,
        fetcher=fail_fetcher,
        cache_root=tmp_path,
        now_utc=now_utc,
        fresh_ttl_seconds=180,
        stale_fallback_seconds=600,
    )

    assert result.source == "rate-limited-skip"
    assert result.snapshots == []
    assert not result.is_stale


def test_fetch_cached_opensky_states_can_bypass_global_rate_limit_for_export(
    tmp_path,
) -> None:
    bbox = _sample_bbox()
    now_utc = datetime(2026, 3, 22, 6, 5, tzinfo=timezone.utc)
    save_aircraft_rate_limit(now_utc - timedelta(seconds=60), cache_root=tmp_path)
    new_snapshot = AircraftSnapshot(
        icao24="export1",
        callsign="EXP001",
        latitude=35.2,
        longitude=133.2,
        baro_altitude_m=9000.0,
        velocity_mps=220.0,
        heading_deg=120.0,
        vertical_rate_mps=0.0,
        on_ground=False,
        last_contact_unix=1710000400,
    )

    def fetcher(_bbox, *, timeout_s):
        return [new_snapshot]

    result = fetch_cached_opensky_states(
        bbox,
        fetcher=fetcher,
        cache_root=tmp_path,
        now_utc=now_utc,
        fresh_ttl_seconds=300,
        stale_fallback_seconds=600,
        enforce_global_rate_limit=False,
    )

    assert result.source == "opensky"
    assert result.snapshots == [new_snapshot]
    assert load_aircraft_rate_limit(cache_root=tmp_path) == now_utc


def test_fetch_cached_opensky_states_skips_gui_update_on_lock_timeout(
    tmp_path,
    monkeypatch,
) -> None:
    bbox = _sample_bbox()

    class FailingLock:
        def __init__(self, *_args, **_kwargs) -> None:
            pass

        def acquire(self) -> None:
            raise Timeout(str(aircraft_fetch_lock_path(cache_root=tmp_path)))

    monkeypatch.setattr("zstarview.aircraft.cache.FileLock", FailingLock)

    result = fetch_cached_opensky_states(
        bbox,
        fetcher=lambda *_args, **_kwargs: [],
        cache_root=tmp_path,
    )

    assert result.source == "rate-limited-skip"
    assert result.snapshots == []


def test_fetch_cached_opensky_states_propagates_export_lock_timeout(
    tmp_path,
    monkeypatch,
) -> None:
    bbox = _sample_bbox()

    class FailingLock:
        def __init__(self, *_args, **_kwargs) -> None:
            pass

        def acquire(self) -> None:
            raise Timeout(str(aircraft_fetch_lock_path(cache_root=tmp_path)))

    monkeypatch.setattr("zstarview.aircraft.cache.FileLock", FailingLock)

    with pytest.raises(Timeout):
        fetch_cached_opensky_states(
            bbox,
            fetcher=lambda *_args, **_kwargs: [],
            cache_root=tmp_path,
            enforce_global_rate_limit=False,
        )


def test_cleanup_aircraft_cache_removes_old_files(tmp_path) -> None:
    bbox = _sample_bbox()
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_aircraft_cache(bbox, [_sample_snapshot()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    cleanup_aircraft_cache(
        cache_root=tmp_path,
        now_utc=fetched_at + timedelta(seconds=601),
        max_age_seconds=600,
    )

    assert not aircraft_cache_path(bbox, cache_root=tmp_path).exists()


def test_cleanup_aircraft_cache_keeps_rate_limit_metadata(tmp_path) -> None:
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_aircraft_rate_limit(fetched_at, cache_root=tmp_path)

    cleanup_aircraft_cache(
        cache_root=tmp_path,
        now_utc=fetched_at + timedelta(seconds=601),
        max_age_seconds=600,
    )

    assert aircraft_rate_limit_path(cache_root=tmp_path).is_file()


def test_fetch_cached_opensky_states_raises_when_stale_cache_is_too_old(tmp_path) -> None:
    bbox = _sample_bbox()
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_aircraft_cache(bbox, [_sample_snapshot()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    def failing_fetcher(*args, **kwargs):
        raise RuntimeError("opensky unavailable")

    with pytest.raises(RuntimeError, match="opensky unavailable"):
        fetch_cached_opensky_states(
            bbox,
            fetcher=failing_fetcher,
            cache_root=tmp_path,
            now_utc=fetched_at + timedelta(seconds=601),
            fresh_ttl_seconds=180,
            stale_fallback_seconds=600,
        )
