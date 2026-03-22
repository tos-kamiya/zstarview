from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from zstarview.satellites.cache import (
    cleanup_satellite_cache,
    fetch_cached_satellite_elements,
    load_satellite_cache,
    satellite_group_cache_path,
    save_satellite_cache,
)


def _sample_record() -> dict[str, object]:
    return {
        "OBJECT_NAME": "ISS (ZARYA)",
        "OBJECT_ID": "1998-067A",
        "EPOCH": "2026-03-22T12:00:00.000000",
        "MEAN_MOTION": "15.50000000",
        "ECCENTRICITY": "0.0006703",
        "INCLINATION": "51.6400",
        "RA_OF_ASC_NODE": "257.7243",
        "ARG_OF_PERICENTER": "130.5360",
        "MEAN_ANOMALY": "325.0288",
        "EPHEMERIS_TYPE": "0",
        "CLASSIFICATION_TYPE": "U",
        "NORAD_CAT_ID": "25544",
        "ELEMENT_SET_NO": "999",
        "REV_AT_EPOCH": "12345",
        "BSTAR": "0.00010234",
        "MEAN_MOTION_DOT": "0.00002182",
        "MEAN_MOTION_DDOT": "0.00000000",
    }


def test_save_and_load_satellite_cache_round_trip(tmp_path) -> None:
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_satellite_cache("iss", [_sample_record()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    cached = load_satellite_cache("iss", cache_root=tmp_path)

    assert cached is not None
    assert cached.group_key == "iss"
    assert cached.fetched_at_utc == fetched_at
    assert len(cached.records) == 1
    assert cached.records[0]["OBJECT_NAME"] == "ISS (ZARYA)"


def test_fetch_cached_satellite_elements_uses_fresh_cache_without_network(tmp_path) -> None:
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_satellite_cache("iss", [_sample_record()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    called = False

    def fail_fetcher(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("network fetch should not run")

    result = fetch_cached_satellite_elements(
        "iss",
        fetcher=fail_fetcher,
        cache_root=tmp_path,
        now_utc=fetched_at + timedelta(hours=1),
        fresh_ttl_seconds=24 * 60 * 60,
    )

    assert not called
    assert result.source == "cache-fresh"
    assert len(result.records) == 1


def test_fetch_cached_satellite_elements_refreshes_and_overwrites_cache(tmp_path) -> None:
    original_fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_satellite_cache("iss", [_sample_record()], fetched_at_utc=original_fetched_at, cache_root=tmp_path)
    new_record = dict(_sample_record())
    new_record["OBJECT_NAME"] = "ISS NEW"

    def fetcher(group_key, *, timeout_s):
        assert group_key == "iss"
        assert timeout_s == 12.0
        return [new_record]

    now_utc = original_fetched_at + timedelta(days=2)
    result = fetch_cached_satellite_elements(
        "iss",
        fetcher=fetcher,
        timeout_s=12.0,
        cache_root=tmp_path,
        now_utc=now_utc,
        fresh_ttl_seconds=24 * 60 * 60,
    )

    assert result.source == "celestrak"
    assert result.records[0]["OBJECT_NAME"] == "ISS NEW"
    cached = load_satellite_cache("iss", cache_root=tmp_path)
    assert cached is not None
    assert cached.records[0]["OBJECT_NAME"] == "ISS NEW"
    assert cached.fetched_at_utc == now_utc


def test_cleanup_satellite_cache_removes_old_files(tmp_path) -> None:
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_satellite_cache("iss", [_sample_record()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    cleanup_satellite_cache(
        cache_root=tmp_path,
        now_utc=fetched_at + timedelta(days=2),
        max_age_seconds=24 * 60 * 60,
    )

    assert not satellite_group_cache_path("iss", cache_root=tmp_path).exists()


def test_fetch_cached_satellite_elements_raises_when_cache_is_stale_and_refetch_fails(tmp_path) -> None:
    fetched_at = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_satellite_cache("iss", [_sample_record()], fetched_at_utc=fetched_at, cache_root=tmp_path)

    def failing_fetcher(*args, **kwargs):
        raise RuntimeError("celestrak unavailable")

    with pytest.raises(RuntimeError, match="celestrak unavailable"):
        fetch_cached_satellite_elements(
            "iss",
            fetcher=failing_fetcher,
            cache_root=tmp_path,
            now_utc=fetched_at + timedelta(days=2),
            fresh_ttl_seconds=24 * 60 * 60,
        )
