from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json

import pytest

from zstarview.satellites.cache import (
    fetch_cached_satellite_elements,
    load_satellite_cache,
    satellite_cache_scope_key,
    resolve_satellite_elements_for_time,
    save_satellite_fetch_failure,
    satellite_group_cache_path,
    save_satellite_cache,
)
from zstarview.satellite_constants import (
    SATELLITE_CACHE_FORMAT_VERSION,
    SATELLITE_GROUP_VALIDITY_SECONDS,
)


def _sample_record(epoch: str = "2026-03-22T12:00:00.000000") -> dict[str, object]:
    return {
        "OBJECT_NAME": "ISS (ZARYA)",
        "OBJECT_ID": "1998-067A",
        "EPOCH": epoch,
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
    element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "iss",
        [_sample_record()],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )

    cached = load_satellite_cache("iss", cache_root=tmp_path)

    assert cached is not None
    assert cached.group_key == "iss"
    assert cached.element_epoch_utc == element_epoch
    assert cached.fetched_at_utc == fetched_at
    assert cached.last_fetch_attempt_utc == fetched_at
    assert cached.last_fetch_failed is False
    assert len(cached.records) == 1
    assert cached.records[0]["OBJECT_NAME"] == "ISS (ZARYA)"


def test_save_satellite_cache_embeds_format_version(tmp_path) -> None:
    element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    path = save_satellite_cache(
        "iss",
        [_sample_record()],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["cache_format_version"] == SATELLITE_CACHE_FORMAT_VERSION


def test_save_and_load_satellite_cache_round_trip_with_scope_key(tmp_path) -> None:
    element_epoch = datetime(2026, 4, 17, 12, 0, tzinfo=timezone.utc)
    fetched_at = datetime(2026, 4, 17, 12, 5, tzinfo=timezone.utc)
    scope_key = satellite_cache_scope_key(
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=50.0,
    )
    save_satellite_cache(
        "horizons",
        [
            {
                "OBJECT_NAME": "JWST",
                "EPOCH": element_epoch.isoformat(),
                "ALT_DEG": 12.5,
                "AZ_DEG": 220.0,
                "_SOURCE": "horizons",
            }
        ],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
        cache_scope_key=scope_key,
        source="horizons",
    )

    cached = load_satellite_cache("horizons", cache_root=tmp_path, cache_scope_key=scope_key)

    assert cached is not None
    assert cached.group_key == "horizons"
    assert cached.records[0]["OBJECT_NAME"] == "JWST"
    assert satellite_group_cache_path("horizons", cache_root=tmp_path, cache_scope_key=scope_key).is_file()


def test_load_satellite_cache_rejects_old_format(tmp_path) -> None:
    path = satellite_group_cache_path("iss", cache_root=tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "group_key": "iss",
                "element_epoch_utc": "2026-03-22T12:00:00+00:00",
                "fetched_at_utc": "2026-03-22T12:05:00+00:00",
                "source": "celestrak",
                "records": [_sample_record()],
            },
            separators=(",", ":"),
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    cached = load_satellite_cache("iss", cache_root=tmp_path)

    assert cached is None


def test_fetch_cached_satellite_elements_uses_fresh_cache_without_network(tmp_path) -> None:
    element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "iss",
        [_sample_record()],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )

    called = False

    def fail_fetcher(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("network fetch should not run")

    result = fetch_cached_satellite_elements(
        "iss",
        fetcher=fail_fetcher,
        cache_root=tmp_path,
        now_utc=element_epoch + timedelta(hours=1),
        fresh_ttl_seconds=6 * 60 * 60,
    )

    assert not called
    assert result.source == "cache-fresh"
    assert result.element_epoch_utc == element_epoch
    assert result.fetched_at_utc == fetched_at
    assert result.last_fetch_attempt_utc == fetched_at
    assert len(result.records) == 1


def test_fetch_cached_satellite_elements_refreshes_current_cache(tmp_path) -> None:
    original_element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    original_fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "iss",
        [_sample_record()],
        element_epoch_utc=original_element_epoch,
        fetched_at_utc=original_fetched_at,
        cache_root=tmp_path,
    )
    new_record = dict(_sample_record(epoch="2026-03-22T18:00:00.000000"))
    new_record["OBJECT_NAME"] = "ISS NEW"
    new_record["_SOURCE"] = "wheretheiss"

    def fetcher(group_key, *, timeout_s):
        assert group_key == "iss"
        assert timeout_s == 12.0
        return [new_record]

    now_utc = original_element_epoch + timedelta(hours=7)
    result = fetch_cached_satellite_elements(
        "iss",
        fetcher=fetcher,
        timeout_s=12.0,
        cache_root=tmp_path,
        now_utc=now_utc,
        fresh_ttl_seconds=6 * 60 * 60,
    )

    assert result.source == "wheretheiss"
    assert result.records[0]["OBJECT_NAME"] == "ISS NEW"
    cached = load_satellite_cache("iss", cache_root=tmp_path)
    assert cached is not None
    assert cached.records[0]["OBJECT_NAME"] == "ISS NEW"
    assert cached.element_epoch_utc == datetime(2026, 3, 22, 18, 0, tzinfo=timezone.utc)
    assert cached.fetched_at_utc == now_utc
    assert cached.last_fetch_failed is False
    assert cached.last_fetch_attempt_utc == now_utc


def test_fetch_cached_satellite_elements_uses_iss_ttl(tmp_path) -> None:
    element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "iss",
        [_sample_record()],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )

    called = False

    def fail_fetcher(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("iss fetch should not run")

    now_utc = element_epoch + timedelta(hours=4)
    result = fetch_cached_satellite_elements(
        "iss",
        fetcher=fail_fetcher,
        cache_root=tmp_path,
        now_utc=now_utc,
    )

    assert SATELLITE_GROUP_VALIDITY_SECONDS["iss"] == 24 * 60 * 60
    assert SATELLITE_GROUP_VALIDITY_SECONDS["horizons"] == 24 * 60 * 60
    assert called is False
    assert result.source == "cache-fresh"


def test_resolve_satellite_elements_for_non_present_rejects_time_shifted_views(tmp_path) -> None:
    save_satellite_cache(
        "iss",
        [_sample_record()],
        element_epoch_utc=datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc),
        fetched_at_utc=datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc),
        cache_root=tmp_path,
    )

    with pytest.raises(RuntimeError, match="time-shifted view is not supported"):
        resolve_satellite_elements_for_time(
            "iss",
            target_time_utc=datetime(2026, 3, 21, 0, 0, tzinfo=timezone.utc),
            time_mode="past",
            cache_root=tmp_path,
            validity_seconds=6 * 60 * 60,
        )

    with pytest.raises(RuntimeError, match="time-shifted view is not supported"):
        resolve_satellite_elements_for_time(
            "iss",
            target_time_utc=datetime(2026, 3, 23, 0, 0, tzinfo=timezone.utc),
            time_mode="future",
            cache_root=tmp_path,
            validity_seconds=6 * 60 * 60,
        )

def test_failed_fetch_persists_backoff_and_reuses_stale_cache(tmp_path) -> None:
    element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "iss",
        [_sample_record()],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )
    failed_at = element_epoch + timedelta(hours=25)

    with pytest.raises(RuntimeError, match="timed out"):
        fetch_cached_satellite_elements(
            "iss",
            fetcher=lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("timed out")),
            cache_root=tmp_path,
            now_utc=failed_at,
        )

    cached = load_satellite_cache("iss", cache_root=tmp_path)
    assert cached is not None
    assert cached.last_fetch_failed is True
    assert cached.last_fetch_error == "timed out"
    assert cached.last_fetch_failure_utc == failed_at
    assert cached.failure_backoff_until_utc == failed_at + timedelta(hours=2)

    called = False

    def fail_if_called(*_args, **_kwargs):
        nonlocal called
        called = True
        raise AssertionError("network fetch should not run during backoff")

    reused = fetch_cached_satellite_elements(
        "iss",
        fetcher=fail_if_called,
        cache_root=tmp_path,
        now_utc=failed_at + timedelta(minutes=30),
    )

    assert called is False
    assert reused.source == "cache-backoff"
    assert reused.last_fetch_failed is True
    assert reused.last_fetch_error == "timed out"


def test_save_satellite_fetch_failure_creates_metadata_only_payload(tmp_path) -> None:
    failed_at = datetime(2026, 3, 23, 1, 0, tzinfo=timezone.utc)
    save_satellite_fetch_failure(
        "iss",
        attempted_at_utc=failed_at,
        error_text="timed out",
        cache_root=tmp_path,
    )

    payload = json.loads(satellite_group_cache_path("iss", cache_root=tmp_path).read_text(encoding="utf-8"))

    assert payload["group_key"] == "iss"
    assert payload["last_fetch_failed"] is True
    assert payload["last_fetch_error"] == "timed out"
    assert payload["last_fetch_attempt_utc"] == failed_at.isoformat()
