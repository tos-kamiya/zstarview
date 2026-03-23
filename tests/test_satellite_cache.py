from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json

import pytest

from zstarview.satellites.cache import (
    fetch_cached_satellite_elements,
    load_satellite_cache,
    resolve_satellite_elements_for_time,
    satellite_group_cache_path,
    save_satellite_cache,
)
from zstarview.satellite_constants import SATELLITE_GROUP_VALIDITY_SECONDS


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
        "station",
        [_sample_record()],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )

    cached = load_satellite_cache("station", cache_root=tmp_path)

    assert cached is not None
    assert cached.group_key == "station"
    assert cached.element_epoch_utc == element_epoch
    assert cached.fetched_at_utc == fetched_at
    assert len(cached.records) == 1
    assert cached.records[0]["OBJECT_NAME"] == "ISS (ZARYA)"


def test_fetch_cached_satellite_elements_uses_fresh_cache_without_network(tmp_path) -> None:
    element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "station",
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
        "station",
        fetcher=fail_fetcher,
        cache_root=tmp_path,
        now_utc=element_epoch + timedelta(hours=1),
        fresh_ttl_seconds=6 * 60 * 60,
    )

    assert not called
    assert result.source == "cache-fresh"
    assert result.element_epoch_utc == element_epoch
    assert result.fetched_at_utc == fetched_at
    assert len(result.records) == 1


def test_fetch_cached_satellite_elements_refreshes_current_cache(tmp_path) -> None:
    original_element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    original_fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "station",
        [_sample_record()],
        element_epoch_utc=original_element_epoch,
        fetched_at_utc=original_fetched_at,
        cache_root=tmp_path,
    )
    new_record = dict(_sample_record(epoch="2026-03-22T18:00:00.000000"))
    new_record["OBJECT_NAME"] = "ISS NEW"

    def fetcher(group_key, *, timeout_s):
        assert group_key == "station"
        assert timeout_s == 12.0
        return [new_record]

    now_utc = original_element_epoch + timedelta(hours=7)
    result = fetch_cached_satellite_elements(
        "station",
        fetcher=fetcher,
        timeout_s=12.0,
        cache_root=tmp_path,
        now_utc=now_utc,
        fresh_ttl_seconds=6 * 60 * 60,
    )

    assert result.source == "celestrak"
    assert result.records[0]["OBJECT_NAME"] == "ISS NEW"
    cached = load_satellite_cache("station", cache_root=tmp_path)
    assert cached is not None
    assert cached.records[0]["OBJECT_NAME"] == "ISS NEW"
    assert cached.element_epoch_utc == datetime(2026, 3, 22, 18, 0, tzinfo=timezone.utc)
    assert cached.fetched_at_utc == now_utc


def test_fetch_cached_satellite_elements_uses_group_specific_ttl(tmp_path) -> None:
    element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "station",
        [_sample_record()],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )
    save_satellite_cache(
        "starlink",
        [_sample_record()],
        element_epoch_utc=element_epoch,
        fetched_at_utc=fetched_at,
        cache_root=tmp_path,
    )

    station_called = False
    starlink_called = False

    def fail_station(*args, **kwargs):
        nonlocal station_called
        station_called = True
        raise AssertionError("station fetch should not run")

    def refresh_starlink(*args, **kwargs):
        nonlocal starlink_called
        starlink_called = True
        return [_sample_record(epoch="2026-03-22T16:00:00.000000")]

    now_utc = element_epoch + timedelta(hours=4)
    station_result = fetch_cached_satellite_elements(
        "station",
        fetcher=fail_station,
        cache_root=tmp_path,
        now_utc=now_utc,
    )
    starlink_result = fetch_cached_satellite_elements(
        "starlink",
        fetcher=refresh_starlink,
        cache_root=tmp_path,
        now_utc=now_utc,
    )

    assert SATELLITE_GROUP_VALIDITY_SECONDS["station"] == 24 * 60 * 60
    assert SATELLITE_GROUP_VALIDITY_SECONDS["starlink"] == 3 * 60 * 60
    assert station_called is False
    assert starlink_called is True
    assert station_result.source == "cache-fresh"
    assert starlink_result.source == "celestrak"


def test_resolve_satellite_elements_for_non_present_rejects_time_shifted_views(tmp_path) -> None:
    save_satellite_cache(
        "station",
        [_sample_record()],
        element_epoch_utc=datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc),
        fetched_at_utc=datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc),
        cache_root=tmp_path,
    )

    with pytest.raises(RuntimeError, match="time-shifted view is not supported"):
        resolve_satellite_elements_for_time(
            "station",
            target_time_utc=datetime(2026, 3, 21, 0, 0, tzinfo=timezone.utc),
            time_mode="past",
            cache_root=tmp_path,
            validity_seconds=6 * 60 * 60,
        )

    with pytest.raises(RuntimeError, match="time-shifted view is not supported"):
        resolve_satellite_elements_for_time(
            "station",
            target_time_utc=datetime(2026, 3, 23, 0, 0, tzinfo=timezone.utc),
            time_mode="future",
            cache_root=tmp_path,
            validity_seconds=6 * 60 * 60,
        )


def test_load_satellite_cache_derives_element_time_from_existing_payload(tmp_path) -> None:
    css_tianhe = dict(_sample_record(epoch="2026-03-22T12:30:00.000000"))
    css_tianhe["OBJECT_NAME"] = "CSS (TIANHE)"
    css_tianhe["NORAD_CAT_ID"] = "48274"
    crew_dragon = dict(_sample_record(epoch="2026-03-22T12:15:00.000000"))
    crew_dragon["OBJECT_NAME"] = "CREW DRAGON 12"
    crew_dragon["NORAD_CAT_ID"] = "99999"
    satellite_group_cache_path("station", cache_root=tmp_path).write_text(
        json.dumps(
            {
                "group_key": "station",
                "fetched_at_utc": "2026-03-22T12:05:00+00:00",
                "source": "cache",
                "records": [_sample_record(), css_tianhe, crew_dragon],
            }
        ),
        encoding="utf-8",
    )

    cached = load_satellite_cache("station", cache_root=tmp_path)

    assert cached is not None
    assert cached.element_epoch_utc == datetime(2026, 3, 22, 12, 30, tzinfo=timezone.utc)
    assert cached.fetched_at_utc == datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    assert [record["NORAD_CAT_ID"] for record in cached.records] == ["25544", "48274"]
