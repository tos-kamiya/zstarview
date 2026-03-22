from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json

import pytest

from zstarview.satellites.cache import (
    cleanup_satellite_cache,
    fetch_cached_satellite_elements,
    load_satellite_cache,
    resolve_satellite_elements_for_time,
    satellite_archive_cache_path,
    satellite_group_cache_path,
    save_satellite_cache,
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


def test_fetch_cached_satellite_elements_refreshes_and_archives_previous_current(tmp_path) -> None:
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
    archived_path = satellite_archive_cache_path("station", original_element_epoch, cache_root=tmp_path)
    assert archived_path.exists()


def test_cleanup_satellite_cache_removes_old_archive_files_only(tmp_path) -> None:
    current_element_epoch = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
    current_fetched_at = datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc)
    old_archive_epoch = datetime(2026, 3, 18, 12, 0, tzinfo=timezone.utc)
    old_archive_fetched_at = datetime(2026, 3, 18, 12, 5, tzinfo=timezone.utc)
    save_satellite_cache(
        "station",
        [_sample_record()],
        element_epoch_utc=current_element_epoch,
        fetched_at_utc=current_fetched_at,
        cache_root=tmp_path,
    )
    archive_path = satellite_archive_cache_path("station", old_archive_epoch, cache_root=tmp_path)
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    archive_path.write_text(
        json.dumps(
            {
                "group_key": "station",
                "element_epoch_utc": old_archive_epoch.isoformat(),
                "fetched_at_utc": old_archive_fetched_at.isoformat(),
                "source": "archive",
                "records": [_sample_record(epoch="2026-03-18T12:00:00.000000")],
            }
        ),
        encoding="utf-8",
    )

    cleanup_satellite_cache(
        cache_root=tmp_path,
        now_utc=current_element_epoch,
        max_age_seconds=3 * 24 * 60 * 60,
    )

    assert satellite_group_cache_path("station", cache_root=tmp_path).exists()
    assert not archive_path.exists()


def test_resolve_satellite_elements_for_past_uses_nearest_archive_snapshot(tmp_path) -> None:
    current_element_epoch = datetime(2026, 3, 22, 18, 0, tzinfo=timezone.utc)
    archive_element_epoch = datetime(2026, 3, 22, 6, 0, tzinfo=timezone.utc)
    save_satellite_cache(
        "station",
        [_sample_record(epoch="2026-03-22T18:00:00.000000")],
        element_epoch_utc=current_element_epoch,
        fetched_at_utc=current_element_epoch + timedelta(minutes=5),
        cache_root=tmp_path,
    )
    archive_path = satellite_archive_cache_path("station", archive_element_epoch, cache_root=tmp_path)
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    archive_path.write_text(
        json.dumps(
            {
                "group_key": "station",
                "element_epoch_utc": archive_element_epoch.isoformat(),
                "fetched_at_utc": (archive_element_epoch + timedelta(minutes=5)).isoformat(),
                "source": "archive",
                "records": [_sample_record(epoch="2026-03-22T06:00:00.000000")],
            }
        ),
        encoding="utf-8",
    )

    result = resolve_satellite_elements_for_time(
        "station",
        target_time_utc=datetime(2026, 3, 22, 7, 0, tzinfo=timezone.utc),
        time_mode="past",
        cache_root=tmp_path,
        validity_seconds=6 * 60 * 60,
    )

    assert result.source == "cache-archive"
    assert result.element_epoch_utc == archive_element_epoch


def test_resolve_satellite_elements_for_past_rejects_cache_outside_validity_window(tmp_path) -> None:
    save_satellite_cache(
        "station",
        [_sample_record()],
        element_epoch_utc=datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc),
        fetched_at_utc=datetime(2026, 3, 22, 12, 5, tzinfo=timezone.utc),
        cache_root=tmp_path,
    )

    with pytest.raises(RuntimeError, match="validity window"):
        resolve_satellite_elements_for_time(
            "station",
            target_time_utc=datetime(2026, 3, 21, 0, 0, tzinfo=timezone.utc),
            time_mode="past",
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
