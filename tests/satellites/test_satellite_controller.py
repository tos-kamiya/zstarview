from __future__ import annotations

import logging
from datetime import datetime, timezone
from urllib.error import URLError

import astropy.time

from zstarview.gui.satellite_controller import SatelliteController
from zstarview.satellites import CachedSatelliteElementSet


def test_expected_cache_miss_is_logged_without_warning(caplog) -> None:
    failures: list[dict[str, object]] = []

    def fetcher(group_key: str, *, target_time_utc: datetime, time_mode: str) -> CachedSatelliteElementSet:
        assert time_mode == "present"
        raise RuntimeError("Satellites: time-shifted view is not supported")

    controller = SatelliteController(
        fetcher=fetcher,
        projector=lambda *args, **kwargs: [],
    )
    controller.satellite_failed.connect(failures.append)

    with caplog.at_level(logging.INFO):
        controller._run_update(
            observer_lat=35.0,
            observer_lon=139.0,
            observer_height_m=0.0,
            time_obj=astropy.time.Time(datetime.now(timezone.utc)),
            enabled_groups=("iss",),
            reason="test",
            request_id=0,
        )

    assert failures == [{"banner": "Satellites: time-shifted view is not supported"}]
    assert "Satellite element fetch unavailable for iss" in caplog.text
    assert "Satellite element fetch failed for iss" not in caplog.text


def test_all_expected_cache_misses_surface_specific_failure_banner() -> None:
    failures: list[dict[str, object]] = []

    def fetcher(group_key: str, *, target_time_utc: datetime, time_mode: str) -> CachedSatelliteElementSet:
        raise RuntimeError("Satellites: time-shifted view is not supported")

    controller = SatelliteController(
        fetcher=fetcher,
        projector=lambda *args, **kwargs: [],
    )
    controller.satellite_failed.connect(failures.append)

    controller._run_update(
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
        time_obj=astropy.time.Time(datetime.now(timezone.utc)),
        enabled_groups=("iss",),
        reason="test",
        request_id=0,
    )

    assert failures == [{"banner": "Satellites: time-shifted view is not supported"}]


def test_timeout_urlerror_is_logged_without_traceback(caplog) -> None:
    failures: list[dict[str, object]] = []
    calls: list[str] = []

    def fetcher(group_key: str, *, target_time_utc: datetime, time_mode: str) -> CachedSatelliteElementSet:
        calls.append(group_key)
        raise URLError(TimeoutError("timed out"))

    controller = SatelliteController(
        fetcher=fetcher,
        projector=lambda *args, **kwargs: [],
    )
    controller.satellite_failed.connect(failures.append)

    with caplog.at_level(logging.WARNING):
        controller._run_update(
            observer_lat=35.0,
            observer_lon=139.0,
            observer_height_m=0.0,
            time_obj=astropy.time.Time(datetime.now(timezone.utc)),
            enabled_groups=("iss",),
            reason="test",
            request_id=0,
        )

    assert calls == ["iss"]
    assert "Satellite element fetch timed out for iss" in caplog.text
    assert "Traceback" not in caplog.text
    assert failures == [{"banner": "Satellites: <urlopen error timed out>"}]


def test_satellite_controller_does_not_project_in_fetch_stage() -> None:
    projector_calls: list[object] = []

    def fetcher(group_key: str, *, target_time_utc: datetime, time_mode: str) -> CachedSatelliteElementSet:
        return CachedSatelliteElementSet(
            group_key=group_key,
            element_epoch_utc=datetime.now(timezone.utc),
            records=[],
            fetched_at_utc=datetime.now(timezone.utc),
            source="test",
        )

    def projector(*args, **kwargs):
        projector_calls.append((args, kwargs))
        raise AssertionError("projector should not run during fetch stage")

    controller = SatelliteController(fetcher=fetcher, projector=projector)
    controller._run_update(
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
        time_obj=astropy.time.Time(datetime.now(timezone.utc)),
        enabled_groups=("iss",),
        reason="test",
        request_id=0,
    )

    assert projector_calls == []
