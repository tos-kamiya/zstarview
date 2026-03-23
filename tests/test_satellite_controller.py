from __future__ import annotations

from datetime import datetime, timezone
import logging
from urllib.error import URLError

import astropy.time

from zstarview.gui.satellite_controller import SatelliteController
from zstarview.satellites import CachedSatelliteElementSet


def test_expected_group_cache_miss_is_logged_without_warning(caplog) -> None:
    payloads: list[dict[str, object]] = []

    def fetcher(group_key: str, *, target_time_utc: datetime, time_mode: str) -> CachedSatelliteElementSet:
        assert time_mode == "present"
        if group_key == "station":
            raise RuntimeError("Satellites: time-shifted view is not supported")
        return CachedSatelliteElementSet(
            group_key=group_key,
            element_epoch_utc=datetime(2026, 3, 23, 5, 30, tzinfo=timezone.utc),
            fetched_at_utc=datetime(2026, 3, 23, 5, 31, tzinfo=timezone.utc),
            source="cache-fresh",
            records=[{"OBJECT_NAME": "STARLINK TEST"}],
        )

    controller = SatelliteController(
        fetcher=fetcher,
        projector=lambda *args, **kwargs: [],
    )
    controller.satellite_ready.connect(payloads.append)

    with caplog.at_level(logging.INFO):
        controller._run_update(
            observer_lat=35.0,
            observer_lon=139.0,
            observer_height_m=0.0,
            time_obj=astropy.time.Time(datetime.now(timezone.utc)),
            enabled_groups=("station", "starlink"),
            reason="test",
            request_id=0,
        )

    assert payloads
    assert payloads[0]["banner"] == "Satellites: partial (station unavailable)"
    assert isinstance(payloads[0]["refreshed_at_utc"], datetime)
    assert "Satellite element fetch unavailable for station" in caplog.text
    assert "Satellite element fetch failed for station" not in caplog.text


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
        enabled_groups=("station", "starlink"),
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
            enabled_groups=("station", "starlink"),
            reason="test",
            request_id=0,
        )

    assert calls == ["station"]
    assert "Satellite element fetch timed out for station" in caplog.text
    assert "Satellite element fetch timed out for starlink" not in caplog.text
    assert "Traceback" not in caplog.text
    assert failures == [{"banner": "Satellites: <urlopen error timed out>"}]
