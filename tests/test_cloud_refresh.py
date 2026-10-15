from __future__ import annotations

from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from zstarview.cloud_refresh import CloudRefreshSchedule
from zstarview.gui.cloud_state import CloudImageState
from zstarview.gui.geosatellite_state import GeoSatelliteState
from zstarview.gui.window_state import SkyWindowState
from zstarview.gui.window_update_cloud import SkyWindowCloudUpdatesMixin

T = datetime(2026, 10, 16, 0, 30, tzinfo=timezone.utc)


@pytest.mark.parametrize(
    "age,expected", [(0, 600), (1799, 600), (1800, 300), (3600, 300)]
)
def test_refresh_age_threshold(age: int, expected: int) -> None:
    assert (
        CloudRefreshSchedule().next_delay_seconds(
            observed_at=T - timedelta(seconds=age), completed_at=T, source="HIMAWARI"
        )
        == expected
    )


def test_repeated_old_observation_does_not_keep_accelerating() -> None:
    policy = CloudRefreshSchedule()
    for minute, expected in ((0, 300), (5, 600), (15, 600)):
        assert (
            policy.next_delay_seconds(
                observed_at=T - timedelta(minutes=40),
                completed_at=T + timedelta(minutes=minute),
                source="HIMAWARI",
            )
            == expected
        )
    assert (
        policy.next_delay_seconds(
            observed_at=T - timedelta(minutes=30),
            completed_at=T + timedelta(minutes=25),
            source="HIMAWARI",
        )
        == 300
    )
    assert (
        policy.next_delay_seconds(
            observed_at=T - timedelta(minutes=30),
            completed_at=T + timedelta(minutes=25),
            source="G19",
        )
        == 300
    )


def test_missing_future_and_timezone_timestamps() -> None:
    policy = CloudRefreshSchedule()
    assert (
        policy.next_delay_seconds(observed_at=None, completed_at=T, source="x") == 600
    )
    assert (
        policy.next_delay_seconds(
            observed_at=T + timedelta(seconds=1), completed_at=T, source="x"
        )
        == 600
    )
    observed = (T - timedelta(minutes=30)).astimezone(timezone(timedelta(hours=9)))
    assert (
        policy.next_delay_seconds(observed_at=observed, completed_at=T, source="x")
        == 300
    )
    assert (
        policy.next_delay_seconds(
            observed_at=observed.astimezone(timezone.utc).replace(tzinfo=None),
            completed_at=T.replace(tzinfo=None),
            source="x",
        )
        == 600
    )


def test_native_source_success_advances_then_restores_refresh() -> None:
    owner: Any = SimpleNamespace(
        cloud_state=CloudImageState(),
        state=SkyWindowState(render_view_center=(45, 180)),
        reproject_cloud_overlay=lambda **_: None,
        _geo_satellite_mode_active=lambda: False,
    )
    observed = T - timedelta(minutes=40)
    payload = dict(
        satellite="HIMAWARI",
        source=SimpleNamespace(time_utc=observed),
        refreshed_at_utc=T,
    )
    SkyWindowCloudUpdatesMixin._on_cloud_source_ready(owner, payload)
    assert owner.state.cloud_next_refresh_utc == T + timedelta(minutes=5)
    payload["refreshed_at_utc"] = T + timedelta(minutes=5)
    SkyWindowCloudUpdatesMixin._on_cloud_source_ready(owner, payload)
    assert owner.state.cloud_next_refresh_utc == T + timedelta(minutes=15)


def test_geo_success_advances_then_restores_refresh() -> None:
    owner: Any = SimpleNamespace(
        geosatellite_state=GeoSatelliteState(),
        state=SkyWindowState(render_view_center=(45, 180)),
        _disc_generation=0,
        _is_shutting_down=False,
        _compositor=SimpleNamespace(invalidate=lambda: None),
        _safe_request_cloud_repaint=lambda: None,
    )
    payload = dict(
        image=np.zeros((1, 1, 4)),
        az=180,
        captured_at_utc=T - timedelta(minutes=40),
        refreshed_at_utc=T,
    )
    SkyWindowCloudUpdatesMixin._on_geosatellite_ready(owner, payload)
    assert owner.state.cloud_next_refresh_utc == T + timedelta(minutes=5)
    payload["refreshed_at_utc"] = T + timedelta(minutes=5)
    SkyWindowCloudUpdatesMixin._on_geosatellite_ready(owner, payload)
    assert owner.state.cloud_next_refresh_utc == T + timedelta(minutes=15)


def test_failure_keeps_last_frame_and_normal_retry(monkeypatch) -> None:
    import zstarview.gui.window_update_cloud as module

    class Clock(datetime):
        @classmethod
        def now(cls, tz=None):
            return T

    monkeypatch.setattr(module, "datetime", Clock)
    image = np.zeros((1, 1, 4))
    owner: Any = SimpleNamespace(
        cloud_state=CloudImageState(image=image),
        geosatellite_state=GeoSatelliteState(image=image),
        state=SkyWindowState(render_view_center=(45, 180)),
        _safe_request_cloud_repaint=lambda: None,
    )
    SkyWindowCloudUpdatesMixin._on_cloud_failed(owner, {})
    assert owner.state.cloud_next_refresh_utc == T + timedelta(minutes=10)
    assert owner.cloud_state.image is image
    SkyWindowCloudUpdatesMixin._on_geosatellite_failed(owner, {})
    assert owner.state.cloud_next_refresh_utc == T + timedelta(minutes=10)
    assert owner.geosatellite_state.image is image


def test_clearing_cloud_data_resets_refresh_history() -> None:
    for state in (CloudImageState(), GeoSatelliteState()):
        assert (
            state.refresh_schedule.next_delay_seconds(
                observed_at=T - timedelta(minutes=40), completed_at=T, source="x"
            )
            == 300
        )
        state.clear_data()
        assert (
            state.refresh_schedule.next_delay_seconds(
                observed_at=T - timedelta(minutes=40), completed_at=T, source="x"
            )
            == 300
        )
