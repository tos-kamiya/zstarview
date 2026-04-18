from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np

from zstarview.clouddisc.types import CloudMeta, CloudSourceData, SourceKey
from zstarview.gui.cloud_controller import CloudController


class _FakeCloudDisc:
    def __init__(self) -> None:
        self._meta = CloudMeta(
            satellite="HIMAWARI",
            product="ISatSS-B13",
            time_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
            src_paths=[Path("/tmp/fake.nc")],
        )

    def fetch_source(self, *, lat: float, lon: float):  # pragma: no cover - not used in direct tests
        return CloudSourceData(
            source_key=SourceKey(
                satellite="HIMAWARI",
                provider="HIMAWARI",
                timeslot_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
            ),
            data_array=object(),
            satellite="HIMAWARI",
            product="ISatSS-B13",
            time_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
            src_paths=[Path("/tmp/fake.nc")],
        )

    def render_from_source_with_coverage(self, **kwargs):
        radius_px = kwargs.get("radius_px")
        size = int(radius_px) * 2 if isinstance(radius_px, int) else 8
        return (
            np.full((size, size, 4), 255, dtype=np.uint8),
            self._meta,
            np.zeros((size, size), dtype=np.uint8),
            1.0,
        )


def test_cloud_render_discards_stale_request_id() -> None:
    controller = CloudController(_FakeCloudDisc())
    source = CloudSourceData(
        source_key=SourceKey(
            satellite="HIMAWARI",
            provider="HIMAWARI",
            timeslot_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
        ),
        data_array=object(),
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
        src_paths=[Path("/tmp/fake.nc")],
    )
    controller._latest_source = source
    controller._latest_request_id = 2
    emitted = []
    controller.cloud_ready.connect(lambda payload: emitted.append(payload))

    controller._run_render_update(
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=128,
        content_fov_deg=90.0,
        reason="manual",
        request_id=1,
    )

    assert emitted == []


def test_cloud_update_keeps_latest_pending_render_request() -> None:
    controller = CloudController(_FakeCloudDisc())
    controller._latest_source = object()
    controller._render_is_running = True
    controller._source_is_running = True

    controller.update(
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="manual",
    )
    controller.update(
        lat=35.0,
        lon=139.0,
        alt=50.0,
        az=200.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="manual",
    )

    assert controller._pending_render_request is not None
    assert controller._pending_render_request["alt"] == 50.0
    assert controller._pending_render_request["az"] == 200.0
