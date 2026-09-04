from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np

from zstarview.clouddisc.altaz_grid import CloudAltAzGrid
from zstarview.clouddisc.types import CloudMeta, CloudSourceData, SourceKey
from zstarview.gui.cloud_controller import (
    ActiveCloudRenderRequest,
    CloudController,
    CloudRenderRequest,
)


def _make_grid(source_key, lat, lon):
    return CloudAltAzGrid(
        amount=np.zeros((4, 4), dtype=np.float32),
        missing_mask=np.zeros((4, 4), dtype=np.uint8),
        alt_min_deg=0.0,
        alt_max_deg=90.0,
        az_min_deg=0.0,
        az_max_deg=360.0,
        observer_lat=lat,
        observer_lon=lon,
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
        shells_km=(6374.0, 6376.0, 6378.0),
        source_key=source_key,
        coverage_ratio=1.0,
    )


class _FakeCloudDisc:
    def __init__(self) -> None:
        self.cfg = object()
        self._meta = CloudMeta(
            satellite="HIMAWARI",
            product="ISatSS-B13",
            time_utc=dt.datetime(2026, 3, 4, 12, 30, tzinfo=dt.timezone.utc),
            src_paths=[Path("/tmp/fake.nc")],
        )

    def fetch_source(
        self, *, lat: float, lon: float, abort_event=None
    ):  # pragma: no cover - not used in direct tests
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

    def build_altaz_grid_from_source(self, **kwargs):
        return _make_grid(kwargs["source"].source_key, kwargs["lat"], kwargs["lon"])


def _stub_render_altaz(width, height, **kwargs):
    return np.full((height, width, 4), 255, dtype=np.uint8)


def _stub_render_missing(width, height, **kwargs):
    return np.zeros((height, width), dtype=np.uint8)


def test_cloud_render_discards_stale_request_id(monkeypatch) -> None:
    disc = _FakeCloudDisc()
    controller = CloudController(disc)
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
    source.altaz_grid = disc.build_altaz_grid_from_source(
        source=source, lat=35.0, lon=139.0
    )
    controller._latest_source = source
    controller._latest_render_request_id = 2
    emitted = []
    controller.cloud_ready.connect(lambda payload: emitted.append(payload))
    monkeypatch.setattr(
        "zstarview.gui.cloud_controller.render_altaz_grid_circles",
        _stub_render_altaz,
    )
    monkeypatch.setattr(
        "zstarview.gui.cloud_controller.render_altaz_missing_mask",
        _stub_render_missing,
    )

    controller._run_render_update(
        request=ActiveCloudRenderRequest(
            request=CloudRenderRequest(
                alt=45.0,
                az=180.0,
                radius_px=128,
                content_fov_deg=90.0,
                reason="manual",
                render_generation=0,
            ),
            request_id=1,
            request_key=(),
            source_id=1,
            source_key=None,
        )
    )

    assert emitted == []


def test_cloud_update_keeps_latest_pending_render_request() -> None:
    controller = CloudController(_FakeCloudDisc())
    controller._latest_source = object()
    controller._render_is_running = True

    controller.update_render(
        alt=45.0,
        az=180.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="manual",
    )
    controller.update_render(
        alt=50.0,
        az=200.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="manual",
    )

    assert controller._pending_render_request is not None
    assert controller._pending_render_request.alt == 50.0
    assert controller._pending_render_request.az == 200.0
