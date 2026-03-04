from __future__ import annotations

from zstarview.ui.cloud_controller import CloudController


class _DummyCloudDisc:
    def fetch_source(self, *, lat: float, lon: float):  # pragma: no cover - not used
        raise RuntimeError("unused")

    def render_from_source(self, **kwargs):  # pragma: no cover - not used
        raise RuntimeError("unused")


def test_cloud_update_keeps_latest_pending_source_request() -> None:
    controller = CloudController(_DummyCloudDisc())
    controller._source_is_running = True
    controller._render_is_running = True
    controller._latest_source = object()

    controller.update(lat=35.0, lon=139.0, alt=45.0, az=180.0, radius_px=256, reason="timer")
    controller.update(lat=36.0, lon=140.0, alt=50.0, az=200.0, radius_px=256, reason="timer")

    assert controller._pending_source_request is not None
    assert controller._pending_source_request["lat"] == 36.0
    assert controller._pending_source_request["lon"] == 140.0
