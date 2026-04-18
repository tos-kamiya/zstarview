from __future__ import annotations

from zstarview.gui.cloud_controller import CloudController


class _DummyCloudDisc:
    def fetch_source(self, *, lat: float, lon: float):  # pragma: no cover - not used
        raise RuntimeError("unused")

    def render_from_source_with_coverage(self, **kwargs):  # pragma: no cover - not used
        raise RuntimeError("unused")


def test_cloud_update_keeps_latest_pending_source_request() -> None:
    controller = CloudController(_DummyCloudDisc())
    controller._source_is_running = True
    controller._render_is_running = True
    controller._latest_source = object()

    controller.update(
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="timer",
    )
    controller.update(
        lat=36.0,
        lon=140.0,
        alt=50.0,
        az=200.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="timer",
    )

    assert controller._pending_source_request is not None
    assert controller._pending_source_request["lat"] == 36.0
    assert controller._pending_source_request["lon"] == 140.0


def test_source_completion_queues_rerender_when_render_is_running() -> None:
    class _SourceOnlyCloudDisc(_DummyCloudDisc):
        def fetch_source(self, *, lat: float, lon: float):
            return object()

    controller = CloudController(_SourceOnlyCloudDisc())
    controller._render_is_running = True
    controller._last_render_request = {
        "lat": 35.0,
        "lon": 139.0,
        "alt": 45.0,
        "az": 180.0,
        "radius_px": 256,
        "content_fov_deg": 90.0,
        "reason": "manual",
        "request_id": 10,
    }

    controller._run_source_update(lat=35.0, lon=139.0, reason="manual")

    assert controller._pending_render_request is not None
    assert controller._pending_render_request["request_id"] == 10
