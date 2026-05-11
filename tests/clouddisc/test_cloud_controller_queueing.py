from __future__ import annotations

import threading
import time

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


def test_source_completion_keeps_pending_render_queued_when_render_is_running() -> None:
    class _SourceOnlyCloudDisc(_DummyCloudDisc):
        def fetch_source(self, *, lat: float, lon: float):
            return object()

    controller = CloudController(_SourceOnlyCloudDisc())
    controller._render_is_running = True
    controller._pending_render_request = {
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


def test_cloud_update_defers_render_until_source_is_ready() -> None:
    class _SourceAndRenderCloudDisc(_DummyCloudDisc):
        def fetch_source(self, *, lat: float, lon: float):
            return object()

        def render_from_source_with_coverage(self, **kwargs):
            raise AssertionError("render should not start before source finishes")

    controller = CloudController(_SourceAndRenderCloudDisc())
    controller._cleanup_counter = 1
    controller._latest_source = None
    calls: list[tuple[str, dict[str, object]]] = []

    def fake_spawn_worker(*, target, kwargs, label):
        calls.append((label, dict(kwargs)))

    controller._spawn_worker = fake_spawn_worker  # type: ignore[method-assign]

    controller.update(
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="initial",
    )

    assert [label for label, _kwargs in calls] == ["source"]
    assert controller._pending_render_request is not None
    assert controller._pending_render_request["request_id"] == 1


def test_cloud_shutdown_waits_for_active_worker_threads(monkeypatch) -> None:
    controller = CloudController(_DummyCloudDisc())
    controller._latest_source = object()
    controller._render_is_running = True

    worker_started = threading.Event()
    worker_release = threading.Event()

    def fake_run_source_update(**kwargs):
        worker_started.set()
        worker_release.wait(timeout=2.0)

    monkeypatch.setattr(controller, "_run_source_update", fake_run_source_update)

    controller.update(
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="manual",
    )

    assert worker_started.wait(timeout=1.0)

    shutdown_done = threading.Event()

    def run_shutdown() -> None:
        controller.shutdown(wait_timeout_s=1.0)
        shutdown_done.set()

    shutdown_thread = threading.Thread(target=run_shutdown)
    shutdown_thread.start()

    time.sleep(0.05)
    assert not shutdown_done.is_set()

    worker_release.set()
    assert shutdown_done.wait(timeout=1.0)
    shutdown_thread.join(timeout=1.0)

    assert not controller._active_workers
