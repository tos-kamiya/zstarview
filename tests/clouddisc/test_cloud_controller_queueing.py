from __future__ import annotations

import threading
import time

from zstarview.clouddisc.types import DownloadCancelledError
from zstarview.gui.cloud_controller import CloudController


class _DummyCloudDisc:
    def __init__(self) -> None:
        self.cfg = object()

    def fetch_source(
        self, *, lat: float, lon: float, abort_event=None
    ):  # pragma: no cover - not used
        raise RuntimeError("unused")


def test_cloud_update_keeps_latest_pending_source_request() -> None:
    controller = CloudController(_DummyCloudDisc())
    controller._source_is_running = True

    controller.update(
        lat=35.0,
        lon=139.0,
        reason="timer",
    )
    controller.update(
        lat=36.0,
        lon=140.0,
        reason="timer",
    )

    assert controller._pending_source_request is not None
    assert controller._pending_source_request["lat"] == 36.0
    assert controller._pending_source_request["lon"] == 140.0


def test_cloud_update_skips_duplicate_pending_source_request() -> None:
    controller = CloudController(_DummyCloudDisc())
    controller._source_is_running = True
    controller._active_source_request_key = controller._source_request_key(
        lat=35.0,
        lon=139.0,
    )

    started = controller.update_source(
        lat=35.0,
        lon=139.0,
        reason="scheduler",
    )

    assert started is False
    assert controller._latest_source_request_id == 0
    assert controller._pending_source_request is None


def test_source_completion_keeps_pending_render_queued_when_render_is_running(
    monkeypatch,
) -> None:
    controller = CloudController(_DummyCloudDisc())
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
    controller._latest_source_request_id = 1

    monkeypatch.setattr(
        "zstarview.gui.cloud_controller.run_cloud_source_worker_process",
        lambda *args, **kwargs: object(),
    )

    controller._run_source_update(lat=35.0, lon=139.0, reason="manual", request_id=1)

    assert controller._pending_render_request is not None
    assert controller._pending_render_request["request_id"] == 10


def test_cloud_update_defers_render_until_source_is_ready() -> None:
    class _SourceAndRenderCloudDisc(_DummyCloudDisc):
        def fetch_source(self, *, lat: float, lon: float, abort_event=None):
            return object()

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
        reason="initial",
    )

    assert [label for label, _kwargs in calls] == ["source"]
    assert controller._pending_render_request is None


def test_cloud_update_skips_duplicate_pending_render_request() -> None:
    controller = CloudController(_DummyCloudDisc())
    source = type("Source", (), {"source_key": object()})()
    controller._latest_source = source
    controller._render_is_running = True
    controller._active_render_request_key = controller._render_request_key(
        source_key=source.source_key,
        alt=45.0,
        az=180.0,
        radius_px=256,
        content_fov_deg=90.0,
        render_generation=0,
    )

    started = controller.update_render(
        lat=35.0,
        lon=139.0,
        alt=45.0,
        az=180.0,
        radius_px=256,
        content_fov_deg=90.0,
        reason="scheduler",
    )

    assert started is False
    assert controller._latest_render_request_id == 0
    assert controller._pending_render_request is None


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


def test_cloud_shutdown_cancels_active_source_download(monkeypatch) -> None:
    started = threading.Event()

    def fake_run_cloud_source_worker_process(*_args, abort_event=None, **_kwargs):
        assert abort_event is not None
        started.set()
        while not abort_event.is_set():
            time.sleep(0.01)
        raise DownloadCancelledError("Cancelled while downloading")

    controller = CloudController(_DummyCloudDisc())
    controller._latest_source = None
    monkeypatch.setattr(
        "zstarview.gui.cloud_controller.run_cloud_source_worker_process",
        fake_run_cloud_source_worker_process,
    )

    controller.update(
        lat=35.0,
        lon=139.0,
        reason="initial",
    )

    assert started.wait(timeout=1.0)

    shutdown_done = threading.Event()

    def run_shutdown() -> None:
        controller.shutdown(wait_timeout_s=1.0)
        shutdown_done.set()

    shutdown_thread = threading.Thread(target=run_shutdown)
    shutdown_thread.start()

    assert shutdown_done.wait(timeout=1.0)
    shutdown_thread.join(timeout=1.0)
