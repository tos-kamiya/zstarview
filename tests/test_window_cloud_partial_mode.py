from __future__ import annotations

from datetime import datetime, timezone
from types import SimpleNamespace

import numpy as np

from zstarview.gui.cloud_state import CloudImageState
from zstarview.gui.window import SkyWindow
from zstarview.gui.window_state import SkyWindowState
from zstarview.render.qt_image import np_rgba_to_qimage


class _DummyCompositor:
    def __init__(self) -> None:
        self.invalidated = False

    def invalidate(self) -> None:
        self.invalidated = True


def test_on_cloud_ready_sets_partial_fields() -> None:
    dummy = SimpleNamespace()
    dummy.cloud_state = CloudImageState(banner_text="Clouds: downloading…")
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy._compositor = _DummyCompositor()
    dummy._disc_generation = 0
    dummy._is_shutting_down = False
    repaint_calls: list[str] = []
    dummy._safe_request_cloud_repaint = lambda: repaint_calls.append("repaint")
    dummy.start_background_cloud_update = lambda **_kwargs: repaint_calls.append("restart")

    image = np.zeros((8, 8, 4), dtype=np.uint8)
    image[..., :3] = 255
    image[..., 3] = 100
    missing = np.zeros((8, 8, 4), dtype=np.uint8)
    missing[2:4, 2:4, 3] = 255
    meta = SimpleNamespace(
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=datetime(2026, 3, 5, 1, 30, tzinfo=timezone.utc),
    )
    payload = {
        "image": np_rgba_to_qimage(image),
        "meta": meta,
        "az": 180.0,
        "time_utc": datetime(2026, 3, 5, 1, 31, tzinfo=timezone.utc),
        "stripe_density": None,
        "missing_mask": np_rgba_to_qimage(missing),
        "coverage_ratio": 0.75,
        "source_key": None,
        "request_id": 123,
        "render_generation": 0,
    }

    SkyWindow._on_cloud_ready(dummy, payload)

    assert dummy._compositor.invalidated is True
    assert dummy.cloud_state.missing_mask is not None
    assert dummy.cloud_state.coverage_ratio == 0.75
    assert repaint_calls == ["repaint"]


def test_on_cloud_ready_discards_stale_generation_and_restarts_render() -> None:
    dummy = SimpleNamespace()
    dummy.cloud_state = CloudImageState(banner_text="Clouds: downloading…")
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy._compositor = _DummyCompositor()
    dummy._disc_generation = 4
    dummy._is_shutting_down = False
    calls: list[str] = []
    dummy._safe_request_cloud_repaint = lambda: calls.append("repaint")
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(str(kwargs.get("reason")))

    payload = {
        "image": np_rgba_to_qimage(np.zeros((8, 8, 4), dtype=np.uint8)),
        "meta": None,
        "az": 180.0,
        "time_utc": datetime(2026, 3, 5, 1, 31, tzinfo=timezone.utc),
        "stripe_density": None,
        "missing_mask": None,
        "coverage_ratio": 1.0,
        "source_key": None,
        "request_id": 123,
        "render_generation": 3,
    }

    SkyWindow._on_cloud_ready(dummy, payload)

    assert dummy._compositor.invalidated is False
    assert dummy.cloud_state.image is None
    assert calls == ["stale-render"]
