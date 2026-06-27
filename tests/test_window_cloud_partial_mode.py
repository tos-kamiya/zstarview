from __future__ import annotations

from datetime import datetime, timezone
from types import SimpleNamespace

import numpy as np

from zstarview.gui.cloud_state import CloudImageState
from zstarview.gui.geosatellite_state import GeoSatelliteState
from zstarview.gui.window import SkyWindow
from zstarview.gui.window_state import SkyWindowState


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
    missing = np.zeros((8, 8), dtype=np.uint8)
    missing[2:4, 2:4] = 255
    meta = SimpleNamespace(
        satellite="HIMAWARI",
        product="ISatSS-B13",
        time_utc=datetime(2026, 3, 5, 1, 30, tzinfo=timezone.utc),
    )
    payload = {
        "image": image,
        "meta": meta,
        "az": 180.0,
        "time_utc": datetime(2026, 3, 5, 1, 31, tzinfo=timezone.utc),
        "cloud_amount_field": None,
        "missing_mask": missing,
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
        "image": np.zeros((8, 8, 4), dtype=np.uint8),
        "meta": None,
        "az": 180.0,
        "time_utc": datetime(2026, 3, 5, 1, 31, tzinfo=timezone.utc),
        "cloud_amount_field": None,
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


def test_on_cloud_source_ready_schedules_projection_without_immediate_repaint() -> None:
    dummy = SimpleNamespace()
    dummy.cloud_state = CloudImageState()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy._disc_generation = 0
    dummy._is_shutting_down = False
    dummy._geo_satellite_enabled = False
    dummy._geosatellite_controller = None
    dummy._cloud_controller = SimpleNamespace()
    dummy.viewer_data = SimpleNamespace(location=(35.0, 139.0))
    dummy._cloud_layer_enabled = lambda: True
    dummy.request_client_update = lambda: (_ for _ in ()).throw(
        AssertionError("request_client_update should not be called")
    )
    projection_calls: list[str] = []
    dummy.reproject_cloud_overlay = lambda **kwargs: projection_calls.append(
        str(kwargs.get("reason"))
    )
    dummy.state.cloud_next_refresh_utc = None
    dummy.state.cloud_projection_next_refresh_utc = None
    dummy.state.interaction_mode = False
    dummy.state.cloud_repaint_deferred = False

    refreshed_at = datetime(2026, 3, 5, 1, 30, tzinfo=timezone.utc)
    payload = {
        "satellite": "HIMAWARI",
        "source_key": object(),
        "refreshed_at_utc": refreshed_at,
        "banner": "",
        "altaz_grid": object(),
    }

    SkyWindow._on_cloud_source_ready(dummy, payload)

    assert dummy.cloud_state.source_refreshed_at_utc == refreshed_at
    assert dummy.cloud_state.current_satellite == "HIMAWARI"
    assert dummy.cloud_state.altaz_grid is payload["altaz_grid"]
    assert projection_calls == ["source-ready"]
    assert dummy.state.cloud_projection_next_refresh_utc is not None


def test_on_geosatellite_source_ready_schedules_projection_without_immediate_repaint() -> None:
    dummy = SimpleNamespace()
    dummy.geosatellite_state = GeoSatelliteState()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy._disc_generation = 0
    dummy._is_shutting_down = False
    dummy._geo_satellite_enabled = True
    dummy._geosatellite_controller = object()
    dummy.viewer_data = SimpleNamespace(location=(35.0, 139.0))
    dummy._cloud_controller = SimpleNamespace()
    dummy.request_client_update = lambda: (_ for _ in ()).throw(
        AssertionError("request_client_update should not be called")
    )
    dummy.state.cloud_projection_next_refresh_utc = datetime(2026, 3, 5, 1, 29, tzinfo=timezone.utc)
    dummy.state.cloud_next_refresh_utc = None
    dummy.state.interaction_mode = False
    dummy.state.cloud_repaint_deferred = False

    refreshed_at = datetime(2026, 3, 5, 1, 30, tzinfo=timezone.utc)
    payload = {
        "refreshed_at_utc": refreshed_at,
        "banner": "",
    }

    SkyWindow._on_geosatellite_source_ready(dummy, payload)

    assert dummy.geosatellite_state.source_refreshed_at_utc == refreshed_at
    assert dummy.geosatellite_state.current_source == "Geo-sat"
    assert dummy.geosatellite_state.banner_text == "Geo-sat + Projecting"
    assert dummy.state.cloud_projection_next_refresh_utc is None
