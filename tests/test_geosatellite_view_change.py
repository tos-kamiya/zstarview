from __future__ import annotations

from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
from typing import Any

import numpy as np

from zstarview.gui.geosatellite_state import GeoSatelliteState
from zstarview.gui.window_render_cache import SkyWindowRenderCacheMixin
from zstarview.gui.window_state import SkyWindowState
from zstarview.gui.window_update_cloud import SkyWindowCloudUpdatesMixin

T = datetime(2026, 11, 21, 0, 0, tzinfo=timezone.utc)


def _owner(view_center: tuple[float, float]) -> Any:
    owner: Any = SimpleNamespace(
        geosatellite_state=GeoSatelliteState(
            image=np.zeros((1, 1, 4)),
            image_hidden_until_reprojected=True,
        ),
        state=SkyWindowState(render_view_center=view_center),
        viewer_data=SimpleNamespace(view_center=view_center),
        _disc_generation=0,
        _is_shutting_down=False,
        _compositor=SimpleNamespace(invalidate=lambda: None),
        _safe_request_cloud_repaint=lambda: None,
    )
    owner._geo_satellite_payload_matches_view = lambda payload: (
        SkyWindowCloudUpdatesMixin._geo_satellite_payload_matches_view(owner, payload)
    )
    return owner


def _payload(alt: float, az: float) -> dict[str, object]:
    return {
        "image": np.ones((1, 1, 4)),
        "alt": alt,
        "az": az,
        "captured_at_utc": T - timedelta(minutes=10),
        "refreshed_at_utc": T,
    }


def test_geo_result_for_current_view_shows_image() -> None:
    owner = _owner((20.0, 0.0))

    SkyWindowCloudUpdatesMixin._on_geosatellite_ready(owner, _payload(20.0, 360.0))

    assert owner.geosatellite_state.image_hidden_until_reprojected is False


def test_geo_result_for_previous_view_stays_hidden() -> None:
    owner = _owner((20.0, 90.0))

    SkyWindowCloudUpdatesMixin._on_geosatellite_ready(owner, _payload(20.0, 85.0))

    assert owner.geosatellite_state.image_hidden_until_reprojected is True


def test_geo_failure_falls_back_to_previous_image() -> None:
    owner = _owner((20.0, 90.0))
    image = owner.geosatellite_state.image

    SkyWindowCloudUpdatesMixin._on_geosatellite_failed(owner, {"banner": "x"})

    assert owner.geosatellite_state.image_hidden_until_reprojected is False
    assert owner.geosatellite_state.image is image


def test_hidden_geo_image_is_not_rendered() -> None:
    geo_state = GeoSatelliteState(image=np.zeros((1, 1, 4)))
    owner: Any = SimpleNamespace(_render_cloud_state=lambda: geo_state)

    assert SkyWindowRenderCacheMixin._render_cloud_image(owner) is geo_state.image
    geo_state.image_hidden_until_reprojected = True
    assert SkyWindowRenderCacheMixin._render_cloud_image(owner) is None
    geo_state.clear_data()
    assert geo_state.image_hidden_until_reprojected is False
