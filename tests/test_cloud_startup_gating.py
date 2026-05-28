from __future__ import annotations

from datetime import timedelta
from types import SimpleNamespace

from zstarview.gui.window import SkyWindow


def test_apply_startup_delta_t_keeps_explicit_cloud_disable() -> None:
    dummy = SimpleNamespace()
    dummy.delta_t = timedelta(0)
    dummy._cloud_toggle_supported = True
    dummy._cloud_requested_enabled = False
    dummy._cloud_alpha_when_enabled = 0.2
    dummy.cloud_disc_alpha = 0.0
    dummy._geo_satellite_enabled = True
    dummy._satellite_toggle_supported = True
    dummy._satellite_requested_enabled = True
    dummy._satellite_opacity_when_enabled = 0.7
    dummy.satellite_opacity = 0.7
    dummy._aircraft_toggle_supported = True
    dummy._aircraft_requested_enabled = True
    dummy._aircraft_opacity_when_enabled = 1.0
    dummy.aircraft_opacity = 1.0
    dummy._action_toggle_clouds = None
    dummy._action_toggle_satellites = None
    dummy._action_toggle_aircraft = None
    dummy._cloud_gui_allowed = True
    dummy._satellite_gui_allowed = True
    dummy._aircraft_gui_allowed = True

    SkyWindow.apply_startup_delta_t(dummy, timedelta(0))

    assert dummy.cloud_disc_alpha == 0.0
