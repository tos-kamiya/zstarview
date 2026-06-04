from __future__ import annotations

from datetime import timedelta
from types import SimpleNamespace

from zstarview.gui.window import SkyWindow


class _ActionStub:
    def __init__(self) -> None:
        self.enabled = True
        self.checked = False

    def setEnabled(self, value: bool) -> None:  # noqa: N802 - Qt naming
        self.enabled = bool(value)

    def setChecked(self, value: bool) -> None:  # noqa: N802 - Qt naming
        self.checked = bool(value)


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


def test_geo_satellite_action_is_disabled_outside_supported_band() -> None:
    action = _ActionStub()
    dummy = SimpleNamespace(
        viewer_data=SimpleNamespace(location=(10.0, 0.0)),
        _geo_satellite_enabled=True,
        _geo_satellite_location_resolved=True,
        _geo_satellite_toggle_supported=lambda: False,
        _action_toggle_geo_satellite=action,
    )

    SkyWindow._sync_geo_satellite_action_state(dummy)

    assert action.enabled is False
    assert action.checked is False


def test_geo_satellite_action_is_enabled_inside_supported_band() -> None:
    action = _ActionStub()
    dummy = SimpleNamespace(
        viewer_data=SimpleNamespace(location=(51.5, -0.1)),
        _geo_satellite_enabled=True,
        _geo_satellite_location_resolved=True,
        _geo_satellite_toggle_supported=lambda: True,
        _action_toggle_geo_satellite=action,
    )

    SkyWindow._sync_geo_satellite_action_state(dummy)

    assert action.enabled is True
    assert action.checked is True


def test_geo_satellite_action_stays_enabled_before_location_is_resolved() -> None:
    action = _ActionStub()
    dummy = SimpleNamespace(
        viewer_data=SimpleNamespace(location=(0.0, 0.0)),
        _geo_satellite_enabled=True,
        _geo_satellite_location_resolved=False,
        _geo_satellite_toggle_supported=lambda: True,
        _action_toggle_geo_satellite=action,
    )

    SkyWindow._sync_geo_satellite_action_state(dummy)

    assert action.enabled is True
    assert action.checked is True
