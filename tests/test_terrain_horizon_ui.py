from __future__ import annotations

from types import SimpleNamespace

from zstarview.ui.window import SkyWindow
from zstarview.ui.window_inputs import prepare_window_user_options
from zstarview.ui.window_updates import SkyWindowUpdatesMixin


class _DummyAction:
    def __init__(self, checked: bool) -> None:
        self._checked = checked

    def isChecked(self) -> bool:  # noqa: N802 - Qt naming
        return self._checked

    def setChecked(self, checked: bool) -> None:  # noqa: N802 - Qt naming
        self._checked = checked


def test_prepare_window_user_options_normalizes_terrain_horizon_fields() -> None:
    options = prepare_window_user_options(
        terrain_horizon_opacity=1.5,
        terrain_horizon_gui_allowed=False,
    )

    assert options.terrain_horizon_opacity == 1.0
    assert options.terrain_horizon_gui_allowed is False


def test_status_line_message_combines_cloud_and_terrain_segments() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_status_line = lambda: "Clouds [AUTO]: downloading"
    dummy._terrain_horizon_status_line = lambda: "Terrain horizon: loading DEM..."

    got = SkyWindowUpdatesMixin._status_line_message(dummy)

    assert got == "Clouds [AUTO]: downloading | Terrain horizon: loading DEM..."


def test_toggle_terrain_horizon_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._terrain_horizon_gui_allowed = False
    dummy.terrain_horizon_opacity = 0.0
    dummy._terrain_horizon_opacity_when_enabled = 0.25
    dummy._action_toggle_terrain_horizon = _DummyAction(False)
    dummy.start_background_terrain_horizon_update = lambda **_kwargs: (_ for _ in ()).throw(AssertionError("should not start"))
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint"))

    SkyWindow.toggle_terrain_horizon(dummy)

    assert dummy.terrain_horizon_opacity == 0.0
    assert dummy._action_toggle_terrain_horizon.isChecked() is False


def test_toggle_terrain_horizon_enables_opacity_and_requests_background_update() -> None:
    dummy = SimpleNamespace()
    dummy._terrain_horizon_gui_allowed = True
    dummy.terrain_horizon_opacity = 0.0
    dummy._terrain_horizon_opacity_when_enabled = 0.25
    dummy._action_toggle_terrain_horizon = _DummyAction(False)
    calls: list[str] = []
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(str(kwargs.get("reason")))
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_terrain_horizon(dummy)

    assert dummy.terrain_horizon_opacity == 0.25
    assert dummy._action_toggle_terrain_horizon.isChecked() is True
    assert calls == ["toggle-on", "update"]
