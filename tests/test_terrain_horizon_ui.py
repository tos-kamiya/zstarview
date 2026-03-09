from __future__ import annotations

from types import SimpleNamespace

import zstarview.ui.window as window_module
from zstarview.ui.window import SkyWindow
from zstarview.ui.window_inputs import prepare_window_user_options
from zstarview.ui.window_updates import SkyWindowUpdatesMixin


class _DummyAction:
    def __init__(self, checked: bool) -> None:
        self._checked = checked
        self._enabled = True

    def isChecked(self) -> bool:  # noqa: N802 - Qt naming
        return self._checked

    def setChecked(self, checked: bool) -> None:  # noqa: N802 - Qt naming
        self._checked = checked

    def isEnabled(self) -> bool:  # noqa: N802 - Qt naming
        return self._enabled

    def setEnabled(self, enabled: bool) -> None:  # noqa: N802 - Qt naming
        self._enabled = enabled


def test_prepare_window_user_options_normalizes_terrain_horizon_fields() -> None:
    options = prepare_window_user_options(
        terrain_horizon_opacity=1.5,
        ground_tint_opacity=1.5,
        terrain_horizon_gui_allowed=False,
    )

    assert options.terrain_horizon_opacity == 1.0
    assert options.ground_tint_opacity == 1.0
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
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(str(kwargs.get("reason")))
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_terrain_horizon(dummy)

    assert dummy.terrain_horizon_opacity == 0.25
    assert dummy._action_toggle_terrain_horizon.isChecked() is True
    assert calls == ["invalidate", "toggle-on", "update"]


def test_toggle_sky_disc_enables_gradient_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.sky_disc_alpha = 0.0
    dummy._sky_disc_alpha_when_enabled = 0.3
    dummy._action_toggle_sky_disc = _DummyAction(False)
    calls: list[str] = []
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.request_sky_data_update = lambda: calls.append("request")
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_sky_disc(dummy)

    assert dummy.sky_disc_alpha == 0.3
    assert dummy._action_toggle_sky_disc.isChecked() is True
    assert calls == ["invalidate", "request", "update"]


def test_toggle_sky_disc_switches_to_flat_disc_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.sky_disc_alpha = 0.3
    dummy._sky_disc_alpha_when_enabled = 0.3
    dummy._action_toggle_sky_disc = _DummyAction(True)
    calls: list[str] = []
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.request_sky_data_update = lambda: calls.append("request")
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_sky_disc(dummy)

    assert dummy.sky_disc_alpha == 0.0
    assert dummy._action_toggle_sky_disc.isChecked() is False
    assert calls == ["invalidate", "request", "update"]


def test_sync_view_altitude_actions_disables_raise_at_zenith() -> None:
    dummy = SimpleNamespace()
    dummy.viewer_data = SimpleNamespace(view_center=(90.0, 180.0))
    dummy._action_raise_view = _DummyAction(False)
    dummy._action_lower_view = _DummyAction(False)

    SkyWindow._sync_view_altitude_actions(dummy)

    assert dummy._action_raise_view.isEnabled() is False
    assert dummy._action_lower_view.isEnabled() is True


def test_sync_view_altitude_actions_disables_lower_at_horizon() -> None:
    dummy = SimpleNamespace()
    dummy.viewer_data = SimpleNamespace(view_center=(0.0, 180.0))
    dummy._action_raise_view = _DummyAction(False)
    dummy._action_lower_view = _DummyAction(False)

    SkyWindow._sync_view_altitude_actions(dummy)

    assert dummy._action_raise_view.isEnabled() is True
    assert dummy._action_lower_view.isEnabled() is False


def test_jump_to_search_target_keeps_negative_target_alt_for_highlight(monkeypatch) -> None:
    monkeypatch.setattr(window_module, "radec_to_altaz", lambda *_args, **_kwargs: (-12.5, 210.0))

    dummy = SimpleNamespace()
    dummy.viewer_data = SimpleNamespace(location=(35.0, 139.0), view_center=(20.0, 30.0))
    dummy.state = SimpleNamespace(
        jump_highlight_name=None,
        jump_highlight_altaz=None,
        jump_highlight_until_ms=0.0,
    )
    sync_calls: list[str] = []
    dummy._sync_view_altitude_actions = lambda: sync_calls.append("sync")
    dummy._current_time_obj = lambda: object()
    dummy._begin_interaction_mode = lambda: sync_calls.append("begin")
    dummy.request_sky_data_update = lambda: sync_calls.append("request")
    dummy.update = lambda: sync_calls.append("update")

    target = SimpleNamespace(label="Circlet", ra_hours=1.0, dec_deg=2.0)
    SkyWindow._jump_to_search_target(dummy, target)

    assert dummy.viewer_data.view_center == (0.0, 210.0)
    assert dummy.state.jump_highlight_name == "Circlet"
    assert dummy.state.jump_highlight_altaz == (-12.5, 210.0)
    assert dummy.state.jump_highlight_until_ms > 0.0
    assert sync_calls == ["sync", "begin", "request", "update"]


def test_rotate_view_in_orientation_mode_updates_render_center_without_full_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.viewer_data = SimpleNamespace(view_center=(20.0, 30.0))
    dummy.state = SimpleNamespace(render_view_center=(20.0, 30.0))
    calls: list[str] = []
    dummy._begin_orientation_interaction_mode = lambda: calls.append("begin-orientation")
    dummy._begin_interaction_mode = lambda: calls.append("begin")
    dummy._sync_view_altitude_actions = lambda: calls.append("sync")
    dummy.request_sky_data_update = lambda: calls.append("request")
    dummy.update = lambda: calls.append("update")

    SkyWindow._rotate_view(dummy, d_alt=5.0, d_az=15.0, interactive_orientation=True)

    assert dummy.viewer_data.view_center == (25.0, 45.0)
    assert dummy.state.render_view_center == (25.0, 45.0)
    assert calls == ["begin-orientation", "sync", "update"]


def test_end_orientation_interaction_mode_requests_full_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(orientation_interaction_mode=True)
    calls: list[str] = []
    dummy.request_sky_data_update = lambda: calls.append("sky")
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(str(kwargs.get("reason")))
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(str(kwargs.get("reason")))
    dummy.update = lambda: calls.append("update")

    SkyWindow._end_orientation_interaction_mode(dummy)

    assert dummy.state.orientation_interaction_mode is False
    assert calls == ["sky", "view-change-idle", "view-change-idle", "update"]
