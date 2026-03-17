from __future__ import annotations

from datetime import datetime, timedelta, timezone
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
        urban_outline_opacity=1.5,
        aircraft_opacity=1.5,
        ground_tint_opacity=1.5,
        sky_disc_gui_allowed=False,
        cloud_gui_allowed=False,
        aircraft_gui_allowed=False,
        terrain_horizon_gui_allowed=False,
        urban_outline_gui_allowed=False,
    )

    assert options.terrain_horizon_opacity == 1.0
    assert options.urban_outline_opacity == 1.0
    assert options.aircraft_opacity == 1.0
    assert options.ground_tint_opacity == 1.0
    assert options.sky_disc_gui_allowed is False
    assert options.cloud_gui_allowed is False
    assert options.aircraft_gui_allowed is False
    assert options.terrain_horizon_gui_allowed is False
    assert options.urban_outline_gui_allowed is False


def test_toggle_clouds_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_toggle_supported = True
    dummy._cloud_gui_allowed = False
    dummy._clouddisc = object()
    dummy.cloud_disc_alpha = 0.0
    dummy._cloud_alpha_when_enabled = 0.2
    dummy._action_toggle_clouds = _DummyAction(False)
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint"))

    SkyWindow.toggle_clouds(dummy)

    assert dummy.cloud_disc_alpha == 0.0
    assert dummy._action_toggle_clouds.isChecked() is False


def test_toggle_sky_disc_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._sky_disc_gui_allowed = False
    dummy.sky_disc_alpha = 0.0
    dummy._sky_disc_alpha_when_enabled = 0.3
    dummy._action_toggle_sky_disc = _DummyAction(False)
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint"))

    SkyWindow.toggle_sky_disc(dummy)

    assert dummy.sky_disc_alpha == 0.0
    assert dummy._action_toggle_sky_disc.isChecked() is False


def test_status_line_message_combines_cloud_and_terrain_segments() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_status_line = lambda: "Clouds [AUTO]: downloading"
    dummy._terrain_horizon_status_line = lambda: "Terrain horizon: loading DEM..."
    dummy._urban_outline_status_line = lambda: "Urban outline: downloading..."

    got = SkyWindowUpdatesMixin._status_line_message(dummy)

    assert got == "Clouds [AUTO]: downloading | Terrain horizon: loading DEM... | Urban outline: downloading..."


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


def test_toggle_aircraft_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._aircraft_toggle_supported = True
    dummy._aircraft_gui_allowed = False
    dummy.aircraft_opacity = 0.0
    dummy._aircraft_opacity_when_enabled = 1.0
    dummy._action_toggle_aircraft = _DummyAction(False)
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint"))

    SkyWindow.toggle_aircraft(dummy)

    assert dummy.aircraft_opacity == 0.0
    assert dummy._action_toggle_aircraft.isChecked() is False


def test_toggle_aircraft_uses_cached_state_without_fetch() -> None:
    now = datetime.now(timezone.utc)
    dummy = SimpleNamespace()
    dummy._aircraft_toggle_supported = True
    dummy._aircraft_gui_allowed = True
    dummy.aircraft_opacity = 0.0
    dummy._aircraft_opacity_when_enabled = 1.0
    dummy._action_toggle_aircraft = _DummyAction(False)
    dummy.aircraft_state = SimpleNamespace(
        snapshots=[object()],
        last_success_utc=now - timedelta(seconds=30),
    )
    calls: list[str] = []
    dummy._enable_aircraft_layer = lambda **kwargs: calls.append(str(kwargs.get("reason")))
    dummy._stop_aircraft_timers = lambda: calls.append("stop")
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_aircraft(dummy)

    assert dummy.aircraft_opacity == 1.0
    assert dummy._action_toggle_aircraft.isChecked() is True
    assert calls == ["toggle-on", "update"]


def test_start_background_aircraft_update_skips_when_layer_hidden() -> None:
    controller_calls: list[str] = []
    dummy = SimpleNamespace()
    dummy._is_shutting_down = False
    dummy.aircraft_opacity = 0.0
    dummy._aircraft_controller = SimpleNamespace(
        update=lambda **_kwargs: controller_calls.append("update") or True,
    )
    dummy.viewer_data = SimpleNamespace(location=(35.0, 135.0), observer_height_m=1.7)
    dummy._current_time_obj = lambda: "time-obj"

    started = SkyWindowUpdatesMixin.start_background_aircraft_update(dummy, reason="manual")

    assert started is False
    assert controller_calls == []


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


def test_toggle_urban_outline_enables_opacity_and_requests_background_update() -> None:
    dummy = SimpleNamespace()
    dummy._urban_outline_gui_allowed = True
    dummy.urban_outline_opacity = 0.0
    dummy._urban_outline_opacity_when_enabled = 0.2
    dummy.show_urban_outline_layer = False
    dummy._action_toggle_urban_outline = _DummyAction(False)
    calls: list[str] = []
    dummy.start_background_urban_outline_update = lambda **kwargs: calls.append(str(kwargs.get("reason")))
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_urban_outline(dummy)

    assert dummy.urban_outline_opacity == 0.2
    assert dummy.show_urban_outline_layer is True
    assert dummy._action_toggle_urban_outline.isChecked() is True
    assert calls == ["toggle-on", "update"]


def test_toggle_urban_outline_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._urban_outline_gui_allowed = False
    dummy.urban_outline_opacity = 0.0
    dummy._urban_outline_opacity_when_enabled = 0.2
    dummy.show_urban_outline_layer = False
    dummy._action_toggle_urban_outline = _DummyAction(False)
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint"))

    SkyWindow.toggle_urban_outline(dummy)

    assert dummy.urban_outline_opacity == 0.0
    assert dummy.show_urban_outline_layer is False
    assert dummy._action_toggle_urban_outline.isChecked() is False


def test_toggle_sky_disc_enables_gradient_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy._sky_disc_gui_allowed = True
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
    dummy._sky_disc_gui_allowed = True
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
    dummy._begin_viewport_interaction_mode = lambda: calls.append("begin-viewport")
    dummy._begin_interaction_mode = lambda: calls.append("begin")
    dummy._sync_view_altitude_actions = lambda: calls.append("sync")
    dummy._update_viewport_interaction_stars = lambda: calls.append("bright-stars")
    dummy.request_sky_data_update = lambda: calls.append("request")
    dummy.update = lambda: calls.append("update")

    SkyWindow._rotate_view(dummy, d_alt=5.0, d_az=15.0, interactive_viewport=True)

    assert dummy.viewer_data.view_center == (25.0, 45.0)
    assert dummy.state.render_view_center == (25.0, 45.0)
    assert calls == ["begin-viewport", "sync", "bright-stars", "update"]


def test_end_viewport_interaction_mode_requests_full_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(viewport_interaction_mode=True, viewport_interaction_stars=object())
    calls: list[str] = []
    dummy.request_sky_data_update = lambda: calls.append("sky")
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(str(kwargs.get("reason")))
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(str(kwargs.get("reason")))
    dummy.update = lambda: calls.append("update")

    SkyWindow._end_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    assert calls == ["sky", "view-change-idle", "view-change-idle", "update"]


def test_show_menu_syncs_actions_before_opening_menu() -> None:
    calls: list[str] = []

    class _DummyButton:
        def height(self) -> int:
            return 30

        def mapToGlobal(self, point) -> tuple[int, int]:
            calls.append(f"map:{point.x()},{point.y()}")
            return (100, 200)

    class _DummyMenu:
        def exec(self, pos) -> None:
            calls.append(f"exec:{pos}")

    dummy = SimpleNamespace()
    dummy.menu_button = _DummyButton()
    dummy.menu = _DummyMenu()
    dummy._sync_view_altitude_actions = lambda: calls.append("sync")

    SkyWindow.show_menu(dummy)

    assert calls == [
        "sync",
        "map:0,30",
        "exec:(100, 200)",
    ]


def test_menu_button_style_sheet_uses_opaque_background_for_night_preset() -> None:
    dummy = SimpleNamespace(visual_preset="night")

    style = SkyWindow._menu_button_style_sheet(dummy)

    assert "background-color: #0a0c10;" in style
    assert "border-radius: 6px;" in style


def test_menu_button_style_sheet_uses_light_background_for_day_preset() -> None:
    dummy = SimpleNamespace(visual_preset="day")

    style = SkyWindow._menu_button_style_sheet(dummy)

    assert "background-color: #f0f8ff;" in style


def test_update_viewport_interaction_stars_uses_bright_limit(monkeypatch) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        window_module,
        "calculate_visible_stars",
        lambda catalog, lat, lon, observer_height_m, time_obj, view_center, max_vmag: (
            captured.update(
                {
                    "catalog": catalog,
                    "lat": lat,
                    "lon": lon,
                    "observer_height_m": observer_height_m,
                    "time_obj": time_obj,
                    "view_center": view_center,
                    "max_vmag": max_vmag,
                }
            )
            or {"name": []},
            object(),
        ),
    )

    dummy = SimpleNamespace()
    dummy.star_catalog_lod6_np = object()
    dummy.viewer_data = SimpleNamespace(location=(35.0, 139.0), observer_height_m=12.0)
    dummy.state = SimpleNamespace(
        celestial_data=object(),
        render_view_center=(55.0, 210.0),
        viewport_interaction_stars=None,
    )
    dummy._current_time_obj = lambda: "time"

    SkyWindow._update_viewport_interaction_stars(dummy)

    assert dummy.state.viewport_interaction_stars == {"name": []}
    assert captured == {
        "catalog": dummy.star_catalog_lod6_np,
        "lat": 35.0,
        "lon": 139.0,
        "observer_height_m": 12.0,
        "time_obj": "time",
        "view_center": (55.0, 210.0),
        "max_vmag": 4.0,
    }


def test_compute_star_render_upscale_factor_matches_downsampled_star_surface() -> None:
    factor = SkyWindow.compute_star_render_upscale_factor(
        disc_width_px=2400,
        expected_width_px=600,
    )

    assert factor == 2.0
