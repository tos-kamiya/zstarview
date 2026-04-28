from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import zstarview.gui.terrain_controller as terrain_controller_module
import zstarview.gui.window as window_module
from zstarview.terrain.dem import COPERNICUS_DEM_BUCKET
from zstarview.gui.window import SkyWindow
from zstarview.gui.terrain_controller import TerrainHorizonController
from zstarview.gui.window_inputs import prepare_window_user_options
from zstarview.gui.window_updates import SkyWindowUpdatesMixin


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


class _DummySignal:
    def connect(self, _callback) -> None:
        return None


class _DummyMenuAction(_DummyAction):
    def __init__(self, text: str, _parent=None, *, separator: bool = False) -> None:
        super().__init__(False)
        self.text = text
        self.separator = separator
        self.triggered = _DummySignal()
        self.shortcut = None
        self.shortcut_context = None
        self.checkable = False

    def setShortcut(self, shortcut) -> None:  # noqa: N802 - Qt naming
        self.shortcut = shortcut

    def setShortcutContext(self, context) -> None:  # noqa: N802 - Qt naming
        self.shortcut_context = context

    def setCheckable(self, checkable: bool) -> None:  # noqa: N802 - Qt naming
        self.checkable = checkable


class _DummyMenu:
    def __init__(self, title: str, _parent=None) -> None:
        self.title = title
        self.entries: list[object] = []

    def addMenu(self, menu):  # noqa: N802 - Qt naming
        self.entries.append(menu)
        return menu

    def addAction(self, action):  # noqa: N802 - Qt naming
        if isinstance(action, str):
            action = _DummyMenuAction(action)
        self.entries.append(action)
        return action

    def addSeparator(self):  # noqa: N802 - Qt naming
        action = _DummyMenuAction("", separator=True)
        self.entries.append(action)
        return action


class _DummyAircraftState:
    def __init__(self) -> None:
        self.banner_text = None
        self.last_result = None

    def set_result(self, snapshots, *, overlay_points, bbox, refreshed_at_utc) -> None:
        self.last_result = {
            "snapshots": snapshots,
            "overlay_points": overlay_points,
            "bbox": bbox,
            "refreshed_at_utc": refreshed_at_utc,
        }

    def set_banner(self, text: str) -> None:
        self.banner_text = text


class _DummyImage:
    def __init__(self) -> None:
        self.saved_paths: list[Path] = []
        self.saved_formats: list[str | None] = []

    def save(self, path: str, image_format: str | None = None) -> bool:
        output_path = Path(path)
        output_path.write_bytes(b"debug-png")
        self.saved_paths.append(output_path)
        self.saved_formats.append(image_format)
        return True


def _noop_request_client_update() -> None:
    return None


class _DummyPaintEvent:
    def __init__(self) -> None:
        self.accepted = False

    def accept(self) -> None:
        self.accepted = True


def _install_menu_action_helpers(dummy: SimpleNamespace, added_actions: list[object]) -> None:
    def _add_menu_action(menu, text, *, shortcut=None, enabled=True, triggered=None):
        action = window_module.QAction(text, dummy)
        menu.addAction(action)
        if shortcut is not None:
            action.setShortcut(shortcut)
            action.setShortcutContext(window_module.Qt.ShortcutContext.WindowShortcut)
            dummy.addAction(action)
        action.setEnabled(enabled)
        if triggered is not None:
            action.triggered.connect(triggered)
        return action

    def _add_checkable_menu_action(
        menu, text, *, checked=False, shortcut=None, enabled=True, triggered=None
    ):
        action = window_module.QAction(text, dummy)
        menu.addAction(action)
        action.setCheckable(True)
        action.setChecked(checked)
        if shortcut is not None:
            action.setShortcut(shortcut)
            action.setShortcutContext(window_module.Qt.ShortcutContext.WindowShortcut)
            dummy.addAction(action)
        action.setEnabled(enabled)
        if triggered is not None:
            action.triggered.connect(triggered)
        return action

    dummy._add_menu_action = _add_menu_action
    dummy._add_checkable_menu_action = _add_checkable_menu_action


def test_paint_event_does_not_delegate_to_qmainwindow(monkeypatch) -> None:
    called = []

    def _fail(*_args, **_kwargs):
        called.append("called")
        raise AssertionError("QMainWindow.paintEvent should not be invoked")

    monkeypatch.setattr(window_module.QMainWindow, "paintEvent", _fail)

    event = _DummyPaintEvent()
    SkyWindow.paintEvent(SimpleNamespace(), event)

    assert event.accepted is True
    assert called == []


def test_prepare_window_user_options_normalizes_terrain_horizon_fields() -> None:
    options = prepare_window_user_options(
        sky_disc_alpha=0.17,
        cloud_disc_alpha=0.075,
        satellite_opacity=0.7,
        terrain_horizon_opacity=1.5,
        earth_guide_opacity=1.5,
        urban_outline_opacity=1.5,
        aircraft_opacity=0.4,
        ground_tint_opacity=0.2,
        enlarge_moon=False,
        bright_bodies_mode="outline",
        star_base_radius=4.0,
        vmag_limit=6.0,
        visual_preset="night",
        star_visibility_boost=1.0,
        visibility_boost=2.0,
        show_dso_initial=None,
        show_asterisms_initial=None,
        show_guidelines_initial=None,
        observation_info_mode=None,
        sky_disc_gui_allowed=False,
        cloud_gui_allowed=False,
        satellite_gui_allowed=True,
        aircraft_gui_allowed=False,
        terrain_horizon_gui_allowed=False,
        earth_guide_gui_allowed=False,
        urban_outline_gui_allowed=False,
    )

    assert options.terrain_horizon_opacity == 1.0
    assert options.earth_guide_opacity == 1.0
    assert options.urban_outline_opacity == 1.0
    assert options.sky_disc_alpha == pytest.approx(0.2125)
    assert options.cloud_disc_alpha == pytest.approx(0.09375)
    assert options.satellite_opacity == 0.875
    assert options.aircraft_opacity == 0.5
    assert options.ground_tint_opacity == 0.25
    assert options.sky_disc_gui_allowed is False
    assert options.cloud_gui_allowed is False
    assert options.aircraft_gui_allowed is False
    assert options.terrain_horizon_gui_allowed is False
    assert options.earth_guide_gui_allowed is False
    assert options.urban_outline_gui_allowed is False
    assert options.observation_info_mode is None


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
    dummy._sky_disc_alpha_when_enabled = 0.17
    dummy._action_toggle_sky_disc = _DummyAction(False)
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint"))

    SkyWindow.toggle_sky_disc(dummy)

    assert dummy.sky_disc_alpha == 0.0
    assert dummy._action_toggle_sky_disc.isChecked() is False


def test_build_window_menu_flattens_file_actions_for_frameless(monkeypatch) -> None:
    monkeypatch.setattr(window_module, "QMenu", _DummyMenu)
    monkeypatch.setattr(window_module, "QAction", _DummyMenuAction)

    added_actions: list[object] = []
    dummy = SimpleNamespace(
        _frameless_window=True,
        state=SimpleNamespace(rotation_step=5.0),
        enlarge_moon=False,
        show_dso=False,
        dso_catalog_np=None,
        show_asterisms=False,
        show_guidelines=True,
        show_observation_info=True,
        observation_info_mode="auto",
        observation_info_pinned=False,
        sky_disc_alpha=0.2,
        cloud_disc_alpha=0.2,
        satellite_opacity=0.5,
        aircraft_opacity=0.5,
        terrain_horizon_opacity=0.1,
        earth_guide_opacity=0.1,
        urban_outline_opacity=0.2,
        vmag_limit=6.0,
        _rotate_view=lambda **_kwargs: None,
        _open_named_star_jump_dialog=lambda: None,
        _open_named_star_search_dialog=lambda: None,
        _open_place_search_dialog=lambda: None,
        toggle_enlarge_moon=lambda: None,
        toggle_dso=lambda: None,
        toggle_asterisms=lambda: None,
        toggle_guidelines=lambda: None,
        toggle_observation_info=lambda: None,
        toggle_sky_disc=lambda: None,
        toggle_clouds=lambda: None,
        toggle_satellites=lambda: None,
        toggle_aircraft=lambda: None,
        toggle_terrain_horizon=lambda: None,
        toggle_earth_guide=lambda: None,
        toggle_urban_outline=lambda: None,
        toggle_fullscreen=lambda: None,
        square_client_area=lambda: None,
        addAction=lambda action: added_actions.append(action),
        _vmag_limit_menu_text=lambda: "Vmag limit 6.0",
    )
    _install_menu_action_helpers(dummy, added_actions)

    SkyWindow._build_window_menu(dummy)

    root_titles = [
        entry.title for entry in dummy.menu.entries if isinstance(entry, _DummyMenu)
    ]
    root_actions = [
        entry.text
        for entry in dummy.menu.entries
        if isinstance(entry, _DummyMenuAction) and not entry.separator
    ]

    assert root_titles == ["Search", "Layers", "View Direction", "Help"]
    assert root_actions[-3:] == ["Square Client Area", "Fullscreen", "Exit"]


def test_build_window_menu_keeps_file_submenu_for_standard_window(monkeypatch) -> None:
    monkeypatch.setattr(window_module, "QMenu", _DummyMenu)
    monkeypatch.setattr(window_module, "QAction", _DummyMenuAction)

    dummy = SimpleNamespace(
        _frameless_window=False,
        state=SimpleNamespace(rotation_step=5.0),
        enlarge_moon=False,
        show_dso=False,
        dso_catalog_np=None,
        show_asterisms=False,
        show_guidelines=True,
        show_observation_info=True,
        observation_info_mode="auto",
        observation_info_pinned=False,
        sky_disc_alpha=0.2,
        cloud_disc_alpha=0.2,
        satellite_opacity=0.5,
        aircraft_opacity=0.5,
        terrain_horizon_opacity=0.1,
        earth_guide_opacity=0.1,
        urban_outline_opacity=0.2,
        vmag_limit=6.0,
        _rotate_view=lambda **_kwargs: None,
        _open_named_star_jump_dialog=lambda: None,
        _open_named_star_search_dialog=lambda: None,
        _open_place_search_dialog=lambda: None,
        toggle_enlarge_moon=lambda: None,
        toggle_dso=lambda: None,
        toggle_asterisms=lambda: None,
        toggle_guidelines=lambda: None,
        toggle_observation_info=lambda: None,
        toggle_sky_disc=lambda: None,
        toggle_clouds=lambda: None,
        toggle_satellites=lambda: None,
        toggle_aircraft=lambda: None,
        toggle_terrain_horizon=lambda: None,
        toggle_earth_guide=lambda: None,
        toggle_urban_outline=lambda: None,
        toggle_fullscreen=lambda: None,
        square_client_area=lambda: None,
        addAction=lambda _action: None,
        _vmag_limit_menu_text=lambda: "Vmag limit 6.0",
    )
    _install_menu_action_helpers(dummy, [])

    SkyWindow._build_window_menu(dummy)

    root_titles = [
        entry.title for entry in dummy.menu.entries if isinstance(entry, _DummyMenu)
    ]
    root_actions = [
        entry.text
        for entry in dummy.menu.entries
        if isinstance(entry, _DummyMenuAction) and not entry.separator
    ]
    file_actions = [
        entry.text
        for entry in dummy.file_menu.entries
        if isinstance(entry, _DummyMenuAction) and not entry.separator
    ]

    assert root_titles == ["File", "Search", "Layers", "View Direction"]
    assert root_actions == []
    assert file_actions[:3] == ["Square Client Area", "Fullscreen", "Exit"]


def test_toggle_guidelines_disables_and_restores_opacity() -> None:
    dummy = SimpleNamespace()
    dummy.show_guidelines = False
    dummy._action_toggle_guidelines = _DummyAction(False)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    updates: list[str] = []
    dummy.update = lambda: updates.append("update")

    SkyWindow.toggle_guidelines(dummy)

    assert dummy.show_guidelines is True
    assert dummy._action_toggle_guidelines.isChecked() is True

    SkyWindow.toggle_guidelines(dummy)

    assert dummy.show_guidelines is False
    assert dummy._action_toggle_guidelines.isChecked() is False
    assert calls == ["request", "request"]


def test_toggle_observation_info_updates_check_state() -> None:
    dummy = SimpleNamespace()
    dummy.show_observation_info = False
    dummy.observation_info_mode = "auto"
    dummy.observation_info_pinned = False
    dummy._action_toggle_observation_info = _DummyAction(False)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    updates: list[str] = []
    dummy.update = lambda: updates.append("update")

    SkyWindow.toggle_observation_info(dummy)

    assert dummy.show_observation_info is True
    assert dummy._action_toggle_observation_info.isChecked() is True

    SkyWindow.toggle_observation_info(dummy)

    assert dummy.show_observation_info is False
    assert dummy._action_toggle_observation_info.isChecked() is False
    assert calls == ["request", "request"]


def test_vmag_limit_menu_text_formats_current_limit() -> None:
    dummy = SimpleNamespace(vmag_limit=6.5)

    assert SkyWindow._vmag_limit_menu_text(dummy) == "Vmag limit 6.5"


def test_square_client_area_resizes_height_to_match_width() -> None:
    calls: list[tuple[int, int]] = []
    dummy = SimpleNamespace()
    dummy._target_client_size = (640, 480)
    dummy.isFullScreen = lambda: False
    dummy.isMaximized = lambda: False
    dummy.client_width = lambda: 960
    dummy.client_height = lambda: 540
    dummy.width = lambda: 1000
    dummy.height = lambda: 600
    dummy.resize = lambda width, height: calls.append((width, height))

    SkyWindow._resize_client_area(dummy, 960, 960)

    assert dummy._target_client_size == (960, 960)
    assert calls == [(1000, 1020)]


def test_status_line_message_combines_cloud_and_terrain_segments() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_status_line = lambda: "Clouds [AUTO]: downloading"
    dummy._satellite_status_line = lambda: ""
    dummy._aircraft_status_line = lambda: ""
    dummy._jpl_small_body_status_line = lambda: ""
    dummy._terrain_horizon_status_line = lambda: "Terrain horizon: loading DEM..."
    dummy._urban_outline_status_line = lambda: "Urban outline: downloading..."

    got = SkyWindowUpdatesMixin._status_line_message(dummy)

    assert (
        got
        == "Clouds [AUTO]: downloading | Terrain horizon: loading DEM... | Urban outline: downloading..."
    )


def test_status_line_message_keeps_placeholder_icons_for_hidden_layers() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_status_line = lambda: "☁ ---"
    dummy._satellite_status_line = lambda: "🛰 ---"
    dummy._aircraft_status_line = lambda: "✈ ---"
    dummy._jpl_small_body_status_line = lambda: ""
    dummy._terrain_horizon_status_line = lambda: "▲ ---"
    dummy._urban_outline_status_line = lambda: "🂓 ---"

    got = SkyWindowUpdatesMixin._status_line_message(dummy)

    assert got == "☁ --- | 🛰 --- | ✈ --- | ▲ --- | 🂓 ---"


def test_jpl_small_body_status_line_includes_altaz() -> None:
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(
        persistent_search_target=SimpleNamespace(
            label="Voyager 1",
            alt_deg=23.3005,
            az_deg=91.2564,
            jpl_group="mb",
        ),
        persistent_search_next_refresh_utc=None,
        persistent_search_last_error="",
    )
    dummy._target_altaz_suffix = lambda target: SkyWindowUpdatesMixin._target_altaz_suffix(
        dummy,
        target,
    )

    altaz_suffix = SkyWindowUpdatesMixin._target_altaz_suffix(
        dummy,
        dummy.state.persistent_search_target,
    )
    got = SkyWindowUpdatesMixin._jpl_small_body_status_line(dummy)

    assert altaz_suffix == " [alt=23.3 az=91.3]"
    assert got == "JPL [Voyager 1]: held [alt=23.3 az=91.3]"


def test_jpl_small_body_status_line_omits_literal_none_error() -> None:
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(
        persistent_search_target=SimpleNamespace(
            label="Voyager 1",
            alt_deg=48.6,
            az_deg=245.6,
            jpl_group="sb",
        ),
        persistent_search_next_refresh_utc=None,
        persistent_search_last_error=None,
    )
    dummy._target_altaz_suffix = lambda target: SkyWindowUpdatesMixin._target_altaz_suffix(
        dummy,
        target,
    )

    got = SkyWindowUpdatesMixin._jpl_small_body_status_line(dummy)

    assert got == "JPL [Voyager 1]: retry pending [alt=48.6 az=245.6]"


def test_toggle_terrain_horizon_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._terrain_horizon_gui_allowed = False
    dummy.terrain_horizon_opacity = 0.0
    dummy._terrain_horizon_opacity_when_enabled = 0.25
    dummy._action_toggle_terrain_horizon = _DummyAction(False)
    dummy.start_background_terrain_horizon_update = lambda **_kwargs: (
        _ for _ in ()
    ).throw(AssertionError("should not start"))
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
    dummy.request_client_update = lambda: calls.append("request")
    dummy._enable_aircraft_layer = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy._stop_aircraft_timers = lambda: calls.append("stop")
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_aircraft(dummy)

    assert dummy.aircraft_opacity == 1.0
    assert dummy._action_toggle_aircraft.isChecked() is True
    assert calls == ["toggle-on", "request"]


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

    started = SkyWindowUpdatesMixin.start_background_aircraft_update(
        dummy, reason="manual"
    )

    assert started is False
    assert controller_calls == []


def test_on_aircraft_ready_saves_debug_snapshot_when_enabled(
    monkeypatch, tmp_path: Path
) -> None:
    refreshed_at = datetime(2026, 3, 24, 12, 34, 56, tzinfo=timezone.utc)
    dummy_image = _DummyImage()
    dummy = SimpleNamespace()
    dummy.aircraft_state = _DummyAircraftState()
    dummy.aircraft_opacity = 1.0
    dummy.state = SimpleNamespace(aircraft_overlay_points=None)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._schedule_next_aircraft_refresh = lambda: calls.append("schedule")
    dummy.update = lambda: calls.append("update")
    dummy.render_current_image = lambda **kwargs: (
        calls.append(f"render:{kwargs.get('include_hud')}") or dummy_image
    )
    dummy._maybe_save_aircraft_debug_snapshot = lambda payload: (
        SkyWindowUpdatesMixin._maybe_save_aircraft_debug_snapshot(
            dummy,
            payload,
        )
    )
    dummy._resolve_aircraft_debug_snapshot_dir = lambda: (
        SkyWindowUpdatesMixin._resolve_aircraft_debug_snapshot_dir(dummy)
    )
    monkeypatch.setenv("ZSTARVIEW_DEBUG_SAVE_AIRCRAFT_READY_FRAME", str(tmp_path))

    SkyWindowUpdatesMixin._on_aircraft_ready(
        dummy,
        {
            "snapshots": ["s1"],
            "overlay_points": ["p1"],
            "bbox": "bbox",
            "refreshed_at_utc": refreshed_at,
            "source": "OpenSky cache",
        },
    )

    assert dummy.state.aircraft_overlay_points == ["p1"]
    assert calls == ["schedule", "request", "render:True"]
    assert len(dummy_image.saved_paths) == 1
    assert (
        dummy_image.saved_paths[0].name
        == "aircraft-ready-20260324T123456Z-opensky-cache.png"
    )
    assert dummy_image.saved_paths[0].parent == tmp_path
    assert dummy_image.saved_formats == ["PNG"]


def test_on_aircraft_ready_skips_debug_snapshot_when_disabled(monkeypatch) -> None:
    refreshed_at = datetime(2026, 3, 24, 12, 34, 56, tzinfo=timezone.utc)
    dummy = SimpleNamespace()
    dummy.aircraft_state = _DummyAircraftState()
    dummy.aircraft_opacity = 1.0
    dummy.state = SimpleNamespace(aircraft_overlay_points=None)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._schedule_next_aircraft_refresh = lambda: calls.append("schedule")
    dummy.update = lambda: calls.append("update")
    dummy.render_current_image = lambda **kwargs: (_ for _ in ()).throw(
        AssertionError("should not render")
    )
    dummy._maybe_save_aircraft_debug_snapshot = lambda payload: (
        SkyWindowUpdatesMixin._maybe_save_aircraft_debug_snapshot(
            dummy,
            payload,
        )
    )
    dummy._resolve_aircraft_debug_snapshot_dir = lambda: (
        SkyWindowUpdatesMixin._resolve_aircraft_debug_snapshot_dir(dummy)
    )
    monkeypatch.delenv("ZSTARVIEW_DEBUG_SAVE_AIRCRAFT_READY_FRAME", raising=False)

    SkyWindowUpdatesMixin._on_aircraft_ready(
        dummy,
        {
            "snapshots": ["s1"],
            "overlay_points": ["p1"],
            "bbox": "bbox",
            "refreshed_at_utc": refreshed_at,
            "source": "opensky",
        },
    )

    assert dummy.state.aircraft_overlay_points == ["p1"]
    assert calls == ["schedule", "request"]


def test_toggle_terrain_horizon_enables_opacity_and_requests_background_update() -> (
    None
):
    dummy = SimpleNamespace()
    dummy._terrain_horizon_gui_allowed = True
    dummy.terrain_horizon_opacity = 0.0
    dummy._terrain_horizon_opacity_when_enabled = 0.25
    dummy._action_toggle_terrain_horizon = _DummyAction(False)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_terrain_horizon(dummy)

    assert dummy.terrain_horizon_opacity == 0.25
    assert dummy._action_toggle_terrain_horizon.isChecked() is True
    assert calls == ["invalidate", "toggle-on", "request"]


def test_toggle_earth_guide_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._earth_guide_gui_allowed = False
    dummy.earth_guide_opacity = 0.0
    dummy._earth_guide_opacity_when_enabled = 0.25
    dummy._action_toggle_earth_guide = _DummyAction(False)
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: (_ for _ in ()).throw(AssertionError("should not invalidate"))
    )
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint"))

    SkyWindow.toggle_earth_guide(dummy)

    assert dummy.earth_guide_opacity == 0.0
    assert dummy._action_toggle_earth_guide.isChecked() is False


def test_toggle_earth_guide_enables_opacity_and_invalidates_compositor() -> None:
    dummy = SimpleNamespace()
    dummy._earth_guide_gui_allowed = True
    dummy.earth_guide_opacity = 0.0
    dummy._earth_guide_opacity_when_enabled = 0.25
    dummy._action_toggle_earth_guide = _DummyAction(False)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_earth_guide(dummy)

    assert dummy.earth_guide_opacity == 0.25
    assert dummy._action_toggle_earth_guide.isChecked() is True
    assert calls == ["invalidate", "request"]


def test_toggle_urban_outline_enables_opacity_and_requests_background_update() -> None:
    dummy = SimpleNamespace()
    dummy._urban_outline_gui_allowed = True
    dummy.urban_outline_opacity = 0.0
    dummy._urban_outline_opacity_when_enabled = 0.2
    dummy.show_urban_outline_layer = False
    dummy._action_toggle_urban_outline = _DummyAction(False)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy.start_background_urban_outline_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_urban_outline(dummy)

    assert dummy.urban_outline_opacity == 0.2
    assert dummy.show_urban_outline_layer is True
    assert dummy._action_toggle_urban_outline.isChecked() is True
    assert calls == ["toggle-on", "request"]


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


def test_terrain_controller_treats_missing_ocean_tiles_as_empty_profile(
    tmp_path, monkeypatch
) -> None:
    controller = TerrainHorizonController(cache_dir=tmp_path)
    ready_payloads: list[object] = []
    failed_payloads: list[object] = []
    controller.terrain_ready.connect(ready_payloads.append)
    controller.terrain_failed.connect(failed_payloads.append)

    def _raise_no_tiles(**_kwargs):
        raise RuntimeError(
            "No Copernicus DEM tiles were downloaded for the requested area."
        )

    monkeypatch.setattr(
        terrain_controller_module, "fetch_copernicus_dem", _raise_no_tiles
    )

    controller._run_update(lat=20.0, lon=-30.0, observer_height_m=1.7, reason="initial")

    assert failed_payloads == []
    assert ready_payloads == [
        {
            "profile_altaz": [],
            "profile_distances_m": [],
            "secondary_profile_altaz_layers": [],
            "secondary_profile_distances_m_layers": [],
            "source": f"{COPERNICUS_DEM_BUCKET}:ocean",
        }
    ]


def test_toggle_sky_disc_enables_gradient_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy._sky_disc_gui_allowed = True
    dummy.sky_disc_alpha = 0.0
    dummy._sky_disc_alpha_when_enabled = 0.17
    dummy._action_toggle_sky_disc = _DummyAction(False)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request-client")
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.request_sky_data_update = lambda: calls.append("request")
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_sky_disc(dummy)

    assert dummy.sky_disc_alpha == 0.17
    assert dummy._action_toggle_sky_disc.isChecked() is True
    assert calls == ["invalidate", "request", "request-client"]


def test_toggle_sky_disc_switches_to_flat_disc_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy._sky_disc_gui_allowed = True
    dummy.sky_disc_alpha = 0.17
    dummy._sky_disc_alpha_when_enabled = 0.17
    dummy._action_toggle_sky_disc = _DummyAction(True)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request-client")
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.request_sky_data_update = lambda: calls.append("request")
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_sky_disc(dummy)

    assert dummy.sky_disc_alpha == 0.0
    assert dummy._action_toggle_sky_disc.isChecked() is False
    assert calls == ["invalidate", "request", "request-client"]


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
    assert dummy._action_lower_view.isEnabled() is True


def test_jump_to_search_target_keeps_negative_target_alt_for_highlight(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        window_module, "radec_to_altaz", lambda *_args, **_kwargs: (-12.5, 210.0)
    )

    dummy = SimpleNamespace()
    dummy.viewer_data = SimpleNamespace(
        location=(35.0, 139.0), view_center=(20.0, 30.0), observer_height_m=1.7
    )
    dummy.state = SimpleNamespace(
        jump_highlight_name=None,
        jump_highlight_altaz=None,
        jump_highlight_until_ms=0.0,
    )
    sync_calls: list[str] = []
    dummy.request_client_update = lambda: sync_calls.append("request-client")
    dummy._sync_view_altitude_actions = lambda: sync_calls.append("sync")
    dummy._current_time_obj = lambda: object()
    dummy._begin_interaction_mode = lambda: sync_calls.append("begin")
    dummy.request_sky_data_update = lambda: sync_calls.append("request")
    dummy.update = lambda: sync_calls.append("update")

    target = SimpleNamespace(label="Circlet", ra_hours=1.0, dec_deg=2.0, kind="star")
    SkyWindow._jump_to_search_target(dummy, target)

    assert dummy.viewer_data.view_center == (-5.0, 210.0)
    assert dummy.state.jump_highlight_name == "Circlet"
    assert dummy.state.jump_highlight_altaz == (-12.5, 210.0)
    assert dummy.state.jump_highlight_until_ms > 0.0
    assert sync_calls == ["sync", "begin", "request", "request-client"]


def test_rotate_view_in_orientation_mode_updates_render_center_without_full_refresh() -> (
    None
):
    dummy = SimpleNamespace()
    dummy.viewer_data = SimpleNamespace(view_center=(20.0, 30.0))
    dummy.state = SimpleNamespace(render_view_center=(20.0, 30.0))
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request-client")
    dummy._begin_viewport_interaction_mode = lambda *args, **kwargs: calls.append(
        "begin-viewport"
    )
    dummy._begin_interaction_mode = lambda: calls.append("begin")
    dummy._sync_view_altitude_actions = lambda: calls.append("sync")
    dummy._update_viewport_interaction_stars = lambda: calls.append("bright-stars")
    dummy.request_sky_data_update = lambda: calls.append("request")
    dummy.update = lambda: calls.append("update")

    SkyWindow._rotate_view(dummy, d_alt=5.0, d_az=15.0, interactive_viewport=True)

    assert dummy.viewer_data.view_center == (25.0, 45.0)
    assert dummy.state.render_view_center == (25.0, 45.0)
    assert calls == ["begin-viewport", "sync", "bright-stars", "request-client"]


def test_end_viewport_interaction_mode_requests_full_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(
        viewport_interaction_mode=True, viewport_interaction_stars=object()
    )
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request-client")
    dummy.request_sky_data_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy.update = lambda: calls.append("update")

    SkyWindow._end_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    assert calls == [
        "viewport-interaction-idle",
        "view-change-idle",
        "view-change-idle",
        "request-client",
    ]


def test_begin_viewport_interaction_mode_clears_cloud_buffers_and_invalidates_old_render() -> (
    None
):
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(viewport_interaction_mode=False)
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="render-key",
        request_id=42,
        missing_mask_key=99,
    )
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: calls.append("invalidate-compositor")
    )
    dummy._cloud_controller = SimpleNamespace(
        invalidate_pending_render_results=lambda: calls.append("invalidate-cloud")
    )
    dummy._viewport_interaction_idle_timer = SimpleNamespace(
        start=lambda: calls.append("start-timer")
    )

    SkyWindow._begin_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is True
    assert dummy.cloud_state.image is None
    assert dummy.cloud_state.missing_mask is None
    assert dummy.cloud_state.cloud_amount_field is None
    assert dummy.cloud_state.render_key is None
    assert dummy.cloud_state.request_id is None
    assert dummy.cloud_state.missing_mask_key is None
    assert calls == ["invalidate-compositor", "invalidate-cloud", "start-timer"]


def test_begin_viewport_interaction_mode_clears_cloud_buffers_even_while_cloud_update_is_running() -> (
    None
):
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(viewport_interaction_mode=False)
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="render-key",
        request_id=42,
        missing_mask_key=99,
    )
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: calls.append("invalidate-compositor")
    )
    dummy._cloud_controller = SimpleNamespace(
        has_in_flight_update=lambda: True,
        invalidate_pending_render_results=lambda: calls.append("invalidate-cloud"),
    )
    dummy._viewport_interaction_idle_timer = SimpleNamespace(
        start=lambda: calls.append("start-timer")
    )

    SkyWindow._begin_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is True
    assert dummy.cloud_state.image is None
    assert dummy.cloud_state.missing_mask is None
    assert dummy.cloud_state.cloud_amount_field is None
    assert dummy.cloud_state.render_key is None
    assert dummy.cloud_state.request_id is None
    assert dummy.cloud_state.missing_mask_key is None
    assert calls == ["invalidate-compositor", "invalidate-cloud", "start-timer"]


def test_handle_client_resize_preserves_visible_cloud_buffers() -> None:
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy._disc_generation = 0
    dummy._frameless_frame = None
    dummy.menu_button = None
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="render-key",
        request_id=42,
        missing_mask_key=99,
    )
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.request_sky_data_update = lambda: calls.append("sky")
    dummy.request_client_update = lambda: calls.append("client")
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy._begin_viewport_interaction_mode = lambda *args, **kwargs: (
        calls.append(f"begin:{kwargs.get('preserve_cloud_buffers', False)}")
    )

    SkyWindow._handle_client_resize(dummy, SimpleNamespace())

    assert dummy._disc_generation == 1
    assert dummy.cloud_state.image is not None
    assert dummy.cloud_state.missing_mask is not None
    assert dummy.cloud_state.cloud_amount_field is not None
    assert dummy.cloud_state.render_key == "render-key"
    assert dummy.cloud_state.request_id == 42
    assert dummy.cloud_state.missing_mask_key == 99
    assert calls == ["begin:True", "invalidate", "sky", "client", "resize"]


def test_cloud_failed_preserves_last_visible_cloud_frame() -> None:
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(interaction_mode=False)
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        banner_text=None,
        set_error_banner=lambda text: calls.append(f"banner:{text}"),
    )
    dummy._safe_request_cloud_repaint = lambda: calls.append("repaint")

    SkyWindowUpdatesMixin._on_cloud_failed(
        dummy, {"banner": "Clouds: temporary failure"}
    )

    assert dummy.cloud_state.image is not None
    assert dummy.cloud_state.missing_mask is not None
    assert dummy.cloud_state.cloud_amount_field is not None
    assert calls == ["banner:Clouds: temporary failure", "repaint"]


def test_cloud_failed_repaints_status_line_during_interaction() -> None:
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(interaction_mode=True)
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        banner_text=None,
        set_error_banner=lambda text: calls.append(f"banner:{text}"),
    )
    dummy._safe_request_cloud_repaint = lambda: calls.append("repaint")

    SkyWindowUpdatesMixin._on_cloud_failed(dummy, {"banner": "Clouds: timed out"})

    assert calls == ["banner:Clouds: timed out", "repaint"]


def test_discard_stale_disc_images_clears_cached_sky_and_cloud_buffers() -> None:
    compositor_calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(sky_disc_image=object())
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="render-key",
        request_id=42,
        missing_mask_key=99,
    )
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: compositor_calls.append("invalidate")
    )

    SkyWindow._discard_stale_disc_images(dummy)

    assert dummy.state.sky_disc_image is None
    assert dummy.cloud_state.image is None
    assert dummy.cloud_state.missing_mask is None
    assert dummy.cloud_state.cloud_amount_field is None
    assert dummy.cloud_state.render_key is None
    assert dummy.cloud_state.request_id is None
    assert dummy.cloud_state.missing_mask_key is None
    assert compositor_calls == ["invalidate"]


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


def test_menu_button_style_sheet_uses_translucent_background_for_night_preset() -> None:
    dummy = SimpleNamespace(visual_preset="night")

    style = SkyWindow._menu_button_style_sheet(dummy)

    assert "background: transparent;" in style
    assert "background-color: rgba(255, 255, 255, 26);" in style
    assert "border-radius:" not in style


def test_menu_button_style_sheet_uses_light_background_for_day_preset() -> None:
    dummy = SimpleNamespace(visual_preset="day")

    style = SkyWindow._menu_button_style_sheet(dummy)

    assert "background: transparent;" in style
    assert "background-color: rgba(255, 255, 255, 26);" in style


def test_size_grip_style_sheet_is_transparent() -> None:
    dummy = SimpleNamespace(visual_preset="night")

    style = SkyWindow._size_grip_style_sheet(dummy)

    assert "background: transparent;" in style


def test_update_viewport_interaction_stars_uses_bright_limit(monkeypatch) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        window_module,
        "calculate_visible_stars",
        lambda catalog, lat, lon, observer_height_m, time_obj, view_center, max_vmag, subset_indices=None: (
            captured.update(
                {
                    "catalog": catalog,
                    "lat": lat,
                    "lon": lon,
                    "observer_height_m": observer_height_m,
                    "time_obj": time_obj,
                    "view_center": view_center,
                    "max_vmag": max_vmag,
                    "subset_indices": subset_indices,
                }
            )
            or {"name": []},
            object(),
        ),
    )

    dummy = SimpleNamespace()
    dummy.star_catalog_np = object()
    dummy.star_catalog_lod6_indices = np.array([0, 2, 4], dtype=np.int32)
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
        "catalog": dummy.star_catalog_np,
        "lat": 35.0,
        "lon": 139.0,
        "observer_height_m": 12.0,
        "time_obj": "time",
        "view_center": (55.0, 210.0),
        "max_vmag": 4.0,
        "subset_indices": dummy.star_catalog_lod6_indices,
    }


def test_compute_star_render_upscale_factor_matches_downsampled_star_surface() -> None:
    factor = SkyWindow.compute_star_render_upscale_factor(
        disc_width_px=2400,
        expected_width_px=600,
    )

    assert factor == 2.0
