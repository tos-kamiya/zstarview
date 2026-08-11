from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from PySide6.QtCore import QRect, Qt
from PySide6.QtGui import QImage, QPainter

import zstarview.gui.terrain_controller as terrain_controller_module
import zstarview.gui.window as window_module
import zstarview.gui.window_widgets as window_widgets_module
from zstarview.__about__ import __version__
from zstarview.cli.args import SKY_OPACITY_DEFAULT
from zstarview.gui.terrain_controller import TerrainHorizonController
from zstarview.gui.window import SkyWindow, SkyWindowCoreMixin
from zstarview.gui.window_inputs import (
    SkyWindowUserOptions,
    prepare_window_user_options,
)
from zstarview.gui.window_updates import SkyWindowUpdatesMixin
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.render import geometry as render_geometry
from zstarview.render.qt_image import qimage_to_np_rgba
from zstarview.search.models import SearchJumpTarget
from zstarview.simplified_view import resolve_simplified_view_mode
from zstarview.terrain import DEFAULT_TERRAIN_DISTANCE_SAMPLE_STEP_M
from zstarview.types import ViewerData


class _DummyAction:
    def __init__(self, checked: bool) -> None:
        self._checked = checked
        self._enabled = True

    def isChecked(self) -> bool:
        return self._checked

    def setChecked(self, checked: bool) -> None:
        self._checked = checked

    def isEnabled(self) -> bool:
        return self._enabled

    def setEnabled(self, enabled: bool) -> None:
        self._enabled = enabled


class _DummySignal:
    def connect(self, _callback) -> None:
        return None


class _FixedDateTime(datetime):
    @classmethod
    def now(cls, tz=None):
        value = cls(2026, 3, 24, 12, 34, 56, tzinfo=timezone.utc)
        return value.astimezone(tz) if tz is not None else value


class _DummyMenuAction(_DummyAction):
    def __init__(self, text: str, _parent=None, *, separator: bool = False) -> None:
        super().__init__(False)
        self.text = text
        self.separator = separator
        self.triggered = _DummySignal()
        self.shortcut = None
        self.shortcut_context = None
        self.checkable = False

    def setShortcut(self, shortcut) -> None:
        self.shortcut = shortcut

    def setShortcutContext(self, context) -> None:
        self.shortcut_context = context

    def setCheckable(self, checkable: bool) -> None:
        self.checkable = checkable


class _DummyMenu:
    def __init__(self, title: str, _parent=None) -> None:
        self.title = title
        self.entries: list[object] = []

    def addMenu(self, menu):
        self.entries.append(menu)
        return menu

    def addAction(self, action):
        if isinstance(action, str):
            action = _DummyMenuAction(action)
        self.entries.append(action)
        return action

    def addSeparator(self):
        action = _DummyMenuAction("", separator=True)
        self.entries.append(action)
        return action


class _DummyAircraftState:
    def __init__(self) -> None:
        self.banner_text = None
        self.last_result = None

    def set_result(self, snapshots, *, bbox, refreshed_at_utc) -> None:
        self.last_result = {
            "snapshots": snapshots,
            "bbox": bbox,
            "refreshed_at_utc": refreshed_at_utc,
        }

    def set_banner(self, text: str) -> None:
        self.banner_text = text


def _noop_request_client_update() -> None:
    return None


class _DummyPaintEvent:
    def __init__(self) -> None:
        self.accepted = False

    def accept(self) -> None:
        self.accepted = True


class _DummyKeyEvent:
    def __init__(
        self,
        key: int,
        *,
        auto_repeat: bool = False,
        modifiers=Qt.KeyboardModifier.NoModifier,
    ) -> None:
        self._key = key
        self._auto_repeat = auto_repeat
        self._modifiers = modifiers
        self.accepted = False

    def key(self) -> int:
        return self._key

    def modifiers(self):
        return self._modifiers

    def isAutoRepeat(self) -> bool:
        return self._auto_repeat

    def accept(self) -> None:
        self.accepted = True


def _install_menu_action_helpers(
    dummy: SimpleNamespace, added_actions: list[object]
) -> None:
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


def test_prepare_window_user_options_normalizes_terrain_horizon_fields() -> None:
    options = prepare_window_user_options(
        sky_disc_alpha=SKY_OPACITY_DEFAULT,
        sky_disc_style="smooth",
        sky_disc_altaz_rings="dimalt",
        sky_disc_altaz_rings_hover="altaz",
        cloud_disc_alpha=0.075,
        satellite_opacity=0.7,
        tropical_cyclone_opacity=0.4,
        terrain_horizon_opacity=1.5,
        earth_guide_opacity=1.5,
        urban_outline_opacity=1.5,
        aircraft_opacity=0.4,
        ground_tint_opacity=0.2,
        water_overlay_opacity=0.4,
        overlay_font_size=20.5,
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
        tropical_cyclone_gui_allowed=False,
        terrain_horizon_gui_allowed=False,
        earth_guide_gui_allowed=False,
        urban_outline_gui_allowed=False,
    )

    assert options.terrain_horizon_opacity == 1.0
    assert options.earth_guide_opacity == 1.0
    assert options.earth_guide_visibility_boost == 2.0
    assert options.urban_outline_opacity == 1.0
    assert options.tropical_cyclone_opacity == pytest.approx(0.8)
    assert options.sky_disc_alpha == pytest.approx(0.32)
    assert options.sky_disc_style == "smooth"
    assert options.sky_disc_altaz_rings == "dimalt"
    assert options.sky_disc_altaz_rings_hover == "altaz"
    assert options.cloud_disc_alpha == pytest.approx(0.15)
    assert options.satellite_opacity == 1.0
    assert options.aircraft_opacity == 0.8
    assert options.ground_tint_opacity == 0.4
    assert options.asterism_visibility_boost == pytest.approx(2.0)
    assert options.sky_disc_gui_allowed is False
    assert options.cloud_gui_allowed is False


def test_sky_window_user_options_defaults_ridge_glow_opacity() -> None:
    options = SkyWindowUserOptions()

    assert options.ridge_glow_opacity == 0.04


def test_toggle_clouds_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_toggle_supported = True
    dummy._cloud_gui_allowed = False
    dummy._clouddisc = object()
    dummy.cloud_disc_alpha = 0.0
    dummy._cloud_alpha_when_enabled = 0.2
    dummy._action_toggle_clouds = _DummyAction(False)
    dummy._sync_cloud_action_state = lambda: SkyWindow._sync_cloud_action_state(dummy)
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint"))

    SkyWindow.toggle_clouds(dummy)

    assert dummy.cloud_disc_alpha == 0.0
    assert dummy._action_toggle_clouds.isChecked() is False


def test_apply_startup_delta_t_syncs_cloud_action_with_enabled_state() -> None:
    dummy = SimpleNamespace()
    dummy.delta_t = timedelta(0)
    dummy._clouddisc = object()
    dummy._geo_satellite_enabled = False
    dummy._cloud_toggle_supported = True
    dummy._cloud_gui_allowed = True
    dummy._cloud_requested_enabled = True
    dummy._cloud_alpha_when_enabled = 0.2
    dummy.cloud_disc_alpha = 0.0
    dummy._satellite_toggle_supported = True
    dummy._satellite_requested_enabled = False
    dummy._satellite_opacity_when_enabled = 0.5
    dummy._satellite_gui_allowed = True
    dummy.satellite_opacity = 0.0
    dummy._aircraft_toggle_supported = True
    dummy._aircraft_requested_enabled = False
    dummy._aircraft_opacity_when_enabled = 0.5
    dummy._aircraft_gui_allowed = True
    dummy.aircraft_opacity = 0.0
    dummy._action_toggle_clouds = _DummyAction(False)
    dummy._action_toggle_satellites = _DummyAction(False)
    dummy._action_toggle_aircraft = _DummyAction(False)
    dummy._action_toggle_tropical_cyclone = None
    dummy._tropical_cyclone_controller = None
    dummy._sync_cloud_action_state = lambda: SkyWindow._sync_cloud_action_state(dummy)

    SkyWindow.apply_startup_delta_t(dummy, timedelta(0))

    assert dummy.cloud_disc_alpha == pytest.approx(0.2)
    assert dummy._action_toggle_clouds.isChecked() is True
    assert dummy._action_toggle_clouds.isEnabled() is True


def test_toggle_sky_disc_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._sky_disc_gui_allowed = False
    dummy.sky_disc_alpha = 0.0
    dummy._sky_disc_alpha_when_enabled = SKY_OPACITY_DEFAULT
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
        night_light_opacity=0.022,
        road_night_lights_opacity=0.0,
        ridge_glow_opacity=0.25,
        cloud_disc_alpha=0.2,
        _geo_satellite_enabled=False,
        _geo_satellite_toggle_supported=lambda: True,
        _ridge_glow_toggle_supported=True,
        water_overlay_opacity=0.4,
        satellite_opacity=0.5,
        aircraft_opacity=0.5,
        terrain_horizon_opacity=0.1,
        earth_guide_opacity=0.1,
        urban_outline_opacity=0.2,
        show_tropical_cyclone_overlay=False,
        _tropical_cyclone_toggle_supported=False,
        _night_light_toggle_supported=True,
        _road_light_toggle_supported=True,
        _water_overlay_gui_allowed=True,
        _terrain_horizon_gui_allowed=True,
        _water_overlay_action_enabled=lambda: True,
        vmag_limit=6.0,
        _rotate_view=lambda **_kwargs: None,
        _open_named_star_jump_dialog=lambda: None,
        _open_named_star_search_dialog=lambda: None,
        _open_place_search_dialog=lambda: None,
        _open_view_direction_dialog=lambda: None,
        toggle_enlarge_moon=lambda: None,
        toggle_dso=lambda: None,
        toggle_asterisms=lambda: None,
        toggle_guidelines=lambda: None,
        toggle_observation_info=lambda: None,
        toggle_sky_disc=lambda: None,
        toggle_clouds=lambda: None,
        toggle_geo_satellite=lambda: None,
        toggle_satellites=lambda: None,
        toggle_aircraft=lambda: None,
        toggle_terrain_horizon=lambda: None,
        toggle_water_overlay=lambda: None,
        toggle_earth_guide=lambda: None,
        toggle_urban_outline=lambda: None,
        toggle_tropical_cyclone_overlay=lambda: None,
        toggle_night_lights=lambda: None,
        toggle_road_lights=lambda: None,
        toggle_fullscreen=lambda: None,
        square_client_area=lambda: None,
        toggle_square_window=lambda: None,
        _restore_default_window_size=lambda: None,
        _fit_client_area_to_screen=lambda: None,
        _request_application_quit=lambda: None,
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
    help_actions = [
        entry.text
        for entry in dummy.help_menu.entries
        if isinstance(entry, _DummyMenuAction) and not entry.separator
    ]
    help_enabled = [
        entry.isEnabled()
        for entry in dummy.help_menu.entries
        if isinstance(entry, _DummyMenuAction) and not entry.separator
    ]

    assert root_titles == ["Search", "Layers", "View Direction", "Help"]
    assert root_actions == [
        "Square Window",
        "Default Window Size",
        "Fit to Screen",
        "Fullscreen",
        "Exit",
    ]
    assert help_actions == [
        "Open-Meteo Terms...",
        "Licenses and Data Sources...",
        f"Version {__version__}",
    ]
    assert help_enabled == [True, True, False]


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
        night_light_opacity=0.022,
        road_night_lights_opacity=0.0,
        ridge_glow_opacity=0.25,
        cloud_disc_alpha=0.2,
        _geo_satellite_enabled=False,
        _geo_satellite_toggle_supported=lambda: True,
        _ridge_glow_toggle_supported=True,
        water_overlay_opacity=0.4,
        satellite_opacity=0.5,
        aircraft_opacity=0.5,
        terrain_horizon_opacity=0.1,
        earth_guide_opacity=0.1,
        urban_outline_opacity=0.2,
        show_tropical_cyclone_overlay=False,
        _tropical_cyclone_toggle_supported=False,
        _night_light_toggle_supported=True,
        _road_light_toggle_supported=True,
        _water_overlay_gui_allowed=True,
        _terrain_horizon_gui_allowed=True,
        _water_overlay_action_enabled=lambda: True,
        vmag_limit=6.0,
        _rotate_view=lambda **_kwargs: None,
        _open_named_star_jump_dialog=lambda: None,
        _open_named_star_search_dialog=lambda: None,
        _open_place_search_dialog=lambda: None,
        _open_view_direction_dialog=lambda: None,
        toggle_enlarge_moon=lambda: None,
        toggle_dso=lambda: None,
        toggle_asterisms=lambda: None,
        toggle_guidelines=lambda: None,
        toggle_observation_info=lambda: None,
        toggle_sky_disc=lambda: None,
        toggle_clouds=lambda: None,
        toggle_geo_satellite=lambda: None,
        toggle_satellites=lambda: None,
        toggle_aircraft=lambda: None,
        toggle_terrain_horizon=lambda: None,
        toggle_water_overlay=lambda: None,
        toggle_earth_guide=lambda: None,
        toggle_urban_outline=lambda: None,
        toggle_tropical_cyclone_overlay=lambda: None,
        toggle_night_lights=lambda: None,
        toggle_road_lights=lambda: None,
        toggle_fullscreen=lambda: None,
        square_client_area=lambda: None,
        toggle_square_window=lambda: None,
        _restore_default_window_size=lambda: None,
        _fit_client_area_to_screen=lambda: None,
        _request_application_quit=lambda: None,
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
    help_actions = [
        entry.text
        for entry in dummy.help_menu.entries
        if isinstance(entry, _DummyMenuAction) and not entry.separator
    ]
    help_enabled = [
        entry.isEnabled()
        for entry in dummy.help_menu.entries
        if isinstance(entry, _DummyMenuAction) and not entry.separator
    ]

    assert root_titles == ["File", "Search", "Layers", "View Direction", "Help"]
    assert root_actions == []
    assert file_actions[:5] == [
        "Square Window",
        "Default Window Size",
        "Fit to Screen",
        "Fullscreen",
        "Exit",
    ]
    assert help_actions == [
        "Open-Meteo Terms...",
        "Licenses and Data Sources...",
        f"Version {__version__}",
    ]
    assert help_enabled == [True, True, False]


def test_build_window_menu_groups_layers_by_sky_and_ground(monkeypatch) -> None:
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
        _geo_satellite_enabled=False,
        _geo_satellite_toggle_supported=lambda: True,
        ridge_glow_opacity=0.25,
        _ridge_glow_toggle_supported=True,
        water_overlay_opacity=0.4,
        satellite_opacity=0.5,
        aircraft_opacity=0.5,
        terrain_horizon_opacity=0.1,
        earth_guide_opacity=0.1,
        night_light_opacity=0.022,
        road_night_lights_opacity=0.0,
        _night_light_toggle_supported=True,
        _road_light_toggle_supported=True,
        urban_outline_opacity=0.2,
        show_tropical_cyclone_overlay=False,
        _tropical_cyclone_toggle_supported=False,
        _water_overlay_action_enabled=lambda: True,
        vmag_limit=6.0,
        toggle_night_lights=lambda: None,
        toggle_road_lights=lambda: None,
        toggle_urban_outline=lambda: None,
        toggle_tropical_cyclone_overlay=lambda: None,
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
        toggle_geo_satellite=lambda: None,
        toggle_satellites=lambda: None,
        toggle_aircraft=lambda: None,
        toggle_terrain_horizon=lambda: None,
        toggle_water_overlay=lambda: None,
        toggle_earth_guide=lambda: None,
        toggle_fullscreen=lambda: None,
        square_client_area=lambda: None,
        toggle_square_window=lambda: None,
        _restore_default_window_size=lambda: None,
        _open_view_direction_dialog=lambda: None,
        _fit_client_area_to_screen=lambda: None,
        _request_application_quit=lambda: None,
        addAction=lambda action: None,
        _vmag_limit_menu_text=lambda: "Vmag limit 6.0",
    )
    _install_menu_action_helpers(dummy, [])

    SkyWindow._build_window_menu(dummy)

    layer_entries = [
        entry
        for entry in dummy.display_menu.entries
        if isinstance(entry, _DummyMenuAction)
    ]
    layer_labels = [entry.text for entry in layer_entries if not entry.separator]
    separator_indexes = [
        index
        for index, entry in enumerate(dummy.display_menu.entries)
        if getattr(entry, "separator", False)
    ]

    assert layer_labels == [
        "Enlarge Moon",
        "DSO",
        "Asterisms",
        "Sky Guides",
        "Observation Info",
        "Sky Color",
        "Clouds",
        "Geo-satellite",
        "Satellites",
        "Aircraft",
        "Typhoon / Cyclone",
        "Forecast Precipitation",
        "Night Lights",
        "Road Lights",
        "Urban Outline",
        "Terrain Horizon",
        "Water Surface",
        "Earth Guide",
        "Vmag limit 6.0",
    ]
    assert len(separator_indexes) == 4


def test_build_window_menu_disables_water_surface_when_terrain_horizon_off(
    monkeypatch,
) -> None:
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
        night_light_opacity=0.022,
        road_night_lights_opacity=0.0,
        ridge_glow_opacity=0.25,
        cloud_disc_alpha=0.2,
        _geo_satellite_enabled=False,
        _geo_satellite_toggle_supported=lambda: True,
        _ridge_glow_toggle_supported=True,
        water_overlay_opacity=0.4,
        satellite_opacity=0.5,
        aircraft_opacity=0.5,
        terrain_horizon_opacity=0.0,
        earth_guide_opacity=0.1,
        urban_outline_opacity=0.2,
        show_tropical_cyclone_overlay=False,
        _tropical_cyclone_toggle_supported=False,
        _night_light_toggle_supported=True,
        _road_light_toggle_supported=True,
        _water_overlay_gui_allowed=True,
        _terrain_horizon_gui_allowed=True,
        _water_overlay_action_enabled=lambda: False,
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
        toggle_geo_satellite=lambda: None,
        toggle_satellites=lambda: None,
        toggle_aircraft=lambda: None,
        toggle_terrain_horizon=lambda: None,
        toggle_water_overlay=lambda: None,
        toggle_earth_guide=lambda: None,
        toggle_urban_outline=lambda: None,
        toggle_tropical_cyclone_overlay=lambda: None,
        toggle_night_lights=lambda: None,
        toggle_road_lights=lambda: None,
        toggle_fullscreen=lambda: None,
        square_client_area=lambda: None,
        toggle_square_window=lambda: None,
        _restore_default_window_size=lambda: None,
        _open_view_direction_dialog=lambda: None,
        _fit_client_area_to_screen=lambda: None,
        _request_application_quit=lambda: None,
        addAction=lambda _action: None,
        _vmag_limit_menu_text=lambda: "Vmag limit 6.0",
    )
    _install_menu_action_helpers(dummy, [])

    SkyWindow._build_window_menu(dummy)

    water_action = next(
        entry
        for entry in dummy.display_menu.entries
        if isinstance(entry, _DummyMenuAction) and entry.text == "Water Surface"
    )

    assert water_action.isEnabled() is False


def test_toggle_terrain_horizon_disables_and_restores_water_surface_action() -> None:
    dummy = SimpleNamespace()
    dummy._terrain_horizon_gui_allowed = True
    dummy._water_overlay_gui_allowed = True
    dummy.terrain_horizon_opacity = 0.25
    dummy._terrain_horizon_opacity_when_enabled = 0.25
    dummy._action_toggle_terrain_horizon = _DummyAction(True)
    dummy._action_toggle_water_overlay = _DummyAction(True)
    dummy.terrain_horizon_state = SimpleNamespace(ground_elevation_m=42.0)
    dummy._refresh_water_overlay_active_dots = lambda: None
    dummy._sync_water_overlay_action_enabled = lambda: None
    dummy._water_overlay_action_enabled = lambda: (
        SkyWindowUpdatesMixin._water_overlay_action_enabled(dummy)
    )
    dummy._sync_water_overlay_action_enabled = lambda: (
        SkyWindowUpdatesMixin._sync_water_overlay_action_enabled(dummy)
    )
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_terrain_horizon(dummy)

    assert dummy.terrain_horizon_opacity == 0.0
    assert dummy._action_toggle_terrain_horizon.isChecked() is False
    assert dummy._action_toggle_water_overlay.isEnabled() is False
    assert calls == ["invalidate", "request"]

    SkyWindow.toggle_terrain_horizon(dummy)

    assert dummy.terrain_horizon_opacity == 0.25
    assert dummy._action_toggle_terrain_horizon.isChecked() is True
    assert dummy._action_toggle_water_overlay.isEnabled() is True
    assert calls == ["invalidate", "request", "invalidate", "toggle-on", "request"]


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


def test_square_window_action_resizes_client_area_once() -> None:
    calls: list[tuple[int, int]] = []
    request_calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.square_client_area = lambda: calls.append((640, 640))
    dummy.request_client_update = lambda: request_calls.append("request")
    dummy._action_square_window = _DummyAction(False)

    SkyWindow.toggle_square_window(dummy)

    assert calls == [(640, 640)]
    assert request_calls == ["request"]
    assert dummy._action_square_window.isEnabled() is True


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


def test_restore_default_window_size_uses_default_client_dimensions() -> None:
    calls: list[tuple[int, int]] = []
    dummy = SimpleNamespace()
    dummy.isFullScreen = lambda: False
    dummy.isMaximized = lambda: False
    dummy._resize_client_area = lambda width, height: calls.append((width, height))

    SkyWindow._restore_default_window_size(dummy)

    assert calls == [(window_module.WINDOW_WIDTH, window_module.WINDOW_HEIGHT)]


def test_status_line_message_combines_cloud_and_terrain_segments() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_status_line = lambda: "Clouds [AUTO]: downloading"
    dummy._satellite_status_line = lambda: ""
    dummy._aircraft_status_line = lambda: ""
    dummy._tropical_cyclone_status_line = lambda: ""
    dummy._jpl_small_body_status_line = lambda: ""
    dummy._terrain_horizon_status_line = lambda: "△ loading DEM..."
    dummy._water_overlay_status_line = lambda: ""
    dummy._urban_outline_status_line = lambda: "🂓 downloading..."
    dummy._effective_simplified_view_mode = lambda: "normal"

    got = SkyWindowUpdatesMixin._status_line_message(dummy)

    assert got == "⎮ Clouds [AUTO]: downloading ⎮ △ loading DEM... ⎮ 🂓 downloading... ⎮"


def test_status_line_message_keeps_placeholder_icons_for_hidden_layers() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_status_line = lambda: "☁ ---"
    dummy._satellite_status_line = lambda: "🛰 ---"
    dummy._aircraft_status_line = lambda: "✈ ---"
    dummy._tropical_cyclone_status_line = lambda: ""
    dummy._jpl_small_body_status_line = lambda: ""
    dummy._terrain_horizon_status_line = lambda: "△ ---"
    dummy._water_overlay_status_line = lambda: ""
    dummy._urban_outline_status_line = lambda: "🂓 ---"
    dummy._effective_simplified_view_mode = lambda: "normal"

    got = SkyWindowUpdatesMixin._status_line_message(dummy)

    assert got == "⎮ ☁ --- ⎮ 🛰 --- ⎮ ✈ --- ⎮ △ --- ⎮ 🂓 --- ⎮"


def test_status_line_message_orders_cyclone_before_satellite_and_aircraft() -> None:
    dummy = SimpleNamespace()
    dummy._cloud_status_line = lambda: "☁ cloud"
    dummy._satellite_status_line = lambda: "🛰 sat"
    dummy._aircraft_status_line = lambda: "✈ aircraft"
    dummy._tropical_cyclone_status_line = lambda: "TC cyclone"
    dummy._jpl_small_body_status_line = lambda: ""
    dummy._terrain_horizon_status_line = lambda: "△ terrain"
    dummy._water_overlay_status_line = lambda: "W water"
    dummy._urban_outline_status_line = lambda: "🂓 urban"
    dummy._effective_simplified_view_mode = lambda: "normal"

    got = SkyWindowUpdatesMixin._status_line_message(dummy)

    assert (
        got
        == "⎮ ☁ cloud ⎮ TC cyclone ⎮ 🛰 sat ⎮ ✈ aircraft ⎮ △ terrain ⎮ W water ⎮ 🂓 urban ⎮"
    )


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
    dummy._target_altaz_suffix = lambda target: (
        SkyWindowUpdatesMixin._target_altaz_suffix(
            dummy,
            target,
        )
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
    dummy._target_altaz_suffix = lambda target: (
        SkyWindowUpdatesMixin._target_altaz_suffix(
            dummy,
            target,
        )
    )

    got = SkyWindowUpdatesMixin._jpl_small_body_status_line(dummy)

    assert got == "JPL [Voyager 1]: retry pending [alt=48.6 az=245.6]"


def test_toggle_terrain_horizon_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._terrain_horizon_gui_allowed = False
    dummy.terrain_horizon_opacity = 0.0
    dummy._terrain_horizon_opacity_when_enabled = 0.25
    dummy._action_toggle_terrain_horizon = _DummyAction(False)
    dummy._action_toggle_water_overlay = _DummyAction(False)
    dummy._sync_water_overlay_action_enabled = lambda: None
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
    dummy._viewport_interaction_active = lambda: False

    started = SkyWindowUpdatesMixin.start_background_aircraft_update(
        dummy, reason="manual"
    )

    assert started is False
    assert controller_calls == []


def test_periodic_debug_snapshot_timer_queues_snapshot_when_enabled(
    monkeypatch, tmp_path: Path
) -> None:
    dummy = SimpleNamespace()
    dummy.aircraft_state = _DummyAircraftState()
    dummy.aircraft_state.snapshots = ["s1"]
    dummy.aircraft_opacity = 1.0
    dummy.state = SimpleNamespace()
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._resolve_periodic_debug_snapshot_dir = (
        SkyWindowCoreMixin._resolve_periodic_debug_snapshot_dir
    )
    dummy._queue_periodic_debug_snapshot = lambda payload: (
        SkyWindowCoreMixin._queue_periodic_debug_snapshot(dummy, payload)
    )
    monkeypatch.setenv("ZSTARVIEW_DEBUG_SAVE_PERIODIC_FRAME", str(tmp_path))

    monkeypatch.setattr(
        "zstarview.gui.window_updates.datetime",
        _FixedDateTime,
    )
    SkyWindowUpdatesMixin._on_periodic_debug_snapshot_timer(dummy)

    assert calls == ["request"]
    assert dummy._pending_periodic_debug_snapshot_path == (
        tmp_path / "periodic-20260324T123456Z.png"
    )
    assert list(tmp_path.iterdir()) == []


def test_periodic_debug_snapshot_timer_queues_before_aircraft_data(
    monkeypatch, tmp_path: Path
) -> None:
    dummy = SimpleNamespace()
    dummy.aircraft_state = _DummyAircraftState()
    dummy.aircraft_state.snapshots = None
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._resolve_periodic_debug_snapshot_dir = (
        SkyWindowCoreMixin._resolve_periodic_debug_snapshot_dir
    )
    dummy._queue_periodic_debug_snapshot = lambda payload: (
        SkyWindowCoreMixin._queue_periodic_debug_snapshot(dummy, payload)
    )
    monkeypatch.setenv("ZSTARVIEW_DEBUG_SAVE_PERIODIC_FRAME", str(tmp_path))

    monkeypatch.setattr(
        "zstarview.gui.window_updates.datetime",
        _FixedDateTime,
    )
    SkyWindowUpdatesMixin._on_periodic_debug_snapshot_timer(dummy)

    assert calls == ["request"]
    assert dummy._pending_periodic_debug_snapshot_path == (
        tmp_path / "periodic-20260324T123456Z.png"
    )


def test_periodic_debug_snapshot_is_saved_during_paint(
    tmp_path: Path,
) -> None:
    dummy = SimpleNamespace()
    dummy._pending_periodic_debug_snapshot_path = (
        tmp_path / "periodic-20260324T123456Z.png"
    )
    dummy._save_periodic_debug_snapshot_image = lambda image, output_path: (
        SkyWindowCoreMixin._save_periodic_debug_snapshot_image(image, output_path)
    )

    frame = QImage(8, 8, QImage.Format.Format_ARGB32_Premultiplied)
    frame.fill(0)

    SkyWindowCoreMixin._flush_periodic_debug_snapshot_save(dummy, frame)

    saved = list(tmp_path.iterdir())
    assert len(saved) == 1
    assert saved[0].name == "periodic-20260324T123456Z.png"
    assert dummy._pending_periodic_debug_snapshot_path is None


def test_periodic_debug_snapshot_saves_each_queued_snapshot(
    monkeypatch, tmp_path: Path
) -> None:
    refreshed_at_1 = datetime(2026, 3, 24, 12, 34, 56, tzinfo=timezone.utc)
    refreshed_at_2 = datetime(2026, 3, 24, 12, 34, 57, tzinfo=timezone.utc)
    dummy = SimpleNamespace()
    monkeypatch.setenv("ZSTARVIEW_DEBUG_SAVE_PERIODIC_FRAME", str(tmp_path))
    dummy._save_periodic_debug_snapshot_image = lambda image, output_path: (
        SkyWindowCoreMixin._save_periodic_debug_snapshot_image(image, output_path)
    )

    frame = QImage(8, 8, QImage.Format.Format_ARGB32_Premultiplied)
    frame.fill(0)

    SkyWindowCoreMixin._queue_periodic_debug_snapshot(
        dummy,
        {"refreshed_at_utc": refreshed_at_1, "source": "OpenSky cache"},
    )
    SkyWindowCoreMixin._flush_periodic_debug_snapshot_save(
        dummy,
        frame,
    )
    SkyWindowCoreMixin._queue_periodic_debug_snapshot(
        dummy,
        {"refreshed_at_utc": refreshed_at_2, "source": "OpenSky cache"},
    )
    SkyWindowCoreMixin._flush_periodic_debug_snapshot_save(
        dummy,
        frame,
    )

    saved = sorted(path.name for path in tmp_path.iterdir())
    assert saved == [
        "periodic-20260324T123456Z.png",
        "periodic-20260324T123457Z.png",
    ]


def test_on_aircraft_ready_skips_debug_snapshot_when_disabled(monkeypatch) -> None:
    refreshed_at = datetime(2026, 3, 24, 12, 34, 56, tzinfo=timezone.utc)
    dummy = SimpleNamespace()
    dummy.aircraft_state = _DummyAircraftState()
    dummy.aircraft_opacity = 1.0
    dummy.state = SimpleNamespace()
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._schedule_next_aircraft_refresh = lambda: calls.append("schedule")
    dummy.reproject_aircraft_overlay = lambda: calls.append("reproject")
    dummy.update = lambda: calls.append("update")
    dummy.render_current_image = lambda **kwargs: (_ for _ in ()).throw(
        AssertionError("should not render")
    )

    SkyWindowUpdatesMixin._on_aircraft_ready(
        dummy,
        {
            "snapshots": ["s1"],
            "bbox": "bbox",
            "refreshed_at_utc": refreshed_at,
            "source": "opensky",
        },
    )

    assert dummy.state is not None
    assert calls == ["schedule", "reproject"]


def test_on_aircraft_ready_skips_debug_snapshot_for_cache_fresh(
    monkeypatch, tmp_path
) -> None:
    refreshed_at = datetime(2026, 3, 24, 12, 34, 56, tzinfo=timezone.utc)
    dummy = SimpleNamespace()
    dummy.aircraft_state = _DummyAircraftState()
    dummy.aircraft_opacity = 1.0
    dummy.state = SimpleNamespace()
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._schedule_next_aircraft_refresh = lambda: calls.append("schedule")
    dummy.reproject_aircraft_overlay = lambda: calls.append("reproject")
    dummy.update = lambda: calls.append("update")
    dummy.render_current_image = lambda **kwargs: (_ for _ in ()).throw(
        AssertionError("should not render")
    )
    SkyWindowUpdatesMixin._on_aircraft_ready(
        dummy,
        {
            "snapshots": ["s1"],
            "bbox": "bbox",
            "refreshed_at_utc": refreshed_at,
            "source": "cache-fresh",
        },
    )

    assert dummy.state is not None
    assert calls == ["schedule", "reproject"]
    assert not hasattr(dummy, "_pending_periodic_debug_snapshot_path")
    assert list(tmp_path.iterdir()) == []


def test_toggle_terrain_horizon_enables_opacity_and_requests_background_update() -> (
    None
):
    dummy = SimpleNamespace()
    dummy._terrain_horizon_gui_allowed = True
    dummy.terrain_horizon_opacity = 0.0
    dummy._terrain_horizon_opacity_when_enabled = 0.25
    dummy._action_toggle_terrain_horizon = _DummyAction(False)
    dummy.terrain_horizon_state = SimpleNamespace(ground_elevation_m=17.5)
    dummy._refresh_water_overlay_active_dots = lambda: None
    dummy._sync_water_overlay_action_enabled = lambda: None
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


def test_toggle_terrain_horizon_off_keeps_retained_ground_elevation() -> None:
    dummy = SimpleNamespace()
    dummy._terrain_horizon_gui_allowed = True
    dummy.terrain_horizon_opacity = 0.25
    dummy._terrain_horizon_opacity_when_enabled = 0.25
    dummy._action_toggle_terrain_horizon = _DummyAction(True)
    dummy.terrain_horizon_state = SimpleNamespace(ground_elevation_m=42.0)
    dummy._refresh_water_overlay_active_dots = lambda: None
    dummy._sync_water_overlay_action_enabled = lambda: None
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_terrain_horizon(dummy)

    assert dummy.terrain_horizon_opacity == 0.0
    assert dummy.terrain_horizon_state.ground_elevation_m == 42.0
    assert dummy._action_toggle_terrain_horizon.isChecked() is False
    assert calls == ["invalidate", "request"]


def test_water_overlay_ground_elevation_prefers_retained_terrain_ground() -> None:
    dummy = SimpleNamespace(
        terrain_horizon_state=SimpleNamespace(ground_elevation_m=123.4),
        viewer_data=SimpleNamespace(ground_elevation_m=17.5),
    )

    got = SkyWindowUpdatesMixin._water_overlay_ground_elevation_m(dummy)

    assert got == 123.4


def test_terrain_horizon_failed_keeps_retained_ground_elevation() -> None:
    dummy = SimpleNamespace()
    dummy.terrain_horizon_state = SimpleNamespace(
        ground_elevation_m=58.0,
        clear_profile=lambda: None,
        set_error_banner=lambda _text: None,
    )
    dummy.state = SimpleNamespace(
        terrain_horizon_profile=[(1.0, 2.0)],
        terrain_horizon_profile_distances_m=[1.0],
        terrain_secondary_ridges_altaz_layers=[[(1.0, 2.0)]],
        terrain_secondary_ridges_distances_m_layers=[[1.0]],
    )
    dummy._refresh_water_overlay_active_dots = lambda: None
    dummy._is_shutting_down = True
    dummy._startup_initial_load_started = False
    dummy._startup_initial_data_loaded = False
    dummy._sync_water_overlay_action_enabled = lambda: None
    dummy.start_background_water_overlay_update = lambda **_kwargs: (
        _ for _ in ()
    ).throw(AssertionError("should not restart water overlay"))
    dummy._compositor = SimpleNamespace(invalidate=lambda: None)
    dummy.request_client_update = lambda: None
    dummy._continue_initial_data_load = lambda: None

    SkyWindowUpdatesMixin._on_terrain_horizon_failed(
        dummy, {"banner": "Terrain horizon: unavailable"}
    )

    assert dummy.terrain_horizon_state.ground_elevation_m == 58.0


def test_water_overlay_ready_invalidates_and_requests_refresh() -> None:
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.water_overlay_state = SimpleNamespace(
        dots=None,
        sea_level_dots=None,
        dem_dots=None,
        banner_text="",
        set_dem_dots_result=lambda dots, source=None: calls.append(
            f"dem:{len(dots)}:{source}"
        ),
        set_sea_level_dots_result=lambda dots, source=None: calls.append(
            f"sea:{len(dots)}:{source}"
        ),
        set_error_banner=lambda text: calls.append(f"error:{text}"),
    )
    dummy._refresh_water_overlay_active_dots = lambda: calls.append("refresh")
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.request_client_update = lambda: calls.append("request")
    dummy._startup_initial_load_started = False
    dummy._startup_initial_data_loaded = False
    dummy._continue_initial_data_load = lambda: None

    SkyWindowUpdatesMixin._on_water_overlay_ready(
        dummy,
        {
            "dots": [object()],
            "sea_dots": [object()],
            "dem_dots": None,
            "mode": "sea",
            "source": "ready",
        },
    )

    assert calls == ["sea:1:ready", "refresh", "invalidate", "request"]


def test_start_initial_data_load_only_kicks_off_sky() -> None:
    dummy = SimpleNamespace()
    dummy._startup_initial_load_started = False
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.0
    dummy.terrain_horizon_opacity = 0.2
    dummy.water_overlay_opacity = 0.2
    dummy.urban_outline_opacity = 0.0
    dummy._satellite_layer_enabled = lambda: False
    dummy._aircraft_layer_enabled = lambda: False
    calls: list[str] = []
    dummy.start_background_sky_data_update = lambda **kwargs: calls.append(
        f"sky:{kwargs.get('is_initial_load')}"
    )
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(
        f"cloud:{kwargs.get('reason')}"
    )
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(
        f"terrain:{kwargs.get('reason')}"
    )
    dummy.start_background_urban_outline_update = lambda **kwargs: calls.append(
        f"urban:{kwargs.get('reason')}"
    )

    SkyWindow.start_initial_data_load(dummy)

    assert calls == ["sky:True"]


def test_start_initial_data_load_skips_water_when_terrain_disabled() -> None:
    dummy = SimpleNamespace()
    dummy._startup_initial_load_started = False
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.0
    dummy.terrain_horizon_opacity = 0.0
    dummy.water_overlay_opacity = 0.2
    dummy.urban_outline_opacity = 0.0
    dummy._satellite_layer_enabled = lambda: False
    dummy._aircraft_layer_enabled = lambda: False
    calls: list[str] = []
    dummy.start_background_sky_data_update = lambda **kwargs: calls.append(
        f"sky:{kwargs.get('is_initial_load')}"
    )
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(
        f"cloud:{kwargs.get('reason')}"
    )
    dummy.start_background_terrain_horizon_update = lambda **kwargs: calls.append(
        f"terrain:{kwargs.get('reason')}"
    )
    dummy.start_background_urban_outline_update = lambda **kwargs: calls.append(
        f"urban:{kwargs.get('reason')}"
    )

    SkyWindow.start_initial_data_load(dummy)

    assert calls == ["sky:True"]


def test_initial_data_load_advances_through_terrain_and_urban() -> None:
    dummy = SimpleNamespace()
    dummy._is_shutting_down = False
    dummy._startup_initial_load_started = True
    dummy._startup_initial_data_loaded = False
    dummy._startup_initial_sky_loaded = True
    dummy._startup_initial_terrain_loaded = False
    dummy._startup_initial_water_loaded = False
    dummy._startup_initial_urban_loaded = False
    dummy._startup_initial_night_light_loaded = False
    dummy.terrain_horizon_opacity = 0.2
    dummy.water_overlay_opacity = 0.2
    dummy.urban_outline_opacity = 0.2
    dummy.night_light_opacity = 0.0
    dummy.terrain_horizon_state = SimpleNamespace(profile_altaz=None)
    calls: list[str] = []
    dummy.start_background_terrain_horizon_update = lambda **kwargs: (
        calls.append(f"terrain:{kwargs.get('reason')}") or True
    )
    dummy.start_background_urban_outline_update = lambda **kwargs: (
        calls.append(f"urban:{kwargs.get('reason')}") or True
    )
    dummy._finish_initial_data_load = lambda: calls.append("finish")

    SkyWindowUpdatesMixin._continue_initial_data_load(dummy)
    assert calls == ["terrain:initial"]

    dummy._startup_initial_terrain_loaded = True
    dummy.terrain_horizon_state.profile_altaz = [(0.0, 0.0)]
    SkyWindowUpdatesMixin._continue_initial_data_load(dummy)
    assert calls == ["terrain:initial", "urban:initial"]

    dummy._startup_initial_urban_loaded = True
    SkyWindowUpdatesMixin._continue_initial_data_load(dummy)
    assert calls == ["terrain:initial", "urban:initial", "finish"]


def test_initial_data_load_does_not_wait_for_night_light_before_finish() -> None:
    dummy = SimpleNamespace()
    dummy._is_shutting_down = False
    dummy._startup_initial_load_started = True
    dummy._startup_initial_data_loaded = False
    dummy._startup_initial_sky_loaded = True
    dummy._startup_initial_terrain_loaded = True
    dummy._startup_initial_water_loaded = True
    dummy._startup_initial_urban_loaded = True
    dummy._startup_initial_night_light_loaded = False
    dummy.terrain_horizon_opacity = 0.2
    dummy.water_overlay_opacity = 0.2
    dummy.urban_outline_opacity = 0.2
    dummy.night_light_opacity = 0.2
    dummy.ridge_glow_opacity = 0.0
    calls: list[str] = []
    dummy.start_background_terrain_horizon_update = lambda **kwargs: (
        calls.append(f"terrain:{kwargs.get('reason')}") or True
    )
    dummy.start_background_urban_outline_update = lambda **kwargs: (
        calls.append(f"urban:{kwargs.get('reason')}") or True
    )
    dummy._finish_initial_data_load = lambda: calls.append("finish")

    SkyWindowUpdatesMixin._continue_initial_data_load(dummy)
    assert calls == ["finish"]


def test_initial_data_load_does_not_wait_for_ridge_glow_before_finish() -> None:
    dummy = SimpleNamespace()
    dummy._is_shutting_down = False
    dummy._startup_initial_load_started = True
    dummy._startup_initial_data_loaded = False
    dummy._startup_initial_sky_loaded = True
    dummy._startup_initial_terrain_loaded = True
    dummy._startup_initial_water_loaded = True
    dummy._startup_initial_urban_loaded = True
    dummy._startup_initial_night_light_loaded = False
    dummy.terrain_horizon_opacity = 0.2
    dummy.water_overlay_opacity = 0.2
    dummy.urban_outline_opacity = 0.2
    dummy.night_light_opacity = 0.0
    dummy.ridge_glow_opacity = 0.2
    calls: list[str] = []
    dummy.start_background_terrain_horizon_update = lambda **kwargs: (
        calls.append(f"terrain:{kwargs.get('reason')}") or True
    )
    dummy.start_background_urban_outline_update = lambda **kwargs: (
        calls.append(f"urban:{kwargs.get('reason')}") or True
    )
    dummy._finish_initial_data_load = lambda: calls.append("finish")

    SkyWindowUpdatesMixin._continue_initial_data_load(dummy)
    assert calls == ["finish"]


def test_sky_data_ready_marks_startup_night_light_loaded_at_night() -> None:
    dummy = SimpleNamespace()
    dummy._is_shutting_down = False
    dummy._startup_initial_load_started = True
    dummy._startup_initial_data_loaded = False
    dummy._startup_initial_sky_loaded = False
    dummy._startup_initial_terrain_loaded = True
    dummy._startup_initial_water_loaded = True
    dummy._startup_initial_urban_loaded = True
    dummy._startup_initial_night_light_loaded = False
    dummy.terrain_horizon_opacity = 0.2
    dummy.night_light_opacity = 0.2
    dummy.ridge_glow_opacity = 0.0
    dummy.sky_update_interval = 60
    dummy.viewer_data = ViewerData(
        location=(0.0, 0.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(0.0, 0.0),
    )
    dummy.state = SimpleNamespace(
        viewport_interaction_release_pending=False,
        viewport_interaction_completion_reason=None,
        viewport_interaction_mode=False,
        viewport_interaction_stars=None,
        render_view_center=(0.0, 0.0),
        night_light_glow_profile=None,
        sky_update_pending=False,
        pending_star_vmag_limit=None,
        cloud_repaint_deferred=False,
        sky_next_refresh_utc=None,
        sky_disc_next_refresh_utc=None,
    )
    dummy._compositor = SimpleNamespace(invalidate=lambda: None)
    dummy.request_client_update = lambda: None
    dummy._safe_request_cloud_repaint = lambda: None
    dummy._continue_initial_data_load = lambda: None
    dummy.request_sky_data_update = lambda **kwargs: None
    dummy._sky_disc_update_interval = lambda: 60
    dummy.client_width = lambda: 640
    dummy.client_height = lambda: 480
    payload = {
        "celestial": SimpleNamespace(planets=[SimpleNamespace(name="sun", alt=-5.0)]),
        "sky_disc": object(),
        "night_light_glow_profile": object(),
        "view_center": (0.0, 0.0),
        "geometry": render_geometry.get_screen_geometry(
            640, 480, dummy.viewer_data.view_alt_deg
        ),
        "render_generation": 0,
    }

    SkyWindowUpdatesMixin._on_sky_data_calculated(dummy, payload)

    assert dummy._startup_initial_night_light_loaded is True


def test_sky_data_ready_marks_startup_ridge_glow_loaded_at_night() -> None:
    dummy = SimpleNamespace()
    dummy._is_shutting_down = False
    dummy._startup_initial_load_started = True
    dummy._startup_initial_data_loaded = False
    dummy._startup_initial_sky_loaded = False
    dummy._startup_initial_terrain_loaded = True
    dummy._startup_initial_water_loaded = True
    dummy._startup_initial_urban_loaded = True
    dummy._startup_initial_night_light_loaded = False
    dummy.terrain_horizon_opacity = 0.2
    dummy.night_light_opacity = 0.0
    dummy.ridge_glow_opacity = 0.2
    dummy.sky_update_interval = 60
    dummy.viewer_data = ViewerData(
        location=(0.0, 0.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(0.0, 0.0),
    )
    dummy.state = SimpleNamespace(
        viewport_interaction_release_pending=False,
        viewport_interaction_completion_reason=None,
        viewport_interaction_mode=False,
        viewport_interaction_stars=None,
        render_view_center=(0.0, 0.0),
        night_light_glow_profile=None,
        sky_update_pending=False,
        pending_star_vmag_limit=None,
        cloud_repaint_deferred=False,
        sky_next_refresh_utc=None,
        sky_disc_next_refresh_utc=None,
    )
    dummy._compositor = SimpleNamespace(invalidate=lambda: None)
    dummy.request_client_update = lambda: None
    dummy._safe_request_cloud_repaint = lambda: None
    dummy._continue_initial_data_load = lambda: None
    dummy.request_sky_data_update = lambda **kwargs: None
    dummy._sky_disc_update_interval = lambda: 60
    dummy.client_width = lambda: 640
    dummy.client_height = lambda: 480
    payload = {
        "celestial": SimpleNamespace(planets=[SimpleNamespace(name="sun", alt=-5.0)]),
        "sky_disc": object(),
        "night_light_glow_profile": object(),
        "view_center": (0.0, 0.0),
        "geometry": render_geometry.get_screen_geometry(
            640, 480, dummy.viewer_data.view_alt_deg
        ),
        "render_generation": 0,
    }

    SkyWindowUpdatesMixin._on_sky_data_calculated(dummy, payload)

    assert dummy._startup_initial_night_light_loaded is True


def test_terrain_horizon_ready_triggers_water_overlay_update() -> None:
    dummy = SimpleNamespace()
    dummy._is_shutting_down = False
    dummy._startup_initial_load_started = True
    dummy._startup_initial_data_loaded = False
    dummy._startup_initial_sky_loaded = True
    dummy._startup_initial_terrain_loaded = False
    dummy._startup_initial_water_loaded = False
    dummy._startup_initial_urban_loaded = False
    dummy.terrain_horizon_opacity = 0.2
    dummy.water_overlay_opacity = 0.2
    dummy.terrain_horizon_state = SimpleNamespace(
        ground_elevation_m=None,
        set_result=lambda *args, **kwargs: None,
    )
    dummy.viewer_data = ViewerData(
        location=(0.0, 0.0), timezone_name="UTC", city_name="Test"
    )
    dummy.state = SimpleNamespace(
        terrain_horizon_profile=None,
        terrain_horizon_profile_distances_m=None,
        terrain_secondary_ridges_altaz_layers=None,
        terrain_secondary_ridges_distances_m_layers=None,
    )
    dummy._refresh_water_overlay_active_dots = lambda: None
    dummy._sync_water_overlay_action_enabled = lambda: None
    dummy._compositor = SimpleNamespace(invalidate=lambda: None)
    dummy.request_client_update = lambda: None
    dummy.request_sky_data_update = lambda **_kwargs: None
    calls: list[str] = []
    dummy.start_background_water_overlay_update = lambda **kwargs: (
        calls.append(str(kwargs.get("reason"))) or True
    )
    dummy._continue_initial_data_load = lambda: calls.append("continue")

    SkyWindowUpdatesMixin._on_terrain_horizon_ready(
        dummy,
        {
            "profile_altaz": [(0.0, 0.0)],
            "profile_distances_m": [1.0],
            "secondary_ridges_altaz_layers": [],
            "secondary_ridges_distances_m_layers": [],
            "ground_elevation_m": 42.0,
            "source": "test",
        },
    )

    assert calls == ["initial"]


def test_post_startup_background_updates_start_cloud_immediately() -> None:
    dummy = SimpleNamespace()
    dummy._is_shutting_down = False
    dummy._post_startup_background_updates_started = False
    dummy._startup_initial_data_loaded = True
    dummy._clouddisc = object()
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 600
    dummy.state = SimpleNamespace(
        sky_next_refresh_utc=None,
        cloud_next_refresh_utc=None,
        satellite_next_refresh_utc=None,
        aircraft_next_refresh_utc=None,
    )
    dummy._tropical_cyclone_controller = None
    timer_calls: list[int] = []
    dummy._scheduler_tick_timer = SimpleNamespace(
        isActive=lambda: False,
        start=lambda ms=None: timer_calls.append(0 if ms is None else int(ms)),
    )
    dummy._satellite_layer_enabled = lambda: False
    dummy._aircraft_layer_enabled = lambda: False
    calls: list[str] = []
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy._on_scheduler_tick = lambda: calls.append("tick")

    SkyWindow._start_post_startup_background_updates(dummy)

    assert calls == ["initial", "tick"]
    assert timer_calls == [0]
    assert dummy.state.sky_next_refresh_utc is not None
    assert dummy.state.cloud_next_refresh_utc is None


def test_toggle_earth_guide_respects_cli_lockout() -> None:
    dummy = SimpleNamespace()
    dummy._earth_guide_gui_allowed = False
    dummy.earth_guide_opacity = 0.0
    dummy._earth_guide_opacity_when_enabled = 0.25
    dummy._action_toggle_earth_guide = _DummyAction(False)
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: (_ for _ in ()).throw(
            AssertionError("should not invalidate")
        )
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


def test_toggle_simplified_view_flips_state_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(
        simplified_view_enabled=False,
        simplified_view_labels_enabled=True,
    )
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request")
    dummy._simplified_view_enabled = lambda: bool(
        getattr(dummy.state, "simplified_view_enabled", False)
    )
    dummy._simplified_view_labels_enabled = lambda: bool(
        getattr(dummy.state, "simplified_view_labels_enabled", True)
    )

    SkyWindow.toggle_simplified_view(dummy)

    assert dummy.state.simplified_view_enabled is True
    assert dummy.state.simplified_view_labels_enabled is False
    assert calls == ["request"]

    SkyWindow.toggle_simplified_view(dummy)

    assert dummy.state.simplified_view_enabled is True
    assert dummy.state.simplified_view_labels_enabled is True
    assert calls == ["request", "request"]

    SkyWindow.toggle_simplified_view(dummy)

    assert dummy.state.simplified_view_enabled is False
    assert dummy.state.simplified_view_labels_enabled is True
    assert calls == ["request", "request", "request"]


def test_resolve_simplified_view_mode_matrix() -> None:
    assert (
        resolve_simplified_view_mode(
            base_enabled=False,
            labels_enabled=True,
        )
        == "normal"
    )
    assert (
        resolve_simplified_view_mode(
            base_enabled=True,
            labels_enabled=False,
        )
        == "nolabels"
    )
    assert (
        resolve_simplified_view_mode(
            base_enabled=True,
            labels_enabled=True,
        )
        == "labels"
    )


def test_handle_client_key_press_triggers_simplified_view_toggle_for_space() -> None:
    dummy = SimpleNamespace()
    dummy._startup_input_blocked = lambda: False
    dummy.state = SimpleNamespace(simplified_view_enabled=False)
    calls: list[str] = []
    dummy.toggle_simplified_view = lambda: calls.append("simplified")

    event = _DummyKeyEvent(window_module.Qt.Key.Key_Space)

    SkyWindow._handle_client_key_press(dummy, event)

    assert calls == ["simplified"]
    assert event.accepted is True


def test_status_line_message_returns_simplified_label_when_enabled() -> None:
    dummy = SimpleNamespace()
    dummy._effective_simplified_view_mode = lambda: "labels"

    assert (
        SkyWindowUpdatesMixin._status_line_message(dummy)
        == "Simplified: with labels [Space]"
    )


def test_status_line_message_returns_simplified_no_labels_message() -> None:
    dummy = SimpleNamespace()
    dummy._effective_simplified_view_mode = lambda: "nolabels"

    assert (
        SkyWindowUpdatesMixin._status_line_message(dummy)
        == "Simplified: no labels [Space]"
    )


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


def test_terrain_controller_reports_unavailable_when_dem_is_missing(
    tmp_path, monkeypatch, caplog
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

    with caplog.at_level("WARNING", logger=terrain_controller_module.__name__):
        controller._run_update(
            lat=20.0, lon=-30.0, observer_height_m=1.7, reason="initial"
        )

    assert ready_payloads == []
    assert failed_payloads == [{"banner": "Terrain horizon: unavailable"}]
    assert "Traceback" not in caplog.text


def test_terrain_controller_uses_shared_dem_scan_step(tmp_path) -> None:
    controller = TerrainHorizonController(cache_dir=tmp_path)

    assert controller._sample_step_m == pytest.approx(
        DEFAULT_TERRAIN_DISTANCE_SAMPLE_STEP_M
    )


def test_toggle_sky_disc_enables_gradient_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy._sky_disc_gui_allowed = True
    dummy.sky_disc_alpha = 0.0
    dummy._sky_disc_alpha_when_enabled = SKY_OPACITY_DEFAULT
    dummy._action_toggle_sky_disc = _DummyAction(False)
    calls: list[str] = []
    dummy.request_client_update = lambda: calls.append("request-client")
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy.request_sky_data_update = lambda: calls.append("request")
    dummy.update = lambda: calls.append("update")

    SkyWindow.toggle_sky_disc(dummy)

    assert dummy.sky_disc_alpha == SKY_OPACITY_DEFAULT
    assert dummy._action_toggle_sky_disc.isChecked() is True
    assert calls == ["invalidate", "request", "request-client"]


def test_toggle_sky_disc_switches_to_flat_disc_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy._sky_disc_gui_allowed = True
    dummy.sky_disc_alpha = SKY_OPACITY_DEFAULT
    dummy._sky_disc_alpha_when_enabled = SKY_OPACITY_DEFAULT
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
    monkeypatch.setattr(
        window_module.QTimer,
        "singleShot",
        lambda ms, func: func(),
    )

    dummy = SimpleNamespace()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SimpleNamespace(
        jump_highlight_name=None,
        jump_highlight_altaz=None,
        jump_highlight_until_ms=0.0,
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=False,
        persistent_search_target=None,
        persistent_search_reference_time_utc=None,
        persistent_search_next_refresh_utc=None,
        persistent_search_last_refresh_utc=None,
        persistent_search_last_error=None,
    )
    sync_calls: list[str] = []
    dummy.request_client_update = lambda: sync_calls.append("request-client")
    dummy._sync_view_altitude_actions = lambda: sync_calls.append("sync")
    dummy._current_time_obj = lambda: object()

    def _begin_viewport_interaction_mode(*args, **kwargs) -> None:
        sync_calls.append("begin-viewport")
        dummy.state.viewport_interaction_mode = True

    dummy._begin_viewport_interaction_mode = _begin_viewport_interaction_mode
    dummy._update_viewport_interaction_stars = lambda: sync_calls.append("update-stars")

    def _end_viewport_interaction_mode(*args, **kwargs) -> None:
        sync_calls.append("request")
        dummy.state.viewport_interaction_mode = False

    dummy._end_viewport_interaction_mode = _end_viewport_interaction_mode
    dummy._finalize_view_direction_change = lambda: (
        SkyWindow._finalize_view_direction_change(dummy)
    )
    dummy.request_sky_data_update = lambda: sync_calls.append("request")
    dummy.update = lambda: sync_calls.append("update")
    dummy._clear_persistent_search = lambda: sync_calls.append("clear")
    dummy._search_view_center_base = (20.0, 30.0)
    dummy._search_view_center_alt_specified = False
    dummy._search_view_center_az_specified = False
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    target = SearchJumpTarget(
        label="Circlet",
        kind="star",
        sort_key=(0.0, "circlet"),
        ra_hours=1.0,
        dec_deg=2.0,
    )
    SkyWindow._jump_to_search_target(dummy, target)

    assert dummy.viewer_data.view_center == (-12.5, 210.0)
    assert dummy.state.jump_highlight_name == "Circlet"
    assert dummy.state.jump_highlight_altaz == (-12.5, 210.0)
    assert dummy.state.jump_highlight_until_ms > 0.0
    assert sync_calls == [
        "begin-viewport",
        "sync",
        "update-stars",
        "request-client",
        "clear",
        "request",
    ]


def test_jump_to_search_target_can_keep_marker_for_local_star(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        window_module, "radec_to_altaz", lambda *_args, **_kwargs: (14.25, 87.0)
    )
    monkeypatch.setattr(
        window_module.QTimer,
        "singleShot",
        lambda ms, func: func(),
    )

    dummy = SimpleNamespace()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SimpleNamespace(
        jump_highlight_name=None,
        jump_highlight_altaz=None,
        jump_highlight_until_ms=0.0,
        persistent_search_target=None,
        persistent_search_reference_time_utc=None,
        persistent_search_next_refresh_utc=None,
        persistent_search_last_refresh_utc=None,
        persistent_search_last_error=None,
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=False,
    )
    sync_calls: list[str] = []
    dummy.request_client_update = lambda: sync_calls.append("request-client")
    dummy._sync_view_altitude_actions = lambda: sync_calls.append("sync")
    dummy._current_time_obj = lambda: object()

    def _begin_viewport_interaction_mode(*args, **kwargs) -> None:
        sync_calls.append("begin-viewport")
        dummy.state.viewport_interaction_mode = True

    dummy._begin_viewport_interaction_mode = _begin_viewport_interaction_mode
    dummy._update_viewport_interaction_stars = lambda: sync_calls.append("update-stars")

    def _end_viewport_interaction_mode(*args, **kwargs) -> None:
        sync_calls.append("request")
        dummy.state.viewport_interaction_mode = False

    dummy._end_viewport_interaction_mode = _end_viewport_interaction_mode
    dummy._finalize_view_direction_change = lambda: (
        SkyWindow._finalize_view_direction_change(dummy)
    )
    dummy.request_sky_data_update = lambda: sync_calls.append("request")
    dummy._clear_persistent_search = lambda: sync_calls.append("clear")
    dummy._log_persistent_search_target_update = lambda **_kwargs: sync_calls.append(
        "log"
    )
    dummy._search_view_center_base = (20.0, 30.0)
    dummy._search_view_center_alt_specified = False
    dummy._search_view_center_az_specified = False
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    target = SearchJumpTarget(
        label="Sirius",
        kind="star",
        sort_key=(0.0, "sirius"),
        ra_hours=6.75,
        dec_deg=-16.7,
        persistent_keep_marker=True,
    )
    SkyWindow._jump_to_search_target(dummy, target)

    assert dummy.viewer_data.view_center == (14.25, 87.0)
    assert dummy.state.render_view_center == (14.25, 87.0)
    assert dummy.state.jump_highlight_name == "Sirius"
    assert dummy.state.jump_highlight_altaz == (14.25, 87.0)
    assert dummy.state.jump_highlight_until_ms > 0.0
    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.label == "Sirius"
    assert dummy.state.persistent_search_target.persistent_keep_marker is True
    assert dummy.state.persistent_search_next_refresh_utc is None
    assert sync_calls == [
        "begin-viewport",
        "sync",
        "update-stars",
        "request-client",
        "log",
        "request",
    ]


def test_rotate_view_in_orientation_mode_updates_render_center_without_full_refresh() -> (
    None
):
    dummy = SimpleNamespace()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
    )
    dummy.state = SimpleNamespace(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=False,
        interaction_mode=False,
    )
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
    dummy._viewport_rotation_keys_down = set()
    dummy._startup_input_blocked = lambda: False

    SkyWindow._rotate_view(dummy, d_alt=5.0, d_az=15.0, interactive_viewport=True)

    assert dummy.viewer_data.view_center == (25.0, 45.0)
    assert dummy.state.render_view_center == (25.0, 45.0)
    assert calls == ["begin-viewport", "sync", "bright-stars", "request-client"]


def test_end_viewport_interaction_mode_requests_full_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(
        viewport_interaction_mode=True, viewport_interaction_stars=object()
    )
    dummy.menu_button = SimpleNamespace(setVisible=lambda *_args, **_kwargs: None)
    dummy._startup_initial_load_started = True
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
    dummy.request_cloud_projection_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy.reproject_tropical_cyclone_overlay = lambda **kwargs: calls.append("cyclone")
    dummy.update = lambda: calls.append("update")

    SkyWindow._end_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    assert calls == [
        "viewport-interaction-idle",
        "view-change-idle",
        "view-change-idle",
        "cyclone",
        "request-client",
    ]


def test_begin_viewport_interaction_mode_clears_cloud_buffers_and_invalidates_old_render() -> (
    None
):
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(viewport_interaction_mode=False, sky_disc_image=None)
    dummy.menu_button = None
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="render-key",
        request_id=42,
        missing_mask_key=99,
    )
    dummy.geosatellite_state = SimpleNamespace(
        image=None,
        missing_mask=None,
        cloud_amount_field=None,
        altaz_grid=None,
        render_key=None,
        request_id=None,
        missing_mask_key=None,
    )
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: calls.append("invalidate-compositor")
    )
    dummy._cloud_controller = SimpleNamespace(
        invalidate_pending_render_results=lambda: calls.append("invalidate-cloud")
    )
    dummy._viewport_interaction_idle_timer = SimpleNamespace(
        isActive=lambda: False,
        stop=lambda: calls.append("stop-timer"),
        start=lambda: calls.append("start-timer"),
    )
    dummy._startup_initial_load_started = True

    SkyWindow._begin_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is True
    assert dummy.cloud_state.image is None
    assert dummy.cloud_state.missing_mask is None
    assert dummy.cloud_state.cloud_amount_field is None
    assert dummy.cloud_state.render_key is None
    assert dummy.cloud_state.request_id is None
    assert dummy.cloud_state.missing_mask_key is None
    assert calls == ["invalidate-compositor", "invalidate-cloud", "start-timer"]


def test_begin_viewport_interaction_mode_preserves_geo_satellite_cloud_buffers() -> None:
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(viewport_interaction_mode=False)
    dummy.menu_button = None
    geo_image = object()
    geo_mask = object()
    geo_field = object()
    geo_grid = object()
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="cloud-render",
        request_id=1,
        missing_mask_key=2,
    )
    dummy.geosatellite_state = SimpleNamespace(
        image=geo_image,
        missing_mask=geo_mask,
        cloud_amount_field=geo_field,
        altaz_grid=geo_grid,
        render_key="geo-render",
        request_id=2,
        missing_mask_key=3,
    )
    dummy._geo_satellite_mode_active = lambda: True
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: calls.append("invalidate-compositor")
    )
    dummy._cloud_controller = SimpleNamespace(
        invalidate_pending_render_results=lambda: calls.append("invalidate-cloud")
    )
    dummy._viewport_interaction_idle_timer = SimpleNamespace(
        isActive=lambda: False,
        stop=lambda: calls.append("stop-timer"),
        start=lambda: calls.append("start-timer"),
    )
    dummy._startup_initial_load_started = True

    SkyWindow._begin_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is True
    assert dummy.geosatellite_state.image is geo_image
    assert dummy.geosatellite_state.missing_mask is geo_mask
    assert dummy.geosatellite_state.cloud_amount_field is geo_field
    assert dummy.geosatellite_state.altaz_grid is geo_grid
    assert dummy.cloud_state.image is None
    assert calls == ["invalidate-compositor", "invalidate-cloud", "start-timer"]


def test_begin_viewport_interaction_mode_clears_cloud_buffers_even_while_cloud_update_is_running() -> (
    None
):
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(viewport_interaction_mode=False)
    dummy.menu_button = None
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="render-key",
        request_id=42,
        missing_mask_key=99,
    )
    dummy.geosatellite_state = SimpleNamespace(
        image=None,
        missing_mask=None,
        cloud_amount_field=None,
        altaz_grid=None,
        render_key=None,
        request_id=None,
        missing_mask_key=None,
    )
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: calls.append("invalidate-compositor")
    )
    dummy._cloud_controller = SimpleNamespace(
        has_in_flight_update=lambda: True,
        invalidate_pending_render_results=lambda: calls.append("invalidate-cloud"),
    )
    dummy._viewport_interaction_idle_timer = SimpleNamespace(
        isActive=lambda: False,
        stop=lambda: calls.append("stop-timer"),
        start=lambda: calls.append("start-timer"),
    )
    dummy._startup_initial_load_started = True

    SkyWindow._begin_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is True
    assert dummy.cloud_state.image is None
    assert dummy.cloud_state.missing_mask is None
    assert dummy.cloud_state.cloud_amount_field is None
    assert dummy.cloud_state.render_key is None
    assert dummy.cloud_state.request_id is None
    assert dummy.cloud_state.missing_mask_key is None
    assert calls == ["invalidate-compositor", "invalidate-cloud", "start-timer"]


def test_clear_cloud_render_buffers_preserves_buffers_when_requested() -> None:
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(cloud_projection_next_refresh_utc=object())
    cloud_image = object()
    cloud_mask = object()
    cloud_field = object()
    geo_image = object()
    geo_mask = object()
    geo_field = object()
    geo_grid = object()
    dummy.cloud_state = SimpleNamespace(
        image=cloud_image,
        missing_mask=cloud_mask,
        cloud_amount_field=cloud_field,
        render_key="cloud-render",
        request_id=1,
        missing_mask_key="cloud-mask",
    )
    dummy.geosatellite_state = SimpleNamespace(
        image=geo_image,
        missing_mask=geo_mask,
        cloud_amount_field=geo_field,
        altaz_grid=geo_grid,
        render_key="geo-render",
        request_id=2,
        missing_mask_key="geo-mask",
    )

    cleared = SkyWindow._clear_cloud_render_buffers(
        dummy,
        preserve_cloud_buffers=True,
    )

    assert cleared is False
    assert dummy.cloud_state.image is cloud_image
    assert dummy.cloud_state.missing_mask is cloud_mask
    assert dummy.cloud_state.cloud_amount_field is cloud_field
    assert dummy.cloud_state.render_key == "cloud-render"
    assert dummy.cloud_state.request_id == 1
    assert dummy.cloud_state.missing_mask_key == "cloud-mask"
    assert dummy.geosatellite_state.image is geo_image
    assert dummy.geosatellite_state.missing_mask is geo_mask
    assert dummy.geosatellite_state.cloud_amount_field is geo_field
    assert dummy.geosatellite_state.altaz_grid is geo_grid
    assert dummy.geosatellite_state.render_key == "geo-render"
    assert dummy.geosatellite_state.request_id == 2
    assert dummy.geosatellite_state.missing_mask_key == "geo-mask"
    assert dummy.state.cloud_projection_next_refresh_utc is not None


def test_begin_viewport_interaction_mode_restarts_idle_timer_when_already_active() -> (
    None
):
    calls: list[str] = []
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(viewport_interaction_mode=True)
    dummy.menu_button = None
    dummy.cloud_state = SimpleNamespace(
        image=None,
        missing_mask=None,
        cloud_amount_field=None,
        render_key=None,
        request_id=None,
        missing_mask_key=None,
    )
    dummy.geosatellite_state = SimpleNamespace(
        image=None,
        missing_mask=None,
        cloud_amount_field=None,
        altaz_grid=None,
        render_key=None,
        request_id=None,
        missing_mask_key=None,
    )
    dummy._compositor = SimpleNamespace(
        invalidate=lambda: calls.append("invalidate-compositor")
    )
    dummy._cloud_controller = SimpleNamespace(
        invalidate_pending_render_results=lambda: calls.append("invalidate-cloud")
    )

    class _Timer:
        def __init__(self) -> None:
            self.active = True

        def isActive(self) -> bool:
            calls.append("check-active")
            return self.active

        def stop(self) -> None:
            calls.append("stop-timer")
            self.active = False

        def start(self) -> None:
            calls.append("start-timer")
            self.active = True

    dummy._viewport_interaction_idle_timer = _Timer()
    dummy._startup_initial_load_started = True

    SkyWindow._begin_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is True
    assert calls == [
        "invalidate-cloud",
        "check-active",
        "stop-timer",
        "start-timer",
    ]


def test_handle_client_resize_discards_cached_sky_disc_and_requests_refresh() -> None:
    calls: list[str] = []
    sky_disc_image = object()
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(
        viewport_interaction_mode=False,
        sky_disc_image=sky_disc_image,
    )
    dummy.menu_button = None
    dummy._disc_generation = 0
    dummy._frameless_frame = object()
    dummy.menu_button = SimpleNamespace(
        raise_=lambda: calls.append("raise-menu"),
        setVisible=lambda *_args, **_kwargs: None,
    )
    dummy.size_grip = SimpleNamespace(raise_=lambda: calls.append("raise-grip"))
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="render-key",
        request_id=42,
        missing_mask_key=99,
    )
    dummy.geosatellite_state = SimpleNamespace(
        image=None,
        missing_mask=None,
        cloud_amount_field=None,
        altaz_grid=None,
        render_key=None,
        request_id=None,
        missing_mask_key=None,
    )
    dummy._compositor = SimpleNamespace(invalidate=lambda: calls.append("invalidate"))
    dummy._cloud_controller = SimpleNamespace(
        invalidate_pending_render_results=lambda: calls.append("invalidate-cloud")
    )
    dummy._viewport_interaction_idle_timer = SimpleNamespace(
        isActive=lambda: False,
        stop=lambda: calls.append("stop-timer"),
        start=lambda: calls.append("start-timer"),
    )
    dummy._startup_initial_load_started = True
    dummy.water_overlay_opacity = 0.0
    dummy.request_sky_data_update = lambda **kwargs: calls.append(
        f"sky:{kwargs.get('reason')}:{kwargs.get('allow_during_viewport_interaction')}"
    )
    dummy.request_client_update = lambda: calls.append("client")
    dummy.start_background_cloud_update = lambda **kwargs: calls.append(
        str(kwargs.get("reason"))
    )
    dummy._raise_overlay_widgets = lambda: (
        calls.append("raise-menu"),
        calls.append("raise-grip"),
    )
    dummy._discard_stale_disc_images = lambda: SkyWindow._discard_stale_disc_images(  # type: ignore[attr-defined]
        dummy
    )
    dummy._begin_viewport_interaction_mode = lambda *args, **kwargs: (
        SkyWindow._begin_viewport_interaction_mode(dummy, *args, **kwargs)
    )

    SkyWindow._handle_client_resize(dummy, SimpleNamespace())

    assert dummy._disc_generation == 1
    assert dummy.state.sky_disc_image is None
    assert dummy.cloud_state.image is None
    assert dummy.cloud_state.missing_mask is None
    assert dummy.cloud_state.cloud_amount_field is None
    assert dummy.cloud_state.render_key is None
    assert dummy.cloud_state.request_id is None
    assert dummy.cloud_state.missing_mask_key is None
    assert calls == [
        "invalidate",
        "invalidate-cloud",
        "start-timer",
        "invalidate",
        "invalidate",
        "sky:resize:True",
        "client",
        "raise-menu",
        "raise-grip",
    ]


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


def test_discard_stale_disc_images_clears_cached_cloud_buffers() -> None:
    compositor_calls: list[str] = []
    sky_disc_image = object()
    dummy = SimpleNamespace()
    dummy.state = SimpleNamespace(sky_disc_image=sky_disc_image)
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=object(),
        render_key="render-key",
        request_id=42,
        missing_mask_key=99,
    )
    dummy.geosatellite_state = SimpleNamespace(
        image=None,
        missing_mask=None,
        cloud_amount_field=None,
        altaz_grid=None,
        render_key=None,
        request_id=None,
        missing_mask_key=None,
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


def test_size_grip_style_sheet_is_transparent() -> None:
    dummy = SimpleNamespace(theme=THEME_STYLES_BY_PRESET["night"])

    style = SkyWindow._size_grip_style_sheet(dummy)

    assert "background: transparent;" in style


def test_resize_grip_widget_paints_a_visible_marker() -> None:
    img = QImage(30, 30, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    window_widgets_module._draw_resize_grip_marker(painter, QRect(0, 0, 30, 30))
    painter.end()

    arr = qimage_to_np_rgba(img)
    assert int(arr[9, 21, 3]) > 0
    assert int(arr[14, 16, 3]) > 0
    assert int(arr[19, 11, 3]) > 0


def test_update_viewport_interaction_stars_uses_bright_limit(monkeypatch) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        window_module,
        "calculate_visible_stars",
        lambda catalog, lat, lon, observer_height_m, time_obj, view_center, max_vmag, subset_indices=None, star_data_policy="scenic_view_scoped": (
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
                    "star_data_policy": star_data_policy,
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
        "star_data_policy": "scenic_view_scoped",
    }


def test_compute_star_render_upscale_factor_matches_downsampled_star_surface() -> None:
    factor = SkyWindow.compute_star_render_upscale_factor(
        disc_width_px=2400,
        expected_width_px=600,
    )

    assert factor == 2.0
