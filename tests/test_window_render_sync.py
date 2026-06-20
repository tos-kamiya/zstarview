from __future__ import annotations

import math
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from dataclasses import fields
from types import SimpleNamespace
from unittest.mock import Mock

import astropy.time
import numpy as np
import pytest
from PySide6.QtCore import QPoint, QPointF, QRect, QSize, Qt
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QResizeEvent
from PySide6.QtWidgets import QApplication, QWidget

import zstarview.gui.window as window_module
from zstarview.gui.draggable_window import DraggableWindow
import zstarview.gui.window_render as window_render_module
import zstarview.gui.window_updates as window_updates_module
import zstarview.render.guides as render_guides_module
import zstarview.render.overlay_info as render_overlay_info_module
import zstarview.render.pipeline as pipeline_module
import zstarview.render.geometry as render_geometry
import zstarview.render.solar_system as render_solar_system_module
import zstarview.render.terrain as render_terrain_module
import zstarview.render.text as render_text_module
from zstarview.gui.famous_star_shortcuts import SearchJumpTarget
from zstarview.gui.window import SkyWindow
from zstarview.gui.window_state import SkyWindowState
from zstarview.location_resolver import PlaceTargetProjection
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.satellites.types import SatelliteOverlayPoint
from zstarview.types import (
    CelestialData,
    PlanetBody,
    StarCatalogMeta,
    UrbanOutlinePolyline,
    ViewerData,
)

_app = QApplication.instance() or QApplication([])


def _viewer(
    view_center: tuple[float, float] = (45.0, 180.0),
    *,
    edge_fov_deg: float = 95.0,
    content_fov_deg: float = 110.0,
) -> SimpleNamespace:
    return SimpleNamespace(
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )


def _terrain_horizon_spec(
    *,
    opacity: float,
    line_width_scale: float = 1.0,
    fast_mode: bool,
) -> object:
    return render_terrain_module.TerrainHorizonRenderSpec(
        opacity=opacity,
        base_width=render_terrain_module.TERRAIN_HORIZON_FAST_WIDTH,
        far_base_width=render_terrain_module.TERRAIN_HORIZON_FAR_BASE_WIDTH,
        fg_alpha=render_terrain_module.terrain_horizon_line_alpha(opacity),
        line_width_scale=line_width_scale,
        color_rgb=render_terrain_module.TERRAIN_HORIZON_LINE_COLOR,
        fast_mode=fast_mode,
        distance_widths=True,
    )


class _DummyTimer:
    def __init__(self, active: bool) -> None:
        self._active = active
        self.started_with: list[int] = []

    def isActive(self) -> bool:  # noqa: N802 - Qt naming
        return self._active

    def start(self, ms: int | None = None) -> None:
        self._active = True
        self.started_with.append(0 if ms is None else ms)


class _DummySignal:
    def __init__(self) -> None:
        self.calls = 0

    def emit(self) -> None:
        self.calls += 1


class _DummyCompositor:
    def __init__(self) -> None:
        self.invalidated = False

    def invalidate(self) -> None:
        self.invalidated = True


class _WindowStub:
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        values = self.__dict__
        self._frameless_window = values.get("_frameless_window", False)
        self.observation_info_mode = values.get("observation_info_mode", "bottom")
        self.observation_info_pinned = values.get("observation_info_pinned", False)
        self.show_observation_info = values.get("show_observation_info", True)
        self.show_dso = values.get("show_dso", False)
        self.show_asterisms = values.get("show_asterisms", False)
        self.show_guidelines = values.get("show_guidelines", True)
        self.enlarge_moon = values.get("enlarge_moon", False)
        self.star_base_radius = values.get("star_base_radius", 4.0)
        self.star_visibility_boost = values.get("star_visibility_boost", 1.0)
        self.asterism_visibility_boost = values.get("asterism_visibility_boost", 1.0)
        self.earth_guide_visibility_boost = values.get(
            "earth_guide_visibility_boost", 1.0
        )
        self.vmag_limit = values.get("vmag_limit", 6.0)
        self.cloud_disc_alpha = values.get("cloud_disc_alpha", 0.0)
        self.satellite_opacity = values.get("satellite_opacity", 0.0)
        self.aircraft_opacity = values.get("aircraft_opacity", 0.0)
        self.tropical_cyclone_opacity = values.get("tropical_cyclone_opacity", 0.25)
        self.terrain_horizon_opacity = values.get("terrain_horizon_opacity", 0.25)
        self.earth_guide_opacity = values.get("earth_guide_opacity", 0.25)
        self.urban_outline_opacity = values.get("urban_outline_opacity", 0.2)
        self.show_urban_outline_layer = values.get("show_urban_outline_layer", True)
        self.show_tropical_cyclone_overlay = values.get(
            "show_tropical_cyclone_overlay",
            False,
        )
        self._night_light_toggle_supported = values.get(
            "_night_light_toggle_supported", True
        )
        self._water_overlay_gui_allowed = values.get("_water_overlay_gui_allowed", True)
        self._tropical_cyclone_toggle_supported = values.get(
            "_tropical_cyclone_toggle_supported", True
        )
        self._terrain_horizon_gui_allowed = values.get(
            "_terrain_horizon_gui_allowed", True
        )
        self.status_line_font = values.get("status_line_font", object())
        self.theme = values.get("theme", THEME_STYLES_BY_PRESET["night"])
        self.bright_bodies_mode = values.get("bright_bodies_mode", "outline")
        self.sky_disc_style = values.get("sky_disc_style", "smooth")
        self.sky_disc_altaz_rings = values.get("sky_disc_altaz_rings", "dimalt")
        self.sky_disc_altaz_rings_hover = values.get(
            "sky_disc_altaz_rings_hover", "altaz"
        )
        self.night_light_opacity = values.get("night_light_opacity", 0.0)
        self.ridge_glow_opacity = values.get("ridge_glow_opacity", 0.03)
        self.water_overlay_opacity = values.get("water_overlay_opacity", 0.4)
        self._star_render_expected_width = values.get(
            "_star_render_expected_width", 600
        )
        self._enabled_satellite_groups = values.get(
            "_enabled_satellite_groups", ("iss",)
        )
        self._disc_generation = values.get("_disc_generation", 0)
        self._client_widget = values.get("_client_widget", None)
        self.menu_button = values.get("menu_button", None)
        self.size_grip = values.get("size_grip", None)
        self.cloud_state = values.get(
            "cloud_state",
            SimpleNamespace(image=None, missing_mask=None, cloud_amount_field=None),
        )
        self.geosatellite_state = values.get(
            "geosatellite_state",
            SimpleNamespace(image=None, missing_mask=None, cloud_amount_field=None),
        )
        self.satellite_state = values.get(
            "satellite_state",
            SimpleNamespace(element_epoch_utc=None, records_by_group=None),
        )
        self.aircraft_state = values.get(
            "aircraft_state",
            SimpleNamespace(snapshots=None),
        )
        self.viewer_data = values.get("viewer_data", None)
        self.delta_t = values.get("delta_t", timedelta(0))
        self._cloud_controller = values.get("_cloud_controller", None)
        self._geosatellite_controller = values.get("_geosatellite_controller", None)
        self._geo_satellite_enabled = values.get("_geo_satellite_enabled", False)
        self._sky_worker = values.get("_sky_worker", None)
        self._satellite_controller = values.get("_satellite_controller", None)
        self._aircraft_controller = values.get("_aircraft_controller", None)
        self._jpl_small_body_controller = values.get("_jpl_small_body_controller", None)
        self._terrain_horizon_controller = values.get(
            "_terrain_horizon_controller", None
        )
        self._water_overlay_controller = values.get("_water_overlay_controller", None)
        self._urban_outline_controller = values.get("_urban_outline_controller", None)
        self._ephemeris = values.get("_ephemeris", None)
        self._sky_disc_alpha_when_enabled = values.get(
            "_sky_disc_alpha_when_enabled", 0.1
        )
        self._terrain_horizon_opacity_when_enabled = values.get(
            "_terrain_horizon_opacity_when_enabled", 0.25
        )
        self._water_overlay_opacity_when_enabled = values.get(
            "_water_overlay_opacity_when_enabled", 0.12
        )
        self.start_background_water_overlay_update = values.get(
            "start_background_water_overlay_update",
            lambda **_kwargs: False,
        )
        self._night_light_opacity_when_enabled = values.get(
            "_night_light_opacity_when_enabled", 0.06
        )
        self._urban_outline_opacity_when_enabled = values.get(
            "_urban_outline_opacity_when_enabled", 0.2
        )
        self._action_toggle_water_overlay = values.get(
            "_action_toggle_water_overlay", None
        )
        self._startup_initial_load_started = values.get(
            "_startup_initial_load_started", True
        )
        self._startup_initial_data_loaded = values.get(
            "_startup_initial_data_loaded", False
        )
        self._startup_input_blocked_state = values.get(
            "_startup_input_blocked_state", False
        )
        self._search_view_center_base = values.get(
            "_search_view_center_base", (20.0, 30.0)
        )
        self._search_view_center_alt_specified = values.get(
            "_search_view_center_alt_specified", False
        )
        self._search_view_center_az_specified = values.get(
            "_search_view_center_az_specified", False
        )
        self._viewport_rotation_keys_down = values.get(
            "_viewport_rotation_keys_down", set()
        )
        self._frame_cache_key = values.get("_frame_cache_key", None)
        self._frame_cache_image = values.get("_frame_cache_image", None)
        self._present_frame_cache_key = values.get("_present_frame_cache_key", None)
        self._present_frame_cache_image = values.get("_present_frame_cache_image", None)
        self._cached_base_label_candidates = values.get(
            "_cached_base_label_candidates", []
        )
        self.terrain_horizon_state = values.get(
            "terrain_horizon_state",
            SimpleNamespace(gound_elevation_m=None, ground_elevation_m=None),
        )
        self.state = values.get("state", None)
        self.tropical_cyclone_state = values.get(
            "tropical_cyclone_state",
            SimpleNamespace(
                snapshots=None,
                snapshot_collection=None,
                banner_text=None,
                next_check_utc=None,
                next_refresh_utc=None,
            ),
        )

    def _geo_satellite_mode_active(self) -> bool:
        return bool(
            self._geo_satellite_enabled and self._geosatellite_controller is not None
        )

    def _simplified_view_enabled(self) -> bool:
        state = self.state
        return (
            bool(getattr(state, "simplified_view_enabled", False))
            if state is not None
            else False
        )

    def _simplified_view_labels_enabled(self) -> bool:
        state = self.state
        return (
            bool(getattr(state, "simplified_view_labels_enabled", True))
            if state is not None
            else True
        )

    def _render_cache_stamp(self, value):
        if value is None:
            return None
        if hasattr(value, "cacheKey"):
            return int(value.cacheKey())
        return id(value)

    def _render_cloud_state(self):
        if self._geo_satellite_mode_active():
            return self.geosatellite_state
        return self.cloud_state

    def client_width(self) -> int:
        width = self.width
        return int(width() if callable(width) else width)

    def client_height(self) -> int:
        height = self.height
        return int(height() if callable(height) else height)

    def client_size(self):
        size = self.size
        return size() if callable(size) else size

    def client_rect(self):
        rect = self.rect
        return rect() if callable(rect) else rect

    def request_client_update(self) -> None:
        update = self.update
        if callable(update):
            update()

    def request_cloud_projection_update(self, *_args, **_kwargs) -> None:
        return None

    def reproject_tropical_cyclone_overlay(
        self,
        *_args,
        **_kwargs,
    ) -> None:
        return None

    def start_background_tropical_cyclone_update(self, *_args, **_kwargs) -> bool:
        return False

    def _begin_viewport_interaction_mode(self, *args, **kwargs) -> None:
        begin = self.__dict__.get("_begin_viewport_interaction_mode")
        if callable(begin):
            begin(*args, **kwargs)
            return
        state = self.__dict__.get("state")
        if state is not None:
            setattr(state, "viewport_interaction_mode", True)

    def _update_viewport_interaction_stars(self) -> None:
        update = self.__dict__.get("_update_viewport_interaction_stars")
        if callable(update):
            update()

    def _end_viewport_interaction_mode(self, *args, **kwargs) -> None:
        end = self.__dict__.get("_end_viewport_interaction_mode")
        if callable(end):
            end(*args, **kwargs)
            return
        window_module.SkyWindow._end_viewport_interaction_mode(self, *args, **kwargs)

    def _sync_viewport_interaction_chrome_visibility(self) -> None:
        sync = self.__dict__.get("_sync_viewport_interaction_chrome_visibility")
        if callable(sync):
            sync()
            return
        menu_button = self.__dict__.get("menu_button")
        state = self.__dict__.get("state")
        if menu_button is not None and state is not None:
            menu_button.setVisible(not bool(state.viewport_interaction_mode))

    def _target_time_utc(self):
        target_time_utc = self.__dict__.get("_target_time_utc")
        if callable(target_time_utc):
            return target_time_utc()
        return datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    def _current_time_obj(self):
        value = self.__dict__.get("_current_time_obj")
        if callable(value):
            return value()
        return astropy.time.Time("2026-04-18T12:00:00", scale="utc")

    def _build_horizons_command(self, target: dict[str, object], *, group: str) -> str:
        spkid = str(target.get("spkid", "")).strip()
        if spkid:
            if group == "sb":
                return f"DES={spkid};"
            return spkid
        pdes = str(target.get("pdes", "")).strip()
        if pdes:
            if group == "sb":
                return f"DES={pdes};"
            return pdes
        name = str(target.get("name", "")).strip()
        if name:
            if group == "sb":
                return name if name.endswith(";") else f"{name};"
            return name
        return ""

    def _extract_horizons_altaz(
        self, rows: list[list[str]]
    ) -> tuple[float, float] | None:
        for row in rows:
            numeric_values: list[float] = []
            for value in row:
                try:
                    numeric_values.append(float(str(value).strip()))
                except (TypeError, ValueError):
                    continue
            if len(numeric_values) >= 2:
                return numeric_values[1], numeric_values[0]
        return None

    def _clear_persistent_search(self) -> None:
        if self.state is None:
            return
        self.state.persistent_search_target = None
        self.state.persistent_search_reference_time_utc = None
        self.state.persistent_search_next_refresh_utc = None
        self.state.persistent_search_last_refresh_utc = None
        self.state.persistent_search_last_error = None

    def _schedule_persistent_search_refresh(self) -> None:
        return None

    def _log_persistent_search_target_update(self, **kwargs) -> None:
        window_module.SkyWindowCoreMixin._log_persistent_search_target_update(
            self,
            **kwargs,
        )

    def _viewport_interaction_active(self) -> bool:
        state = self.state
        if state is None:
            return False
        return bool(state.viewport_interaction_mode or state.interaction_mode)

    def _startup_input_blocked(self) -> bool:
        return bool(self._startup_input_blocked_state)

    def _continue_initial_data_load(self) -> None:
        return None

    def _startup_splash_visible(self) -> bool:
        overlay = self.__dict__.get("_startup_log_overlay")
        return bool(overlay is not None and overlay.isVisible())

    def _water_overlay_action_enabled(self) -> bool:
        return (
            bool(self._water_overlay_gui_allowed) and self.terrain_horizon_opacity > 0.0
        )

    def _sync_water_overlay_action_enabled(self) -> None:
        if self._action_toggle_water_overlay is not None:
            self._action_toggle_water_overlay.setEnabled(
                self._water_overlay_action_enabled()
            )


def _make_scene(
    *,
    viewer: ViewerData | None = None,
    celestial_data: CelestialData | object | None = None,
    terrain_horizon_profile: object | None = None,
    terrain_horizon_profile_distances_m: object | None = None,
    terrain_secondary_ridges_altaz_layers: object | None = None,
    terrain_secondary_ridges_distances_m_layers: object | None = None,
    urban_outlines: object | None = None,
) -> pipeline_module.RenderSceneData:
    if viewer is None:
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        )
    if celestial_data is None:
        celestial_data = CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[],
            stars={
                "name": [],
                "source_id": [],
                "alt": [],
                "az": [],
                "vmag": [],
                "bv": [],
                "size_factor": [],
                "color_factor_base": [],
            },
            deep_sky_objects={},
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
        )
    return pipeline_module.RenderSceneData(
        viewer=viewer,
        celestial_data=celestial_data,
        sky_disc_image=None,
        cloud_image=None,
        cloud_missing_mask=None,
        cloud_amount_field=None,
        terrain_horizon_profile=terrain_horizon_profile,
        terrain_horizon_profile_distances_m=terrain_horizon_profile_distances_m,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
        urban_outlines=urban_outlines,
        satellite_records_by_group=None,
        aircraft_snapshots=None,
    )


def _make_frame(
    scene: pipeline_module.RenderSceneData,
    geometry,
    viewport_rect,
) -> pipeline_module.FrameContext:
    return pipeline_module.FrameContext(
        viewer=scene.viewer,
        time_obj=scene.time_obj,
        geometry=geometry,
        viewport_rect=viewport_rect,
    )


def _make_style(**overrides) -> pipeline_module.RenderStyle:
    values = {
        "visual_preset": "night",
        "text_font": object(),
        "status_line_font": object(),
        "show_background_gradient": True,
        "show_custom_window_frame": False,
        "show_observation_info": True,
        "show_dso": False,
        "show_asterisms": False,
        "show_guidelines": True,
        "enlarge_moon": False,
        "bright_bodies_mode": "outline",
        "star_base_radius": 4.0,
        "star_visibility_boost": 1.0,
        "asterism_visibility_boost": 1.0,
        "earth_guide_visibility_boost": 1.0,
        "vmag_limit": 6.0,
        "sky_disc_altaz_rings": "dimalt",
        "sky_disc_altaz_rings_hover": "altaz",
        "cloud_disc_alpha": 0.0,
        "satellite_opacity": 0.0,
        "terrain_horizon_opacity": 0.25,
        "earth_guide_opacity": 0.25,
        "urban_outline_opacity": 0.2,
        "show_urban_outline_layer": True,
        "aircraft_opacity": 0.0,
        "tropical_cyclone_opacity": 0.4,
        "show_tropical_cyclone_overlay": True,
        "star_render_expected_width": 600,
    }
    allowed = {field.name for field in fields(pipeline_module.RenderStyle)}
    values.update({key: value for key, value in overrides.items() if key in allowed})
    if "theme" not in overrides:
        values["theme"] = pipeline_module.THEME_STYLES_BY_PRESET.get(
            values["visual_preset"],
            pipeline_module.THEME_STYLES_BY_PRESET["night"],
        )
    return pipeline_module.RenderStyle(**values)


def _make_hud(**overrides) -> pipeline_module.RenderHudState:
    values = {
        "mouse_pos": None,
        "overlay_info_bottom_left": False,
        "viewport_interaction_mode": False,
        "viewport_interaction_stars": None,
        "simplified_view_enabled": False,
        "simplified_view_labels_enabled": True,
        "status_message": None,
    }
    values.update(overrides)
    return pipeline_module.RenderHudState(**values)


def test_viewer_data_for_render_uses_render_view_center() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=123.0,
        height_add_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(60.0, 210.0))

    got = SkyWindow._viewer_data_for_render(dummy)
    assert got.view_center == (60.0, 210.0)
    assert got.location == (35.0, 139.0)
    assert got.observer_height_m == 123.0
    assert got.height_add_m == 1.7
    assert got.location_height_label is None
    assert got.location_height_m == 0.0


def test_render_style_uses_window_observation_info_toggle() -> None:
    dummy = _WindowStub()
    dummy.visual_preset = "night"
    dummy.text_font = object()
    dummy.status_line_font = object()
    dummy._frameless_window = False
    dummy.show_observation_info = False
    dummy.show_dso = True
    dummy.show_asterisms = False
    dummy.show_guidelines = True
    dummy.enlarge_moon = False
    dummy.star_base_radius = 4.0
    dummy.star_visibility_boost = 1.0
    dummy.vmag_limit = 6.0
    dummy.cloud_disc_alpha = 0.0
    dummy.satellite_opacity = 0.0
    dummy.terrain_horizon_opacity = 0.25
    dummy.urban_outline_opacity = 0.2
    dummy.ridge_glow_opacity = 0.03
    dummy.show_urban_outline_layer = True
    dummy.aircraft_opacity = 0.0
    dummy._star_render_expected_width = 600

    style = SkyWindow._render_style(dummy)

    assert style.show_custom_window_frame is False
    assert style.show_observation_info is False


def test_render_hud_state_uses_upper_third_to_switch_overlay_to_bottom_left() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = QPoint(10, 20)
    dummy.state.overlay_info_bottom_left = False
    dummy.height = lambda: 300
    dummy._status_line_message = lambda: None

    hud = SkyWindow._render_hud_state(dummy)

    assert hud.overlay_info_bottom_left is True
    assert dummy.state.overlay_info_bottom_left is True


def test_render_hud_state_keeps_overlay_anchor_in_middle_third() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = QPoint(10, 150)
    dummy.state.overlay_info_bottom_left = True
    dummy.height = lambda: 300
    dummy._status_line_message = lambda: None

    hud = SkyWindow._render_hud_state(dummy)

    assert hud.overlay_info_bottom_left is True
    assert dummy.state.overlay_info_bottom_left is True


def test_render_hud_state_uses_lower_third_to_switch_overlay_to_top_left() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = QPoint(10, 280)
    dummy.state.overlay_info_bottom_left = True
    dummy.height = lambda: 300
    dummy._status_line_message = lambda: None

    hud = SkyWindow._render_hud_state(dummy)

    assert hud.overlay_info_bottom_left is False
    assert dummy.state.overlay_info_bottom_left is False


def test_status_line_text_always_uses_night_style(monkeypatch) -> None:
    class DummyFontMetrics:
        def lineSpacing(self) -> int:  # noqa: N802 - Qt naming
            return 12

    class DummyPainter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setFont(self, _font) -> None:
            pass

        def fontMetrics(self):
            return DummyFontMetrics()

    calls: list[tuple[str, bool]] = []

    def fake_get_text_style(theme, *, status_line: bool = False):
        calls.append((theme, status_line))
        return (object(), object())

    monkeypatch.setattr(render_text_module, "get_text_style", fake_get_text_style)
    monkeypatch.setattr(
        render_text_module, "draw_outlined_text", lambda *_args, **_kwargs: None
    )

    render_text_module._draw_status_line_text(
        painter=DummyPainter(),
        message="loading",
        status_line_font=QFont(),
        viewport_rect=SimpleNamespace(bottom=lambda: 100),
        theme=THEME_STYLES_BY_PRESET["day"],
    )

    assert calls == [(THEME_STYLES_BY_PRESET["day"], True)]


def test_on_sky_data_calculated_updates_render_snapshot_once() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 0
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.request_sky_data_update = lambda *_args, **_kwargs: None
    dummy._safe_request_cloud_repaint = lambda: None

    update_calls: list[str] = []
    dummy.update = lambda: update_calls.append("update")

    celestial = object()
    sky_disc = object()
    payload = {
        "celestial": celestial,
        "sky_disc": sky_disc,
        "view_center": (20.0, 30.0),
        "geometry": render_geometry.get_screen_geometry(
            640, 480, dummy.viewer_data.view_alt_deg
        ),
        "render_generation": 0,
    }
    SkyWindow._on_sky_data_calculated(dummy, payload)

    assert dummy.state.render_view_center == (20.0, 30.0)
    assert dummy.state.celestial_data is celestial
    assert dummy.state.sky_disc_image is sky_disc
    assert dummy._compositor.invalidated is True
    assert update_calls == ["update"]


def test_on_sky_data_calculated_preserves_render_center_during_viewport_interaction() -> (
    None
):
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(40.0, 150.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(40.0, 150.0),
        viewport_interaction_mode=True,
    )
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._startup_initial_load_started = False
    dummy._startup_initial_data_loaded = True
    dummy._disc_generation = 0
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.request_sky_data_update = lambda *_args, **_kwargs: None
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.update = lambda: None

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (40.0, 150.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.render_view_center == (40.0, 150.0)


def test_on_sky_data_calculated_triggers_release_followup_updates() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(40.0, 150.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(40.0, 150.0),
        viewport_interaction_release_pending=True,
    )
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 0
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.request_sky_data_update = lambda *_args, **_kwargs: None
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.request_client_update = lambda: None
    dummy.reproject_tropical_cyclone_overlay = Mock()

    class _MenuButton:
        def __init__(self) -> None:
            self.visible = False

        def setVisible(self, value: bool) -> None:
            self.visible = value

    dummy.menu_button = _MenuButton()
    cloud_calls: list[str] = []
    dummy.start_background_cloud_update = lambda **kwargs: cloud_calls.append(
        str(kwargs.get("reason"))
    )
    dummy.start_background_terrain_horizon_update = lambda **kwargs: cloud_calls.append(
        str(kwargs.get("reason"))
    )

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (40.0, 150.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.viewport_interaction_release_pending is False
    assert dummy.state.viewport_interaction_mode is False
    assert dummy.menu_button.visible is True
    assert cloud_calls == ["view-change-release", "view-change-release"]
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with()


def test_on_sky_data_calculated_uses_idle_completion_reason() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(40.0, 150.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(40.0, 150.0),
        viewport_interaction_mode=True,
        viewport_interaction_release_pending=True,
        viewport_interaction_completion_reason="view-change-idle",
    )
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 0
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.request_sky_data_update = lambda *_args, **_kwargs: None
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.request_client_update = lambda: None
    dummy.reproject_tropical_cyclone_overlay = Mock()

    class _MenuButton:
        def __init__(self) -> None:
            self.visible = False

        def setVisible(self, value: bool) -> None:
            self.visible = value

    dummy.menu_button = _MenuButton()
    cloud_calls: list[str] = []
    dummy.start_background_cloud_update = lambda **kwargs: cloud_calls.append(
        str(kwargs.get("reason"))
    )
    dummy.start_background_terrain_horizon_update = lambda **kwargs: cloud_calls.append(
        str(kwargs.get("reason"))
    )

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (40.0, 150.0),
            "geometry": render_geometry.get_screen_geometry(
                640,
                480,
                dummy.viewer_data.view_alt_deg,
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.viewport_interaction_release_pending is False
    assert dummy.state.viewport_interaction_completion_reason is None
    assert dummy.state.viewport_interaction_mode is False
    assert dummy.menu_button.visible is True
    assert cloud_calls == ["view-change-idle", "view-change-idle"]
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with()


def test_on_sky_data_calculated_keeps_existing_cloud_refresh_deadline() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    cloud_deadline = datetime(2026, 5, 25, 0, 10, tzinfo=timezone.utc)
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        cloud_next_refresh_utc=cloud_deadline,
    )
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = object()
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 0
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.request_sky_data_update = lambda *_args, **_kwargs: None
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.update = lambda: None

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (20.0, 30.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.cloud_next_refresh_utc == cloud_deadline


def test_apply_startup_delta_t_disables_tropical_cyclone_layer_when_time_shifted() -> (
    None
):
    dummy = _WindowStub(
        tropical_cyclone_opacity=0.4,
        show_tropical_cyclone_overlay=True,
        _tropical_cyclone_controller=object(),
        _tropical_cyclone_requested_enabled=True,
        _tropical_cyclone_opacity_when_enabled=0.4,
        _action_toggle_satellites=Mock(),
        _action_toggle_aircraft=Mock(),
        _action_toggle_tropical_cyclone=Mock(),
    )

    window_module.SkyWindow.apply_startup_delta_t(dummy, timedelta(hours=-10))

    assert dummy.tropical_cyclone_opacity == 0.0
    assert dummy.show_tropical_cyclone_overlay is True
    dummy._action_toggle_tropical_cyclone.setEnabled.assert_called_once_with(False)
    dummy._action_toggle_tropical_cyclone.setChecked.assert_not_called()


def test_on_sky_data_calculated_discards_stale_view_center_after_jump() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(14.25, 87.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(14.25, 87.0),
        viewport_interaction_mode=False,
    )
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 0
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    retry_calls: list[str] = []
    dummy.request_sky_data_update = lambda *_args, **kwargs: retry_calls.append(
        str(kwargs.get("reason"))
    )
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.update = lambda: None

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (0.0, 180.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.render_view_center == (14.25, 87.0)
    assert retry_calls == ["stale-view-center"]


def test_schedule_satellite_retry_after_failure_uses_two_hour_backoff() -> None:
    dummy = _WindowStub()
    dummy.satellite_opacity = 0.5
    dummy._is_shutting_down = False
    dummy.state = SimpleNamespace(satellite_next_refresh_utc=None)
    dummy._satellite_layer_enabled = lambda: True

    SkyWindow._schedule_satellite_retry_after_failure(dummy)

    assert dummy.state.satellite_next_refresh_utc is not None
    assert dummy.state.satellite_next_refresh_utc > datetime.now(timezone.utc)


def test_satellite_validity_remaining_ms_uses_refresh_time() -> None:
    dummy = _WindowStub()
    dummy.satellite_state = SimpleNamespace(
        refreshed_at_utc=datetime.now(timezone.utc),
        element_epoch_utc=datetime(2020, 1, 1, tzinfo=timezone.utc),
    )

    remaining_ms = SkyWindow._satellite_validity_remaining_ms(dummy)

    assert remaining_ms is not None
    assert remaining_ms > 0


def test_on_satellite_failed_schedules_failure_backoff() -> None:
    dummy = _WindowStub()
    dummy.satellite_opacity = 0.5
    dummy.satellite_state = SimpleNamespace(set_error_banner=Mock())
    retry_calls: list[str] = []
    dummy._schedule_satellite_retry_after_failure = lambda: retry_calls.append("retry")
    dummy.update = Mock()

    SkyWindow._on_satellite_failed(dummy, {"banner": "Satellites: timed out"})

    dummy.satellite_state.set_error_banner.assert_called_once_with(
        "Satellites: timed out"
    )
    assert retry_calls == ["retry"]
    dummy.update.assert_called_once()


def test_jump_to_satellite_target_uses_cached_satellite_records_below_horizon(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        window_module, "find_satellite_altaz", lambda *_args, **_kwargs: (-12.0, 123.0)
    )
    monkeypatch.setattr(
        window_module.QTimer,
        "singleShot",
        lambda _ms, func: func(),
    )

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy.satellite_state = SimpleNamespace(
        records_by_group={"iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]},
        set_banner=Mock(),
    )
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()
    dummy._current_time_obj = lambda: astropy.time.Time("2026-03-22T12:00:00Z")
    dummy._find_satellite_jump_altaz = lambda key: SkyWindow._find_satellite_jump_altaz(
        dummy, key
    )

    SkyWindow._jump_to_search_target(
        dummy,
        SearchJumpTarget(
            label="ISS",
            ra_hours=0.0,
            dec_deg=0.0,
            kind="satellite",
            sort_key=(99.0, "iss"),
            subtitle="Satellite",
            object_key="ISS",
        ),
    )

    assert dummy.viewer_data.view_center == (-12.0, 123.0)
    assert dummy.state.jump_highlight_name == "ISS"
    assert dummy.state.jump_highlight_altaz == (-12.0, 123.0)
    dummy.request_sky_data_update.assert_called_once()


def test_find_satellite_jump_altaz_falls_back_to_disk_cache(monkeypatch) -> None:
    monkeypatch.setattr(
        window_module,
        "find_satellite_altaz",
        lambda records_by_group, **_kwargs: (
            (-40.0, 151.0) if records_by_group.get("iss") else None
        ),
    )
    monkeypatch.setattr(
        window_module,
        "load_satellite_cache",
        lambda group_key: (
            SimpleNamespace(records=[{"OBJECT_NAME": "ISS (ZARYA)"}])
            if group_key == "iss"
            else None
        ),
    )

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(40.7128, -74.0060),
        timezone_name="America/New_York",
        city_name="New York City",
        view_center=(20.0, 30.0),
        observer_height_m=10.0,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={})
    dummy._enabled_satellite_groups = ("iss",)
    dummy._current_time_obj = lambda: astropy.time.Time("2026-03-23T12:13:24Z")
    dummy._load_cached_satellite_records = lambda groups: (
        SkyWindow._load_cached_satellite_records(dummy, groups)
    )

    altaz = SkyWindow._find_satellite_jump_altaz(dummy, "ISS")

    assert altaz == (-40.0, 151.0)


def test_jump_to_satellite_target_sets_banner_when_not_available() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy.update = Mock()
    dummy._load_cached_satellite_records = lambda _groups: {}
    dummy._find_satellite_jump_altaz = lambda key: SkyWindow._find_satellite_jump_altaz(
        dummy, key
    )

    SkyWindow._jump_to_search_target(
        dummy,
        SearchJumpTarget(
            label="ISS",
            ra_hours=0.0,
            dec_deg=0.0,
            kind="satellite",
            sort_key=(99.0, "iss"),
            subtitle="Satellite",
            object_key="ISS",
        ),
    )

    dummy.satellite_state.set_banner.assert_called_once_with(
        "Satellites: ISS not available"
    )
    dummy.update.assert_called_once()


def test_jump_to_place_target_uses_projected_altaz(monkeypatch) -> None:
    monkeypatch.setattr(
        window_module,
        "_project_place_targets_to_altaz",
        lambda **kwargs: (
            PlaceTargetProjection(
                alt_deg=-3.5,
                az_deg=145.0,
                distance_km=12.0,
                target_latitude_deg=float(kwargs["target_latitude_deg"][0]),
                target_longitude_deg=float(kwargs["target_longitude_deg"][0]),
                target_height_m=float(kwargs["target_height_m"][0]),
            ),
        ),
    )
    monkeypatch.setattr(
        window_module.QTimer,
        "singleShot",
        lambda _ms, func: func(),
    )

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    SkyWindow._jump_to_search_target(
        dummy,
        SearchJumpTarget(
            label="Tokyo Station",
            kind="place",
            sort_key=(0.0, "tokyo station"),
            subtitle="Place / railway / station",
            latitude_deg=35.681236,
            longitude_deg=139.767125,
        ),
    )

    # Allow the view center to go below the horizon (>= -45.0°) per policy.
    assert dummy.viewer_data.view_center == (-3.5, 145.0)
    assert dummy.state.jump_highlight_name == "Tokyo Station"
    assert dummy.state.jump_highlight_altaz == (-3.5, 145.0)
    dummy.request_sky_data_update.assert_called_once()


def test_jump_to_jpl_small_body_target_can_set_persistent_overlay(caplog) -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda _target, **_kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda _target, **_kwargs: (12.5, 220.0),
        )
        with caplog.at_level("INFO"):
            SkyWindow._jump_to_search_target(
                dummy,
                SearchJumpTarget(
                    label="Ceres",
                    kind="jpl_small_body",
                    sort_key=(0.0, "ceres"),
                    subtitle="Asteroid / 1 Ceres",
                    object_key="20000001",
                    command="DES=20000001;",
                    persistent_keep_marker=True,
                ),
            )

    assert dummy.viewer_data.view_center == (12.5, 220.0)
    assert dummy.state.jump_highlight_name == "Ceres"
    assert dummy.state.jump_highlight_altaz == (12.5, 220.0)
    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.label == "Ceres"
    assert dummy.state.persistent_search_target.persistent_keep_marker is True
    assert dummy.state.persistent_search_next_refresh_utc == datetime(
        2026, 4, 18, 13, 0, tzinfo=timezone.utc
    )
    assert (
        "JPL persistent target set: label=Ceres kind=jpl_small_body group=<none>"
        in caplog.text
    )
    assert "target_time_utc=2026-04-18T12:00:00+00:00" in caplog.text
    assert "alt=12.5 az=220.0 command=DES=20000001;" in caplog.text


def test_jump_to_jpl_small_body_target_uses_state_vector_when_present(
    monkeypatch,
) -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()
    dummy._current_time_obj = lambda: astropy.time.Time("2026-04-18T12:00:00Z")

    monkeypatch.setattr(
        window_module,
        "resolve_jpl_target_state_vector",
        lambda _target, **_kwargs: (
            datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
            (1.0, 2.0, 3.0),
            (0.1, 0.2, 0.3),
        ),
    )
    monkeypatch.setattr(
        window_module,
        "project_jpl_target_altaz_from_state_vector",
        lambda _target, **_kwargs: (11.5, 221.0),
    )

    SkyWindow._jump_to_search_target(
        dummy,
        SearchJumpTarget(
            label="Ceres",
            kind="jpl_small_body",
            sort_key=(0.0, "ceres"),
            subtitle="Asteroid / 1 Ceres",
            object_key="20000001",
            command="DES=20000001;",
            persistent_keep_marker=True,
        ),
    )

    assert dummy.viewer_data.view_center == (11.5, 221.0)
    assert dummy.state.jump_highlight_altaz == (11.5, 221.0)
    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.horizons_position_km == (1.0, 2.0, 3.0)
    assert dummy.state.persistent_search_target.horizons_velocity_km_s == (
        0.1,
        0.2,
        0.3,
    )


def test_jump_to_jpl_small_body_target_honors_fixed_search_axes() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(5.0, 210.0),
        observer_height_m=1.7,
    )
    dummy._search_view_center_base = (5.0, 210.0)
    dummy._search_view_center_alt_specified = True
    dummy._search_view_center_az_specified = False
    dummy.state = SkyWindowState(
        render_view_center=(5.0, 210.0),
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda _target, **_kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda _target, **_kwargs: (12.5, 220.0),
        )
        SkyWindow._jump_to_search_target(
            dummy,
            SearchJumpTarget(
                label="Ceres",
                kind="jpl_small_body",
                sort_key=(0.0, "ceres"),
                subtitle="Asteroid / 1 Ceres",
                object_key="20000001",
                command="DES=20000001;",
            ),
        )

    assert dummy.viewer_data.view_center == (5.0, 220.0)


def test_jump_to_jpl_small_body_target_can_disable_fixed_search_axes() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(5.0, 210.0),
        observer_height_m=1.7,
    )
    dummy._search_view_center_base = (5.0, 210.0)
    dummy._search_view_center_alt_specified = True
    dummy._search_view_center_az_specified = True
    dummy.state = SkyWindowState(
        render_view_center=(5.0, 210.0),
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda _target, **_kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda _target, **_kwargs: (12.5, 220.0),
        )
        SkyWindow._jump_to_search_target(
            dummy,
            SearchJumpTarget(
                label="Ceres",
                kind="jpl_small_body",
                sort_key=(0.0, "ceres"),
                subtitle="Asteroid / 1 Ceres",
                object_key="20000001",
                command="DES=20000001;",
                preserve_cli_view_center=False,
            ),
        )

    assert dummy.viewer_data.view_center == (12.5, 220.0)


def test_jump_to_jpl_small_body_target_without_keep_flags_clears_overlay() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=SearchJumpTarget(
            label="Old",
            kind="jpl_small_body",
            sort_key=(0.0, "old"),
            alt_deg=10.0,
            az_deg=30.0,
            persistent_keep_marker=True,
        ),
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda _target, **_kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda _target, **_kwargs: (12.5, 220.0),
        )
        SkyWindow._jump_to_search_target(
            dummy,
            SearchJumpTarget(
                label="Ceres",
                kind="jpl_small_body",
                sort_key=(0.0, "ceres"),
                subtitle="Asteroid / 1 Ceres",
                object_key="20000001",
                command="DES=20000001;",
            ),
        )

    assert dummy.state.persistent_search_target is None


def test_jpl_small_body_failure_reschedules_one_hour_later() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    current_target = SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "ceres"),
        subtitle="Asteroid / 1 Ceres",
        object_key="20000001",
        command="DES=20000001;",
        alt_deg=12.5,
        az_deg=220.0,
        target_time_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
        persistent_keep_marker=True,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=current_target,
        persistent_search_next_refresh_utc=datetime(
            2026, 4, 18, 13, 0, tzinfo=timezone.utc
        ),
        persistent_search_reference_time_utc=datetime(
            2026, 4, 18, 12, 0, tzinfo=timezone.utc
        ),
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy.request_client_update = Mock()
    dummy._schedule_persistent_search_refresh = Mock()

    SkyWindow._on_jpl_failed(
        dummy,
        {
            "target": current_target,
            "target_time_utc": datetime(2026, 4, 18, 13, 0, tzinfo=timezone.utc),
            "refreshed_at_utc": datetime(2026, 4, 18, 13, 2, tzinfo=timezone.utc),
            "banner": "JPL: timed out",
            "error": "timed out",
            "reason": "timer",
        },
    )

    assert dummy.state.persistent_search_last_error == "timed out"
    assert dummy.state.persistent_search_last_refresh_utc == datetime(
        2026, 4, 18, 13, 2, tzinfo=timezone.utc
    )
    assert dummy.state.persistent_search_next_refresh_utc == datetime(
        2026, 4, 18, 14, 2, tzinfo=timezone.utc
    )
    dummy._schedule_persistent_search_refresh.assert_called_once()
    dummy.request_client_update.assert_called_once()


def test_on_jpl_ready_logs_refreshed_persistent_target(caplog) -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._current_time_obj = lambda: astropy.time.Time("2026-04-18T13:00:00Z")
    current_target = SearchJumpTarget(
        label="Voyager 1",
        kind="jpl_small_body",
        sort_key=(0.0, "voyager 1"),
        subtitle="spacecraft",
        object_key="-31",
        command="DES=-31;",
        alt_deg=48.6,
        az_deg=245.6,
        horizons_epoch_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
        horizons_position_km=(1.0, 2.0, 3.0),
        horizons_velocity_km_s=(0.1, 0.2, 0.3),
        target_time_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
        jpl_group="sb",
        persistent_keep_marker=True,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=current_target,
        persistent_search_next_refresh_utc=datetime(
            2026, 4, 18, 13, 0, tzinfo=timezone.utc
        ),
        persistent_search_reference_time_utc=datetime(
            2026, 4, 18, 12, 0, tzinfo=timezone.utc
        ),
    )
    dummy.request_client_update = Mock()
    dummy._schedule_persistent_search_refresh = Mock()
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda _target, **_kwargs: (49.1, 244.7),
        )

        with caplog.at_level("INFO"):
            SkyWindow._on_jpl_ready(
                dummy,
                {
                    "target": current_target,
                    "target_time_utc": datetime(
                        2026, 4, 18, 13, 0, tzinfo=timezone.utc
                    ),
                    "refreshed_at_utc": datetime(
                        2026, 4, 18, 13, 2, tzinfo=timezone.utc
                    ),
                    "horizons_epoch_utc": datetime(
                        2026, 4, 18, 13, 0, tzinfo=timezone.utc
                    ),
                    "horizons_position_km": (4.0, 5.0, 6.0),
                    "horizons_velocity_km_s": (0.4, 0.5, 0.6),
                    "rows": [["2026-Apr-18 13:00:00", "*", "m", "244.7", "49.1"]],
                    "reason": "timer",
                },
            )

    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.alt_deg == 49.1
    assert dummy.state.persistent_search_target.az_deg == 244.7
    assert dummy.state.persistent_search_target.horizons_position_km == (4.0, 5.0, 6.0)
    assert dummy.state.persistent_search_target.horizons_velocity_km_s == (
        0.4,
        0.5,
        0.6,
    )
    assert (
        "JPL persistent target refreshed: label=Voyager 1 kind=jpl_small_body group=sb"
        in caplog.text
    )
    assert "target_time_utc=2026-04-18T13:00:00+00:00" in caplog.text
    assert "alt=49.1 az=244.7 command=DES=-31;" in caplog.text


def test_refresh_projected_persistent_search_target_reprojects_state_vector(
    monkeypatch,
) -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._current_time_obj = lambda: astropy.time.Time("2026-04-18T13:00:00Z")
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=SearchJumpTarget(
            label="Voyager 1",
            kind="jpl_small_body",
            sort_key=(0.0, "voyager 1"),
            command="DES=-31;",
            alt_deg=48.6,
            az_deg=245.6,
            horizons_epoch_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
            horizons_position_km=(1.0, 2.0, 3.0),
            horizons_velocity_km_s=(0.1, 0.2, 0.3),
            persistent_keep_marker=True,
        ),
    )
    dummy.request_client_update = Mock()
    monkeypatch.setattr(
        window_updates_module,
        "project_jpl_target_altaz_from_state_vector",
        lambda _target, **_kwargs: (12.5, 220.0),
    )

    window_updates_module.SkyWindowUpdatesMixin.refresh_projected_persistent_search_target(
        dummy
    )

    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.alt_deg == 12.5
    assert dummy.state.persistent_search_target.az_deg == 220.0
    assert dummy.request_client_update.called


def test_handle_client_resize_discards_cached_sky_disc_and_requests_refresh() -> None:
    dummy = _WindowStub()
    sky_disc_image = QImage(4, 4, QImage.Format.Format_ARGB32_Premultiplied)
    dummy._frameless_frame = None
    dummy.menu_button = None
    dummy.size_grip = None
    dummy._compositor = _DummyCompositor()
    dummy._begin_viewport_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.request_client_update = Mock()
    dummy.start_background_cloud_update = Mock()
    dummy._raise_overlay_widgets = Mock()
    dummy._discard_stale_disc_images = lambda: (
        window_module.SkyWindowCoreMixin._discard_stale_disc_images(  # type: ignore[attr-defined]
            dummy
        )
    )
    dummy.cloud_controller = None
    dummy._cloud_controller = None
    dummy.cloud_state = SimpleNamespace(
        image=np.zeros((2, 2, 4), dtype=np.uint8),
        missing_mask=np.ones((2, 2), dtype=np.uint8),
        cloud_amount_field=SimpleNamespace(),
        render_key=object(),
        request_id=1,
        missing_mask_key=2,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        sky_disc_image=sky_disc_image,
    )
    dummy.width = lambda: 200
    dummy.height = lambda: 100
    dummy.client_width = lambda: 200
    dummy.client_height = lambda: 100

    event = QResizeEvent(QSize(220, 120), QSize(200, 100))

    SkyWindow._handle_client_resize(dummy, event)

    assert dummy._disc_generation == 1
    assert dummy.state.sky_disc_image is None
    assert dummy.cloud_state.image is None
    assert dummy.cloud_state.missing_mask is None
    assert dummy.cloud_state.cloud_amount_field is None
    assert dummy.cloud_state.render_key is None
    assert dummy.cloud_state.request_id is None
    assert dummy.cloud_state.missing_mask_key is None
    assert dummy._compositor.invalidated is True
    dummy._begin_viewport_interaction_mode.assert_called_once()
    dummy.request_sky_data_update.assert_called_once_with(
        reason="resize",
        allow_during_viewport_interaction=True,
    )
    dummy.request_client_update.assert_called_once()


def test_search_satellite_targets_resolves_known_artificial_satellites() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    targets = SkyWindow._search_satellite_targets(dummy, "ISS")

    assert len(targets) == 1
    assert targets[0].label == "ISS"
    assert targets[0].kind == "satellite"
    assert targets[0].alt_deg is None
    assert targets[0].az_deg is None


def test_search_jpl_targets_skips_solar_system_bodies(monkeypatch) -> None:
    lookup_calls: list[str] = []

    def fake_lookup(_query: str, *, group: str):
        lookup_calls.append(group)
        return {"count": 0, "result": []}

    monkeypatch.setattr(window_module, "fetch_horizons_lookup", fake_lookup)

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    targets = SkyWindow._search_jpl_targets(dummy, "Mars")

    assert targets == []
    assert lookup_calls == []


def test_search_jpl_targets_limits_candidates_to_500(monkeypatch) -> None:
    lookup_calls: list[str] = []

    def fake_lookup(_query: str, *, group: str):
        lookup_calls.append(group)
        if group == "mb":
            return {
                "count": 600,
                "result": [
                    {
                        "name": f"PANSTARRS-{idx}",
                        "pdes": str(idx),
                        "spkid": str(1000 + idx),
                        "type": "asteroid",
                    }
                    for idx in range(600)
                ],
            }
        return {"count": 0, "result": []}

    monkeypatch.setattr(window_module, "fetch_horizons_lookup", fake_lookup)

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    targets = SkyWindow._search_jpl_targets(dummy, "PANSTARRS")

    assert lookup_calls == ["sct", "mb", "sb"]
    assert len(targets) == 500
    assert targets[0].label == "PANSTARRS-0"
    assert targets[-1].label == "PANSTARRS-499"


def test_search_jpl_targets_skips_solar_system_bodies_directly() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )

    assert SkyWindow._search_jpl_targets(dummy, "Sun") == []
    assert SkyWindow._search_jpl_targets(dummy, "Moon") == []
    assert SkyWindow._search_jpl_targets(dummy, "Mars") == []


def test_handle_client_key_press_rotates_view_immediately() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy._sync_view_altitude_actions = Mock()
    dummy.request_client_update = Mock()
    dummy._viewport_interaction_idle_timer = _DummyTimer(active=False)
    dummy._rotate_view = lambda *args, **kwargs: SkyWindow._rotate_view(
        dummy, *args, **kwargs
    )
    dummy._begin_viewport_interaction_mode = lambda *args, **kwargs: (
        SkyWindow._begin_viewport_interaction_mode(dummy, *args, **kwargs)
    )
    dummy._update_viewport_interaction_stars = lambda: (
        SkyWindow._update_viewport_interaction_stars(dummy)
    )
    dummy._viewport_rotation_keys = lambda: SkyWindow._viewport_rotation_keys(dummy)

    event = SimpleNamespace(
        key=lambda: Qt.Key.Key_Left,
        modifiers=lambda: Qt.KeyboardModifier.NoModifier,
        isAutoRepeat=lambda: False,
        accept=Mock(),
    )

    SkyWindow._handle_client_key_press(dummy, event)

    assert dummy.viewer_data.view_center == (20.0, 25.0)
    assert dummy.state.render_view_center == (20.0, 25.0)
    assert dummy.state.viewport_interaction_mode is True
    assert dummy._viewport_interaction_idle_timer.started_with == []
    assert dummy._viewport_rotation_keys_down == {Qt.Key.Key_Left}
    dummy._sync_view_altitude_actions.assert_called_once()
    dummy.request_client_update.assert_called_once()
    event.accept.assert_called_once()


def test_set_view_center_leaves_viewport_fast_mode_after_dialog_change() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
        viewport_interaction_stars=object(),
    )
    dummy._sync_view_altitude_actions = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.request_client_update = Mock()
    dummy.start_background_cloud_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy._end_viewport_interaction_mode = lambda *args, **kwargs: (
        SkyWindow._end_viewport_interaction_mode(dummy, *args, **kwargs)
    )
    dummy._finalize_view_direction_change = lambda: (
        SkyWindow._finalize_view_direction_change(dummy)
    )
    dummy.request_cloud_projection_update = Mock()

    SkyWindow._set_view_center(
        dummy,
        25.0,
        45.0,
        interactive_viewport=False,
        start_viewport_idle_timer=False,
    )

    assert dummy.viewer_data.view_center == (25.0, 45.0)
    assert dummy.state.render_view_center == (25.0, 45.0)
    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    dummy.request_sky_data_update.assert_called_once_with(
        reason="view-change-idle",
        allow_during_viewport_interaction=True,
    )
    dummy.request_cloud_projection_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.start_background_terrain_horizon_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.request_client_update.assert_called_once()
    dummy._sync_view_altitude_actions.assert_called_once()


def test_open_view_direction_dialog_shows_fast_frame_before_release(
    monkeypatch,
) -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=False,
    )
    set_calls: list[tuple[tuple[float, float], dict[str, object]]] = []
    end_calls: list[str] = []

    class _FakeDialog:
        def __init__(self, view_center: tuple[float, float], parent) -> None:
            self._view_center = view_center
            self._parent = parent

        def exec(self) -> int:
            return 1

        def selected_view_center(self) -> tuple[float, float]:
            return (25.0, 45.0)

    monkeypatch.setattr(window_module, "ViewDirectionDialog", _FakeDialog)
    monkeypatch.setattr(
        window_module.QTimer,
        "singleShot",
        lambda ms, func: end_calls.append(f"timer:{ms}") or func(),
    )
    dummy._set_view_center = lambda *args, **kwargs: set_calls.append(
        (args, dict(kwargs))
    )
    dummy._end_viewport_interaction_mode = lambda *_args, **kwargs: end_calls.append(
        str(kwargs.get("reason"))
    )
    dummy._finalize_view_direction_change = lambda: (
        SkyWindow._finalize_view_direction_change(dummy)
    )
    dummy._finalize_view_direction_dialog_change = lambda: (
        SkyWindow._finalize_view_direction_dialog_change(dummy)
    )

    SkyWindow._open_view_direction_dialog(dummy)

    assert set_calls == [
        (
            (25.0, 45.0),
            {
                "interactive_viewport": True,
                "start_viewport_idle_timer": False,
            },
        )
    ]
    assert end_calls == ["timer:0", "view-change-release"]


def test_handle_client_key_release_ends_viewport_interaction_mode() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy._viewport_rotation_keys_down = {Qt.Key.Key_Left}
    dummy._viewport_rotation_keys = lambda: SkyWindow._viewport_rotation_keys(dummy)
    dummy._end_viewport_interaction_mode = lambda *args, **kwargs: (
        SkyWindow._end_viewport_interaction_mode(dummy, *args, **kwargs)
    )
    dummy.request_sky_data_update = Mock()
    dummy.start_background_cloud_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.reproject_tropical_cyclone_overlay = Mock()
    dummy.request_client_update = Mock()

    event = SimpleNamespace(
        key=lambda: Qt.Key.Key_Left,
        isAutoRepeat=lambda: False,
        accept=Mock(),
    )

    SkyWindow._handle_client_key_release(dummy, event)

    assert dummy._viewport_rotation_keys_down == set()
    assert dummy.state.viewport_interaction_release_pending is True
    assert dummy.state.viewport_interaction_mode is True
    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-release",
        allow_during_viewport_interaction=True,
    )
    dummy.start_background_cloud_update.assert_not_called()
    dummy.start_background_terrain_horizon_update.assert_not_called()
    dummy.request_client_update.assert_not_called()
    event.accept.assert_called_once()


def test_handle_client_mouse_move_is_ignored_during_startup_block() -> None:
    dummy = _WindowStub()
    dummy._startup_input_blocked = lambda: True
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        mouse_pos=None,
    )
    dummy.request_client_update = Mock()

    event = SimpleNamespace(pos=lambda: QPoint(10, 20), accept=Mock())

    SkyWindow._handle_client_mouse_move(dummy, event)

    assert dummy.state.mouse_pos is None
    dummy.request_client_update.assert_not_called()
    event.accept.assert_called_once()


def test_render_hud_state_ignores_mouse_position_during_startup_block() -> None:
    dummy = _WindowStub()
    dummy._startup_input_blocked = lambda: True
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        mouse_pos=QPoint(10, 20),
        overlay_info_bottom_left=False,
    )
    dummy.observation_info_pinned = False
    dummy.client_height = lambda: 300
    dummy._status_line_message = lambda: "status"

    hud = window_render_module.SkyWindowRenderMixin._render_hud_state(dummy)

    assert hud.mouse_pos is None
    assert hud.overlay_info_bottom_left is False


def test_end_viewport_interaction_mode_marks_idle_reason() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy.request_sky_data_update = Mock()
    dummy.request_cloud_projection_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.reproject_tropical_cyclone_overlay = Mock()
    dummy.request_client_update = Mock()

    SkyWindow._end_viewport_interaction_mode(dummy)

    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    assert dummy.state.viewport_interaction_release_pending is False
    assert dummy.state.viewport_interaction_completion_reason is None
    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-idle",
        allow_during_viewport_interaction=True,
    )
    dummy.request_cloud_projection_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.start_background_terrain_horizon_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with()
    dummy.request_client_update.assert_called_once()


def test_end_viewport_interaction_mode_release_reprojects_tropical_cyclone() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy.request_sky_data_update = Mock()
    dummy.request_cloud_projection_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.reproject_tropical_cyclone_overlay = Mock()
    dummy.request_client_update = Mock()

    SkyWindow._end_viewport_interaction_mode(
        dummy,
        reason="viewport-interaction-release",
    )

    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-release",
        allow_during_viewport_interaction=True,
    )
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with(
        allow_during_viewport_interaction=True,
    )
    dummy.request_cloud_projection_update.assert_not_called()
    dummy.start_background_terrain_horizon_update.assert_not_called()
    dummy.request_client_update.assert_not_called()


def test_end_viewport_interaction_mode_release_clears_interaction_when_sky_update_is_busy() -> (
    None
):
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy.request_sky_data_update = Mock(return_value=False)
    dummy.request_cloud_projection_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.reproject_tropical_cyclone_overlay = Mock()
    dummy.request_client_update = Mock()

    SkyWindow._end_viewport_interaction_mode(
        dummy,
        reason="viewport-interaction-release",
    )

    assert dummy.state.viewport_interaction_mode is False
    assert dummy.state.viewport_interaction_stars is None
    assert dummy.state.viewport_interaction_release_pending is False
    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-release",
        allow_during_viewport_interaction=True,
    )
    dummy.request_cloud_projection_update.assert_not_called()
    dummy.start_background_terrain_horizon_update.assert_not_called()
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with()
    dummy.request_client_update.assert_called_once()


def test_tropical_cyclone_layer_is_disabled_when_opacity_is_zero() -> None:
    dummy = _WindowStub(
        show_tropical_cyclone_overlay=True,
        tropical_cyclone_opacity=0.0,
        _tropical_cyclone_controller=object(),
    )

    assert (
        window_updates_module.SkyWindowUpdatesMixin._tropical_cyclone_layer_enabled(
            dummy
        )
        is False
    )


def test_jump_to_jpl_major_body_target_keeps_overlay_without_refresh() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda _target, **_kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda _target, **_kwargs: (15.0, 123.0),
        )
        SkyWindow._jump_to_search_target(
            dummy,
            SearchJumpTarget(
                label="Mars",
                kind="jpl_body",
                sort_key=(0.0, "mars"),
                subtitle="major body",
                object_key="499",
                command="499",
                jpl_group="mb",
                persistent_keep_marker=True,
            ),
        )

    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.label == "Mars"
    assert dummy.state.persistent_search_next_refresh_utc == datetime(
        2026, 4, 18, 13, 0, tzinfo=timezone.utc
    )


def test_reproject_satellite_overlay_falls_back_to_disk_cache() -> None:
    dummy = _WindowStub()
    dummy.satellite_opacity = 1.0
    dummy.viewer_data = ViewerData(
        location=(40.7128, -74.0060),
        timezone_name="America/New_York",
        city_name="New York City",
        view_center=(0.0, 151.0),
        observer_height_m=10.0,
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={"iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]},
        element_epoch_utc=datetime.now(timezone.utc),
    )
    dummy.state = SkyWindowState(
        render_view_center=(0.0, 151.0),
        satellite_projection_next_refresh_utc=None,
    )
    dummy._enabled_satellite_groups = ("iss",)
    dummy._satellite_validity_remaining_ms = lambda: 1000
    dummy._current_time_obj = lambda: astropy.time.Time("2026-03-23T12:13:24Z")
    dummy.update = Mock()

    SkyWindow.reproject_satellite_overlay(dummy)

    assert dummy.state.satellite_projection_next_refresh_utc is not None
    dummy.update.assert_called_once()


def test_on_sky_data_calculated_discards_stale_generation_and_requests_refresh() -> (
    None
):
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 3
    dummy.width = lambda: 800
    dummy.height = lambda: 500
    requests: list[str] = []
    dummy.request_sky_data_update = lambda *_args, **_kwargs: requests.append("refresh")
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.update = lambda: (_ for _ in ()).throw(
        AssertionError("should not repaint stale sky payload")
    )

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (15.0, 120.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 2,
        },
    )

    assert dummy.state.sky_disc_image is None
    assert dummy._compositor.invalidated is False
    assert requests == ["refresh"]


def test_draw_viewport_interaction_layers_limits_stars_to_bright_subset(
    monkeypatch,
) -> None:
    calls: list[tuple[str, object]] = []

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: calls.append(("reference", None)),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: calls.append(("terrain", None)),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: calls.append(("urban", None)),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: calls.append(("direction", None)),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: calls.append(("zenith", None)),
    )

    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **kwargs: calls.append(("stars", kwargs.get("draw_vmag_limit"))),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **kwargs: calls.append(("planets", kwargs.get("draw_labels"))),
    )
    pipeline_module._draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(celestial_data=object()),
        style=_make_style(),
        hud=_make_hud(),
    )

    assert calls == [
        ("reference", None),
        ("stars", 4.0),
        ("planets", False),
        ("terrain", None),
    ]


def test_render_cached_frame_image_reuses_existing_image() -> None:
    render_calls: list[str] = []

    dummy = _WindowStub()
    dummy._present_frame_cache_key = None
    dummy._present_frame_cache_image = None
    dummy.size = lambda: window_render_module.QImage(
        24, 16, QImage.Format.Format_ARGB32_Premultiplied
    ).size()

    def render_fn(frame_painter) -> None:
        render_calls.append("render")
        frame_painter.fillRect(0, 0, 24, 16, window_render_module.Qt.GlobalColor.black)

    image_a = window_render_module.SkyWindowRenderMixin._render_cached_frame_image(
        dummy,
        frame_key=("same",),
        render_fn=render_fn,
        cache_key_attr="_present_frame_cache_key",
        cache_image_attr="_present_frame_cache_image",
    )
    image_b = window_render_module.SkyWindowRenderMixin._render_cached_frame_image(
        dummy,
        frame_key=("same",),
        render_fn=render_fn,
        cache_key_attr="_present_frame_cache_key",
        cache_image_attr="_present_frame_cache_image",
    )

    assert render_calls == ["render"]
    assert image_a.cacheKey() == image_b.cacheKey()


def test_render_fast_frame_image_downsamples_base_scene(monkeypatch) -> None:
    base_frame_sizes: list[tuple[int, int]] = []
    status_rect_sizes: list[tuple[int, int]] = []
    call_order: list[str] = []

    def _capture_base_scene(*_args, **kwargs) -> None:
        call_order.append("base")
        frame = kwargs["frame"]
        base_frame_sizes.append(
            (
                int(frame.viewport_rect.width()),
                int(frame.viewport_rect.height()),
            )
        )

    def _capture_status_line(*_args, **kwargs) -> None:
        call_order.append("status")
        viewport_rect = kwargs["viewport_rect"]
        status_rect_sizes.append(
            (
                int(viewport_rect.width()),
                int(viewport_rect.height()),
            )
        )

    def _capture_fast_overlays(*_args, **_kwargs) -> None:
        call_order.append("fast-overlays")

    monkeypatch.setattr(
        window_render_module,
        "render_base_scene_into_painter",
        _capture_base_scene,
    )
    monkeypatch.setattr(
        window_render_module,
        "render_fast_overlay_layers_into_painter",
        _capture_fast_overlays,
    )
    monkeypatch.setattr(
        window_render_module,
        "_draw_status_line",
        _capture_status_line,
    )
    monkeypatch.setattr(
        window_render_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: call_order.append("labels"),
    )

    dummy = _WindowStub(
        state=SkyWindowState(
            render_view_center=(45.0, 180.0),
            viewport_interaction_mode=True,
        ),
    )
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy._compositor = _DummyCompositor()
    dummy._fast_frame_base_cache_key = None
    dummy._fast_frame_base_cache_image = None
    dummy._fast_frame_cache_key = None
    dummy._fast_frame_cache_image = None
    dummy.client_rect = lambda: QRect(0, 0, 1600, 900)
    dummy.client_size = lambda: QSize(1600, 900)

    scene = _make_scene(viewer=dummy.viewer_data)
    style = _make_style(show_custom_window_frame=True)
    hud = _make_hud(viewport_interaction_mode=True, status_message="fast")
    frame = window_render_module.FrameContext(
        viewer=scene.viewer,
        time_obj=scene.time_obj,
        geometry=render_geometry.get_screen_geometry(
            1600,
            900,
            scene.viewer.view_alt_deg,
        ),
        viewport_rect=QRect(0, 0, 1600, 900),
    )

    image = window_render_module.SkyWindowRenderMixin._render_fast_frame_image(
        dummy,
        base_frame_key=("frame",),
        frame=frame,
        scene=scene,
        style=style,
        hud=hud,
    )

    assert base_frame_sizes == [(600, 338)]
    assert status_rect_sizes == [(1600, 900)]
    assert call_order == ["base", "fast-overlays", "labels", "status"]
    assert image.size() == QSize(1600, 900)


def test_draw_background_layer_can_skip_menu_button(monkeypatch) -> None:
    border_menu_flags: list[bool] = []

    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_radial_background",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_window_border",
        lambda *_args, **kwargs: border_menu_flags.append(
            bool(kwargs.get("draw_menu_button", True))
        ),
    )

    pipeline_module._draw_background_layer(
        painter=object(),
        geometry=SimpleNamespace(),
        viewport_rect=QRect(0, 0, 1600, 900),
        scene=_make_scene(),
        style=_make_style(show_custom_window_frame=True),
        draw_menu_button=False,
    )

    assert border_menu_flags == [False]


def test_viewport_interaction_hides_menu_button() -> None:
    class _MenuButton:
        def __init__(self) -> None:
            self.visible = True

        def setVisible(self, value: bool) -> None:
            self.visible = value

    dummy = _WindowStub()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.menu_button = _MenuButton()

    SkyWindow._sync_viewport_interaction_chrome_visibility(dummy)
    assert dummy.menu_button.visible is True

    dummy.state.viewport_interaction_mode = True
    SkyWindow._sync_viewport_interaction_chrome_visibility(dummy)
    assert dummy.menu_button.visible is False


def test_paint_event_skips_rendering_while_startup_overlay_visible(
    monkeypatch,
) -> None:
    class _VisibleOverlay:
        def isVisible(self) -> bool:  # noqa: N802 - Qt naming
            return True

    class _FailPainter:
        def __init__(self, *_args, **_kwargs) -> None:
            raise AssertionError(
                "QPainter should not be constructed while startup overlay is visible"
            )

    dummy = _WindowStub(
        _startup_log_overlay=_VisibleOverlay(),
        state=SkyWindowState(render_view_center=(45.0, 180.0)),
    )

    monkeypatch.setattr(window_render_module, "QPainter", _FailPainter)

    window_render_module.SkyWindowRenderMixin.paintEvent(dummy, SimpleNamespace())


def test_render_frame_cache_key_ignores_hover_and_status_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.enlarge_moon = False
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.aircraft_opacity = 0.4
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._status_line_message = lambda: "initial"
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(), missing_mask=object(), cloud_amount_field=None
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.state.water_overlay_dots = [object()]
    dummy.state.mouse_pos = None
    dummy.state.jump_highlight_name = None
    dummy.state.jump_highlight_altaz = None
    dummy.state.jump_highlight_until_ms = 0.0

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    dummy.state.mouse_pos = SimpleNamespace(x=lambda: 10, y=lambda: 20)
    dummy.state.jump_highlight_name = "Vega"
    dummy.state.jump_highlight_altaz = (20.0, 30.0)
    dummy.state.jump_highlight_until_ms = 12345.0
    dummy._status_line_message = lambda: "changed"

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    assert key_a == key_b


def test_render_frame_cache_key_tracks_water_overlay_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.enlarge_moon = False
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.aircraft_opacity = 0.4
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(), missing_mask=object(), cloud_amount_field=None
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.state.water_overlay_dots = [object()]

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    dummy.state.water_overlay_dots = [object(), object()]

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    assert key_a != key_b


def test_render_frame_cache_key_ignores_projected_tropical_cyclone_state_for_base_cache() -> (
    None
):
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.enlarge_moon = False
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.aircraft_opacity = 0.4
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(), missing_mask=object(), cloud_amount_field=None
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.state.water_overlay_dots = [object()]
    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text=None,
    )
    dummy._current_time_obj = lambda: astropy.time.Time(
        "2026-04-18T12:00:00", scale="utc"
    )

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm",
    )
    dummy._current_time_obj = lambda: astropy.time.Time(
        "2026-04-18T12:00:03", scale="utc"
    )

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    assert key_a == key_b


def test_present_frame_cache_key_tracks_hover_and_status_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.enlarge_moon = False
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.satellite_opacity = 0.4
    dummy.aircraft_opacity = 0.3
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=SimpleNamespace(source_cache_key=123),
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]

    base_key = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )
    key_a = SkyWindow._present_frame_cache_key(
        dummy,
        base_frame_key=base_key,
        hud=_make_hud(status_message="initial"),
    )

    dummy.state.mouse_pos = QPoint(10, 20)
    dummy.state.jump_highlight_name = "Vega"
    dummy.state.jump_highlight_altaz = (20.0, 30.0)
    dummy.state.jump_highlight_until_ms = 12345.0

    key_b = SkyWindow._present_frame_cache_key(
        dummy,
        base_frame_key=base_key,
        hud=_make_hud(mouse_pos=QPoint(10, 20), status_message="changed"),
    )

    assert key_a != key_b


def test_render_frame_cache_key_ignores_fast_overlay_state_for_base_cache() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.enlarge_moon = False
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.satellite_opacity = 0.4
    dummy.aircraft_opacity = 0.3
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=SimpleNamespace(source_cache_key=123),
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm-a",
    )

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    dummy.satellite_opacity = 0.9
    dummy.aircraft_opacity = 0.8
    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm-b",
    )

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    assert key_a == key_b


def test_present_frame_cache_key_tracks_projected_tropical_cyclone_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.enlarge_moon = False
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.satellite_opacity = 0.4
    dummy.aircraft_opacity = 0.3
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=SimpleNamespace(source_cache_key=123),
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm-a",
    )
    dummy._current_time_obj = lambda: astropy.time.Time(
        "2026-04-18T12:00:00", scale="utc"
    )

    base_key = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )
    key_a = SkyWindow._present_frame_cache_key(
        dummy,
        base_frame_key=base_key,
        hud=_make_hud(status_message="initial"),
    )

    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm-b",
    )
    key_b = SkyWindow._present_frame_cache_key(
        dummy,
        base_frame_key=base_key,
        hud=_make_hud(status_message="initial"),
    )

    assert key_a != key_b


def test_resolve_hover_targets_keeps_star_and_satellite_candidates_independent(
    monkeypatch,
) -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    mouse_pos = QPoint(12, 34)

    star_hit = ({"name": "Vega"}, QPointF(10.0, 30.0))
    satellite_hit = (
        SimpleNamespace(satellite_name="ISS"),
        QPointF(14.0, 36.0),
    )

    monkeypatch.setattr(
        window_render_module.render_stars,
        "find_highlighted_object",
        lambda *_args, **_kwargs: star_hit,
    )
    monkeypatch.setattr(
        window_render_module.render_deep_sky_objects,
        "find_highlighted_dso",
        lambda *_args, **_kwargs: {"name": "DSO"},
    )
    monkeypatch.setattr(
        window_render_module.render_satellites,
        "find_highlighted_satellite",
        lambda *_args, **_kwargs: satellite_hit,
    )

    (
        highlighted_object,
        highlighted_dso,
        highlighted_satellite,
        highlighted_tropical_cyclone,
    ) = window_render_module._resolve_hover_targets(
        celestial_data=celestial_data,
        render_viewer=viewer,
        mouse_pos=mouse_pos,
        geometry=geometry,
        satellite_records_by_group={"iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]},
        show_dso=True,
    )

    assert highlighted_object == star_hit
    assert highlighted_dso == {"name": "DSO"}
    assert highlighted_satellite == satellite_hit
    assert highlighted_tropical_cyclone is None


def test_draw_viewport_interaction_layers_prefers_interaction_star_subset(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_fast",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_normal",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )

    full_stars = {"name": ["full"]}
    interaction_stars = {"name": ["bright-only"]}
    seen_stars: list[object] = []
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda _p, **kwargs: seen_stars.append(kwargs["scene"].celestial_data.stars),
    )
    pipeline_module._draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(
            celestial_data=CelestialData(
                time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
                planets=[],
                stars=full_stars,
                deep_sky_objects={},
                celestial_equator_points=[],
                ecliptic_points=[],
                horizon_points=[],
            ),
        ),
        style=_make_style(),
        hud=_make_hud(viewport_interaction_stars=interaction_stars),
    )

    assert seen_stars == [interaction_stars]


def test_draw_viewport_interaction_layers_skips_water_when_terrain_horizon_hidden(
    monkeypatch,
) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: calls.append("water"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_fast",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_normal",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: None,
    )

    scene = replace(
        _make_scene(terrain_horizon_profile=[(1.0, 10.0)]),
        water_overlay_dots=[object()],
    )
    pipeline_module._draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=scene,
        style=_make_style(terrain_horizon_opacity=0.0, water_overlay_opacity=0.5),
        hud=_make_hud(),
    )

    assert calls == ["terrain"]


def test_draw_viewport_interaction_layers_prefers_scene_water_overlay_points(
    monkeypatch,
) -> None:
    terrain_calls: list[str] = []
    seen_water_points: list[object] = []
    sentinel_water_points = [object(), object()]

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: terrain_calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda _p, _g, _viewer, water_points, *_args, **_kwargs: (
            seen_water_points.append(water_points)
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_fast",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_normal",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: None,
    )

    scene = replace(
        _make_scene(terrain_horizon_profile=[(1.0, 10.0)]),
        water_overlay_dots=sentinel_water_points,
    )
    pipeline_module._draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=scene,
        style=_make_style(terrain_horizon_opacity=0.25, water_overlay_opacity=0.5),
        hud=_make_hud(),
    )

    assert terrain_calls == ["terrain"]
    assert seen_water_points == [sentinel_water_points]


def test_render_base_scene_skips_water_when_terrain_horizon_hidden(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module, "_clear_background_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_background_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_sky_cloud_layers", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_guide_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_planet_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_satellite_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_aircraft_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_urban_outline_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_star_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_hover_overlay_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_status_line", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_fast",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_normal",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: calls.append("water"),
    )

    scene = replace(
        _make_scene(terrain_horizon_profile=[(1.0, 10.0)]),
        water_overlay_dots=[object()],
    )

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setCompositionMode(self, *_args, **_kwargs) -> None:
            pass

        def fillRect(self, *_args, **_kwargs) -> None:
            pass

    pipeline_module.render_base_scene_into_painter(
        painter=_Painter(),
        frame=_make_frame(
            scene,
            SimpleNamespace(radius=600),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(terrain_horizon_opacity=0.0, water_overlay_opacity=0.5),
        compositor=object(),
        hud=_make_hud(),
    )

    assert calls == []


def test_draw_urban_outline_layer_skips_when_hidden(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: calls.append("urban"),
    )

    pipeline_module._draw_urban_outline_layer(
        painter=object(),
        geometry=object(),
        scene=_make_scene(
            urban_outlines=[
                UrbanOutlinePolyline(points=[(-1.0, 10.0), (-2.0, 12.0)], height_m=50.0)
            ]
        ),
        style=_make_style(show_urban_outline_layer=False),
    )

    assert calls == []


def test_draw_viewport_interaction_layers_draws_terrain_profile(monkeypatch) -> None:
    seen_main_profiles: list[object] = []
    seen_view_centers: list[object] = []
    seen_line_width_scales: list[float] = []
    fast_calls: list[object] = []
    expected_line_width_scale = pipeline_module.compute_star_render_upscale_factor(
        1200, 600
    )

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda _p, _g, viewer, main_profile, _main_distances, **kwargs: (
            seen_main_profiles.append(main_profile),
            seen_view_centers.append(viewer.view_center),
            seen_line_width_scales.append(float(kwargs["spec"].line_width_scale)),
            fast_calls.append(True),
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    terrain_profile = [(1.0, 10.0), (2.0, 20.0)]
    monkeypatch.setattr(
        pipeline_module, "_draw_star_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_planet_layer", lambda *_args, **_kwargs: None
    )
    pipeline_module._draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=terrain_profile,
        ),
        style=_make_style(),
        hud=_make_hud(),
    )

    assert seen_main_profiles == [terrain_profile]
    assert seen_view_centers == [(50.0, 210.0)]
    assert seen_line_width_scales == [expected_line_width_scale]
    assert fast_calls == [True]


def test_draw_viewport_interaction_layers_skips_urban_outlines(monkeypatch) -> None:
    seen_profiles: list[object] = []
    seen_view_centers: list[object] = []

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda _p, _g, profile, view_center, **_kwargs: (
            seen_profiles.append(profile),
            seen_view_centers.append(view_center),
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    urban_outlines = [[(-1.0, 10.0), (-2.0, 20.0)]]
    monkeypatch.setattr(
        pipeline_module, "_draw_star_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_planet_layer", lambda *_args, **_kwargs: None
    )
    pipeline_module._draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            urban_outlines=urban_outlines,
        ),
        style=_make_style(),
        hud=_make_hud(),
    )

    assert seen_profiles == []
    assert seen_view_centers == []


def test_draw_terrain_layers_scales_asterisms_but_keeps_urban_outline_widths_fixed(
    monkeypatch,
) -> None:
    calls: dict[str, list[float]] = {
        "asterisms": [],
        "dso": [],
        "terrain": [],
        "terrain_secondary": [],
        "reference": [],
        "direction": [],
        "zenith": [],
        "urban": [],
    }
    expected_line_width_scale = pipeline_module.compute_star_render_upscale_factor(
        1200, 600
    )

    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_deep_sky_shapes",
        lambda *_args, **kwargs: calls["dso"].append(
            float(kwargs.get("opacity_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_dso_hover_info",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **kwargs: calls["asterisms"].append(
            float(
                kwargs.get("base_line_width_scale", kwargs.get("line_width_scale", 1.0))
            )
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **kwargs: calls["terrain"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **kwargs: calls["terrain_secondary"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **kwargs: calls["reference"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **kwargs: calls["direction"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **kwargs: calls["zenith"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )

    monkeypatch.setattr(
        pipeline_module,
        "_draw_urban_outline_layer",
        lambda *_args, **kwargs: calls["urban"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    pipeline_module._draw_terrain_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=[(1.0, 10.0), (2.0, 20.0)],
            terrain_secondary_ridges_altaz_layers=[[(1.0, 10.0), (2.0, 20.0)]],
            terrain_secondary_ridges_distances_m_layers=[[10_000.0, 12_000.0]],
        ),
        style=_make_style(
            show_dso=True, show_asterisms=True, asterism_visibility_boost=2.0
        ),
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls["dso"] == [1.0]
    assert calls["asterisms"] == [expected_line_width_scale * 2.0]
    assert calls["terrain"] == []
    assert calls["terrain_secondary"] == [expected_line_width_scale]
    assert calls["reference"] == [1.0]
    assert calls["direction"] == []
    assert calls["zenith"] == []
    assert calls["urban"] == [1.0]


def test_draw_terrain_layers_dims_dso_and_asterisms_in_simplified_view(
    monkeypatch,
) -> None:
    calls: dict[str, list[float]] = {
        "asterisms": [],
        "dso": [],
    }

    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_deep_sky_shapes",
        lambda *_args, **kwargs: calls["dso"].append(
            float(kwargs.get("opacity_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **kwargs: calls["asterisms"].append(
            float(kwargs.get("base_line_alpha_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_urban_outline_layer",
        lambda *_args, **_kwargs: None,
    )

    pipeline_module._draw_terrain_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=[(1.0, 10.0), (2.0, 20.0)],
            terrain_secondary_ridges_altaz_layers=[[(1.0, 10.0), (2.0, 20.0)]],
            terrain_secondary_ridges_distances_m_layers=[[10_000.0, 12_000.0]],
        ),
        style=_make_style(
            show_dso=True, show_asterisms=True, asterism_visibility_boost=2.0
        ),
        simplified_view_active=True,
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls["dso"] == [0.4]
    assert calls["asterisms"] == [0.8]


def test_draw_terrain_layers_does_not_draw_dso_hover_info(monkeypatch) -> None:
    dso_hover_calls: list[str] = []

    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_deep_sky_shapes",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_dso_hover_info",
        lambda *_args, **_kwargs: dso_hover_calls.append("dso-hover"),
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    monkeypatch.setattr(
        pipeline_module, "_draw_urban_outline_layer", lambda *_args, **_kwargs: None
    )
    pipeline_module._draw_terrain_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=[(1.0, 10.0), (2.0, 20.0)],
        ),
        style=_make_style(show_dso=True),
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert dso_hover_calls == []


def test_draw_terrain_layers_skips_secondary_layers_while_simplified_view_active(
    monkeypatch,
) -> None:
    calls: list[str] = []

    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_deep_sky_shapes",
        lambda *_args, **_kwargs: calls.append("dso"),
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **_kwargs: calls.append("asterisms"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: calls.append("guides"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: calls.append("secondary"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: calls.append("water"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_urban_outline_layer",
        lambda *_args, **_kwargs: calls.append("urban"),
    )

    scene = replace(
        _make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=[(1.0, 10.0), (2.0, 20.0)],
            terrain_secondary_ridges_altaz_layers=[[(1.0, 10.0), (2.0, 20.0)]],
            terrain_secondary_ridges_distances_m_layers=[[10_000.0, 12_000.0]],
        ),
        water_overlay_dots=[object()],
    )

    pipeline_module._draw_terrain_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=scene,
        style=_make_style(
            show_dso=True,
            show_asterisms=True,
            show_guidelines=True,
            terrain_horizon_opacity=0.25,
            water_overlay_opacity=0.5,
            show_urban_outline_layer=True,
        ),
        simplified_view_active=True,
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls == ["dso", "asterisms", "guides"]


def test_render_scene_draws_dso_hover_immediately_before_overlay(monkeypatch) -> None:
    calls: list[str] = []

    monkeypatch.setattr(
        pipeline_module,
        "_clear_background_layer",
        lambda *_args, **_kwargs: calls.append("clear"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_background_layer",
        lambda *_args, **_kwargs: calls.append("background"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_guide_layer",
        lambda *_args, **_kwargs: calls.append("guide"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **_kwargs: calls.append("sky-cloud"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_terrain_layers",
        lambda *_args, **_kwargs: calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: calls.append("stars"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: calls.append("planets"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_satellite_layer",
        lambda *_args, **_kwargs: calls.append("satellites"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_aircraft_layer",
        lambda *_args, **_kwargs: calls.append("aircraft"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_overlay_layer",
        lambda *_args, **_kwargs: calls.append("overlay"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_hover_overlay_layer",
        lambda *_args, **_kwargs: calls.append("hover"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "render_text",
        pipeline_module.render_text,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: calls.append("labels"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_status_line",
        lambda *_args, **_kwargs: calls.append("status"),
    )

    geometry = SimpleNamespace(center=(100, 100), radius=80)
    viewport_rect = SimpleNamespace(width=lambda: 200, height=lambda: 200)
    scene = pipeline_module.RenderSceneData(
        viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        ),
        celestial_data=CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[],
            stars={
                "name": [],
                "source_id": [],
                "alt": [],
                "az": [],
                "vmag": [],
                "bv": [],
                "size_factor": [],
                "color_factor_base": [],
            },
            deep_sky_objects={},
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
        ),
        sky_disc_image=None,
        cloud_image=None,
        cloud_missing_mask=None,
        cloud_amount_field=None,
        terrain_horizon_profile=None,
        terrain_horizon_profile_distances_m=None,
        terrain_secondary_ridges_altaz_layers=None,
        terrain_secondary_ridges_distances_m_layers=None,
        urban_outlines=None,
        satellite_records_by_group=None,
        aircraft_snapshots=None,
    )
    style = pipeline_module.RenderStyle(
        theme=pipeline_module.THEME_STYLES_BY_PRESET["night"],
        visual_preset="night",
        text_font=object(),
        status_line_font=object(),
        show_background_gradient=True,
        show_custom_window_frame=False,
        show_observation_info=True,
        show_dso=True,
        show_asterisms=False,
        show_guidelines=True,
        enlarge_moon=False,
        bright_bodies_mode="outline",
        star_base_radius=4.0,
        star_visibility_boost=1.0,
        asterism_visibility_boost=1.0,
        earth_guide_visibility_boost=1.0,
        vmag_limit=6.0,
        sky_disc_altaz_rings="off",
        sky_disc_altaz_rings_hover="altaz",
        cloud_disc_alpha=0.0,
        satellite_opacity=0.0,
        terrain_horizon_opacity=0.0,
        earth_guide_opacity=0.0,
        urban_outline_opacity=0.0,
        show_urban_outline_layer=False,
        aircraft_opacity=0.0,
        tropical_cyclone_opacity=0.4,
        show_tropical_cyclone_overlay=True,
        star_render_expected_width=600,
    )
    hud = pipeline_module.RenderHudState(
        mouse_pos=None,
        overlay_info_bottom_left=False,
        viewport_interaction_mode=False,
        viewport_interaction_stars=None,
        status_message=None,
    )

    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(scene, geometry, viewport_rect),
        scene=scene,
        style=style,
        hud=hud,
        compositor=object(),
    )
    pipeline_module.render_hud_overlay_into_painter(
        painter=object(),
        frame=_make_frame(scene, geometry, viewport_rect),
        scene=scene,
        style=style,
        hud=hud,
        highlighted_object=({"name": "Vega"}, object()),
        highlighted_dso=({"name": "M31"}, object()),
    )

    assert calls == [
        "clear",
        "background",
        "sky-cloud",
        "guide",
        "terrain",
        "stars",
        "planets",
        "satellites",
        "aircraft",
        "labels",
        "hover",
        "overlay",
        "status",
    ]


def test_render_scene_reduces_layers_during_simplified_view(monkeypatch) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        pipeline_module, "_clear_background_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_background_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **kwargs: captured.update(
            {
                "cloud_disc_alpha": kwargs["style"].cloud_disc_alpha,
                "earth_guide_opacity": kwargs["style"].earth_guide_opacity,
                "sky_disc_image": kwargs["scene"].sky_disc_image,
            }
        ),
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_guide_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_terrain_layers",
        lambda *_args, **kwargs: captured.update(
            {"simplified_view_active": kwargs["simplified_view_active"]}
        ),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_main_terrain_profile_layer",
        lambda *_args, **_kwargs: captured.update({"main_terrain_profile": True}),
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_star_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_planet_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_satellite_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_aircraft_layer", lambda *_args, **_kwargs: None
    )

    scene = replace(
        _make_scene(),
        sky_disc_image=object(),
        night_light_glow_profile=object(),
    )
    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(100, 100), radius=80),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(cloud_disc_alpha=0.2, earth_guide_opacity=0.25),
        hud=_make_hud(simplified_view_enabled=True),
        compositor=object(),
    )

    assert captured == {
        "cloud_disc_alpha": 0.0,
        "earth_guide_opacity": 0.0,
        "sky_disc_image": scene.sky_disc_image,
        "main_terrain_profile": True,
        "simplified_view_active": True,
    }


def test_draw_sky_cloud_layers_skips_night_lights_while_simplified_view_active(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    class _Compositor:
        def draw(self, *_args, **kwargs) -> None:
            captured.update(
                {
                    "night_light_glow_profile": kwargs["night_light_glow_profile"],
                    "night_light_opacity": kwargs["night_light_opacity"],
                }
            )

    pipeline_module._draw_sky_cloud_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=80),
        scene=replace(_make_scene(), night_light_glow_profile=object()),
        style=_make_style(night_light_opacity=0.12, ridge_glow_opacity=0.34),
        compositor=_Compositor(),
        star_render_surface_size=(200, 200),
        simplified_view_active=True,
    )

    assert captured == {
        "night_light_glow_profile": None,
        "night_light_opacity": 0.0,
    }


def test_background_press_ignores_drag_exclusions() -> None:
    class _ProbeWindow(DraggableWindow):
        pass

    probe = _ProbeWindow()
    root = QWidget()
    root.resize(200, 200)
    root.show()
    excluded = QWidget(root)
    excluded.setGeometry(160, 160, 20, 20)
    excluded.show()
    probe.add_drag_exclusion(excluded)

    event = SimpleNamespace(
        button=lambda: Qt.MouseButton.LeftButton,
        position=lambda: SimpleNamespace(toPoint=lambda: QPoint(170, 170)),
        globalPosition=lambda: SimpleNamespace(toPoint=lambda: QPoint(170, 170)),
        accept=Mock(),
    )

    assert probe._begin_drag(None, event, root=root) is False
    assert probe._drag_press_pending is False


def test_render_hud_overlay_draws_persistent_search_label(monkeypatch) -> None:
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        pipeline_module,
        "_draw_hover_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_status_line",
        lambda *_args, **_kwargs: None,
    )

    def fake_draw_search_target_overlay(*_args, **kwargs) -> None:
        captured.update(kwargs)
        label_candidates = kwargs.get("label_candidates")
        if isinstance(label_candidates, list):
            label_candidates.append({"text": "Dubhe", "priority": 15})

    monkeypatch.setattr(
        pipeline_module.render_search_overlay,
        "draw_search_target_overlay",
        fake_draw_search_target_overlay,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args: captured.update({"labels": _args[1]}),
    )

    scene = replace(
        _make_scene(),
        viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        ),
    )
    pipeline_module.render_hud_overlay_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(100, 100), radius=80),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(star_render_expected_width=600),
        hud=_make_hud(),
        highlighted_object=None,
        highlighted_dso=None,
        label_candidates=[],
        search_overlay_target=SearchJumpTarget(
            label="Dubhe",
            kind="star",
            sort_key=(0.0, "dubhe"),
            alt_deg=50.0,
            az_deg=20.0,
            persistent_keep_marker=True,
        ),
    )

    assert captured["draw_label"] is True
    assert captured["labels"] == [{"text": "Dubhe", "priority": 15}]


def test_collect_visible_named_star_labels_returns_only_named_visible_stars() -> None:
    scene = _make_scene(
        celestial_data=CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[],
            stars={
                "star_index": np.array([11, 12], dtype=np.int32),
                "alt": np.array([45.0, 10.0]),
                "az": np.array([180.0, 10.0]),
                "vmag": np.array([5.96, 1.0]),
                "bv": np.array([0.869, 0.0]),
                "size_factor": np.array([1.0, 1.0]),
                "color_factor_base": np.array([1.0, 1.0]),
            },
            deep_sky_objects={},
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
            star_catalog_meta=StarCatalogMeta(
                name_indices=np.array([11], dtype=np.int32),
                names=np.array(["Dubhe"], dtype=object),
                source_id_indices=np.array([], dtype=np.int32),
                source_ids=np.array([], dtype=object),
            ),
        ),
    )

    labels = pipeline_module.render_stars.collect_visible_named_star_labels(
        scene.celestial_data,
        scene.viewer,
        SimpleNamespace(center=(200, 200), radius=200),
        star_base_radius=4.0,
        draw_vmag_limit=6.0,
        viewport_size=(400, 400),
    )

    assert len(labels) == 1
    assert labels[0][0] == "Dubhe"
    assert labels[0][1].x() == pytest.approx(200.0)
    assert labels[0][1].y() == pytest.approx(200.0)
    assert labels[0][2] == (255, 210, 161)


def test_render_hud_overlay_draws_simplified_named_star_labels_at_fixed_offset(
    monkeypatch,
) -> None:
    captured: list[tuple[str, float, float, int]] = []

    monkeypatch.setattr(
        pipeline_module,
        "_draw_hover_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_status_line",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "collect_visible_named_star_labels",
        lambda *_args, **_kwargs: [
            ("Dubhe", QPointF(120.0, 80.0), (10, 20, 30)),
            ("Merak", QPointF(150.0, 100.0), (40, 50, 60)),
        ],
    )

    def fake_draw_outlined_text(*_args, **kwargs) -> None:
        text = _args[1]
        pos = _args[2]
        captured.append(
            (text, float(pos.x()), float(pos.y()), kwargs["style"].text_color.alpha())
        )
        style = kwargs["style"]
        assert style.font is not None
        assert style.outline_width == 0.0
        assert style.outline_color.alpha() == 0
        if text == "Dubhe":
            expected = pipeline_module.render_text.blend_color_toward_white(
                QColor(10, 20, 30),
                amount=pipeline_module.render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
            )
            assert style.text_color.red() == expected.red()
            assert style.text_color.green() == expected.green()
            assert style.text_color.blue() == expected.blue()
            assert style.text_color.alpha() == int(round(255 * 0.4))
        if text == "Merak":
            expected = pipeline_module.render_text.blend_color_toward_white(
                QColor(40, 50, 60),
                amount=pipeline_module.render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
            )
            assert style.text_color.red() == expected.red()
            assert style.text_color.green() == expected.green()
            assert style.text_color.blue() == expected.blue()
            assert style.text_color.alpha() == int(round(255 * 0.4))

    monkeypatch.setattr(
        pipeline_module.render_text,
        "draw_outlined_text",
        fake_draw_outlined_text,
    )

    scene = _make_scene()
    style = _make_style(text_font=QFont(), status_line_font=QFont())
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
    pipeline_module.render_hud_overlay_into_painter(
        painter=painter,
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(200, 200), radius=200),
            QRect(0, 0, 400, 400),
        ),
        scene=scene,
        style=style,
        hud=_make_hud(simplified_view_enabled=True),
        highlighted_object=None,
        highlighted_dso=None,
        label_candidates=[],
    )
    painter.end()

    assert len(captured) == 2
    expected_positions = {
        "Dubhe": QPointF(120.0, 80.0),
        "Merak": QPointF(150.0, 100.0),
    }
    for text, x, y, alpha in captured:
        assert alpha == int(round(255 * 0.4))
        bounds = pipeline_module.render_text._text_bounds_at_baseline(
            text,
            QFont(),
            QPointF(0.0, 0.0),
        )
        expected = expected_positions[text]
        assert x + float(bounds.left()) == pytest.approx(float(expected.x()))
        assert y + float(bounds.bottom()) == pytest.approx(float(expected.y()))


def test_render_hud_overlay_skips_simplified_labels_when_disabled(monkeypatch) -> None:
    labels_drawn: list[str] = []

    monkeypatch.setattr(
        pipeline_module,
        "_draw_hover_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_status_line",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "collect_visible_named_star_labels",
        lambda *_args, **_kwargs: [
            ("Dubhe", QPointF(120.0, 80.0), (10, 20, 30)),
        ],
    )

    def fake_draw_outlined_text(*_args, **_kwargs) -> None:
        labels_drawn.append(str(_args[1]))

    monkeypatch.setattr(
        pipeline_module.render_text,
        "draw_outlined_text",
        fake_draw_outlined_text,
    )

    scene = _make_scene()
    style = _make_style(text_font=QFont(), status_line_font=QFont())
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
    pipeline_module.render_hud_overlay_into_painter(
        painter=painter,
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(200, 200), radius=200),
            QRect(0, 0, 400, 400),
        ),
        scene=scene,
        style=style,
        hud=_make_hud(
            simplified_view_enabled=True, simplified_view_labels_enabled=False
        ),
        highlighted_object=None,
        highlighted_dso=None,
        label_candidates=[],
    )
    painter.end()

    assert labels_drawn == []


def test_render_fast_overlay_layers_passes_simplified_satellite_labels(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        pipeline_module.render_satellites,
        "draw_satellite_overlay",
        lambda *_args, **kwargs: captured.update(
            {"draw_simplified_labels": kwargs.get("draw_simplified_labels")}
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_aircraft,
        "draw_aircraft_overlay",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_tropical_cyclones,
        "draw_tropical_cyclone_overlay",
        lambda *_args, **_kwargs: None,
    )

    scene = _make_scene()
    style = _make_style(
        satellite_opacity=1.0,
        aircraft_opacity=0.0,
        tropical_cyclone_opacity=0.0,
        text_font=QFont(),
    )
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
    try:
        pipeline_module.render_fast_overlay_layers_into_painter(
            painter=painter,
            frame=_make_frame(
                scene,
                SimpleNamespace(center=(200, 200), radius=200),
                QRect(0, 0, 400, 400),
            ),
            scene=scene,
            style=style,
            draw_labels=False,
            draw_simplified_satellite_labels=True,
        )
    finally:
        painter.end()

    assert captured == {"draw_simplified_labels": True}


def test_render_hud_overlay_forwards_simplified_satellite_label_flag(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        pipeline_module.render_overlay_info,
        "draw_overlay_info",
        lambda *_args, **kwargs: captured.update(
            {
                "draw_simplified_satellite_labels": kwargs[
                    "draw_simplified_satellite_labels"
                ],
                "highlighted_satellite": _args[9],
            }
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_simplified_named_star_labels",
        lambda *_args, **_kwargs: None,
    )

    scene = _make_scene()
    style = _make_style(text_font=QFont())
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
    highlighted_satellite = (
        SatelliteOverlayPoint(
            group_key="iss",
            satellite_name="ISS",
            alt_deg=12.0,
            az_deg=220.0,
            marker_scale=0.42,
        ),
        QPointF(130.0, 95.0),
    )
    try:
        pipeline_module.render_hud_overlay_into_painter(
            painter=painter,
            frame=_make_frame(
                scene,
                SimpleNamespace(center=(200, 200), radius=200),
                QRect(0, 0, 400, 400),
            ),
            scene=scene,
            style=style,
            hud=_make_hud(
                simplified_view_enabled=True,
                simplified_view_labels_enabled=True,
            ),
            highlighted_object=None,
            highlighted_dso=None,
            highlighted_satellite=highlighted_satellite,
            label_candidates=[],
        )
    finally:
        painter.end()

    assert captured["draw_simplified_satellite_labels"] is True
    assert captured["highlighted_satellite"] == highlighted_satellite


def test_render_scene_keeps_sky_bitmap_during_viewport_interaction(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}
    scene = replace(_make_scene(), sky_disc_image=object())

    monkeypatch.setattr(
        pipeline_module, "_clear_background_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_background_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_guide_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **kwargs: captured.update(
            {
                "cloud_disc_alpha": kwargs["style"].cloud_disc_alpha,
                "sky_disc_image": kwargs["scene"].sky_disc_image,
            }
        ),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_viewport_interaction_layers",
        lambda *_args, **_kwargs: None,
    )

    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(100, 100), radius=80),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(cloud_disc_alpha=0.2),
        hud=_make_hud(viewport_interaction_mode=True),
        compositor=object(),
    )

    assert captured == {"cloud_disc_alpha": 0.0, "sky_disc_image": scene.sky_disc_image}


def test_draw_guide_layer_draws_zenith_marker(monkeypatch) -> None:
    calls: list[str] = []

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: calls.append("direction"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: calls.append("zenith"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_celestial_pole_markers",
        lambda *_args, **_kwargs: calls.append("poles"),
    )

    pipeline_module._draw_guide_layer(
        painter=object(),
        geometry=SimpleNamespace(center=(100, 100), radius=80),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(),
        style=_make_style(show_guidelines=True),
    )

    assert calls == ["direction", "zenith", "poles"]


def test_render_base_scene_can_skip_fast_overlays(monkeypatch) -> None:
    calls: list[str] = []
    scene = _make_scene()
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    viewport_rect = SimpleNamespace(width=lambda: 200, height=lambda: 200)

    monkeypatch.setattr(
        pipeline_module,
        "_clear_background_layer",
        lambda *_args, **_kwargs: calls.append("clear"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_background_layer",
        lambda *_args, **_kwargs: calls.append("background"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_guide_layer",
        lambda *_args, **_kwargs: calls.append("guide"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **_kwargs: calls.append("sky-cloud"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_terrain_layers",
        lambda *_args, **_kwargs: calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: calls.append("stars"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: calls.append("planets"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_satellite_layer",
        lambda *_args, **_kwargs: calls.append("satellites"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_aircraft_layer",
        lambda *_args, **_kwargs: calls.append("aircraft"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "render_text",
        pipeline_module.render_text,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: calls.append("labels"),
    )

    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(scene, geometry, viewport_rect),
        scene=scene,
        style=_make_style(),
        hud=_make_hud(),
        compositor=object(),
        draw_fast_overlays=False,
    )

    assert calls == [
        "clear",
        "background",
        "sky-cloud",
        "guide",
        "terrain",
        "stars",
        "planets",
        "labels",
    ]


def test_draw_planet_layer_passes_marker_scale(monkeypatch) -> None:
    seen_marker_scales: list[float] = []

    monkeypatch.setattr(
        pipeline_module.render_solar_system,
        "draw_solar_system_bodies",
        lambda *_args, **kwargs: seen_marker_scales.append(
            float(kwargs.get("marker_scale", 1.0))
        ),
    )

    pipeline_module._draw_planet_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=_make_scene(),
        style=_make_style(star_render_expected_width=600),
        enlarge_moon=False,
        label_candidates=[],
    )

    assert seen_marker_scales == [
        pipeline_module.compute_star_render_upscale_factor(1200, 600)
    ]


def test_draw_hover_overlay_passes_marker_scale_to_moon(monkeypatch) -> None:
    seen_marker_scales: list[float] = []

    monkeypatch.setattr(
        render_solar_system_module,
        "draw_hovered_moon_overlay",
        lambda *_args, **kwargs: seen_marker_scales.append(
            float(kwargs.get("marker_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        render_overlay_info_module,
        "draw_overlay_info",
        lambda *_args, **_kwargs: None,
    )

    scene = _make_scene(
        celestial_data=CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[
                PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True),
                PlanetBody(
                    name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True
                ),
            ],
            stars={
                "star_index": np.array([], dtype=np.int64),
                "alt": np.array([], dtype=np.float64),
                "az": np.array([], dtype=np.float64),
                "vmag": np.array([], dtype=np.float64),
                "bv": np.array([], dtype=np.float64),
                "size_factor": np.array([], dtype=np.float64),
                "color_factor_base": np.array([], dtype=np.float64),
            },
            deep_sky_objects={
                "id": np.array([], dtype=np.int64),
                "name": np.array([], dtype=object),
                "type": np.array([], dtype=object),
                "alt": np.array([], dtype=np.float64),
                "az": np.array([], dtype=np.float64),
                "vmag": np.array([], dtype=np.float64),
                "major_arcmin": np.array([], dtype=np.float64),
                "minor_arcmin": np.array([], dtype=np.float64),
                "pa_deg": np.array([], dtype=np.float64),
            },
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
        )
    )

    pipeline_module._draw_hover_overlay_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=QRect(0, 0, 1200, 1200),
        scene=scene,
        style=_make_style(
            star_render_expected_width=600,
            show_asterisms=False,
            show_observation_info=False,
        ),
        highlighted_object=({"name": "moon"}, QPointF(10.0, 10.0)),
        highlighted_dso=None,
    )

    assert seen_marker_scales == [
        pipeline_module.compute_star_render_upscale_factor(1200, 600)
    ]


def test_draw_hover_overlay_requires_direction_marker_hover(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "resolve_direction_marker_hover",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_altitude_ring_overlay",
        lambda *_args, **_kwargs: calls.append("background-rings"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_grid_overlay",
        lambda *_args, **_kwargs: calls.append("grid"),
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_solar_system,
        "draw_hovered_moon_overlay",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_dso_hover_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        render_overlay_info_module,
        "draw_overlay_info",
        lambda *_args, **_kwargs: None,
    )

    pipeline_module._draw_hover_overlay_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=QRect(0, 0, 1200, 1200),
        scene=_make_scene(celestial_data=object()),
        style=_make_style(
            show_guidelines=True,
            sky_disc_altaz_rings_hover="altaz",
        ),
        highlighted_object=None,
        highlighted_dso=None,
        mouse_pos=QPoint(600, 600),
    )

    assert calls == []


def test_draw_overlay_layer_skips_static_info_when_disabled(monkeypatch) -> None:
    draw_overlay_info = Mock()
    monkeypatch.setattr(
        pipeline_module.render_overlay_info, "draw_overlay_info", draw_overlay_info
    )

    pipeline_module._draw_overlay_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=100),
        viewport_rect=SimpleNamespace(height=lambda: 200),
        scene=_make_scene(),
        style=_make_style(show_observation_info=False),
        mouse_pos=None,
        overlay_info_bottom_left=False,
        highlighted_object=None,
        highlighted_dso=None,
        enlarge_moon=False,
        label_reservations=[],
        label_candidates=[],
    )

    draw_overlay_info.assert_not_called()


def test_draw_background_layer_skips_gradient_when_disabled(monkeypatch) -> None:
    draw_radial_background = Mock()
    draw_window_border = Mock()
    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_radial_background",
        draw_radial_background,
    )
    monkeypatch.setattr(
        pipeline_module.render_background, "draw_window_border", draw_window_border
    )

    pipeline_module._draw_background_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=100),
        viewport_rect=SimpleNamespace(),
        scene=_make_scene(),
        style=_make_style(show_background_gradient=False),
    )

    draw_radial_background.assert_not_called()
    draw_window_border.assert_not_called()


def test_draw_background_layer_skips_custom_frame_when_disabled(monkeypatch) -> None:
    draw_radial_background = Mock()
    draw_window_border = Mock()
    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_radial_background",
        draw_radial_background,
    )
    monkeypatch.setattr(
        pipeline_module.render_background, "draw_window_border", draw_window_border
    )

    pipeline_module._draw_background_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=100),
        viewport_rect=QRect(0, 0, 200, 200),
        scene=_make_scene(),
        style=_make_style(show_custom_window_frame=False),
    )

    draw_radial_background.assert_called_once()
    draw_window_border.assert_not_called()


def test_draw_hover_overlay_layer_enlarges_hovered_moon_by_name(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_solar_system,
        "draw_hovered_moon_overlay",
        lambda *_args, **_kwargs: calls.append("moon-hover"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_dso_hover_layer",
        lambda *_args, **_kwargs: calls.append("dso-hover"),
    )
    monkeypatch.setattr(
        render_overlay_info_module,
        "draw_overlay_info",
        lambda *_args, **_kwargs: calls.append("overlay-info"),
    )

    pipeline_module._draw_hover_overlay_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=QRect(0, 0, 1200, 1200),
        scene=_make_scene(celestial_data=object()),
        style=_make_style(),
        highlighted_object=({"name": "moon"}, object()),
        highlighted_dso=None,
    )

    assert calls == ["moon-hover", "dso-hover", "overlay-info"]


def test_draw_sky_reference_lines_uses_render_view_center_projection() -> None:
    calls: list[tuple[float, float]] = []

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

    render_guides_module.draw_sky_reference_lines(
        _Painter(),
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        celestial_data=SimpleNamespace(
            celestial_equator_points=[(1.0, 2.0)],
            ecliptic_points=[(3.0, 4.0)],
            horizon_points=[(5.0, 6.0)],
        ),
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(55.0, 200.0),
            observer_height_m=10.0,
        ),
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            calls.append(_view_center) or (alt, az)
        ),
    )

    assert len(calls) == 6
    assert calls == [(55.0, 200.0)] * 6


def test_draw_sky_reference_lines_uses_wider_dash_patterns(monkeypatch) -> None:
    dash_patterns: list[list[int]] = []
    pen_styles: list[object] = []
    pen_alpha_values: list[int] = []
    pen_widths: list[float] = []

    class _FakePen:
        def __init__(self, color, _width, style=None) -> None:
            self._dash_pattern: list[int] = []
            self._style = style
            self._alpha = color.alpha() if hasattr(color, "alpha") else None
            self._width = float(_width)

        def setCosmetic(self, *_args, **_kwargs) -> None:
            pass

        def setDashPattern(self, pattern) -> None:
            self._dash_pattern = list(pattern)

        def setCapStyle(self, *_args, **_kwargs) -> None:
            pass

        def setJoinStyle(self, *_args, **_kwargs) -> None:
            pass

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            dash_patterns.append(list(getattr(pen, "_dash_pattern", [])))
            pen_styles.append(getattr(pen, "_style", None))
            pen_alpha_values.append(getattr(pen, "_alpha", None))
            pen_widths.append(getattr(pen, "_width", None))

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

    monkeypatch.setattr(render_guides_module, "QPen", _FakePen)

    render_guides_module.draw_sky_reference_lines(
        _Painter(),
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        celestial_data=SimpleNamespace(
            celestial_equator_points=[(0.10, 0.20), (0.12, 0.22)],
            ecliptic_points=[(0.30, 0.40), (0.32, 0.42)],
            horizon_points=[(0.50, 0.60), (0.52, 0.62)],
        ),
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(55.0, 200.0),
            observer_height_m=10.0,
        ),
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (alt, az),
    )

    assert dash_patterns[0::3] == [[], [], []]
    assert dash_patterns[1::3] == [[], [], []]
    assert dash_patterns[2::3] == [[16, 6], [4, 6], [10, 1]]
    assert pen_styles[0::3] == [
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
    ]
    assert pen_styles[1::3] == [
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
    ]
    assert pen_styles[2::3] == [None, None, None]
    assert pen_alpha_values[0::3] == [18, 18, 18]
    assert pen_alpha_values[1::3] == [30, 30, 30]
    assert pen_alpha_values[2::3] == [255, 255, 255]
    assert [round(width, 3) for width in pen_widths[0::3]] == [1.100, 1.254, 1.100]
    assert [round(width, 3) for width in pen_widths[1::3]] == [0.750, 0.855, 0.750]
    assert [round(width, 3) for width in pen_widths[2::3]] == [0.510, 0.627, 0.550]


def test_draw_direction_labels_uses_horizon_line_color(monkeypatch) -> None:
    seen_colors: list[tuple[int, int, int]] = []

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setFont(self, _font) -> None:
            pass

        def setPen(self, _pen) -> None:
            pass

        def drawLine(self, *_args, **_kwargs) -> None:
            pass

        def viewport(self):
            return QRect(0, 0, 200, 200)

    def _capture(*_args, style=None, **_kwargs):
        seen_colors.append(
            (style.text_color.red(), style.text_color.green(), style.text_color.blue())
        )

    monkeypatch.setattr(render_guides_module, "draw_outlined_text", _capture)
    monkeypatch.setattr(
        render_guides_module, "is_in_fov", lambda *_args, **_kwargs: True
    )
    monkeypatch.setattr(
        render_guides_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_guides_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    render_guides_module.draw_direction_labels(
        _Painter(),
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(30.0, 40.0),
            edge_fov_deg=95.0,
            content_fov_deg=180.0,
        ),
        text_font=QFont(),
        theme=window_module.THEME_STYLES_BY_PRESET["white"],
    )

    assert seen_colors
    assert all(
        color == render_guides_module.HORIZON_LINE_COLOR for color in seen_colors
    )


def test_draw_zenith_marker_uses_horizon_line_color_for_all_themes(
    monkeypatch,
) -> None:
    seen_colors: list[tuple[int, int, int]] = []

    class _FakePen:
        def __init__(self, color, _width) -> None:
            seen_colors.append((color.red(), color.green(), color.blue()))

    class _Painter:
        def setPen(self, _pen) -> None:
            pass

        def drawLine(self, *_args, **_kwargs) -> None:
            pass

    monkeypatch.setattr(render_guides_module, "QPen", _FakePen)
    monkeypatch.setattr(
        render_guides_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_guides_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_guides_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    for theme in THEME_STYLES_BY_PRESET.values():
        render_guides_module.draw_zenith_marker(
            _Painter(),
            geometry=SimpleNamespace(center=(0, 0), radius=1),
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(30.0, 40.0),
                edge_fov_deg=95.0,
                content_fov_deg=180.0,
            ),
            theme=theme,
        )

    assert seen_colors
    assert all(
        color == render_guides_module.HORIZON_LINE_COLOR for color in seen_colors
    )


def test_draw_celestial_pole_markers_uses_celestial_equator_color_for_all_themes(
    monkeypatch,
) -> None:
    seen_colors: list[tuple[int, int, int]] = []
    seen_positions: list[tuple[float, float]] = []

    class _FakePen:
        def __init__(self, color, _width) -> None:
            seen_colors.append((color.red(), color.green(), color.blue()))

    class _Painter:
        def setPen(self, _pen) -> None:
            pass

        def drawLine(self, *_args, **_kwargs) -> None:
            pass

    monkeypatch.setattr(render_guides_module, "QPen", _FakePen)
    monkeypatch.setattr(
        render_guides_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_guides_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (
            seen_positions.append((float(alt), float(az))) or (float(az), float(alt))
        ),
    )
    monkeypatch.setattr(
        render_guides_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    for theme in THEME_STYLES_BY_PRESET.values():
        render_guides_module.draw_celestial_pole_markers(
            _Painter(),
            geometry=SimpleNamespace(center=(0, 0), radius=1),
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(30.0, 40.0),
                edge_fov_deg=95.0,
                content_fov_deg=180.0,
            ),
        )

    assert seen_colors
    assert all(
        color == render_guides_module.CELESTIAL_EQUATOR_COLOR for color in seen_colors
    )
    assert seen_positions[0::2] == [(35.0, 0.0)] * len(THEME_STYLES_BY_PRESET)
    assert seen_positions[1::2] == [(-35.0, 180.0)] * len(THEME_STYLES_BY_PRESET)


def test_draw_urban_outlines_clips_two_point_outline_out_of_view() -> None:
    class _Painter:
        def __init__(self) -> None:
            self.polylines: list[list[tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "color") or not hasattr(pen, "widthF"):
                return

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, poly) -> None:
            self.polylines.append(
                [(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())]
            )

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(points=[(-10.0, 10.0), (-12.0, 10.3)], height_m=50.0)
        ],
        viewer=_viewer(),
        opacity=0.38,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.polylines == []


def test_draw_urban_outlines_uses_fixed_alpha_and_near_underlay(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.alpha_values: list[int] = []
            self.width_values: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "color") or not hasattr(pen, "widthF"):
                return
            self.alpha_values.append(int(pen.color().alpha()))
            self.width_values.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(-10.0, 10.0), (-12.0, 12.0)],
                height_m=0.0,
                distance_km=0.01,
            ),
            UrbanOutlinePolyline(
                points=[(-10.0, 20.0), (-12.0, 22.0)],
                height_m=50.0,
                distance_km=0.5,
            ),
        ],
        viewer=_viewer(),
        opacity=0.2,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.alpha_values == [3, 2, 1, 51, 3, 2, 1, 51]
    assert [round(width, 1) for width in painter.width_values] == [
        9.2,
        7.2,
        4.4,
        2.3,
        2.4,
        2.0,
        2.0,
        1.3,
    ]


def test_draw_urban_outlines_allows_sub_unit_width_scale(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.width_values: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "widthF"):
                return
            self.width_values.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(-10.0, 10.0), (-12.0, 12.0)],
                height_m=50.0,
                distance_km=0.01,
            )
        ],
        viewer=_viewer(),
        opacity=0.38,
        line_width_scale=0.5,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(alt),
            float(az),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert [round(width, 2) for width in painter.width_values[:4]] == [
        4.6,
        3.6,
        2.2,
        1.14,
    ]


def test_draw_urban_outlines_scales_alpha_for_subpixel_widths(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.alpha_values: list[int] = []
            self.width_values: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "color") or not hasattr(pen, "widthF"):
                return
            self.alpha_values.append(int(pen.color().alpha()))
            self.width_values.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(-10.0, 10.0), (-12.0, 12.0)],
                height_m=0.0,
                distance_km=0.01,
            )
        ],
        viewer=_viewer(),
        opacity=0.2,
        line_width_scale=0.1,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.alpha_values == [3, 1, 0, 12]
    assert [round(width, 3) for width in painter.width_values] == [
        0.92,
        0.72,
        0.44,
        0.228,
    ]


def test_draw_urban_outlines_thickens_tall_buildings(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.width_values: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            if not hasattr(pen, "widthF"):
                return
            self.width_values.append(float(pen.widthF()))

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(
                points=[(-10.0, 10.0), (-12.0, 12.0)],
                height_m=600.0,
                distance_km=0.01,
            )
        ],
        viewer=_viewer(),
        opacity=0.38,
        line_width_scale=0.5,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(alt),
            float(az),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert [round(width, 2) for width in painter.width_values[:4]] == [
        9.2,
        7.2,
        4.4,
        1.14,
    ]


def test_draw_terrain_horizon_line_scales_line_widths(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_widths: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.pen_widths.append(float(pen.widthF()))

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    render_terrain_module._draw_terrain_profile_layer(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=_viewer(),
        terrain_profile_altaz=[(0.0, 0.0), (0.1, 0.1)],
        terrain_profile_distances_m=None,
        spec=_terrain_horizon_spec(opacity=0.38, line_width_scale=2.0, fast_mode=True),
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
        split_by_gaps_func=lambda points: [points],
    )

    assert painter.pen_widths == [
        pytest.approx(render_terrain_module.TERRAIN_HORIZON_FAST_WIDTH * 2.0)
    ]


def test_draw_terrain_horizon_line_scales_widths_by_distance(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_widths: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.pen_widths.append(float(pen.widthF()))

        def drawLine(self, *_args) -> None:
            pass

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module._draw_terrain_profile_layer(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=_viewer(),
        terrain_profile_altaz=[(0.0, 0.0), (0.0, 0.1), (0.0, 0.2)],
        terrain_profile_distances_m=[1_000.0, 50_000.0, 120_000.0],
        spec=_terrain_horizon_spec(opacity=0.38, line_width_scale=1.0, fast_mode=False),
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
        split_by_gaps_func=lambda points: [points],
    )

    assert len(painter.pen_widths) == 2
    assert painter.pen_widths[0] > painter.pen_widths[1]


def test_draw_terrain_secondary_ridges_use_fixed_widths(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_widths: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.pen_widths.append(float(pen.widthF()))

        def drawLine(self, *_args) -> None:
            pass

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        ),
        terrain_secondary_ridges_layers=[
            [(0.0, 0.0), (0.0, 0.1), (0.0, 0.2)],
            [(0.1, 0.0), (0.1, 0.1), (0.1, 0.2)],
            [(0.2, 0.0), (0.2, 0.1), (0.2, 0.2)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1_000.0, 2_000.0, 3_000.0],
            [10_000.0, 12_000.0, 15_000.0],
            [50_000.0, 60_000.0, 70_000.0],
        ],
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.pen_widths) == 6
    assert all(width > 0.0 for width in painter.pen_widths)
    assert painter.pen_widths[0] > painter.pen_widths[1]
    assert painter.pen_widths[1] > painter.pen_widths[-1]


def test_draw_terrain_secondary_ridges_swaps_visible_and_occluded_colors(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_rgbs: list[tuple[int, int, int]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            color = QColor(pen.color())
            self.pen_rgbs.append((color.red(), color.green(), color.blue()))

        def drawLine(self, *_args) -> None:
            pass

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        ),
        terrain_secondary_ridges_layers=[
            [(20.0, 0.0), (20.0, 0.1)],
            [(10.0, 0.0), (10.0, 0.1)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1_000.0, 2_000.0],
            [10_000.0, 12_000.0],
        ],
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.pen_rgbs) == 4
    assert all(
        rgb == render_terrain_module.TERRAIN_SECONDARY_RIDGE_OCCLUDED_COLOR_RGB
        for rgb in painter.pen_rgbs
    )


def test_secondary_ridge_alpha_base_is_lower() -> None:
    assert render_terrain_module.terrain_secondary_ridge_line_alpha(0.38) < 0.08


def test_secondary_ridge_overlay_alpha_is_scaled_down(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_terrain_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.alphas: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.alphas.append(float(pen.color().alphaF()))

        def drawLine(self, *_args) -> None:
            pass

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        ),
        terrain_secondary_ridges_layers=[
            [(20.0, 0.0), (20.0, 0.1)],
            [(10.0, 0.0), (10.0, 0.1)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1_000.0, 2_000.0],
            [10_000.0, 12_000.0],
        ],
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.alphas) == 4
    assert all(0.0 <= alpha <= 1.0 for alpha in painter.alphas)
    assert painter.alphas[0] < painter.alphas[1]


def test_draw_terrain_secondary_ridges_bridges_seam_near_zero(monkeypatch) -> None:
    monkeypatch.setattr(
        render_terrain_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "split_by_gaps",
        lambda points: [points],
    )

    class _Painter:
        def __init__(self) -> None:
            self.lines: list[tuple[tuple[float, float], tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawLine(self, start, end) -> None:
            self.lines.append(
                ((float(start.x()), float(start.y())), (float(end.x()), float(end.y())))
            )

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=SimpleNamespace(
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        ),
        terrain_secondary_ridges_layers=[
            [(5.0, 359.0), (5.0, 0.0), (5.0, 1.0)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1_000.0, 2_000.0, 3_000.0],
        ],
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            math.sin(math.radians(float(az))),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
        split_by_gaps_func=lambda points: [points],
    )

    assert painter.lines == []


def test_draw_terrain_horizon_line_uses_edge_fov_for_projection() -> None:
    class _Painter:
        def __init__(self) -> None:
            self.polylines: list[list[tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, poly) -> None:
            self.polylines.append(
                [(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())]
            )

    def _render(edge_fov_deg: float) -> list[tuple[float, float]]:
        painter = _Painter()
        render_terrain_module._draw_terrain_profile_layer(
            painter,
            geometry=SimpleNamespace(center=(0, 0), radius=1),
            viewer=_viewer(edge_fov_deg=edge_fov_deg, content_fov_deg=180.0),
            terrain_profile_altaz=[(0.0, 180.0), (0.0, 190.0)],
            terrain_profile_distances_m=None,
            spec=_terrain_horizon_spec(opacity=0.38, fast_mode=True),
            is_in_fov_func=lambda *_args, **_kwargs: True,
            altaz_to_normalized_xy_func=render_terrain_module.altaz_to_normalized_xy,
            normalized_to_screen_xy_func=lambda nx, ny, _geometry: (
                float(nx),
                float(ny),
            ),
            split_by_gaps_func=lambda points: [points],
        )
        assert painter.polylines
        return painter.polylines[0]

    points_90 = _render(90.0)
    points_120 = _render(120.0)

    assert points_120[0][1] < points_90[0][1]


def test_draw_terrain_horizon_line_rotates_profile_away_from_north_seam() -> None:
    class _Painter:
        def __init__(self) -> None:
            self.polylines: list[list[tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, poly) -> None:
            self.polylines.append(
                [(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())]
            )

    painter = _Painter()
    render_terrain_module._draw_terrain_profile_layer(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer=_viewer(edge_fov_deg=120.0, content_fov_deg=180.0),
        terrain_profile_altaz=[
            (0.0, 0.0),
            (0.0, 10.0),
            (0.0, 180.0),
            (0.0, 190.0),
        ],
        terrain_profile_distances_m=None,
        spec=_terrain_horizon_spec(opacity=0.38, fast_mode=True),
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        split_by_gaps_func=lambda points: [points],
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.polylines
    assert [point[0] for point in painter.polylines[0]] == [0.0, 10.0, 180.0, 190.0]


def test_draw_star_layer_forwards_outline_flag(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def fake_draw_stars(*_args, **kwargs) -> None:
        captured.update(kwargs)

    monkeypatch.setattr(
        pipeline_module.render_stars, "draw_stars_fast", fake_draw_stars
    )

    pipeline_module._draw_star_layer(
        painter=object(),
        geometry=SimpleNamespace(center=(120, 120), radius=80),
        viewport_rect=QRect(0, 0, 240, 240),
        scene=_make_scene(),
        style=_make_style(bright_bodies_mode="outline"),
        star_render_surface_size=(240, 240),
        fast_mode=True,
    )

    assert captured["outline_bright_bodies"] is True
    assert "fast_mode" not in captured
