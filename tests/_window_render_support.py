# ruff: noqa: F401

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
import zstarview.render.instrument_background as render_instrument_background_module
import zstarview.render.atlas_pipeline as atlas_pipeline_module
import zstarview.render.zstarview_pipeline as zstarview_pipeline_module
import zstarview.render.solar_system as render_solar_system_module
import zstarview.render.terrain as render_terrain_module
import zstarview.render.text as render_text_module
from zstarview.gui.famous_star_shortcuts import SearchJumpTarget
from zstarview.gui.window import SkyWindow
from zstarview.gui.window_state import SkyWindowState
from zstarview.location_resolver import PlaceTargetProjection
from zstarview.paths import NIGHT_LIGHT_DEFAULT_OPACITY, THEME_STYLES_BY_PRESET
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

    def stop(self) -> None:
        self._active = False


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
        self.ridge_glow_opacity = values.get("ridge_glow_opacity", 0.04)
        self.water_overlay_opacity = values.get("water_overlay_opacity", 0.4)
        self.ground_tint_opacity = values.get("ground_tint_opacity", 0.025)
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
            SimpleNamespace(
                image=None,
                missing_mask=None,
                cloud_amount_field=None,
                altaz_grid=None,
            ),
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
            "_sky_disc_alpha_when_enabled", 0.15
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
            "_night_light_opacity_when_enabled", NIGHT_LIGHT_DEFAULT_OPACITY
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

    def reproject_geo_satellite_overlay(self, *_args, **_kwargs) -> None:
        return None

    def reproject_cloud_overlay(self, *_args, **kwargs) -> None:
        start = self.__dict__.get("start_background_cloud_update")
        if callable(start):
            start(**kwargs)

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

    def _finalize_view_direction_change(self) -> None:
        finalize = self.__dict__.get("_finalize_view_direction_change")
        if callable(finalize):
            finalize()
            return
        window_module.SkyWindow._finalize_view_direction_change(self)

    def _sync_viewport_interaction_chrome_visibility(self) -> None:
        sync = self.__dict__.get("_sync_viewport_interaction_chrome_visibility")
        if callable(sync):
            sync()
            return
        menu_button = self.__dict__.get("menu_button")
        state = self.__dict__.get("state")
        if menu_button is not None and state is not None:
            menu_button.setVisible(not bool(state.viewport_interaction_mode))

    def _sync_cloud_action_state(self) -> None:
        sync = self.__dict__.get("_sync_cloud_action_state")
        if callable(sync):
            sync()

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
        cloud_missing_mask=None,
        cloud_altaz_grid=None,
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
        values["theme"] = THEME_STYLES_BY_PRESET.get(
            values["visual_preset"],
            THEME_STYLES_BY_PRESET["night"],
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


__all__ = [name for name in globals() if not name.startswith("__")]
