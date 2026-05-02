from __future__ import annotations

import astropy.time
import math
from dataclasses import replace
from datetime import datetime, timezone
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from PySide6.QtCore import QPoint, QPointF, QRect, Qt
from PySide6.QtGui import QColor, QFont, QImage

import zstarview.render.pipeline as pipeline_module
import zstarview.render.guides as render_guides_module
import zstarview.render.overlay_info as render_overlay_info_module
import zstarview.render.search_overlay as search_overlay_module
import zstarview.render.terrain as render_terrain_module
import zstarview.render.text as render_text_module
import zstarview.gui.window as window_module
import zstarview.gui.window_render as window_render_module
import zstarview.gui.window_updates as window_updates_module
import zstarview.render.solar_system as render_solar_system_module
from zstarview.location_resolver import PlaceTargetProjection
from zstarview.satellite_constants import SATELLITE_FAILURE_RETRY_SECONDS
from zstarview.types import CelestialData, PlanetBody, UrbanOutlinePolyline, ViewerData
from zstarview.gui.famous_star_shortcuts import SearchJumpTarget
from zstarview.gui.window import SkyWindow
from zstarview.gui.window_state import SkyWindowState
from PySide6.QtWidgets import QApplication

_app = QApplication.instance() or QApplication([])


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
        self.observation_info_mode = values.get("observation_info_mode", "auto")
        self.observation_info_pinned = values.get("observation_info_pinned", False)
        self.show_observation_info = values.get("show_observation_info", True)
        self.show_dso = values.get("show_dso", False)
        self.show_asterisms = values.get("show_asterisms", False)
        self.show_guidelines = values.get("show_guidelines", True)
        self.enlarge_moon = values.get("enlarge_moon", False)
        self.star_base_radius = values.get("star_base_radius", 4.0)
        self.star_visibility_boost = values.get("star_visibility_boost", 1.0)
        self.vmag_limit = values.get("vmag_limit", 6.0)
        self.cloud_disc_alpha = values.get("cloud_disc_alpha", 0.0)
        self.satellite_opacity = values.get("satellite_opacity", 0.0)
        self.aircraft_opacity = values.get("aircraft_opacity", 0.0)
        self.terrain_horizon_opacity = values.get("terrain_horizon_opacity", 0.25)
        self.earth_guide_opacity = values.get("earth_guide_opacity", 0.25)
        self.urban_outline_opacity = values.get("urban_outline_opacity", 0.2)
        self.show_urban_outline_layer = values.get("show_urban_outline_layer", True)
        self.status_line_font = values.get("status_line_font", object())
        self._star_render_expected_width = values.get(
            "_star_render_expected_width", 600
        )
        self._enabled_satellite_groups = values.get(
            "_enabled_satellite_groups", ("iss",)
        )
        self._disc_generation = values.get("_disc_generation", 0)
        self._client_widget = values.get("_client_widget", None)
        self.cloud_state = values.get("cloud_state", None)
        self.satellite_state = values.get("satellite_state", None)
        self.viewer_data = values.get("viewer_data", None)
        self._cloud_controller = values.get("_cloud_controller", None)
        self.state = values.get("state", None)

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

    def _target_time_utc(self):
        target_time_utc = self.__dict__.get("_target_time_utc")
        if callable(target_time_utc):
            return target_time_utc()
        return datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

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

    def _extract_horizons_altaz(self, rows: list[list[str]]) -> tuple[float, float] | None:
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


def _make_scene(
    *,
    viewer: ViewerData | None = None,
    celestial_data: CelestialData | object | None = None,
    terrain_horizon_profile: object | None = None,
    terrain_horizon_profile_distances_m: object | None = None,
    terrain_horizon_secondary_profile_altaz_layers: object | None = None,
    terrain_horizon_secondary_profile_distances_m_layers: object | None = None,
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
        terrain_horizon_secondary_profile_altaz_layers=terrain_horizon_secondary_profile_altaz_layers,
        terrain_horizon_secondary_profile_distances_m_layers=terrain_horizon_secondary_profile_distances_m_layers,
        urban_outlines=urban_outlines,
        satellite_overlay_points=None,
        aircraft_overlay_points=None,
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
        "vmag_limit": 6.0,
        "cloud_disc_alpha": 0.0,
        "satellite_opacity": 0.0,
        "terrain_horizon_opacity": 0.25,
        "earth_guide_opacity": 0.25,
        "urban_outline_opacity": 0.2,
        "show_urban_outline_layer": True,
        "aircraft_opacity": 0.0,
        "star_render_expected_width": 600,
    }
    values.update(overrides)
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
    )
    dummy.state = SkyWindowState(render_view_center=(60.0, 210.0))

    got = SkyWindow._viewer_data_for_render(dummy)
    assert got.view_center == (60.0, 210.0)
    assert got.location == (35.0, 139.0)
    assert got.observer_height_m == 123.0
    assert got.location_height_label is None
    assert got.location_height_m is None
    assert got.show_observer_height is False


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
        preset="day",
    )

    assert calls == [("night", True)]


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
        "view_center": (15.0, 120.0),
        "render_width_px": 640,
        "render_height_px": 480,
        "render_generation": 0,
    }
    SkyWindow._on_sky_data_calculated(dummy, payload)

    assert dummy.state.render_view_center == (15.0, 120.0)
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
            "view_center": (15.0, 120.0),
            "render_width_px": 640,
            "render_height_px": 480,
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
            "view_center": (15.0, 120.0),
            "render_width_px": 640,
            "render_height_px": 480,
            "render_generation": 0,
        },
    )

    assert dummy.state.viewport_interaction_release_pending is False
    assert dummy.state.viewport_interaction_mode is False
    assert cloud_calls == ["view-change-release", "view-change-release"]


def test_schedule_satellite_retry_after_failure_uses_two_hour_backoff() -> None:
    dummy = _WindowStub()
    dummy.satellite_opacity = 0.5
    dummy._is_shutting_down = False
    timer = _DummyTimer(active=False)
    dummy._satellite_update_timer = timer
    dummy._satellite_layer_enabled = lambda: True

    SkyWindow._schedule_satellite_retry_after_failure(dummy)

    assert timer.started_with == [SATELLITE_FAILURE_RETRY_SECONDS * 1000]


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
        window_module, "find_satellite_altaz", lambda *args, **kwargs: (-12.0, 123.0)
    )

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0), satellite_overlay_points=None
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={"iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]},
        overlay_points=None,
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

    assert dummy.viewer_data.view_center == (-5.0, 123.0)
    assert dummy.state.jump_highlight_name == "ISS"
    assert dummy.state.jump_highlight_altaz == (-12.0, 123.0)
    dummy.request_sky_data_update.assert_called_once()


def test_find_satellite_jump_altaz_falls_back_to_disk_cache(monkeypatch) -> None:
    monkeypatch.setattr(
        window_module,
        "find_satellite_altaz",
        lambda records_by_group, **kwargs: (
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
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0), satellite_overlay_points=None
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
    dummy.update = Mock()
    dummy._load_cached_satellite_records = lambda groups: {}
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
        "project_place_target_to_altaz",
        lambda **kwargs: PlaceTargetProjection(
            alt_deg=-3.5,
            az_deg=145.0,
            distance_km=12.0,
            target_latitude_deg=float(kwargs["target_latitude_deg"]),
            target_longitude_deg=float(kwargs["target_longitude_deg"]),
            target_height_m=0.0,
        ),
    )

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0), satellite_overlay_points=None
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
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

    # Allow the view center to go below the horizon (>= -5.0°) per policy.
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
        satellite_overlay_points=None,
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda target, **kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda target, **kwargs: (12.5, 220.0),
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
    assert "JPL persistent target set: label=Ceres kind=jpl_small_body group=<none>" in caplog.text
    assert "target_time_utc=2026-04-18T12:00:00+00:00" in caplog.text
    assert "alt=12.5 az=220.0 command=DES=20000001;" in caplog.text


def test_jump_to_jpl_small_body_target_uses_state_vector_when_present(monkeypatch) -> None:
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
        satellite_overlay_points=None,
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()
    dummy._current_time_obj = lambda: astropy.time.Time("2026-04-18T12:00:00Z")

    monkeypatch.setattr(
        window_module,
        "resolve_jpl_target_state_vector",
        lambda target, **kwargs: (
            datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
            (1.0, 2.0, 3.0),
            (0.1, 0.2, 0.3),
        ),
    )
    monkeypatch.setattr(
        window_module,
        "project_jpl_target_altaz_from_state_vector",
        lambda target, **kwargs: (11.5, 221.0),
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
    assert dummy.state.persistent_search_target.horizons_velocity_km_s == (0.1, 0.2, 0.3)


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
        satellite_overlay_points=None,
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda target, **kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda target, **kwargs: (12.5, 220.0),
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
        satellite_overlay_points=None,
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda target, **kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda target, **kwargs: (12.5, 220.0),
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
        satellite_overlay_points=None,
        persistent_search_target=SearchJumpTarget(
            label="Old",
            kind="jpl_small_body",
            sort_key=(0.0, "old"),
            alt_deg=10.0,
            az_deg=30.0,
            persistent_keep_marker=True,
        ),
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda target, **kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda target, **kwargs: (12.5, 220.0),
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
        satellite_overlay_points=None,
        persistent_search_target=current_target,
        persistent_search_next_refresh_utc=datetime(2026, 4, 18, 13, 0, tzinfo=timezone.utc),
        persistent_search_reference_time_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
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
        satellite_overlay_points=None,
        persistent_search_target=current_target,
        persistent_search_next_refresh_utc=datetime(2026, 4, 18, 13, 0, tzinfo=timezone.utc),
        persistent_search_reference_time_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
    )
    dummy.request_client_update = Mock()
    dummy._schedule_persistent_search_refresh = Mock()
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda target, **kwargs: (49.1, 244.7),
        )

        with caplog.at_level("INFO"):
            SkyWindow._on_jpl_ready(
                dummy,
                {
                    "target": current_target,
                    "target_time_utc": datetime(2026, 4, 18, 13, 0, tzinfo=timezone.utc),
                    "refreshed_at_utc": datetime(2026, 4, 18, 13, 2, tzinfo=timezone.utc),
                    "horizons_epoch_utc": datetime(2026, 4, 18, 13, 0, tzinfo=timezone.utc),
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
    assert dummy.state.persistent_search_target.horizons_velocity_km_s == (0.4, 0.5, 0.6)
    assert "JPL persistent target refreshed: label=Voyager 1 kind=jpl_small_body group=sb" in caplog.text
    assert "target_time_utc=2026-04-18T13:00:00+00:00" in caplog.text
    assert "alt=49.1 az=244.7 command=DES=-31;" in caplog.text


def test_refresh_projected_persistent_search_target_reprojects_state_vector(monkeypatch) -> None:
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
        lambda target, **kwargs: (12.5, 220.0),
    )

    window_updates_module.SkyWindowUpdatesMixin.refresh_projected_persistent_search_target(
        dummy
    )

    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.alt_deg == 12.5
    assert dummy.state.persistent_search_target.az_deg == 220.0
    assert dummy.request_client_update.called


def test_search_satellite_targets_resolves_known_artificial_satellites(monkeypatch) -> None:
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

    def fake_lookup(query: str, *, group: str):
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

    def fake_lookup(query: str, *, group: str):
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

    assert lookup_calls == ["mb", "sb"]
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
    dummy._begin_viewport_interaction_mode = (
        lambda *args, **kwargs: SkyWindow._begin_viewport_interaction_mode(
            dummy, *args, **kwargs
        )
    )
    dummy._update_viewport_interaction_stars = (
        lambda: SkyWindow._update_viewport_interaction_stars(dummy)
    )
    dummy._viewport_rotation_keys = lambda: SkyWindow._viewport_rotation_keys(dummy)

    event = SimpleNamespace(
        key=lambda: Qt.Key.Key_Left,
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


def test_handle_client_key_release_ends_viewport_interaction_mode() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy._viewport_rotation_keys_down = {Qt.Key.Key_Left}
    dummy._viewport_rotation_keys = lambda: SkyWindow._viewport_rotation_keys(dummy)
    dummy._end_viewport_interaction_mode = (
        lambda *args, **kwargs: SkyWindow._end_viewport_interaction_mode(
            dummy, *args, **kwargs
        )
    )
    dummy.request_sky_data_update = Mock()
    dummy.start_background_cloud_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
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
        reason="viewport-interaction-release"
    )
    dummy.start_background_cloud_update.assert_not_called()
    dummy.start_background_terrain_horizon_update.assert_not_called()
    dummy.request_client_update.assert_not_called()
    event.accept.assert_called_once()


def test_end_viewport_interaction_mode_marks_idle_reason() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        viewport_interaction_mode=True,
    )
    dummy.request_sky_data_update = Mock()
    dummy.start_background_cloud_update = Mock()
    dummy.start_background_terrain_horizon_update = Mock()
    dummy.request_client_update = Mock()

    SkyWindow._end_viewport_interaction_mode(dummy)

    dummy.request_sky_data_update.assert_called_once_with(
        reason="viewport-interaction-idle"
    )
    dummy.start_background_cloud_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.start_background_terrain_horizon_update.assert_called_once_with(
        reason="view-change-idle"
    )
    dummy.request_client_update.assert_called_once()


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
        satellite_overlay_points=None,
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={}, overlay_points=None, set_banner=Mock()
    )
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda target, **kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda target, **kwargs: (15.0, 123.0),
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


def test_draw_persistent_search_overlay_draws_label_when_marker_is_kept(
    monkeypatch,
) -> None:
    draw_calls: list[str] = []

    class _Painter:
        def viewport(self):
            return QRect(0, 0, 200, 200)

    monkeypatch.setattr(
        search_overlay_module.render_guides,
        "draw_gauge_cross",
        lambda *_args, **_kwargs: draw_calls.append("marker"),
    )
    monkeypatch.setattr(
        search_overlay_module.render_text,
        "draw_outlined_text",
        lambda *_args, **_kwargs: draw_calls.append("label"),
    )
    monkeypatch.setattr(search_overlay_module, "is_in_fov", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(
        search_overlay_module,
        "altaz_to_normalized_xy",
        lambda alt, az, view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        search_overlay_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    dummy = _WindowStub()
    dummy.visual_preset = "white"
    dummy.text_font = QFont()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
        content_fov_deg=180.0,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=SearchJumpTarget(
            label="Ceres",
            kind="jpl_small_body",
            sort_key=(0.0, "ceres"),
            alt_deg=12.5,
            az_deg=220.0,
            persistent_keep_marker=True,
        ),
    )

    window_render_module.SkyWindowRenderMixin._draw_persistent_search_overlay(
        dummy,
        _Painter(),
        geometry=SimpleNamespace(center=(100, 100), radius=80),
    )

    assert draw_calls == ["marker", "label"]


def test_draw_persistent_search_overlay_scales_marker_with_window_scale(
    monkeypatch,
) -> None:
    marker_scales: list[float] = []

    class _Painter:
        def viewport(self):
            return QRect(0, 0, 200, 200)

    monkeypatch.setattr(
        window_render_module,
        "compute_star_render_upscale_factor",
        lambda *_args, **_kwargs: 2.5,
    )
    monkeypatch.setattr(
        search_overlay_module.render_guides,
        "draw_gauge_cross",
        lambda *_args, **kwargs: marker_scales.append(
            float(kwargs.get("scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        search_overlay_module.render_text,
        "draw_outlined_text",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(search_overlay_module, "is_in_fov", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(
        search_overlay_module,
        "altaz_to_normalized_xy",
        lambda alt, az, view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        search_overlay_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    dummy = _WindowStub()
    dummy.visual_preset = "white"
    dummy.text_font = QFont()
    dummy._star_render_expected_width = 600
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
        content_fov_deg=180.0,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=SearchJumpTarget(
            label="Ceres",
            kind="jpl_small_body",
            sort_key=(0.0, "ceres"),
            alt_deg=12.5,
            az_deg=220.0,
            persistent_keep_marker=True,
        ),
    )

    window_render_module.SkyWindowRenderMixin._draw_persistent_search_overlay(
        dummy,
        _Painter(),
        geometry=SimpleNamespace(center=(100, 100), radius=80),
    )

    assert marker_scales == [1.05]


def test_refresh_projected_satellite_overlay_falls_back_to_disk_cache(
    monkeypatch,
) -> None:
    projected_points = [
        SimpleNamespace(satellite_name="ISS (ZARYA)", alt_deg=-40.0, az_deg=151.0)
    ]
    monkeypatch.setattr(
        window_updates_module,
        "project_satellite_records",
        lambda *args, **kwargs: projected_points,
    )

    dummy = _WindowStub()
    dummy.satellite_opacity = 1.0
    dummy.viewer_data = ViewerData(
        location=(40.7128, -74.0060),
        timezone_name="America/New_York",
        city_name="New York City",
        view_center=(0.0, 151.0),
        observer_height_m=10.0,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, overlay_points=None)
    dummy.state = SkyWindowState(
        render_view_center=(0.0, 151.0), satellite_overlay_points=None
    )
    dummy._enabled_satellite_groups = ("iss",)
    dummy._satellite_validity_remaining_ms = lambda: 1000
    dummy._load_cached_satellite_records = lambda groups: {
        "iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]
    }
    dummy._current_time_obj = lambda: astropy.time.Time("2026-03-23T12:13:24Z")
    dummy.update = Mock()

    SkyWindow.refresh_projected_satellite_overlay(dummy)

    assert dummy.satellite_state.overlay_points == projected_points
    assert dummy.state.satellite_overlay_points == projected_points
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
            "render_width_px": 640,
            "render_height_px": 480,
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
        "draw_terrain_horizon_line",
        lambda *_args, **_kwargs: calls.append(("terrain", None)),
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
        ("terrain", None),
    ]


def test_draw_cached_frame_reuses_existing_image() -> None:
    draws: list[tuple[int, int]] = []
    render_calls: list[str] = []

    class _Painter:
        def drawImage(self, x: int, y: int, image: QImage) -> None:  # noqa: N802 - Qt naming
            draws.append((x, y))
            assert not image.isNull()

    dummy = _WindowStub()
    dummy._frame_cache_key = None
    dummy._frame_cache_image = None
    dummy.size = lambda: window_render_module.QImage(
        32, 24, QImage.Format.Format_ARGB32_Premultiplied
    ).size()

    def render_fn(frame_painter) -> None:
        render_calls.append("render")
        frame_painter.fillRect(0, 0, 32, 24, window_render_module.Qt.GlobalColor.black)

    painter = _Painter()

    window_render_module.SkyWindowRenderMixin._draw_cached_frame(
        dummy, painter, ("same",), render_fn
    )
    window_render_module.SkyWindowRenderMixin._draw_cached_frame(
        dummy, painter, ("same",), render_fn
    )

    assert render_calls == ["render"]
    assert draws == [(0, 0), (0, 0)]


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


def test_render_frame_cache_key_ignores_hover_and_status_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = object()
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
    dummy.state.aircraft_overlay_points = [object()]
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


def test_present_frame_cache_key_tracks_hover_and_status_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = object()
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
    dummy.state.satellite_overlay_points = [object()]
    dummy.state.aircraft_overlay_points = [object()]

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
    celestial_data = object()
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
    dummy.state.satellite_overlay_points = [object()]
    dummy.state.aircraft_overlay_points = [object()]

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    dummy.state.satellite_overlay_points = [object(), object()]
    dummy.state.aircraft_overlay_points = [object(), object()]
    dummy.satellite_opacity = 0.9
    dummy.aircraft_opacity = 0.8

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    assert key_a == key_b


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
    celestial_data = object()
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

    highlighted_object, highlighted_dso, highlighted_satellite = (
        window_render_module._resolve_hover_targets(
            celestial_data=celestial_data,
            render_viewer=viewer,
            mouse_pos=mouse_pos,
            geometry=geometry,
            satellite_overlay_points=[object()],
            show_dso=True,
        )
    )

    assert highlighted_object == star_hit
    assert highlighted_dso == {"name": "DSO"}
    assert highlighted_satellite == satellite_hit


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
        "draw_terrain_horizon_line",
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
    seen_profiles: list[object] = []
    seen_view_centers: list[object] = []
    seen_line_width_scales: list[float] = []
    secondary_calls: list[object] = []
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
        "draw_terrain_horizon_line",
        lambda _p, _g, profile, distances, view_center, **kwargs: (
            seen_profiles.append(profile),
            seen_view_centers.append(view_center),
            seen_line_width_scales.append(float(kwargs.get("line_width_scale", 1.0))),
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: secondary_calls.append("called"),
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

    assert seen_profiles == [terrain_profile]
    assert seen_view_centers == [(50.0, 210.0)]
    assert seen_line_width_scales == [expected_line_width_scale]
    assert secondary_calls == []


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
        "draw_terrain_horizon_line",
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


def test_draw_terrain_layers_keeps_asterisms_and_urban_outline_widths_fixed(monkeypatch) -> None:
    calls: dict[str, list[float]] = {
        "asterisms": [],
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
        lambda *_args, **_kwargs: None,
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
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_horizon_line",
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
            terrain_horizon_secondary_profile_altaz_layers=[
                [(1.0, 10.0), (2.0, 20.0)]
            ],
            terrain_horizon_secondary_profile_distances_m_layers=[
                [10_000.0, 12_000.0]
            ],
        ),
        style=_make_style(show_asterisms=True),
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls["asterisms"] == [1.0]
    assert calls["terrain"] == []
    assert calls["terrain_secondary"] == [expected_line_width_scale]
    assert calls["reference"] == [1.0]
    assert calls["direction"] == []
    assert calls["zenith"] == []
    assert calls["urban"] == [1.0]


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
        "draw_terrain_horizon_line",
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
        "_draw_label_layer",
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
        terrain_horizon_secondary_profile_altaz_layers=None,
        terrain_horizon_secondary_profile_distances_m_layers=None,
        urban_outlines=None,
        satellite_overlay_points=None,
        aircraft_overlay_points=None,
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
        vmag_limit=6.0,
        cloud_disc_alpha=0.0,
        satellite_opacity=0.0,
        terrain_horizon_opacity=0.0,
        earth_guide_opacity=0.0,
        urban_outline_opacity=0.0,
        show_urban_outline_layer=False,
        aircraft_opacity=0.0,
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
        geometry=geometry,
        viewport_rect=viewport_rect,
        scene=scene,
        style=style,
        hud=hud,
        compositor=object(),
    )
    pipeline_module.render_hud_overlay_into_painter(
        painter=object(),
        geometry=geometry,
        viewport_rect=viewport_rect,
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


def test_render_scene_hides_cloud_bitmap_during_viewport_interaction(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

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
        geometry=SimpleNamespace(center=(100, 100), radius=80),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=replace(_make_scene(), sky_disc_image=object()),
        style=_make_style(cloud_disc_alpha=0.2),
        hud=_make_hud(viewport_interaction_mode=True),
        compositor=object(),
    )

    assert captured == {"cloud_disc_alpha": 0.0, "sky_disc_image": None}


def test_render_base_scene_can_skip_fast_overlays(monkeypatch) -> None:
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
        "_draw_label_layer",
        lambda *_args, **_kwargs: calls.append("labels"),
    )

    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        geometry=SimpleNamespace(center=(100, 100), radius=80),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(),
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
                PlanetBody(name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True),
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
        scene=_make_scene(celestial_data=object()),
        style=_make_style(),
        highlighted_object=({"name": "moon"}, object()),
        highlighted_dso=None,
    )

    assert calls == ["moon-hover", "dso-hover", "overlay-info"]


def test_draw_sky_reference_lines_uses_render_view_center_projection(
    monkeypatch,
) -> None:
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
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, view_center, **_kwargs: (
            calls.append(view_center) or (alt, az)
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
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, view_center, **_kwargs: (alt, az),
    )

    assert dash_patterns[0::3] == [[], [], []]
    assert dash_patterns[1::3] == [[], [], []]
    assert dash_patterns[2::3] == [[12, 6], [4, 6], [10, 1]]
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
    assert [round(width, 3) for width in pen_widths[0::3]] == [1.254, 1.254, 1.100]
    assert [round(width, 3) for width in pen_widths[1::3]] == [0.855, 0.855, 0.750]
    assert [round(width, 3) for width in pen_widths[2::3]] == [0.627, 0.627, 0.550]


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
        lambda alt, az, view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_guides_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    render_guides_module.draw_direction_labels(
        _Painter(),
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        view_center=(30.0, 40.0),
        text_font=QFont(),
        theme=window_module.THEME_STYLES_BY_PRESET["white"],
        content_fov_deg=180.0,
    )

    assert seen_colors
    assert all(
        color == render_guides_module.HORIZON_LINE_COLOR for color in seen_colors
    )


def test_draw_urban_outlines_clips_two_point_outline_out_of_view(
    monkeypatch,
) -> None:
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
    render_terrain_module.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(points=[(-10.0, 10.0), (-12.0, 10.3)], height_m=50.0)
        ],
        view_center=(45.0, 180.0),
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
            self.alpha_values.append(int(pen.color().alpha()))
            self.width_values.append(float(pen.widthF()))

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
        view_center=(45.0, 180.0),
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
            self.width_values.append(float(pen.widthF()))

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
        view_center=(45.0, 180.0),
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
            self.width_values.append(float(pen.widthF()))

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
        view_center=(45.0, 180.0),
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
    render_terrain_module.draw_terrain_horizon_line(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        terrain_profile_altaz=[(0.0, 0.0), (0.1, 0.1)],
        terrain_profile_distances_m=None,
        view_center=(45.0, 180.0),
        opacity=0.38,
        line_width_scale=2.0,
        fast_mode=True,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert painter.pen_widths == [7.2]


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
    render_terrain_module.draw_terrain_horizon_line(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        terrain_profile_altaz=[(0.0, 0.0), (0.0, 0.1), (0.0, 0.2)],
        terrain_profile_distances_m=[1_000.0, 50_000.0, 120_000.0],
        view_center=(45.0, 180.0),
        opacity=0.38,
        line_width_scale=1.0,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
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
        terrain_secondary_profile_layers=[
            [(0.0, 0.0), (0.0, 0.1), (0.0, 0.2)],
            [(0.1, 0.0), (0.1, 0.1), (0.1, 0.2)],
            [(0.2, 0.0), (0.2, 0.1), (0.2, 0.2)],
        ],
        terrain_secondary_profile_distances_m_layers=[
            [1_000.0, 2_000.0, 3_000.0],
            [10_000.0, 12_000.0, 15_000.0],
            [50_000.0, 60_000.0, 70_000.0],
        ],
        view_center=(45.0, 180.0),
        opacity=0.38,
        line_width_scale=1.0,
        fast_mode=False,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.pen_widths) == 24
    for offset in (0, 8, 16):
        chunk = painter.pen_widths[offset : offset + 8]
        assert chunk[0] > chunk[1]
        assert chunk[2] > chunk[3] > chunk[4]
        assert chunk[5] > chunk[6] > chunk[7]
        assert chunk[4] == pytest.approx(chunk[7])
    assert painter.pen_widths[1] > painter.pen_widths[9] > painter.pen_widths[17]
    assert painter.pen_widths[0] > painter.pen_widths[8] > painter.pen_widths[16]


def test_draw_terrain_secondary_ridges_swaps_visible_and_occluded_colors(monkeypatch) -> None:
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
        terrain_secondary_profile_layers=[
            [(20.0, 0.0), (20.0, 0.1)],
            [(10.0, 0.0), (10.0, 0.1)],
        ],
        terrain_secondary_profile_distances_m_layers=[
            [1_000.0, 2_000.0],
            [10_000.0, 12_000.0],
        ],
        view_center=(45.0, 180.0),
        opacity=0.38,
        line_width_scale=1.0,
        fast_mode=False,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.pen_rgbs) == 7
    assert painter.pen_rgbs[0] == render_terrain_module.TERRAIN_SECONDARY_RIDGE_OCCLUDED_COLOR_RGB
    assert painter.pen_rgbs[1] == render_terrain_module.TERRAIN_SECONDARY_RIDGE_OCCLUDED_COLOR_RGB
    assert painter.pen_rgbs[2] == render_terrain_module.TERRAIN_SECONDARY_RIDGE_VISIBLE_COLOR_RGB
    assert painter.pen_rgbs[3] == render_terrain_module.TERRAIN_SECONDARY_RIDGE_VISIBLE_COLOR_RGB
    assert painter.pen_rgbs[4] == render_terrain_module.TERRAIN_SECONDARY_RIDGE_VISIBLE_COLOR_RGB
    assert painter.pen_rgbs[5] == render_terrain_module.TERRAIN_SECONDARY_RIDGE_OCCLUDED_COLOR_RGB
    assert painter.pen_rgbs[6] == render_terrain_module.TERRAIN_SECONDARY_RIDGE_OCCLUDED_COLOR_RGB


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
        terrain_secondary_profile_layers=[
            [(20.0, 0.0), (20.0, 0.1)],
            [(10.0, 0.0), (10.0, 0.1)],
        ],
        terrain_secondary_profile_distances_m_layers=[
            [1_000.0, 2_000.0],
            [10_000.0, 12_000.0],
        ],
        view_center=(45.0, 180.0),
        opacity=0.38,
        line_width_scale=1.0,
        fast_mode=False,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            float(az),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    assert len(painter.alphas) == 7
    assert painter.alphas[2] < painter.alphas[3] < painter.alphas[4]
    assert painter.alphas[4] < painter.alphas[1] * 0.4


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
            self.lines.append(((float(start.x()), float(start.y())), (float(end.x()), float(end.y()))))

        def drawPolyline(self, *_args) -> None:
            pass

    painter = _Painter()
    render_terrain_module.draw_terrain_secondary_ridges(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        terrain_secondary_profile_layers=[
            [(5.0, 359.0), (5.0, 0.0), (5.0, 1.0)],
        ],
        terrain_secondary_profile_distances_m_layers=[
            [1_000.0, 2_000.0, 3_000.0],
        ],
        view_center=(45.0, 180.0),
        opacity=0.38,
        line_width_scale=1.0,
        fast_mode=False,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            math.sin(math.radians(float(az))),
            float(alt),
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
        split_by_gaps_func=lambda points: [points],
    )

    assert len(painter.lines) == 9
    expected_lines = [
        ((0.0, 5.0), (math.sin(math.radians(1.0)), 5.0)),
        ((math.sin(math.radians(1.0)), 5.0), (math.sin(math.radians(359.0)), 5.0)),
        ((0.0, 5.0), (math.sin(math.radians(359.0)), 5.0)),
    ]
    for expected_start, expected_end in expected_lines:
        assert painter.lines.count((expected_start, expected_end)) == 3


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
            self.polylines.append([(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())])

    def _render(edge_fov_deg: float) -> list[tuple[float, float]]:
        painter = _Painter()
        render_terrain_module.draw_terrain_horizon_line(
            painter,
            geometry=SimpleNamespace(center=(0, 0), radius=1),
            terrain_profile_altaz=[(0.0, 180.0), (0.0, 190.0)],
            terrain_profile_distances_m=None,
            view_center=(45.0, 180.0),
            opacity=0.38,
            edge_fov_deg=edge_fov_deg,
            content_fov_deg=180.0,
            fast_mode=True,
            is_in_fov_func=lambda *_args, **_kwargs: True,
            normalized_to_screen_xy_func=lambda nx, ny, _geometry: (float(nx), float(ny)),
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
            self.polylines.append([(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())])

    painter = _Painter()
    render_terrain_module.draw_terrain_horizon_line(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        terrain_profile_altaz=[
            (0.0, 0.0),
            (0.0, 10.0),
            (0.0, 180.0),
            (0.0, 190.0),
        ],
        terrain_profile_distances_m=None,
        view_center=(45.0, 180.0),
        opacity=0.38,
        edge_fov_deg=120.0,
        content_fov_deg=180.0,
        fast_mode=True,
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

    monkeypatch.setattr(pipeline_module.render_stars, "draw_stars", fake_draw_stars)

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
    assert captured["fast_mode"] is True
