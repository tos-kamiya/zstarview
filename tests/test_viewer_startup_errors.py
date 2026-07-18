from __future__ import annotations

import logging
from types import SimpleNamespace
from typing import Callable

import pytest

import zstarview.gui.viewer as viewer
from zstarview.location_resolver import LocationResolveError


class _DummyApp:
    def __init__(self) -> None:
        self.quit_on_last: list[bool] = []
        self.exit_codes: list[int] = []
        self.terrain_result: tuple[str, object] | None = None

    def setQuitOnLastWindowClosed(self, value: bool) -> None:
        self.quit_on_last.append(bool(value))

    def processEvents(self) -> None:
        return None

    def exit(self, code: int = 0) -> None:
        self.exit_codes.append(int(code))

    def exec(self) -> int:
        if self.terrain_result is not None and _DummyWindow.last_instance is not None:
            signal_name, payload = self.terrain_result
            signal = getattr(
                _DummyWindow.last_instance._terrain_horizon_controller, signal_name
            )
            signal.emit(payload)
            self.terrain_result = None
        return self.exit_codes[-1] if self.exit_codes else 0


class _DummyRootLogger:
    def __init__(self) -> None:
        self.added: list[object] = []
        self.removed: list[object] = []
        self._real_root = logging.getLogger()

    def addHandler(self, handler: object) -> None:
        self.added.append(handler)
        self._real_root.addHandler(handler)

    def removeHandler(self, handler: object) -> None:
        self.removed.append(handler)
        self._real_root.removeHandler(handler)


class _DummyOverlay:
    def __init__(self) -> None:
        self.lines: list[tuple[str, int]] = []
        self.shown = 0
        self.hidden = 0

    def append_line(self, line: str, levelno: int) -> None:
        self.lines.append((line, int(levelno)))

    def show(self) -> None:
        self.shown += 1

    def hide(self) -> None:
        self.hidden += 1

    def raise_(self) -> None:
        return None


class _DummySignal:
    def __init__(self) -> None:
        self._callbacks: list[Callable[..., None]] = []

    def connect(self, fn: Callable[..., None]) -> None:
        self._callbacks.append(fn)

    def emit(self, *args, **kwargs) -> None:
        for fn in list(self._callbacks):
            fn(*args, **kwargs)


class _DummyTerrainController:
    def __init__(self) -> None:
        self.terrain_started = _DummySignal()
        self.terrain_ready = _DummySignal()
        self.terrain_failed = _DummySignal()


class _DummyWindow:
    last_instance = None

    def __init__(self, *args, **kwargs) -> None:
        type(self).last_instance = self
        user_options = kwargs.get("user_options")
        if user_options is None and len(args) >= 3:
            user_options = args[2]
        self.overlay = _DummyOverlay()
        self.shown = 0
        self.overlay_raise_calls = 0
        self.applied_delta = None
        self.applied_viewer_data = None
        self.terrain_horizon_opacity = float(
            getattr(user_options, "terrain_horizon_opacity", 0.0)
        )
        self.post_initial_hidden_states: list[int] = []
        self.initial_data_loaded = _DummySignal()
        self._terrain_horizon_controller = _DummyTerrainController()

    def _ensure_startup_log_overlay(self) -> _DummyOverlay:
        return self.overlay

    def show(self) -> None:
        self.shown += 1

    def apply_startup_delta_t(self, delta_t) -> None:
        self.applied_delta = delta_t

    def apply_startup_viewer_data(self, viewer_data) -> None:
        self.applied_viewer_data = viewer_data

    def start_initial_data_load(self) -> None:
        self.initial_data_loaded.emit()
        self.post_initial_hidden_states.append(self.overlay.hidden)

    def _hide_startup_log_overlay(self) -> None:
        self.overlay.hide()

    def _raise_overlay_widgets(self) -> None:
        self.overlay_raise_calls += 1

    def _jump_to_search_target(self, _target) -> None:
        return None


def _make_args(*, close_on_startup_error: bool) -> SimpleNamespace:
    return SimpleNamespace(
        city="Singapole",
        place=None,
        place_countrycode=None,
        place_lang="en",
        theme="night",
        clear_long_lived_cache=False,
        close_on_startup_error=close_on_startup_error,
        vmag_limit=None,
        vmag_brightness_multiplier=1.0,
        cloud_stripe=("width", 50, 0.85),
        edge_fov_deg=95.0,
        content_fov_deg=110.0,
        view_center_alt=90.0,
        view_center_az=180.0,
        sky_opacity=0.1,
        sky_disc_style="smooth",
        sky_disc_altaz_rings="dimalt",
        sky_disc_altaz_rings_hover="altaz",
        night_light_opacity=0.07,
        ridge_glow_opacity=0.36,
        cloud_opacity=0.075,
        geo_satellite=False,
        satellite_opacity=0.7,
        aircraft_opacity=0.4,
        tropical_cyclone_opacity=0.4,
        terrain_horizon_opacity=0.003,
        earth_guide_opacity=0.028,
        urban_outline_opacity=0.2,
        water_surface_opacity=0.4,
        ground_tint_opacity=0.025,
        overlay_font_size=11,
        enlarge_moon=False,
        bright_bodies="outline",
        star_base_radius=4.0,
        visibility_boost=1.0,
        show_dso_initial=None,
        show_asterisms_initial=None,
        show_guidelines_initial=None,
        observation_info="auto",
        sky_update_interval=60,
        urban_outline_radius_km=2.5,
        urban_outline_skyscraper_radius_km=60.0,
        urban_outline_min_height_m=0.0,
        urban_outline_feature_type="both",
        urban_outline_skyscraper_only=False,
        urban_outline_max_candidates=5000,
        urban_outline_download_timeout_seconds=120.0,
        cloud_missing_tint_opacity=0.2,
        expected_render_width=600,
        window_geometry=None,
        window_frame="frameless",
        datetime=None,
        days=0,
        hours=0,
        timezone=None,
        observer_height_m=None,
        use_building_top=False,
        search="",
        search_keep_marker=False,
        view_center_alt_specified=False,
        view_center_az_specified=False,
    )


def _install_common_mocks(
    monkeypatch,
    args: SimpleNamespace,
    *,
    resolve_location_raises: bool = True,
):
    app = _DummyApp()
    root_logger = _DummyRootLogger()
    monkeypatch.setattr(viewer, "parse_args", lambda: args)
    monkeypatch.setattr(viewer, "_handle_dataset_query_cli", lambda _args: None)
    monkeypatch.setattr(viewer, "setup_root_logger", lambda: root_logger)
    monkeypatch.setattr("zstarview.splash.setup_app", lambda _app_name: app)
    monkeypatch.setattr(viewer.QTimer, "singleShot", lambda _msec, callback: callback())
    monkeypatch.setattr(viewer._StartupBootstrap, "start", lambda self: self._run())
    monkeypatch.setattr(viewer, "_load_star_catalog_for_launch", lambda _vmag_limit: object())
    monkeypatch.setattr(viewer, "_load_dso_catalog_for_launch", lambda: None)
    monkeypatch.setattr(viewer, "prepare_window_catalogs", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(viewer, "prepare_window_user_options", lambda **kwargs: SimpleNamespace(**kwargs))
    monkeypatch.setattr(viewer, "prepare_window_runtime_options", lambda **_kwargs: object())
    if resolve_location_raises:
        def _raise_location_error(*_args, **_kwargs):
            raise LocationResolveError()

        monkeypatch.setattr(viewer, "resolve_launch_location", _raise_location_error)
    else:
        monkeypatch.setattr(
            viewer,
            "resolve_launch_location",
            lambda *_args, **_kwargs: SimpleNamespace(
                display_name="Dummy City",
                lat=35.0,
                lon=139.0,
                tz="UTC",
                observer_height_m=3.0,
                ground_elevation_m=42.0,
                location_height_label=None,
                location_height_m=42.0,
            ),
        )
    monkeypatch.setattr("zstarview.gui.window.SkyWindow", _DummyWindow)
    monkeypatch.setattr("zstarview.gui.window.StandardSkyWindow", _DummyWindow)
    return app, root_logger


def test_main_keeps_startup_error_visible_by_default(monkeypatch) -> None:
    args = _make_args(close_on_startup_error=False)
    app, root_logger = _install_common_mocks(monkeypatch, args)

    with pytest.raises(SystemExit) as exc_info:
        viewer.main()

    assert exc_info.value.code == 0
    assert app.quit_on_last == [True]
    assert app.exit_codes == []
    assert len(root_logger.added) == 1
    assert root_logger.removed == root_logger.added
    assert _DummyWindow.last_instance is not None
    assert _DummyWindow.last_instance.overlay.shown == 1
    assert _DummyWindow.last_instance.overlay.hidden == 0
    assert _DummyWindow.last_instance.overlay_raise_calls >= 1
    assert any("Startup failed" in line for line, _level in _DummyWindow.last_instance.overlay.lines)
    assert any(level >= logging.ERROR for _line, level in _DummyWindow.last_instance.overlay.lines)


def test_main_auto_closes_on_startup_error_when_requested(monkeypatch) -> None:
    args = _make_args(close_on_startup_error=True)
    app, root_logger = _install_common_mocks(monkeypatch, args)
    single_shot_calls: list[int] = []

    def _single_shot(_msec: int, callback) -> None:
        single_shot_calls.append(int(_msec))
        callback()

    monkeypatch.setattr(viewer.QTimer, "singleShot", _single_shot)

    with pytest.raises(SystemExit) as exc_info:
        viewer.main()

    assert exc_info.value.code == 1
    assert app.quit_on_last == [True]
    assert app.exit_codes == [1]
    assert single_shot_calls == [0, 0]
    assert len(root_logger.added) == 1
    assert root_logger.removed == root_logger.added
    assert _DummyWindow.last_instance is not None
    assert _DummyWindow.last_instance.overlay.shown == 1
    assert _DummyWindow.last_instance.overlay.hidden == 0
    assert _DummyWindow.last_instance.overlay_raise_calls >= 1
    assert any("Startup failed" in line for line, _level in _DummyWindow.last_instance.overlay.lines)


def test_main_keeps_overlay_visible_until_terrain_resolves(monkeypatch) -> None:
    args = _make_args(close_on_startup_error=False)
    args.terrain_horizon_opacity = 0.25
    app, root_logger = _install_common_mocks(
        monkeypatch,
        args,
        resolve_location_raises=False,
    )
    app.terrain_result = ("terrain_ready", {"ground_elevation_m": 12.0})

    with pytest.raises(SystemExit) as exc_info:
        viewer.main()

    assert exc_info.value.code == 0
    assert app.quit_on_last == [True]
    assert len(root_logger.added) == 1
    assert root_logger.removed == root_logger.added
    assert _DummyWindow.last_instance is not None
    assert _DummyWindow.last_instance.overlay.shown == 1
    assert _DummyWindow.last_instance.post_initial_hidden_states == [0]
    assert _DummyWindow.last_instance.overlay.hidden == 1


def test_main_hides_overlay_without_terrain_when_not_requested(monkeypatch) -> None:
    args = _make_args(close_on_startup_error=False)
    args.terrain_horizon_opacity = 0.0
    app, root_logger = _install_common_mocks(
        monkeypatch,
        args,
        resolve_location_raises=False,
    )

    with pytest.raises(SystemExit) as exc_info:
        viewer.main()

    assert exc_info.value.code == 0
    assert app.quit_on_last == [True]
    assert len(root_logger.added) == 1
    assert root_logger.removed == root_logger.added
    assert _DummyWindow.last_instance is not None
    assert _DummyWindow.last_instance.overlay.shown == 1
    assert _DummyWindow.last_instance.post_initial_hidden_states == [1]
    assert _DummyWindow.last_instance.overlay.hidden == 1


def test_main_propagates_cli_view_center_flags_to_window(monkeypatch) -> None:
    args = _make_args(close_on_startup_error=False)
    args.view_center_alt_specified = True
    args.view_center_az_specified = False
    app, root_logger = _install_common_mocks(
        monkeypatch,
        args,
        resolve_location_raises=False,
    )

    with pytest.raises(SystemExit) as exc_info:
        viewer.main()

    assert exc_info.value.code == 0
    assert app.quit_on_last == [True]
    assert len(root_logger.added) == 1
    assert root_logger.removed == root_logger.added
    assert _DummyWindow.last_instance is not None
    assert _DummyWindow.last_instance._search_view_center_alt_specified is True
    assert _DummyWindow.last_instance._search_view_center_az_specified is False
