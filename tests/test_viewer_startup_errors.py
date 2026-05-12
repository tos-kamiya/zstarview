from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

import zstarview.gui.viewer as viewer
from zstarview.location_resolver import LocationResolveError


class _DummyApp:
    def __init__(self) -> None:
        self.quit_on_last: list[bool] = []
        self.exit_codes: list[int] = []

    def setQuitOnLastWindowClosed(self, value: bool) -> None:
        self.quit_on_last.append(bool(value))

    def processEvents(self) -> None:
        return None

    def exit(self, code: int = 0) -> None:
        self.exit_codes.append(int(code))

    def exec(self) -> int:
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


class _DummyWindow:
    last_instance = None

    def __init__(self, *args, **kwargs) -> None:
        type(self).last_instance = self
        self.overlay = _DummyOverlay()
        self.shown = 0
        self.applied_delta = None
        self.applied_viewer_data = None
        self.initial_data_loaded = SimpleNamespace(connect=lambda _fn: None)

    def _ensure_startup_log_overlay(self) -> _DummyOverlay:
        return self.overlay

    def show(self) -> None:
        self.shown += 1

    def apply_startup_delta_t(self, delta_t) -> None:
        self.applied_delta = delta_t

    def apply_startup_viewer_data(self, viewer_data) -> None:
        self.applied_viewer_data = viewer_data

    def start_initial_data_load(self) -> None:
        return None

    def _hide_startup_log_overlay(self) -> None:
        self.overlay.hide()

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
        night_light_opacity=0.02,
        cloud_opacity=0.075,
        satellite_opacity=0.7,
        aircraft_opacity=0.4,
        terrain_horizon_opacity=0.003,
        earth_guide_opacity=0.028,
        urban_outline_opacity=0.2,
        ground_tint_opacity=0.04,
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


def _install_common_mocks(monkeypatch, args: SimpleNamespace):
    app = _DummyApp()
    root_logger = _DummyRootLogger()
    monkeypatch.setattr(viewer, "parse_args", lambda: args)
    monkeypatch.setattr(viewer, "_handle_dataset_query_cli", lambda _args: None)
    monkeypatch.setattr(viewer, "setup_root_logger", lambda: root_logger)
    monkeypatch.setattr("zstarview.splash.setup_app", lambda _app_name: app)
    monkeypatch.setattr(viewer, "_load_star_catalog_for_launch", lambda _vmag_limit: object())
    monkeypatch.setattr(viewer, "_load_dso_catalog_for_launch", lambda: None)
    monkeypatch.setattr(viewer, "prepare_window_catalogs", lambda *args, **kwargs: object())
    monkeypatch.setattr(viewer, "prepare_window_user_options", lambda **kwargs: object())
    monkeypatch.setattr(viewer, "prepare_window_runtime_options", lambda **kwargs: object())
    def _raise_location_error(*_args, **_kwargs):
        raise LocationResolveError()

    monkeypatch.setattr(viewer, "resolve_launch_location", _raise_location_error)
    monkeypatch.setattr("zstarview.gui.window.FramelessSkyWindow", _DummyWindow)
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
    assert single_shot_calls == [0]
    assert len(root_logger.added) == 1
    assert root_logger.removed == root_logger.added
    assert _DummyWindow.last_instance is not None
    assert _DummyWindow.last_instance.overlay.shown == 1
    assert _DummyWindow.last_instance.overlay.hidden == 0
    assert any("Startup failed" in line for line, _level in _DummyWindow.last_instance.overlay.lines)
