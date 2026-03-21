from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QObject, QTimer, Signal

from .cli_args import parse_export_image_args
from .paths import APP_DISPLAY_NAME
from .splash import setup_app
from .startup import (
    StartupAbortError,
    _startup_load_dso,
    _startup_load_stars,
    _startup_parse_time_arguments,
    _startup_resolve_city,
    setup_root_logger,
)
from .ui.window import SkyWindow
from .ui.window_inputs import (
    PreparedWindowCatalogs,
    SkyWindowRuntimeOptions,
    SkyWindowUserOptions,
    prepare_window_catalogs,
    prepare_window_runtime_options,
    prepare_window_user_options,
    prepare_window_viewer_data,
)

logger = logging.getLogger(__name__)


@dataclass
class _LayerWaitState:
    required: bool
    ready: bool = False
    failed: bool = False
    detail: str | None = None

    @property
    def terminal(self) -> bool:
        return (not self.required) or self.ready or self.failed


def _build_window_inputs_from_args(
    args: object,
) -> tuple[PreparedWindowCatalogs, object, SkyWindowUserOptions, SkyWindowRuntimeOptions]:
    city = _startup_resolve_city(
        getattr(args, "city", ""),
        place_query=getattr(args, "place", None),
        place_countrycode=getattr(args, "place_countrycode", None),
        place_lang=getattr(args, "place_lang", "en"),
    )
    delta_t = _startup_parse_time_arguments(
        getattr(args, "datetime", None),
        getattr(args, "days", 0),
        getattr(args, "hours", 0),
    )
    star_catalog = _startup_load_stars(getattr(args, "vmag_limit", 6.0))
    dso_catalog = _startup_load_dso()

    view_center = (getattr(args, "view_center_alt", 90.0), getattr(args, "view_center_az", 180.0))
    view_center = (min(90.0, max(0.0, view_center[0])), view_center[1] % 360.0)
    cloud_stripe_count, cloud_stripe_width = getattr(args, "cloud_stripe", (50, 0.2))
    visual_preset = getattr(args, "theme", "night")
    star_visibility_boost = 1.12 if visual_preset == "white" else 1.05 if visual_preset == "day" else 1.0
    vmag_brightness_scale = -math.log10(getattr(args, "vmag_brightness_multiplier", 2.5))

    catalogs = prepare_window_catalogs(
        star_catalog,
        dso_catalog=dso_catalog,
        vmag_brightness_scale=vmag_brightness_scale,
    )
    viewer_data = prepare_window_viewer_data(
        city.display_name,
        (city.lat, city.lon, city.tz),
        view_center,
        content_fov_deg=getattr(args, "content_fov_deg", 100.0),
        observer_height_m=(
            city.observer_height_m
            if getattr(args, "observer_height_m", None) is None
            else getattr(args, "observer_height_m")
        ),
        location_height_label=city.location_height_label,
        location_height_m=city.location_height_m,
        show_observer_height=getattr(args, "observer_height_m", None) is not None,
    )
    user_options = prepare_window_user_options(
        sky_disc_alpha=getattr(args, "sky_opacity", 0.15),
        cloud_disc_alpha=(
            0.0
            if cloud_stripe_count == 0 or cloud_stripe_width == 0.0
            else getattr(args, "cloud_opacity", 0.15)
        ),
        aircraft_opacity=getattr(args, "aircraft_opacity", 0.5),
        terrain_horizon_opacity=getattr(args, "terrain_horizon_opacity", 0.05),
        urban_outline_opacity=getattr(args, "urban_outline_opacity", 0.2),
        ground_tint_opacity=getattr(args, "ground_tint_opacity", 0.1),
        enlarge_moon=bool(getattr(args, "enlarge_moon", False)),
        star_base_radius=getattr(args, "star_base_radius", 4.0),
        vmag_limit=getattr(args, "vmag_limit", 6.0),
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
        show_dso_initial=getattr(args, "show_dso_initial", None),
        show_asterisms_initial=getattr(args, "show_asterisms_initial", None),
        sky_disc_gui_allowed=getattr(args, "sky_opacity", 0.15) > 0.0,
        cloud_gui_allowed=getattr(args, "cloud_opacity", 0.15) > 0.0,
        aircraft_gui_allowed=getattr(args, "aircraft_opacity", 0.5) > 0.0,
        terrain_horizon_gui_allowed=getattr(args, "terrain_horizon_opacity", 0.05) > 0.0,
        urban_outline_gui_allowed=getattr(args, "urban_outline_opacity", 0.2) > 0.0,
    )
    runtime_options = prepare_window_runtime_options(
        delta_t=delta_t,
        sky_update_interval=60,
        urban_outline_radius_km=getattr(args, "urban_outline_radius_km", 2.5),
        urban_outline_min_height_m=getattr(args, "urban_outline_min_height_m", 0.0),
        urban_outline_feature_type=getattr(args, "urban_outline_feature_type", "both"),
        urban_outline_skyscraper_only=bool(getattr(args, "urban_outline_skyscraper_only", False)),
        cloud_stripe_style=(cloud_stripe_count, cloud_stripe_width),
        cloud_missing_tint_opacity=getattr(args, "cloud_missing_tint_opacity", 0.0),
        star_render_expected_width=getattr(args, "expected_render_width", 600),
        content_fov_deg=getattr(args, "content_fov_deg", 100.0),
        window_geometry_arg=None,
    )
    return catalogs, viewer_data, user_options, runtime_options


class _ExportImageRunner(QObject):
    finished = Signal(int)

    def __init__(self, *, window: SkyWindow, output_path: Path, timeout_seconds: float, allow_partial_data: bool) -> None:
        super().__init__(window)
        self._window = window
        self._output_path = output_path
        self._allow_partial_data = bool(allow_partial_data)
        self._done = False
        self._states = {
            "sky": _LayerWaitState(required=True),
            "cloud": _LayerWaitState(required=float(window.cloud_disc_alpha) > 0.0),
            "terrain": _LayerWaitState(required=float(window.terrain_horizon_opacity) > 0.0),
            "urban": _LayerWaitState(required=float(window.urban_outline_opacity) > 0.0),
            "aircraft": _LayerWaitState(required=float(window.aircraft_opacity) > 0.0),
        }
        if self._states["cloud"].required and getattr(window, "_clouddisc", None) is None:
            self._mark_failed("cloud", "cloud service unavailable")
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(max(0, int(round(float(timeout_seconds) * 1000.0))))
        self._timer.timeout.connect(self._on_timeout)

        window.initial_data_loaded.connect(self._on_initial_data_loaded)
        if getattr(window, "_cloud_controller", None) is not None:
            window._cloud_controller.cloud_ready.connect(lambda _payload: self._mark_ready("cloud"))
            window._cloud_controller.cloud_failed.connect(
                lambda payload: self._mark_failed("cloud", str(payload.get("banner", "cloud failed")))
            )
        if getattr(window, "_terrain_horizon_controller", None) is not None:
            window._terrain_horizon_controller.terrain_ready.connect(lambda _payload: self._mark_ready("terrain"))
            window._terrain_horizon_controller.terrain_failed.connect(
                lambda payload: self._mark_failed("terrain", str(payload.get("banner", "terrain failed")))
            )
        if getattr(window, "_urban_outline_controller", None) is not None:
            window._urban_outline_controller.urban_ready.connect(lambda _payload: self._mark_ready("urban"))
            window._urban_outline_controller.urban_failed.connect(
                lambda payload: self._mark_failed("urban", str(payload.get("banner", "urban failed")))
            )
        if getattr(window, "_aircraft_controller", None) is not None:
            window._aircraft_controller.aircraft_ready.connect(lambda _payload: self._mark_ready("aircraft"))
            window._aircraft_controller.aircraft_failed.connect(
                lambda payload: self._mark_failed("aircraft", str(payload.get("banner", "aircraft failed")))
            )

    def start(self) -> None:
        if self._done:
            return
        if any(state.required for name, state in self._states.items() if name != "sky"):
            self._timer.start()
        self._check_finished()

    def _on_initial_data_loaded(self) -> None:
        self._mark_ready("sky")

    def _on_timeout(self) -> None:
        for name, state in self._states.items():
            if state.required and not state.terminal:
                self._mark_failed(name, f"{name} timed out")
        self._check_finished()

    def _mark_ready(self, name: str) -> None:
        state = self._states[name]
        if state.terminal:
            return
        state.ready = True
        state.detail = None
        self._check_finished()

    def _mark_failed(self, name: str, detail: str) -> None:
        state = self._states[name]
        if state.terminal:
            return
        state.failed = True
        state.detail = detail
        self._check_finished()

    def _check_finished(self) -> None:
        if self._done:
            return
        if not all(state.terminal for state in self._states.values()):
            return
        self._done = True
        if self._timer.isActive():
            self._timer.stop()

        failed_layers = [name for name, state in self._states.items() if state.failed]
        for name in failed_layers:
            logger.warning("Export layer unavailable: %s (%s)", name, self._states[name].detail)

        if self._states["sky"].failed:
            logger.error("Export aborted because celestial data never became ready.")
            self._shutdown(1)
            return

        failed_optional_layers = [name for name in failed_layers if name != "sky"]
        if failed_optional_layers and not self._allow_partial_data:
            logger.error("Export aborted because partial data is not allowed.")
            self._shutdown(1)
            return

        image = self._window.render_current_image(include_hud=False)
        self._output_path.parent.mkdir(parents=True, exist_ok=True)
        if not image.save(str(self._output_path), "PNG"):
            logger.error("Failed to save image: %s", self._output_path)
            self._shutdown(1)
            return
        logger.info("Saved image: %s", self._output_path)
        self._shutdown(0)

    def _shutdown(self, code: int) -> None:
        self._window._begin_shutdown()
        self.finished.emit(int(code))


def main() -> None:
    args = parse_export_image_args()
    setup_root_logger()
    logger.info("%s export-image starting...", APP_DISPLAY_NAME)

    try:
        catalogs, viewer_data, user_options, runtime_options = _build_window_inputs_from_args(args)
    except StartupAbortError:
        raise SystemExit(1)

    app = setup_app(f"{APP_DISPLAY_NAME} Export Image")
    app.setQuitOnLastWindowClosed(False)

    window = SkyWindow(
        viewer_data,
        catalogs,
        user_options=user_options,
        runtime_options=runtime_options,
    )
    image_width, image_height = getattr(args, "image_size")
    window.resize(int(image_width), int(image_height))

    runner = _ExportImageRunner(
        window=window,
        output_path=Path(getattr(args, "output")).expanduser(),
        timeout_seconds=float(getattr(args, "layer_timeout_seconds", 30.0)),
        allow_partial_data=bool(getattr(args, "allow_partial_data", False)),
    )
    runner.finished.connect(app.exit)
    QTimer.singleShot(0, runner.start)
    raise SystemExit(app.exec())
