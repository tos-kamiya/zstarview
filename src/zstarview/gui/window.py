"""
The main window of the ZStarView application.

This module defines the `SkyWindow` class, which is the primary user interface
for the application. It handles rendering the celestial objects, sky background,
clouds, and all user interactions like rotation, zooming, and object highlighting.
"""

import logging
import os
import time
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Union

import astropy.time
from PySide6.QtCore import (
    QPoint,
    QRect,
    QSize,
    Qt,
    QTimer,
    QUrl,
    Signal,
)
from PySide6.QtGui import (
    QAction,
    QCloseEvent,
    QDesktopServices,
    QFont,
    QFontDatabase,
    QGuiApplication,
    QIcon,
    QPaintEvent,
    QResizeEvent,
)
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QMenu,  # noqa: F401 - preserve monkeypatch/public compatibility
    QWidget,
)

from ..aircraft_constants import AIRCRAFT_REFRESH_INTERVAL_SECONDS
from ..asterisms import ASTERISM_KEYS_BY_SOURCE_ID
from ..astro import (  # noqa: F401 - calculate_visible_stars remains patchable here
    calculate_visible_stars,
    load_ephemeris,
    radec_to_altaz,
)
from ..cli.args import SKY_OPACITY_DEFAULT
from ..clouddisc import (
    CloudDisc,
    CloudDiscConfig,
)
from ..clouddisc.providers.select import pick_satellite
from ..geosatellite.pipeline import is_within_europe_band
from ..location_resolver import (
    project_place_targets_to_altaz as _project_place_targets_to_altaz,
)
from ..location_resolver import (
    search_place_candidates,
)
from ..overlay_time import overlay_availability_for_delta, target_time_utc_from_delta
from ..paths import (
    APP_ICON_FILE,
    CACHE_PATH,
    CLOUD_MISSING_TINT_RGBA,
    NIGHT_LIGHT_DEFAULT_OPACITY,
    OBSERVER_MAX_ALT_DEG,
    OBSERVER_MIN_ALT_DEG,
    OVERTURE_DERIVED_ROOT_DIR,
    PLATEAU_DERIVED_ROOT_DIR,
    ROAD_LIGHT_DEFAULT_OPACITY,
    STATUS_LINE_FONT_SIZE,
    TEXT_FONT_PATH,
    THEME_STYLES_BY_PRESET,
    TROPICAL_CYCLONE_DEFAULT_OPACITY,
    WINDOW_HEIGHT,
    WINDOW_WIDTH,
)
from ..render import geometry as render_geometry
from ..render import stars as render_stars
from ..render.pipeline import (
    compute_star_render_surface_size,
    compute_star_render_upscale_factor,
)
from ..satellite_constants import (
    SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS,
    SATELLITE_FAILURE_RETRY_SECONDS,
    SATELLITE_HORIZONS_CACHE_KEY,
    SATELLITE_ISS_CACHE_KEY,
)
from ..satellites import (
    fetch_horizons_lookup,
    find_satellite_altaz,
    load_satellite_cache,
    satellite_cache_scope_key,
)
from ..search.jpl import (
    project_jpl_target_altaz_from_state_vector,
    resolve_jpl_target_state_vector,
    search_jpl_targets,
)
from ..search.models import SearchJumpTarget
from ..search.satellites import search_satellite_targets
from ..simplified_view import resolve_simplified_view_mode
from .display_mode import (
    DISPLAY_MODE_INVERTED_CITY,
    DISPLAY_MODE_NORMAL,
    DISPLAY_MODE_SIMPLE_LABELS,
    DISPLAY_MODE_SIMPLE_NO_LABELS,
    display_mode_label,
    next_display_mode,
)
from ..types import ViewerData
from .aircraft_controller import AircraftController
from .aircraft_state import AircraftState
from .meteor_controller import MeteorController
from .meteor_state import MeteorState
from .cloud_controller import CloudController
from .cloud_state import CloudImageState
from .composite import SkyCompositorCache
from .draggable_window import DraggableWindow
from .famous_star_dialog import NamedStarJumpDialog
from .famous_star_search_dialog import NamedStarSearchDialog
from .famous_star_shortcuts import (
    NamedStarShortcut,
    build_place_search_jump_targets,
)
from .geosatellite_controller import GeoSatelliteController
from .geosatellite_state import GeoSatelliteState
from .jpl_small_body_controller import JplSmallBodyController
from .place_search_dialog import PlaceSearchDialog
from .precipitation_controller import PrecipitationController
from .road_night_lights_controller import RoadNightLightsController
from .satellite_controller import SatelliteController
from .satellite_state import SatelliteState
from .sky_worker import SkyDataWorker
from .terrain_controller import TerrainHorizonController
from .terrain_state import TerrainHorizonState
from .tropical_cyclone_controller import TropicalCycloneController
from .tropical_cyclone_state import TropicalCycloneState
from .urban_outline_controller import UrbanOutlineController
from .urban_outline_state import UrbanOutlineState
from .view_direction_dialog import ViewDirectionDialog
from .water_overlay_controller import WaterOverlayController
from .water_overlay_state import WaterOverlayState
from .window_actions import (
    OPEN_METEO_TERMS_URL,
    SkyWindowActionsMixin,
)
from .license_dialog import LicenseDialog
from .moon_hover_controller import MoonHoverController
from .solar_hover_controller import SolarHoverController
from .window_input import SkyWindowInputMixin
from .window_inputs import (
    PreparedWindowCatalogs,
    SkyWindowRuntimeOptions,
    SkyWindowUserOptions,
)
from .window_render import SkyWindowRenderMixin
from .window_state import SkyWindowState
from .window_updates import SkyWindowUpdatesMixin
from .window_widgets import (
    FramelessWindowFrame,
    ShutdownMessageOverlay,
    SkyWindowClientWidget,
    StartupLogOverlay,
)
from .worker_pool import shutdown_gui_worker_pool

PERIODIC_DEBUG_SNAPSHOT_INTERVAL_MS = 60_000

logger = logging.getLogger(__name__)


def _calendar_second_delay_ms(timestamp: float | None = None) -> int:
    """Return the delay to the next wall-clock second boundary."""
    now = time.time() if timestamp is None else float(timestamp)
    next_second = int(now) + 1
    return max(1, int((next_second - now) * 1000.0 + 0.999))


def open_code_data_licenses_and_credits() -> None:
    """Preserve the historical public helper location for callers and tests."""
    LicenseDialog().exec()


def open_open_meteo_terms() -> None:
    """Preserve a public helper location shared with other help actions."""
    QDesktopServices.openUrl(QUrl(OPEN_METEO_TERMS_URL))


def _replace_search_jump_target(
    target: object, /, **changes: object
) -> SearchJumpTarget:
    """Return an updated SearchJumpTarget."""
    return replace(target, **changes)  # type: ignore[arg-type]


def _replace_viewer_data(viewer_data: ViewerData, /, **changes: object) -> ViewerData:
    """Return updated viewer data."""
    return replace(viewer_data, **changes)


def _resize_event_size(event: QResizeEvent, attr: str) -> tuple[int, int]:
    value = getattr(event, attr, None)
    if callable(value):
        try:
            size = value()
            return int(size.width()), int(size.height())
        except Exception:
            pass
    return (-1, -1)


WindowGeometryArg = Union[str, tuple[int, int, int, int]]
DEFAULT_CLOUD_ALT_MIN_DEG = 1.0


def _clamp_window_geometry_to_screen(
    x: int,
    y: int,
    width: int,
    height: int,
    *,
    min_width: int,
    min_height: int,
) -> tuple[int, int, int, int]:
    width = max(int(width), int(min_width))
    height = max(int(height), int(min_height))

    screens = QGuiApplication.screens()
    if not screens:
        return int(x), int(y), width, height

    candidate = QRect(int(x), int(y), width, height)
    available_rect = None
    for screen in screens:
        rect = screen.availableGeometry()
        if rect.intersects(candidate):
            available_rect = rect
            break
    if available_rect is None:
        primary = QGuiApplication.primaryScreen()
        available_rect = (
            primary.availableGeometry()
            if primary is not None
            else screens[0].availableGeometry()
        )

    width = min(width, max(1, available_rect.width()))
    height = min(height, max(1, available_rect.height()))

    max_x = available_rect.right() - width + 1
    max_y = available_rect.bottom() - height + 1
    x = min(max(int(x), available_rect.left()), max_x)
    y = min(max(int(y), available_rect.top()), max_y)
    return x, y, width, height


class SkyWindowCoreMixin(
    SkyWindowInputMixin,
    SkyWindowActionsMixin,
    SkyWindowRenderMixin,
    SkyWindowUpdatesMixin,
):
    """
    The main application window, displaying the sky view.

    This class orchestrates the entire application. It manages fetching and
    processing of celestial and cloud data, renders all visual elements, and
    handles user input for interaction.
    """

    # Signal emitted when initial sky data is fully loaded and window can be shown.
    initial_data_loaded = Signal()
    # Signal to request repaint safely from background threads.
    cloud_repaint_requested = Signal()
    compute_star_render_surface_size = staticmethod(compute_star_render_surface_size)
    compute_star_render_upscale_factor = staticmethod(
        compute_star_render_upscale_factor
    )

    FRAMELESS_WINDOW = False

    def paintEvent(self, event: QPaintEvent) -> None:
        """Keep the host window from clearing the client area during repaints."""
        event.accept()

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self._startup_window_shown = True
        self._maybe_release_startup_input_block()
        if self._frameless_window or self._client_geometry_sync_done:
            return
        self._client_geometry_sync_done = True
        target_client_width, target_client_height = self._target_client_size
        self._resize_client_area(target_client_width, target_client_height)

    def __init__(
        self,
        viewer_data: ViewerData,
        catalogs: PreparedWindowCatalogs,
        user_options: SkyWindowUserOptions,
        runtime_options: SkyWindowRuntimeOptions,
        *,
        defer_initial_load: bool = False,
    ) -> None:
        """
        Initializes the SkyWindow.

        Args:
            viewer_data: Prepared observer/view state for the window.
            catalogs: Precomputed catalog data used by rendering and named-star jumps.
            user_options: User-facing display and toggle options.
            runtime_options: Runtime scheduling and window-hosting options.
        """
        super().__init__()
        self.star_catalog_np = catalogs.star_catalog_np
        self.star_catalog_lod6_indices = catalogs.star_catalog_lod6_indices
        self.star_catalog_meta = catalogs.star_catalog_meta
        self.dso_catalog_np = catalogs.dso_catalog_np
        self.show_dso: bool = self.dso_catalog_np is not None
        if user_options.show_dso_initial is not None:
            self.show_dso = (
                bool(user_options.show_dso_initial) and self.dso_catalog_np is not None
            )
        self.show_asterisms: bool = (
            True
            if user_options.show_asterisms_initial is None
            else bool(user_options.show_asterisms_initial)
        )
        self.show_guidelines: bool = (
            True
            if user_options.show_guidelines_initial is None
            else bool(user_options.show_guidelines_initial)
        )
        # Determine observation info startup mode. The user option may be a mode string
        # ('auto','top','bottom','off') or None.
        raw_observation_mode = user_options.observation_info_mode
        if isinstance(raw_observation_mode, str):
            observation_mode = raw_observation_mode
        elif raw_observation_mode is None:
            observation_mode = "auto"
        else:
            observation_mode = "auto" if bool(raw_observation_mode) else "off"
        self.observation_info_mode = observation_mode
        # Pinned when explicitly top or bottom
        self.observation_info_pinned = observation_mode in ("top", "bottom")
        # Visible unless 'off'
        self.show_observation_info: bool = observation_mode != "off"
        # Initialize overlay_info_bottom_left state: True means bottom-left, False means top-left
        # When pinned, set the HUD state immediately so rendering respects the pin.
        if observation_mode == "bottom":
            try:
                self.state.overlay_info_bottom_left = True
            except Exception:
                # state may not yet be fully initialized; ignore and rely on the render path
                pass
        elif observation_mode == "top":
            try:
                self.state.overlay_info_bottom_left = False
            except Exception:
                pass

        self._named_stars_by_band = catalogs.named_stars_by_band
        self._named_stars_search_all = catalogs.named_stars_search_all
        self.delta_t = runtime_options.delta_t
        self.runtime_options = runtime_options
        overlay_availability = overlay_availability_for_delta(self.delta_t)
        self.sky_disc_alpha = user_options.sky_disc_alpha
        self.sky_disc_style = user_options.sky_disc_style
        self.sky_disc_altaz_rings = user_options.sky_disc_altaz_rings
        self.sky_disc_altaz_rings_hover = user_options.sky_disc_altaz_rings_hover
        self.display_tone_curve = user_options.display_tone_curve
        self._sky_disc_alpha_when_enabled = (
            user_options.sky_disc_alpha
            if user_options.sky_disc_alpha > 0.0
            else SKY_OPACITY_DEFAULT
        )
        requested_cloud_alpha = user_options.cloud_disc_alpha
        self._cloud_requested_enabled = requested_cloud_alpha > 0.0
        # Cloud opacity is disabled if we are looking at a time-shifted view,
        # as we can only fetch current cloud data.
        self._cloud_alpha_when_enabled = (
            requested_cloud_alpha if requested_cloud_alpha > 0.0 else 0.2
        )
        self.cloud_disc_alpha: float = requested_cloud_alpha
        requested_satellite_opacity = user_options.satellite_opacity
        self._satellite_toggle_supported = overlay_availability.satellite
        self._satellite_requested_enabled = requested_satellite_opacity > 0.0
        self._satellite_opacity_when_enabled = (
            requested_satellite_opacity if requested_satellite_opacity > 0.0 else 1.0
        )
        self.satellite_opacity = (
            requested_satellite_opacity if self._satellite_toggle_supported else 0.0
        )
        requested_aircraft_opacity = user_options.aircraft_opacity
        self._aircraft_toggle_supported = overlay_availability.aircraft
        self._aircraft_requested_enabled = requested_aircraft_opacity > 0.0
        self._aircraft_opacity_when_enabled = (
            requested_aircraft_opacity if requested_aircraft_opacity > 0.0 else 1.0
        )
        self.aircraft_opacity = (
            requested_aircraft_opacity if self._aircraft_toggle_supported else 0.0
        )
        requested_meteor_opacity = float(user_options.meteor_trails_opacity)
        self._meteor_opacity_when_enabled = (
            requested_meteor_opacity if requested_meteor_opacity > 0.0 else 0.72
        )
        self.meteor_opacity = requested_meteor_opacity
        self.meteor_trails_max_candidates = int(
            user_options.meteor_trails_max_candidates
        )
        self._meteor_gui_allowed = bool(user_options.meteor_trails_gui_allowed)
        self.terrain_horizon_opacity = user_options.terrain_horizon_opacity
        self.earth_guide_opacity = user_options.earth_guide_opacity
        self.water_overlay_opacity = user_options.water_overlay_opacity
        requested_precipitation_opacity = user_options.precipitation_opacity
        self._precipitation_requested_enabled = requested_precipitation_opacity > 0.0
        self._precipitation_toggle_supported = bool(
            self._precipitation_requested_enabled and overlay_availability.precipitation
        )
        self._precipitation_opacity_when_enabled = (
            requested_precipitation_opacity
            if requested_precipitation_opacity > 0.0
            else 1.0
        )
        self.precipitation_opacity = (
            requested_precipitation_opacity
            if self._precipitation_toggle_supported
            else 0.0
        )
        self.precipitation_status = ""
        self.precipitation_forecast_time_utc = None
        self.precipitation_interval_seconds = None
        requested_road_light_opacity = user_options.road_light_opacity
        self._road_light_toggle_supported = bool(user_options.road_light_gui_allowed)
        self._road_light_opacity_when_enabled = (
            requested_road_light_opacity
            if requested_road_light_opacity > 0.0
            else ROAD_LIGHT_DEFAULT_OPACITY
        )
        self.road_night_lights_opacity = (
            requested_road_light_opacity
            if self._road_light_toggle_supported
            else 0.0
        )
        self.road_night_lights_status = ""
        requested_night_light_opacity = user_options.night_light_opacity
        self.ridge_glow_opacity = user_options.ridge_glow_opacity
        self._night_light_toggle_supported = bool(user_options.night_light_gui_allowed)
        self._night_light_opacity_when_enabled = (
            requested_night_light_opacity
            if requested_night_light_opacity > 0.0
            else NIGHT_LIGHT_DEFAULT_OPACITY
        )
        self.night_light_opacity = (
            requested_night_light_opacity if self._night_light_toggle_supported else 0.0
        )
        self.diffuse_sky_source = str(
            getattr(user_options, "diffuse_sky_source", "gaia")
        ).strip().lower()
        requested_akari_ir_bands_opacity = user_options.akari_ir_bands_opacity
        self._akari_ir_bands_toggle_supported = bool(
            user_options.akari_ir_bands_gui_allowed
        )
        self._akari_ir_bands_opacity_when_enabled = (
            requested_akari_ir_bands_opacity
            if requested_akari_ir_bands_opacity > 0.0
            else 0.1
        )
        self.akari_ir_bands_opacity = (
            requested_akari_ir_bands_opacity
            if self._akari_ir_bands_toggle_supported
            else 0.0
        )
        self.urban_outline_opacity = user_options.urban_outline_opacity
        self.ground_tint_opacity = user_options.ground_tint_opacity
        self._terrain_horizon_opacity_when_enabled = (
            user_options.terrain_horizon_opacity
            if user_options.terrain_horizon_opacity > 0.0
            else 0.25
        )
        self._water_overlay_opacity_when_enabled = (
            user_options.water_overlay_opacity
            if user_options.water_overlay_opacity > 0.0
            else 0.12
        )
        self._tropical_cyclone_toggle_supported = bool(
            user_options.tropical_cyclone_gui_allowed
            and overlay_availability.tropical_cyclone
        )
        self._tropical_cyclone_opacity_when_enabled = (
            user_options.tropical_cyclone_opacity
            if user_options.tropical_cyclone_opacity > 0.0
            else TROPICAL_CYCLONE_DEFAULT_OPACITY
        )
        self._tropical_cyclone_requested_enabled = (
            user_options.tropical_cyclone_opacity > 0.0
        )
        self.tropical_cyclone_opacity = (
            user_options.tropical_cyclone_opacity
            if self._tropical_cyclone_toggle_supported
            else 0.0
        )
        self._earth_guide_opacity_when_enabled = (
            user_options.earth_guide_opacity
            if user_options.earth_guide_opacity > 0.0
            else 0.028
        )
        self._urban_outline_opacity_when_enabled = (
            user_options.urban_outline_opacity
            if user_options.urban_outline_opacity > 0.0
            else 0.2
        )
        self._sky_disc_gui_allowed = bool(user_options.sky_disc_gui_allowed)
        self._cloud_gui_allowed = bool(user_options.cloud_gui_allowed)
        self._geo_satellite_enabled = bool(user_options.geo_satellite)
        self._geo_satellite_location_resolved = False
        self._satellite_gui_allowed = bool(user_options.satellite_gui_allowed)
        self._aircraft_gui_allowed = bool(user_options.aircraft_gui_allowed)
        self._tropical_cyclone_gui_allowed = bool(
            user_options.tropical_cyclone_gui_allowed
        )
        self._terrain_horizon_gui_allowed = bool(
            user_options.terrain_horizon_gui_allowed
        )
        self._earth_guide_gui_allowed = bool(user_options.earth_guide_gui_allowed)
        self._urban_outline_gui_allowed = bool(user_options.urban_outline_gui_allowed)
        self._water_overlay_gui_allowed = True
        self._clouddisc: CloudDisc | None = None
        self.show_urban_outline_layer: bool = self.urban_outline_opacity > 0.0
        self.inverted_city_enabled = False
        self.show_water_overlay_layer: bool = self.water_overlay_opacity > 0.0
        self.show_tropical_cyclone_overlay: bool = self.tropical_cyclone_opacity > 0.0
        self.moon_style = str(user_options.moon_style)
        self.moon_scale = int(user_options.moon_scale)
        self._configured_moon_style = self.moon_style
        self._configured_moon_scale = self.moon_scale
        self.bright_bodies_mode = user_options.bright_bodies_mode
        self.star_base_radius = user_options.star_base_radius
        self.vmag_limit = user_options.vmag_limit
        self.twinkle_count = user_options.twinkle_count
        self.twinkle_enabled = self.twinkle_count > 0
        self.light_background_star_outline = user_options.light_background_star_outline
        self.sky_update_interval = runtime_options.sky_update_interval
        self.urban_outline_radius_km = float(runtime_options.urban_outline_radius_km)
        self.urban_outline_skyscraper_radius_km = float(
            runtime_options.urban_outline_skyscraper_radius_km
        )
        self.urban_outline_min_height_m = float(
            runtime_options.urban_outline_min_height_m
        )
        self.urban_outline_max_candidates = int(
            runtime_options.urban_outline_max_candidates
        )
        self.road_light_max_candidates = int(runtime_options.road_light_max_candidates)
        self.urban_outline_feature_type = str(
            runtime_options.urban_outline_feature_type
        )
        self.urban_outline_skyscraper_only = bool(
            runtime_options.urban_outline_skyscraper_only
        )
        self.urban_outline_download_timeout_seconds = float(
            runtime_options.urban_outline_download_timeout_seconds
        )
        self.visual_preset = user_options.visual_preset
        self.presentation_id = (
            str(user_options.presentation_id).strip().lower() or "scenic"
        )
        self.star_data_policy = (
            str(user_options.star_data_policy).strip().lower() or "scenic_view_scoped"
        )
        self.theme = THEME_STYLES_BY_PRESET.get(
            self.visual_preset,
            THEME_STYLES_BY_PRESET["night"],
        )
        self.star_visibility_boost = user_options.star_visibility_boost
        self.asterism_visibility_boost = user_options.asterism_visibility_boost
        self.earth_guide_visibility_boost = user_options.earth_guide_visibility_boost
        self._star_render_expected_width = runtime_options.star_render_expected_width
        self.content_fov_deg = float(runtime_options.content_fov_deg)
        self._cloud_toggle_supported = overlay_availability.cloud and (
            self._clouddisc is not None or self._geo_satellite_enabled
        )
        if not self._cloud_toggle_supported:
            self.cloud_disc_alpha = 0.0

        # --- Viewer and Window Setup ---
        self.viewer_data = viewer_data
        self._search_view_center_base = tuple(self.viewer_data.view_center)
        self._search_view_center_alt_specified = False
        self._search_view_center_az_specified = False
        self.state = SkyWindowState(
            render_view_center=tuple(self.viewer_data.view_center),
            urban_outlines=None,
        )
        urban_outline_available = bool(
            self.urban_outline_opacity > 0.0 and self._urban_outline_gui_allowed
        )
        startup_mode = (
            DISPLAY_MODE_INVERTED_CITY
            if bool(user_options.inverted_city_initial) and urban_outline_available
            else DISPLAY_MODE_NORMAL
        )
        self.state.default_display_mode = startup_mode
        self.state.current_display_mode = startup_mode
        self.inverted_city_enabled = startup_mode == DISPLAY_MODE_INVERTED_CITY
        self._startup_initial_load_started = False
        self._startup_initial_sky_loaded = False
        self._startup_initial_terrain_loaded = False
        self._startup_initial_water_loaded = False
        self._startup_initial_urban_loaded = False
        self._startup_initial_night_light_loaded = False
        self._startup_initial_data_loaded = False
        self._startup_resize_pending = False
        self._post_startup_background_updates_started = False
        self._startup_window_shown = False
        self._startup_input_release_pending = False
        self._startup_input_blocked_state = True
        self._pending_periodic_debug_snapshot_path = None
        self._sky_refresh_due = False
        self._cloud_refresh_due = False
        self._satellite_refresh_due = False
        self._aircraft_refresh_due = False
        self._persistent_search_refresh_due = False
        self._viewport_rotation_keys_down: set[int] = set()
        # Ensure overlay_info_bottom_left reflects the startup mode now that
        # the mutable state object exists. True==bottom-left, False==top-left.
        try:
            if self.observation_info_mode == "bottom":
                self.state.overlay_info_bottom_left = True
            elif self.observation_info_mode == "top":
                self.state.overlay_info_bottom_left = False
        except Exception:
            pass

        self._enabled_satellite_groups: tuple[str, ...] = (
            SATELLITE_ISS_CACHE_KEY,
            SATELLITE_HORIZONS_CACHE_KEY,
        )
        self._disc_generation = 0
        self._frame_cache_key: object | None = None
        self._frame_cache_image = None
        self._present_frame_cache_key: object | None = None
        self._present_frame_cache_image = None
        self._display_frame_cache_key: object | None = None
        self._display_frame_cache_image = None
        self._fast_frame_base_cache_key: object | None = None
        self._fast_frame_base_cache_image = None
        self._fast_frame_cache_key: object | None = None
        self._fast_frame_cache_image = None
        self.setWindowTitle(self.viewer_data.city_name)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setFocus(Qt.FocusReason.ActiveWindowFocusReason)
        self.setWindowIcon(QIcon(APP_ICON_FILE))
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        initial_x = 100
        initial_y = 100
        initial_width = WINDOW_WIDTH
        initial_height = WINDOW_HEIGHT
        min_width = 400
        min_height = 400
        self.setMinimumSize(min_width, min_height)
        self._window_frame_mode = runtime_options.window_frame_mode
        self._frameless_window = bool(self.FRAMELESS_WINDOW)
        requested_geometry: tuple[int, int, int, int] | None = None
        if runtime_options.window_geometry_arg == "restore":
            load_window_geometry = runtime_options.load_last_window_geometry
            if load_window_geometry is not None:
                requested_geometry = load_window_geometry()
        elif isinstance(runtime_options.window_geometry_arg, tuple):
            requested_geometry = runtime_options.window_geometry_arg
        if requested_geometry is not None:
            initial_x, initial_y, initial_width, initial_height = requested_geometry
        initial_x, initial_y, initial_width, initial_height = (
            _clamp_window_geometry_to_screen(
                initial_x,
                initial_y,
                initial_width,
                initial_height,
                min_width=min_width,
                min_height=min_height,
            )
        )
        self._target_client_size = (int(initial_width), int(initial_height))
        self._client_geometry_sync_done = False
        self.setGeometry(initial_x, initial_y, initial_width, initial_height)
        self._client_widget = SkyWindowClientWidget(self)
        self._startup_log_overlay: StartupLogOverlay | None = None
        self._shutdown_overlay: ShutdownMessageOverlay | None = None
        self._frameless_frame: FramelessWindowFrame | None = None
        self.menu_button: QWidget | None = None
        self.size_grip: QWidget | None = None
        self._action_moon_option: QAction | None = None
        self._action_toggle_clouds: QAction | None = None
        self._action_toggle_geo_satellite: QAction | None = None
        self._action_toggle_satellites: QAction | None = None
        self._action_toggle_aircraft: QAction | None = None
        self._action_toggle_meteors: QAction | None = None
        self._action_toggle_terrain_horizon: QAction | None = None
        self._action_toggle_water_overlay: QAction | None = None
        self._action_toggle_earth_guide: QAction | None = None
        self._action_toggle_night_lights: QAction | None = None
        self._action_toggle_road_lights: QAction | None = None
        self._action_toggle_akari_ir_bands: QAction | None = None
        self._action_toggle_urban_outline: QAction | None = None
        self._action_toggle_tropical_cyclone: QAction | None = None
        self._action_toggle_dso: QAction | None = None
        self._action_toggle_asterisms: QAction | None = None
        self._action_toggle_twinkle: QAction | None = None
        self._action_toggle_guidelines: QAction | None = None
        self._action_toggle_observation_info: QAction | None = None
        self._action_square_window: QAction | None = None
        self._action_toggle_sky_disc: QAction | None = None
        self._action_raise_view: QAction | None = None
        self._action_lower_view: QAction | None = None
        self._build_window_menu()
        self._install_window_host()
        app = QApplication.instance()
        if app is not None:
            app.installEventFilter(self)
        self._client_widget.setFocus(Qt.FocusReason.ActiveWindowFocusReason)

        # --- Fonts ---
        text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
        text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
        self.text_font = QFont(text_font_family)
        self.text_font.setPointSizeF(float(user_options.overlay_font_size))
        self.status_line_font = QFont(text_font_family)
        self.status_line_font.setPointSizeF(float(STATUS_LINE_FONT_SIZE))

        # --- Data Update Timers and State ---
        self._is_shutting_down: bool = False
        self._moon_hover_controller = MoonHoverController(self)
        self._moon_hover_controller.image_ready.connect(self._on_moon_hover_image_ready)
        self._solar_hover_controller = SolarHoverController(self)
        self._solar_hover_controller.image_ready.connect(self._on_solar_hover_image_ready)
        self._setup_update_infrastructure()
        self._ephemeris = load_ephemeris()

        # --- Cloud Data State and Cache ---
        self.cloud_state = CloudImageState()
        self.geosatellite_state = GeoSatelliteState()
        self.satellite_state = SatelliteState()
        self.aircraft_state = AircraftState()
        self.meteor_state = MeteorState()
        self.tropical_cyclone_state = TropicalCycloneState()
        self.terrain_horizon_state = TerrainHorizonState()
        self.water_overlay_state = WaterOverlayState()
        self.road_night_light_polylines = None
        self.urban_outline_state = UrbanOutlineState()
        self._cloud_controller: CloudController | None = None
        self._geosatellite_controller: GeoSatelliteController | None = None
        self._satellite_controller: SatelliteController | None = None
        self._aircraft_controller: AircraftController | None = None
        self._meteor_controller: MeteorController | None = None
        self._tropical_cyclone_controller: TropicalCycloneController | None = None
        self._jpl_small_body_controller: JplSmallBodyController | None = None
        self._terrain_horizon_controller: TerrainHorizonController | None = None
        self._water_overlay_controller: WaterOverlayController | None = None
        self._precipitation_controller: PrecipitationController | None = None
        self._road_night_lights_controller: RoadNightLightsController | None = None
        self._urban_outline_controller: UrbanOutlineController | None = None
        # --- CloudDisc Service Initialization ---
        clouddisc_config = CloudDiscConfig(
            cache_dir=CACHE_PATH,
            sat_priority=("AUTO",),
            bt_warm_k=310.0,
            bt_cold_k=190.0,
            alt_min_deg=DEFAULT_CLOUD_ALT_MIN_DEG,
            search_back_minutes=120,
        )
        try:
            self._clouddisc = CloudDisc(clouddisc_config)
        except Exception as e:
            logger.warning(f"CloudDisc init failed: {e}")
        if self._clouddisc is not None:
            self._cloud_controller = CloudController(self._clouddisc, self)
            self._cloud_controller.cloud_started.connect(self._on_cloud_started)
            self._cloud_controller.cloud_source_ready.connect(
                self._on_cloud_source_ready
            )
            self._cloud_controller.cloud_ready.connect(self._on_cloud_ready)
            self._cloud_controller.cloud_failed.connect(self._on_cloud_failed)
        self._geosatellite_controller = GeoSatelliteController(parent=self)
        self._geosatellite_controller.geo_started.connect(self._on_geosatellite_started)
        self._geosatellite_controller.geo_source_ready.connect(
            self._on_geosatellite_source_ready
        )
        self._geosatellite_controller.geo_ready.connect(self._on_geosatellite_ready)
        self._geosatellite_controller.geo_failed.connect(self._on_geosatellite_failed)
        self._satellite_controller = SatelliteController(parent=self)
        self._satellite_controller.satellite_started.connect(self._on_satellite_started)
        self._satellite_controller.satellite_ready.connect(self._on_satellite_ready)
        self._satellite_controller.satellite_failed.connect(self._on_satellite_failed)
        self._aircraft_controller = AircraftController(parent=self)
        self._aircraft_controller.aircraft_started.connect(self._on_aircraft_started)
        self._aircraft_controller.aircraft_ready.connect(self._on_aircraft_ready)
        self._aircraft_controller.aircraft_failed.connect(self._on_aircraft_failed)
        self._meteor_controller = MeteorController(parent=self)
        self._meteor_controller.meteor_started.connect(self._on_meteor_started)
        self._meteor_controller.meteor_ready.connect(self._on_meteor_ready)
        self._meteor_controller.meteor_failed.connect(self._on_meteor_failed)
        if self._tropical_cyclone_toggle_supported:
            self._tropical_cyclone_controller = TropicalCycloneController(parent=self)
            self._tropical_cyclone_controller.cyclone_started.connect(
                self._on_tropical_cyclone_started
            )
            self._tropical_cyclone_controller.cyclone_ready.connect(
                self._on_tropical_cyclone_ready
            )
            self._tropical_cyclone_controller.cyclone_failed.connect(
                self._on_tropical_cyclone_failed
            )
        self._jpl_small_body_controller = JplSmallBodyController(parent=self)
        self._jpl_small_body_controller.jpl_started.connect(self._on_jpl_started)
        self._jpl_small_body_controller.jpl_ready.connect(self._on_jpl_ready)
        self._jpl_small_body_controller.jpl_failed.connect(self._on_jpl_failed)
        if self._action_toggle_clouds is not None:
            self._action_toggle_clouds.setEnabled(
                self._cloud_toggle_supported
                and (self._clouddisc is not None or self._geo_satellite_enabled)
                and self._cloud_gui_allowed
            )
        if self._action_toggle_satellites is not None:
            self._action_toggle_satellites.setEnabled(
                self._satellite_toggle_supported and self._satellite_gui_allowed
            )
        if self._action_toggle_aircraft is not None:
            self._action_toggle_aircraft.setEnabled(
                self._aircraft_toggle_supported and self._aircraft_gui_allowed
            )
        if self._action_toggle_meteors is not None:
            self._action_toggle_meteors.setEnabled(self._meteor_gui_allowed)
        if self._action_toggle_tropical_cyclone is not None:
            self._action_toggle_tropical_cyclone.setEnabled(
                self._tropical_cyclone_toggle_supported
            )

        terrain_cache_dir = Path(CACHE_PATH) / "copernicus-dem"
        self._terrain_horizon_controller = TerrainHorizonController(
            cache_dir=terrain_cache_dir,
            parent=self,
        )
        self._terrain_horizon_controller.terrain_started.connect(
            self._on_terrain_horizon_started
        )
        self._terrain_horizon_controller.terrain_ready.connect(
            self._on_terrain_horizon_ready
        )
        self._terrain_horizon_controller.terrain_failed.connect(
            self._on_terrain_horizon_failed
        )
        self._water_overlay_controller = WaterOverlayController(parent=self)
        self._water_overlay_controller.water_started.connect(
            self._on_water_overlay_started
        )
        self._water_overlay_controller.water_ready.connect(self._on_water_overlay_ready)
        self._water_overlay_controller.water_failed.connect(
            self._on_water_overlay_failed
        )
        if self.precipitation_opacity > 0.0:
            self._precipitation_controller = PrecipitationController(parent=self)
            self._precipitation_controller.precipitation_started.connect(
                self._on_precipitation_started
            )
            self._precipitation_controller.precipitation_ready.connect(
                self._on_precipitation_ready
            )
            self._precipitation_controller.precipitation_failed.connect(
                self._on_precipitation_failed
            )
        self._road_night_lights_controller = RoadNightLightsController(
            max_candidates=self.road_light_max_candidates,
            parent=self,
        )
        self._road_night_lights_controller.road_started.connect(
            self._on_road_night_lights_started
        )
        self._road_night_lights_controller.road_ready.connect(
            self._on_road_night_lights_ready
        )
        self._road_night_lights_controller.road_failed.connect(
            self._on_road_night_lights_failed
        )
        self._urban_outline_controller = UrbanOutlineController(
            derived_root_dir=Path(OVERTURE_DERIVED_ROOT_DIR),
            min_building_height_m=self.urban_outline_min_height_m,
            max_candidates=self.urban_outline_max_candidates,
            radius_km=self.urban_outline_radius_km,
            skyscraper_outer_radius_km=self.urban_outline_skyscraper_radius_km,
            feature_type=self.urban_outline_feature_type,
            skyscraper_only=self.urban_outline_skyscraper_only,
            plateau_root_dir=Path(PLATEAU_DERIVED_ROOT_DIR),
            download_timeout_s=self.urban_outline_download_timeout_seconds,
            parent=self,
        )
        self._urban_outline_controller.urban_started.connect(
            self._on_urban_outline_started
        )
        self._urban_outline_controller.urban_ready.connect(self._on_urban_outline_ready)
        self._urban_outline_controller.urban_failed.connect(
            self._on_urban_outline_failed
        )
        if self._action_toggle_terrain_horizon is not None:
            self._action_toggle_terrain_horizon.setEnabled(
                self._terrain_horizon_gui_allowed
            )
        if self._action_toggle_earth_guide is not None:
            self._action_toggle_earth_guide.setEnabled(self._earth_guide_gui_allowed)
        if self._action_toggle_night_lights is not None:
            self._action_toggle_night_lights.setEnabled(
                self._night_light_toggle_supported
            )
        if self._action_toggle_road_lights is not None:
            self._action_toggle_road_lights.setEnabled(
                self._road_light_toggle_supported
            )
        if self._action_toggle_akari_ir_bands is not None:
            self._action_toggle_akari_ir_bands.setEnabled(
                self._akari_ir_bands_toggle_supported
            )
        if self._action_toggle_sky_disc is not None:
            self._action_toggle_sky_disc.setEnabled(self._sky_disc_gui_allowed)
        if self._action_toggle_urban_outline is not None:
            self._action_toggle_urban_outline.setEnabled(
                self._urban_outline_gui_allowed
            )

        # --- Composition Cache (moved to dedicated class) ---
        target_stripes, width_factor = runtime_options.cloud_stripe_style
        missing_tint_alpha = int(
            round(255.0 * runtime_options.cloud_missing_tint_opacity)
        )
        missing_tint_rgba = (
            int(CLOUD_MISSING_TINT_RGBA[0]),
            int(CLOUD_MISSING_TINT_RGBA[1]),
            int(CLOUD_MISSING_TINT_RGBA[2]),
            missing_tint_alpha,
        )
        self._compositor = SkyCompositorCache(
            cloud_target_stripes=int(target_stripes),
            cloud_stripe_width_factor=float(width_factor),
            cloud_stripe_mode=runtime_options.cloud_stripe_mode,
            missing_tint_rgba=missing_tint_rgba,
            ground_tint_opacity=self.ground_tint_opacity,
        )

        # Cloud error banner is kept inside CloudImageState

        # --- Initial Data Load ---
        if not defer_initial_load:
            self.start_initial_data_load()

    @staticmethod
    def _resolve_periodic_debug_snapshot_dir() -> Path | None:
        raw = os.getenv("ZSTARVIEW_DEBUG_SAVE_PERIODIC_FRAME", "").strip()
        if not raw:
            return None
        lowered = raw.lower()
        if lowered in {"0", "false", "no", "off"}:
            return None
        if lowered in {"1", "true", "yes", "on"}:
            return Path(CACHE_PATH) / "debug" / "periodic"
        return Path(raw).expanduser()

    @staticmethod
    def _resolve_periodic_debug_snapshot_path(payload: dict) -> Path | None:
        output_dir = SkyWindowCoreMixin._resolve_periodic_debug_snapshot_dir()
        if output_dir is None:
            return None
        refreshed_at = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at, datetime):
            refreshed_at = datetime.now(timezone.utc)
        filename = f"periodic-{refreshed_at.strftime('%Y%m%dT%H%M%SZ')}.png"
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir / filename

    def _queue_periodic_debug_snapshot(self, payload: dict) -> None:
        output_path = SkyWindowCoreMixin._resolve_periodic_debug_snapshot_path(payload)
        if output_path is None:
            return
        self._pending_periodic_debug_snapshot_path = output_path

    @staticmethod
    def _save_periodic_debug_snapshot_image(image, output_path: Path) -> None:
        try:
            if not image.save(str(output_path), "PNG"):
                logger.warning(
                    "Failed to save periodic debug snapshot: %s", output_path
                )
                return
            logger.info("Saved periodic debug snapshot: %s", output_path)
        except Exception as exc:
            logger.warning("Periodic debug snapshot failed: %s", exc, exc_info=True)

    def _flush_periodic_debug_snapshot_save(self, present_frame) -> None:
        output_path = self._pending_periodic_debug_snapshot_path
        if not isinstance(output_path, Path):
            return
        self._pending_periodic_debug_snapshot_path = None
        self._save_periodic_debug_snapshot_image(present_frame, output_path)

    def _setup_update_infrastructure(self) -> None:
        """Initialize timers, worker, and signal wiring for background updates."""
        self._sky_worker = SkyDataWorker(self)
        self._sky_worker.data_ready.connect(self._on_sky_data_calculated)
        self._sky_worker.sky_disc_ready.connect(self._on_sky_disc_calculated)
        self._sky_worker.planet_data_ready.connect(self._on_planet_data_calculated)
        self.cloud_repaint_requested.connect(self.request_client_update)
        self.initial_data_loaded.connect(self._on_initial_data_loaded)

        app = QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(self._begin_shutdown)

        self._scheduler_tick_timer = QTimer(self)
        self._scheduler_tick_timer.setSingleShot(True)
        self._scheduler_tick_timer.timeout.connect(
            self._on_calendar_second_timeout
        )

        self._periodic_debug_snapshot_timer = QTimer(self)
        self._periodic_debug_snapshot_timer.setInterval(
            PERIODIC_DEBUG_SNAPSHOT_INTERVAL_MS
        )
        self._periodic_debug_snapshot_timer.timeout.connect(
            self._on_periodic_debug_snapshot_timer
        )
        if self._resolve_periodic_debug_snapshot_dir() is not None:
            self._periodic_debug_snapshot_timer.start()

        self._interaction_idle_timer = QTimer(self)
        self._interaction_idle_timer.setSingleShot(True)
        self._interaction_idle_timer.setInterval(self.state.interaction_idle_ms)
        self._interaction_idle_timer.timeout.connect(self._end_interaction_mode)

        self._viewport_interaction_idle_timer = QTimer(self)
        self._viewport_interaction_idle_timer.setSingleShot(True)
        self._viewport_interaction_idle_timer.setInterval(
            self.state.viewport_interaction_idle_ms
        )
        self._viewport_interaction_idle_timer.timeout.connect(
            self._end_viewport_interaction_mode
        )

    def apply_startup_viewer_data(self, viewer_data: ViewerData) -> None:
        """Replace the temporary startup viewer data with the resolved location."""
        self.viewer_data = viewer_data
        self._geo_satellite_location_resolved = True
        self._search_view_center_base = tuple(self.viewer_data.view_center)
        self.state.render_view_center = tuple(self.viewer_data.view_center)
        self.setWindowTitle(self.viewer_data.city_name)
        self._sync_geo_satellite_action_state()

    def apply_startup_delta_t(self, delta_t: timedelta) -> None:
        """Replace the temporary startup delta with the resolved launch time delta."""
        self.delta_t = delta_t
        overlay_availability = overlay_availability_for_delta(self.delta_t)
        self._satellite_toggle_supported = overlay_availability.satellite
        if not self._satellite_toggle_supported:
            self.satellite_opacity = 0.0
        elif self._satellite_requested_enabled:
            self.satellite_opacity = self._satellite_opacity_when_enabled
        self._aircraft_toggle_supported = overlay_availability.aircraft
        if not self._aircraft_toggle_supported:
            self.aircraft_opacity = 0.0
        elif self._aircraft_requested_enabled:
            self.aircraft_opacity = self._aircraft_opacity_when_enabled
        self._tropical_cyclone_toggle_supported = (
            overlay_availability.tropical_cyclone
            and (self._tropical_cyclone_controller is not None)
        )
        if not self._tropical_cyclone_toggle_supported:
            self.tropical_cyclone_opacity = 0.0
        elif self._tropical_cyclone_requested_enabled:
            self.tropical_cyclone_opacity = self._tropical_cyclone_opacity_when_enabled
        self._cloud_toggle_supported = overlay_availability.cloud and (
            self._clouddisc is not None or self._geo_satellite_enabled
        )
        if not self._cloud_toggle_supported:
            self.cloud_disc_alpha = 0.0
        elif self._cloud_requested_enabled:
            self.cloud_disc_alpha = self._cloud_alpha_when_enabled
        self._precipitation_toggle_supported = bool(
            overlay_availability.precipitation
            and getattr(self, "_precipitation_requested_enabled", False)
            and getattr(self, "_precipitation_controller", None) is not None
        )
        if not self._precipitation_toggle_supported:
            self.precipitation_opacity = 0.0
        elif getattr(self, "_precipitation_requested_enabled", False):
            self.precipitation_opacity = self._precipitation_opacity_when_enabled
        self._sync_cloud_action_state()
        if self._action_toggle_satellites is not None:
            self._action_toggle_satellites.setEnabled(
                self._satellite_toggle_supported and self._satellite_gui_allowed
            )
        if self._action_toggle_aircraft is not None:
            self._action_toggle_aircraft.setEnabled(
                self._aircraft_toggle_supported and self._aircraft_gui_allowed
            )
        if self._action_toggle_tropical_cyclone is not None:
            self._action_toggle_tropical_cyclone.setEnabled(
                self._tropical_cyclone_toggle_supported
                and self._tropical_cyclone_gui_allowed
            )
        precipitation_action = getattr(self, "_action_toggle_precipitation", None)
        if precipitation_action is not None:
            precipitation_action.setEnabled(
                self._precipitation_toggle_supported
            )
            precipitation_action.setChecked(self.precipitation_opacity > 0.0)

    def start_initial_data_load(self) -> None:
        """Kick off the initial background data load once the startup state is ready."""
        if self._startup_initial_load_started:
            return
        self._startup_initial_load_started = True
        self._startup_initial_sky_loaded = False
        self._startup_initial_terrain_loaded = False
        self._startup_initial_water_loaded = False
        self._startup_initial_urban_loaded = False
        self._startup_initial_night_light_loaded = (
            float(getattr(self, "terrain_horizon_opacity", 0.0)) <= 0.0
            or float(getattr(self, "night_light_opacity", 0.0)) <= 0.0
        )
        self._startup_initial_data_loaded = False
        self._startup_resize_pending = False
        self._post_startup_background_updates_started = False
        self._sky_refresh_due = False
        self._cloud_refresh_due = False
        self._satellite_refresh_due = False
        self._aircraft_refresh_due = False
        self._persistent_search_refresh_due = False
        self.start_background_sky_data_update(is_initial_load=True)

    def _resize_client_area(
        self, target_client_width: int, target_client_height: int
    ) -> None:
        """Resize the host so the client widget reaches the requested size."""
        if self.isFullScreen() or self.isMaximized():
            self.showNormal()
        self._target_client_size = (int(target_client_width), int(target_client_height))
        current_client_width = int(self.client_width())
        current_client_height = int(self.client_height())
        delta_width = int(target_client_width) - current_client_width
        delta_height = int(target_client_height) - current_client_height
        if delta_width == 0 and delta_height == 0:
            return
        self.resize(self.width() + delta_width, self.height() + delta_height)

    def _restore_default_window_size(self) -> None:
        """Return the window to the application's default client size."""
        if self.isFullScreen() or self.isMaximized():
            self.showNormal()
        self._resize_client_area(WINDOW_WIDTH, WINDOW_HEIGHT)

    def _fit_client_area_to_screen(self) -> None:
        """Resize the client area to the largest aspect-preserving size on screen."""
        current_client_width = max(1, int(self.client_width()))
        current_client_height = max(1, int(self.client_height()))
        if self.isFullScreen() or self.isMaximized():
            self.showNormal()

        screen = None
        window_handle = self.windowHandle()
        if window_handle is not None:
            screen = window_handle.screen()
        if screen is None:
            screen = self.screen()
        if screen is None:
            screen = QGuiApplication.primaryScreen()
        if screen is None:
            return

        available_geometry = screen.availableGeometry()
        if not available_geometry.isValid():
            return

        frame_width = max(0, int(self.width()) - int(self.client_width()))
        frame_height = max(0, int(self.height()) - int(self.client_height()))
        max_window_width = max(1, int((available_geometry.width() - frame_width) * 0.9))
        max_window_height = max(
            1, int((available_geometry.height() - frame_height) * 0.9)
        )

        target_window_size = QSize(current_client_width, current_client_height).scaled(
            max_window_width, max_window_height, Qt.AspectRatioMode.KeepAspectRatio
        )
        target_client_width = max(1, int(target_window_size.width()))
        target_client_height = max(1, int(target_window_size.height()))
        if (
            target_client_width == current_client_width
            and target_client_height == current_client_height
        ):
            return
        self._resize_client_area(target_client_width, target_client_height)

    def square_client_area(self) -> None:
        """Resize the client area so width and height match."""
        side = max(1, min(int(self.client_width()), int(self.client_height())))
        self._resize_client_area(side, side)

    def _install_window_host(self) -> None:
        """Install the host-specific window chrome around the shared client widget."""
        raise NotImplementedError

    def _ensure_shutdown_overlay(self) -> ShutdownMessageOverlay:
        if self._shutdown_overlay is None:
            self._shutdown_overlay = ShutdownMessageOverlay(self)
        self._shutdown_overlay.setGeometry(self._shutdown_message_geometry())
        return self._shutdown_overlay

    def _ensure_startup_log_overlay(self) -> StartupLogOverlay:
        if self._startup_log_overlay is None:
            theme = THEME_STYLES_BY_PRESET.get(
                self.visual_preset,
                THEME_STYLES_BY_PRESET["night"],
            )
            startup_log_alpha = (
                255
                if str(getattr(self, "presentation_id", "scenic")).strip().lower()
                == "instrument"
                else 180
            )
            self._startup_log_overlay = StartupLogOverlay(
                self._client_widget,
                text_rgb=theme.splash.info_text_rgb,
                background_rgba=(
                    int(theme.window_background.base_rgb[0]),
                    int(theme.window_background.base_rgb[1]),
                    int(theme.window_background.base_rgb[2]),
                    startup_log_alpha,
                ),
            )
        self._layout_startup_log_overlay()
        return self._startup_log_overlay

    def _layout_startup_log_overlay(self) -> None:
        overlay = self._startup_log_overlay
        if overlay is None:
            return
        overlay.setGeometry(self._client_widget.rect())
        if overlay.isVisible():
            overlay.raise_()

    def _shutdown_message_geometry(self) -> QRect:
        overlay = self._shutdown_overlay
        client_width = max(1, int(self.client_width()))
        width = max(1, int(round(client_width * 0.95)))
        width = min(width, client_width)
        if overlay is None:
            height = 64
        else:
            height = max(1, int(overlay.sizeHint().height()))
        x = max(0, (client_width - width) // 2)
        y = max(0, (int(self.height()) - height) // 2)
        return QRect(x, y, width, height)

    def _show_shutdown_message(self) -> None:
        overlay = self._ensure_shutdown_overlay()
        overlay.raise_()
        overlay.show()
        overlay.setGeometry(self._shutdown_message_geometry())
        QApplication.processEvents()

    def _hide_shutdown_message(self) -> None:
        if self._shutdown_overlay is not None:
            self._shutdown_overlay.hide()

    def _hide_startup_log_overlay(self) -> None:
        if self._startup_log_overlay is not None:
            self._startup_log_overlay.hide()

    def _request_application_quit(self) -> None:
        self._show_shutdown_message()
        QApplication.quit()

    def _size_grip_style_sheet(self) -> str:
        return "QWidget { border: none; background: transparent;}"

    def _raise_overlay_widgets(self) -> None:
        if self.menu_button is not None:
            self.menu_button.raise_()
        if self.size_grip is not None:
            self.size_grip.raise_()

    def _sync_viewport_interaction_chrome_visibility(self) -> None:
        if self.menu_button is None:
            return
        self.menu_button.setVisible(not bool(self.state.viewport_interaction_mode))

    def _simplified_view_enabled(self) -> bool:
        return self.state.current_display_mode in {
            DISPLAY_MODE_SIMPLE_NO_LABELS,
            DISPLAY_MODE_SIMPLE_LABELS,
        }

    def _simplified_view_labels_enabled(self) -> bool:
        return self.state.current_display_mode == DISPLAY_MODE_SIMPLE_LABELS

    def _effective_simplified_view_mode(self) -> str:
        return resolve_simplified_view_mode(
            base_enabled=self._simplified_view_enabled(),
            labels_enabled=self._simplified_view_labels_enabled(),
        )

    def _simplified_view_active(self) -> bool:
        return self._effective_simplified_view_mode() != "normal"

    def toggle_simplified_view(self) -> None:
        urban_outline_available = bool(
            self.urban_outline_opacity > 0.0 and self._urban_outline_gui_allowed
        )
        current_display_mode = getattr(self.state, "current_display_mode", None)
        if current_display_mode is None:
            if bool(getattr(self, "inverted_city_enabled", False)):
                current_display_mode = DISPLAY_MODE_INVERTED_CITY
            elif self._simplified_view_enabled():
                current_display_mode = (
                    DISPLAY_MODE_SIMPLE_LABELS
                    if self._simplified_view_labels_enabled()
                    else DISPLAY_MODE_SIMPLE_NO_LABELS
                )
            else:
                current_display_mode = DISPLAY_MODE_NORMAL
        self.state.current_display_mode = next_display_mode(
            current_display_mode,
            urban_outline_available=urban_outline_available,
        )
        self.inverted_city_enabled = (
            self.state.current_display_mode == DISPLAY_MODE_INVERTED_CITY
        )
        self.state.simplified_view_enabled = self.state.current_display_mode in {
            DISPLAY_MODE_SIMPLE_NO_LABELS,
            DISPLAY_MODE_SIMPLE_LABELS,
        }
        self.state.simplified_view_labels_enabled = (
            self.state.current_display_mode
            in {DISPLAY_MODE_NORMAL, DISPLAY_MODE_SIMPLE_LABELS}
        )
        self.request_client_update()

    def _mode_status_line(self) -> str:
        if self.state.current_display_mode == self.state.default_display_mode:
            return ""
        return f"{display_mode_label(self.state.current_display_mode)} [Space]"

    def client_width(self) -> int:
        if self._client_widget is not None:
            return self._client_widget.width()
        return (
            self.centralWidget().width()
            if self.centralWidget() is not None
            else super().width()
        )

    def client_height(self) -> int:
        if self._client_widget is not None:
            return self._client_widget.height()
        return (
            self.centralWidget().height()
            if self.centralWidget() is not None
            else super().height()
        )

    def client_size(self):
        if self._client_widget is not None:
            return self._client_widget.size()
        return (
            self.centralWidget().size()
            if self.centralWidget() is not None
            else super().size()
        )

    def client_rect(self):
        if self._client_widget is not None:
            return self._client_widget.rect()
        return (
            self.centralWidget().rect()
            if self.centralWidget() is not None
            else super().rect()
        )

    def request_client_update(self) -> None:
        if self._client_widget is not None:
            self._client_widget.update()
            return
        central = self.centralWidget()
        if central is not None:
            central.update()
            return
        super().update()

    def _on_moon_hover_image_ready(self, payload: object) -> None:
        if self._is_shutting_down or not isinstance(payload, tuple) or len(payload) != 2:
            return
        key, result = payload
        if result is None or not hasattr(result, "image"):
            return
        if self.state.moon_hover_image_key != key:
            return
        self.state.moon_hover_image = result
        self.request_client_update()

    def _on_solar_hover_image_ready(self, payload: object) -> None:
        if self._is_shutting_down or not isinstance(payload, tuple) or len(payload) != 2:
            return
        key, result = payload
        if result is None or not hasattr(result, "image"):
            return
        if self.state.solar_hover_image_key != key:
            return
        self.state.solar_hover_image = result
        self.request_client_update()

    def _handle_client_resize(self, event: QResizeEvent) -> None:
        _, _ = _resize_event_size(event, "oldSize")
        _, _ = _resize_event_size(event, "size")
        if not self._startup_initial_load_started:
            self._layout_startup_log_overlay()
            return
        if (
            hasattr(self, "_startup_initial_data_loaded")
            and not self._startup_initial_data_loaded
        ):
            self._startup_resize_pending = True
            self._layout_startup_log_overlay()
            self.request_client_update()
            self._raise_overlay_widgets()
            return
        self._begin_viewport_interaction_mode()
        self._disc_generation = int(self._disc_generation) + 1
        self._discard_stale_disc_images()
        if self._frameless_frame is None and self.menu_button is not None:
            button_size = self.menu_button.size()
            self.menu_button.move(self.client_width() - button_size.width(), 0)

        # Invalidate the composition cache since the size has changed
        self._compositor.invalidate()
        self.request_sky_data_update(
            reason="resize",
            allow_during_viewport_interaction=True,
        )
        self.request_client_update()
        self._raise_overlay_widgets()

    def _discard_stale_disc_images(self) -> None:
        discarded = False
        if self.state.sky_disc_image is not None:
            self.state.sky_disc_image = None
            discarded = True
        discarded = SkyWindow._clear_cloud_render_buffers(self) or discarded
        if discarded:
            self._compositor.invalidate()

    def _clear_cloud_render_buffers(
        self,
        *,
        preserve_cloud_buffers: bool = False,
        preserve_geo_satellite_buffers: bool = False,
    ) -> bool:
        if preserve_cloud_buffers:
            return False
        cleared = False
        if not preserve_cloud_buffers:
            cloud_state = self.cloud_state
            if cloud_state is not None:
                if cloud_state.image is not None:
                    cloud_state.image = None
                    cleared = True
                if cloud_state.missing_mask is not None:
                    cloud_state.missing_mask = None
                    cleared = True
                if cloud_state.cloud_amount_field is not None:
                    cloud_state.cloud_amount_field = None
                    cleared = True
                cloud_state.render_key = None
                cloud_state.request_id = None
                cloud_state.missing_mask_key = None
        geo_state = self.geosatellite_state
        if geo_state is not None and not preserve_geo_satellite_buffers:
            if geo_state.image is not None:
                geo_state.image = None
                cleared = True
            if geo_state.missing_mask is not None:
                geo_state.missing_mask = None
                cleared = True
            if geo_state.cloud_amount_field is not None:
                geo_state.cloud_amount_field = None
                cleared = True
            if geo_state.altaz_grid is not None:
                geo_state.altaz_grid = None
                cleared = True
            geo_state.render_key = None
            geo_state.request_id = None
            geo_state.missing_mask_key = None
        if cleared:
            self.state.cloud_projection_next_refresh_utc = None
        return cleared

    def show_menu(self) -> None:
        if self.menu_button is None:
            return
        self._sync_view_altitude_actions()
        menu_pos = self.menu_button.mapToGlobal(QPoint(0, self.menu_button.height()))
        self.menu.exec(menu_pos)

    def _current_time_obj(self) -> astropy.time.Time:
        return astropy.time.Time(self._target_time_utc())

    def _open_named_star_jump_dialog(self) -> None:
        dialog = NamedStarJumpDialog(self._named_stars_by_band, self)
        if dialog.exec() == 0:
            return
        star = dialog.selected_star()
        if star is None:
            return
        self._jump_to_named_star(star)

    def _open_named_star_search_dialog(self) -> None:
        dialog = NamedStarSearchDialog(
            self._named_stars_search_all,
            self,
            cli_view_center_alt_specified=bool(self._search_view_center_alt_specified),
            cli_view_center_az_specified=bool(self._search_view_center_az_specified),
            satellite_search_callback=self._search_satellite_targets,
            jpl_search_callback=self._search_jpl_targets,
        )
        if dialog.exec() == 0:
            return
        if dialog.should_clear_persistent_marker():
            self._clear_persistent_search()
            return
        target = dialog.selected_target()
        if target is None:
            return
        self._jump_to_search_target(target)

    def _open_place_search_dialog(self) -> None:
        dialog = PlaceSearchDialog(
            self._search_place_jump_targets,
            self,
            cli_view_center_alt_specified=bool(self._search_view_center_alt_specified),
            cli_view_center_az_specified=bool(self._search_view_center_az_specified),
        )
        if dialog.exec() == 0:
            return
        target = dialog.selected_target()
        if target is None:
            return
        self._jump_to_search_target(target)

    def _open_view_direction_dialog(self) -> None:
        dialog = ViewDirectionDialog(tuple(self.viewer_data.view_center), self)
        if dialog.exec() == 0:
            return
        alt_deg, az_deg = dialog.selected_view_center()
        self._set_view_center(
            alt_deg,
            az_deg,
            interactive_viewport=True,
            start_viewport_idle_timer=False,
        )
        QTimer.singleShot(0, self._finalize_view_direction_change)

    def _finalize_view_direction_change(self) -> None:
        self._end_viewport_interaction_mode(reason="view-change-release")

    def _finalize_view_direction_dialog_change(self) -> None:
        self._finalize_view_direction_change()

    def _search_jpl_targets(self, query: str) -> list[SearchJumpTarget]:
        target_time_utc = self._target_time_utc()
        return search_jpl_targets(
            query,
            target_time_utc=target_time_utc,
            lookup_fetch=fetch_horizons_lookup,
        )

    def _search_satellite_targets(self, query: str) -> list[SearchJumpTarget]:
        return search_satellite_targets(query, target_time_utc=self._target_time_utc())

    def _find_satellite_jump_altaz(self, object_key: str) -> tuple[float, float] | None:
        records = self.satellite_state.records_by_group or None
        if not records:
            records = self._load_cached_satellite_records(
                tuple(self._enabled_satellite_groups)
            )
        if not records:
            return None
        return find_satellite_altaz(
            records,
            object_key=object_key,
            observer_lat=float(self.viewer_data.location[0]),
            observer_lon=float(self.viewer_data.location[1]),
            observer_height_m=float(self.viewer_data.observer_height_m),
            time_obj=self._current_time_obj(),
        )

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
                # Horizons observer CSV reports azimuth first and elevation second.
                return numeric_values[1], numeric_values[0]
        return None

    def _jpl_small_body_persistent_target(self) -> SearchJumpTarget | None:
        target = self.state.persistent_search_target
        if target is None:
            return None
        if not bool(target.persistent_keep_marker):
            return None
        return target

    def _clear_persistent_search(self) -> None:
        self.state.persistent_search_target = None
        self.state.persistent_search_reference_time_utc = None
        self.state.persistent_search_next_refresh_utc = None
        self.state.persistent_search_last_refresh_utc = None
        self.state.persistent_search_last_error = None
        self._persistent_search_refresh_due = False

    def _schedule_persistent_search_refresh(self) -> None:
        if self._is_shutting_down:
            return
        target = self._jpl_small_body_persistent_target()
        if target is None:
            self._clear_persistent_search()
            return
        next_refresh_utc = self.state.persistent_search_next_refresh_utc
        if next_refresh_utc is None:
            reference_time_utc = target.target_time_utc or self._target_time_utc()
            self.state.persistent_search_reference_time_utc = reference_time_utc
            next_refresh_utc = reference_time_utc + timedelta(hours=1)
            self.state.persistent_search_next_refresh_utc = next_refresh_utc

    def _start_persistent_search_refresh(self, *, reason: str = "timer") -> bool:
        target = self._jpl_small_body_persistent_target()
        if target is None:
            self._clear_persistent_search()
            return False
        query_time_utc = self.state.persistent_search_next_refresh_utc
        if query_time_utc is None:
            query_time_utc = target.target_time_utc or self._target_time_utc()
        if self._jpl_small_body_controller is None:
            return False
        return self._jpl_small_body_controller.update(
            observer_lat=float(self.viewer_data.location[0]),
            observer_lon=float(self.viewer_data.location[1]),
            observer_height_m=float(self.viewer_data.observer_height_m),
            target=target,
            target_time_utc=query_time_utc,
            reason=reason,
        )

    def _on_jpl_started(self, payload: object) -> None:
        banner = ""
        if isinstance(payload, dict):
            banner = str(payload.get("banner", "")).strip()
        if banner:
            self.state.persistent_search_last_error = None
        self.request_client_update()

    def _log_persistent_search_target_update(
        self,
        *,
        action: str,
        target: SearchJumpTarget,
        target_time_utc: datetime,
        alt_deg: float,
        az_deg: float,
    ) -> None:
        logger.info(
            "JPL persistent target %s: label=%s kind=%s group=%s target_time_utc=%s alt=%.1f az=%.1f command=%s",
            action,
            str(target.label).strip() or "<unnamed>",
            str(target.kind).strip() or "<unknown>",
            str(target.jpl_group).strip() or "<none>",
            target_time_utc.astimezone(timezone.utc).isoformat(),
            float(alt_deg),
            float(az_deg) % 360.0,
            str(target.command).strip() or "<missing>",
        )

    def _on_jpl_ready(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        target = payload.get("target")
        if not isinstance(target, SearchJumpTarget):
            return
        current_target = self.state.persistent_search_target
        if current_target is None:
            return
        if (
            current_target.command != target.command
            or current_target.label != target.label
        ):
            return
        target_time_utc = payload.get("target_time_utc")
        if not isinstance(target_time_utc, datetime):
            target_time_utc = current_target.target_time_utc or self._target_time_utc()
        horizons_epoch_utc = payload.get("horizons_epoch_utc")
        horizons_position_km = payload.get("horizons_position_km")
        horizons_velocity_km_s = payload.get("horizons_velocity_km_s")
        viewer_data = self.viewer_data
        projected_altaz = None
        if not (
            viewer_data is not None
            and isinstance(horizons_epoch_utc, datetime)
            and isinstance(horizons_position_km, (list, tuple))
            and isinstance(horizons_velocity_km_s, (list, tuple))
        ):
            return
        vector_target = _replace_search_jump_target(
            current_target,
            horizons_epoch_utc=horizons_epoch_utc,
            horizons_position_km=tuple(horizons_position_km),
            horizons_velocity_km_s=tuple(horizons_velocity_km_s),
        )
        projected_altaz = project_jpl_target_altaz_from_state_vector(
            vector_target,
            observer_lat=float(viewer_data.location[0]),
            observer_lon=float(viewer_data.location[1]),
            observer_height_m=float(viewer_data.observer_height_m),
            time_obj=self._current_time_obj(),
        )
        if projected_altaz is None:
            return
        alt_deg, az_deg = projected_altaz
        updated_target = _replace_search_jump_target(
            vector_target,
            alt_deg=alt_deg,
            az_deg=az_deg,
            target_time_utc=target_time_utc,
        )
        self._log_persistent_search_target_update(
            action="refreshed",
            target=updated_target,
            target_time_utc=target_time_utc,
            alt_deg=alt_deg,
            az_deg=az_deg,
        )
        self.state.persistent_search_target = updated_target
        self.state.persistent_search_reference_time_utc = target_time_utc
        self.state.persistent_search_last_refresh_utc = payload.get("refreshed_at_utc")
        self.state.persistent_search_last_error = None
        self.state.persistent_search_next_refresh_utc = target_time_utc + timedelta(
            hours=1
        )
        self.request_client_update()
        self._schedule_persistent_search_refresh()

    def _on_jpl_failed(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        target = payload.get("target")
        if not isinstance(target, SearchJumpTarget):
            return
        current_target = self.state.persistent_search_target
        if current_target is None:
            return
        if (
            current_target.command != target.command
            or current_target.label != target.label
        ):
            return
        refreshed_at_utc = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at_utc, datetime):
            refreshed_at_utc = datetime.now(timezone.utc)
        error_text = str(payload.get("error", "")).strip()
        banner_text = str(payload.get("banner", "")).strip()
        last_error = error_text or banner_text
        if last_error.casefold() == "none":
            last_error = ""
        self.state.persistent_search_last_error = last_error
        self.state.persistent_search_last_refresh_utc = refreshed_at_utc
        self.state.persistent_search_next_refresh_utc = refreshed_at_utc + timedelta(
            hours=1
        )
        self.request_client_update()
        self._schedule_persistent_search_refresh()

    def _jump_to_search_target(self, target: SearchJumpTarget) -> None:
        target_kind = target.kind
        state_vector_target = target
        current_time = self._target_time_utc()
        current_time_obj = self._current_time_obj()
        if target_kind == "satellite":
            if target.alt_deg is not None and target.az_deg is not None:
                target_altaz = (float(target.alt_deg), float(target.az_deg) % 360.0)
            else:
                target_altaz = self._find_satellite_jump_altaz(
                    target.object_key or target.label
                )
            if target_altaz is None:
                self.satellite_state.set_banner(
                    f"Satellites: {target.label} not available"
                )
                self.request_client_update()
                return
            target_alt = float(target_altaz[0])
            target_az = float(target_altaz[1]) % 360.0
        elif target_kind == "place":
            if target.latitude_deg is None or target.longitude_deg is None:
                return
            projection = _project_place_targets_to_altaz(
                observer_latitude_deg=float(self.viewer_data.location[0]),
                observer_longitude_deg=float(self.viewer_data.location[1]),
                observer_height_m=float(self.viewer_data.observer_height_m),
                target_latitude_deg=[float(target.latitude_deg)],
                target_longitude_deg=[float(target.longitude_deg)],
                target_height_m=[0.0],
            )[0]
            target_alt = float(projection.alt_deg)
            target_az = float(projection.az_deg) % 360.0
        elif target.kind in ("jpl_small_body", "jpl_body"):
            state_vector_target = target
            if (
                target.horizons_epoch_utc is None
                or target.horizons_position_km is None
                or target.horizons_velocity_km_s is None
            ):
                state_vector = resolve_jpl_target_state_vector(
                    target,
                    target_time_utc=target.target_time_utc or current_time,
                )
                if state_vector is None:
                    return
                horizons_epoch_utc, horizons_position_km, horizons_velocity_km_s = (
                    state_vector
                )
                state_vector_target = _replace_search_jump_target(
                    target,
                    horizons_epoch_utc=horizons_epoch_utc,
                    horizons_position_km=horizons_position_km,
                    horizons_velocity_km_s=horizons_velocity_km_s,
                )
            target_altaz = project_jpl_target_altaz_from_state_vector(
                state_vector_target,
                observer_lat=float(self.viewer_data.location[0]),
                observer_lon=float(self.viewer_data.location[1]),
                observer_height_m=float(self.viewer_data.observer_height_m),
                time_obj=current_time_obj,
            )
            if target_altaz is None:
                return
            target_alt = float(target_altaz[0])
            target_az = float(target_altaz[1]) % 360.0
        else:
            target_alt, target_az = radec_to_altaz(
                target.ra_hours,
                target.dec_deg,
                float(self.viewer_data.location[0]),
                float(self.viewer_data.location[1]),
                float(self.viewer_data.observer_height_m),
                current_time,
            )
            target_alt = float(target_alt)
            target_az = float(target_az) % 360.0
        # When jumping to a target, keep the jump highlight altitude as-reported
        # (may be negative), but clamp the actual view center to the horizon (0°)
        # so the view doesn't go below the horizon automatically.
        base_alt, base_az = tuple(self._search_view_center_base)
        preserve_cli_view_center = target.preserve_cli_view_center
        fixed_alt = bool(self._search_view_center_alt_specified)
        fixed_az = bool(self._search_view_center_az_specified)
        if preserve_cli_view_center is False:
            fixed_alt = False
            fixed_az = False
        new_alt = (
            float(base_alt)
            if fixed_alt
            else max(OBSERVER_MIN_ALT_DEG, min(OBSERVER_MAX_ALT_DEG, target_alt))
        )
        new_az = float(base_az) % 360.0 if fixed_az else target_az
        if not bool(self.state.viewport_interaction_mode):
            self._begin_viewport_interaction_mode(start_idle_timer=False)
        self.viewer_data = _replace_viewer_data(
            self.viewer_data, view_center=(new_alt, new_az)
        )
        self.state.render_view_center = (new_alt, new_az)
        self._sync_view_altitude_actions()
        self._update_viewport_interaction_stars()
        self.request_client_update()

        self.state.jump_highlight_name = target.label
        self.state.jump_highlight_altaz = (target_alt, target_az)
        self.state.jump_highlight_until_ms = (time.monotonic() * 1000.0) + 3000.0
        if bool(target.persistent_keep_marker):
            reference_time_utc = target.target_time_utc or current_time
            horizons_epoch_utc = state_vector_target.horizons_epoch_utc
            horizons_position_km = state_vector_target.horizons_position_km
            horizons_velocity_km_s = state_vector_target.horizons_velocity_km_s
            updated_target = _replace_search_jump_target(
                state_vector_target,
                alt_deg=target_alt,
                az_deg=target_az,
                persistent_keep_marker=True,
                horizons_epoch_utc=horizons_epoch_utc,
                horizons_position_km=horizons_position_km,
                horizons_velocity_km_s=horizons_velocity_km_s,
            )
            self._log_persistent_search_target_update(
                action="set",
                target=updated_target,
                target_time_utc=reference_time_utc,
                alt_deg=target_alt,
                az_deg=target_az,
            )
            self.state.persistent_search_target = updated_target
            self.state.persistent_search_reference_time_utc = reference_time_utc
            self.state.persistent_search_last_error = None
            self.state.persistent_search_last_refresh_utc = target.target_time_utc
            if target.kind in {"jpl_small_body", "jpl_body"}:
                self.state.persistent_search_next_refresh_utc = (
                    reference_time_utc + timedelta(hours=1)
                )
                self._schedule_persistent_search_refresh()
            else:
                self.state.persistent_search_next_refresh_utc = None
        else:
            self._clear_persistent_search()
        QTimer.singleShot(0, self._finalize_view_direction_change)

    def _search_place_jump_targets(
        self,
        query: str,
        countrycode: str | None = None,
        language: str = "en",
    ) -> list[SearchJumpTarget]:
        candidates = search_place_candidates(
            query, countrycode=countrycode, language=language
        )
        return build_place_search_jump_targets(candidates)

    def _jump_to_named_star(self, star: NamedStarShortcut) -> None:
        self._jump_to_search_target(
            SearchJumpTarget(
                label=star.name,
                ra_hours=star.ra_hours,
                dec_deg=star.dec_deg,
                kind=star.kind,
                subtitle=star.subtitle or f"Vmag {star.vmag:.2f}",
                sort_key=(star.vmag, star.name.casefold()),
                object_key=star.name if star.kind == "satellite" else "",
            )
        )

    def _load_cached_satellite_records(
        self, enabled_groups: tuple[str, ...]
    ) -> dict[str, list[dict[str, object]]]:
        records_by_group: dict[str, list[dict[str, object]]] = {}
        for group_key in enabled_groups:
            cache_scope = None
            if group_key == SATELLITE_HORIZONS_CACHE_KEY:
                cache_scope = satellite_cache_scope_key(
                    observer_lat=float(self.viewer_data.location[0]),
                    observer_lon=float(self.viewer_data.location[1]),
                    observer_height_m=float(self.viewer_data.observer_height_m),
                )
            if cache_scope is None:
                cached = load_satellite_cache(group_key)
            else:
                cached = load_satellite_cache(group_key, cache_scope_key=cache_scope)
            if cached is None:
                continue
            records_by_group[str(group_key)] = list(cached.records)
        return records_by_group

    def _startup_input_blocked(self) -> bool:
        return bool(self._startup_input_blocked_state)

    def _on_initial_data_loaded(self) -> None:
        self._startup_initial_data_loaded = True
        self._maybe_release_startup_input_block()

    def _maybe_release_startup_input_block(self) -> None:
        if self._is_shutting_down:
            return
        if not self._startup_input_blocked_state:
            return
        if not self._startup_initial_data_loaded:
            return
        if not self._startup_window_shown:
            return
        if self._startup_input_release_pending:
            return
        self._startup_input_release_pending = True
        QTimer.singleShot(0, self._release_startup_input_block)

    def _release_startup_input_block(self) -> None:
        self._startup_input_release_pending = False
        self._startup_input_blocked_state = False
        self._start_post_startup_background_updates()

    def _start_post_startup_background_updates(self) -> None:
        if self._is_shutting_down:
            return
        if self._post_startup_background_updates_started:
            return
        if not self._startup_initial_data_loaded:
            return
        self._post_startup_background_updates_started = True
        now = datetime.now(timezone.utc)
        self.state.sky_next_refresh_utc = now + timedelta(
            seconds=self.sky_update_interval
        )
        if self.cloud_disc_alpha > 0.0:
            # Start the first cloud fetch immediately after startup so the overlay
            # does not sit in the idle state for a full refresh interval.
            self.start_background_cloud_update(reason="initial")
        if self._satellite_layer_enabled():
            if (
                self.satellite_state.records_by_group
                and self.satellite_state.element_epoch_utc is not None
            ):
                self._schedule_next_satellite_refresh()
            else:
                self.state.satellite_next_refresh_utc = now
        if self._aircraft_layer_enabled():
            if (
                self.aircraft_state.snapshots
                and self.aircraft_state.last_success_utc is not None
            ):
                self._schedule_next_aircraft_refresh()
            else:
                self.state.aircraft_next_refresh_utc = now
        if float(getattr(self, "meteor_opacity", 0.0)) > 0.0:
            self.state.meteor_next_refresh_utc = now
        if self._tropical_cyclone_controller is not None:
            cyclone_snapshots = self.tropical_cyclone_state.snapshots
            if (
                cyclone_snapshots
                and self.tropical_cyclone_state.cached_at_utc is not None
            ):
                if any(
                    snapshot.has_projectable_timeline()
                    for snapshot in cyclone_snapshots
                ):
                    self.tropical_cyclone_state.next_check_utc = (
                        self.tropical_cyclone_state.cached_at_utc
                        + timedelta(minutes=90)
                    )
                    self.tropical_cyclone_state.next_refresh_utc = (
                        self.tropical_cyclone_state.cached_at_utc + timedelta(hours=3)
                    )
                else:
                    self.tropical_cyclone_state.next_check_utc = now
                    self.tropical_cyclone_state.next_refresh_utc = now
            else:
                self.tropical_cyclone_state.next_check_utc = now
                self.tropical_cyclone_state.next_refresh_utc = now
        self._arm_scheduler_tick_timer()
        self._on_scheduler_tick()

    def _begin_shutdown(self) -> None:
        """Stop scheduling new background work while the app is closing."""
        if self._is_shutting_down:
            return
        self._is_shutting_down = True
        logger.info("Shutdown requested; closing application.")
        self._show_shutdown_message()
        try:
            self._sky_worker.shutdown()
            if self._cloud_controller is not None:
                self._cloud_controller.shutdown()
            if self._geosatellite_controller is not None:
                self._geosatellite_controller.shutdown()
            if self._satellite_controller is not None:
                self._satellite_controller.shutdown(wait_timeout_s=2.0)
            if self._aircraft_controller is not None:
                self._aircraft_controller.shutdown()
            meteor_controller = getattr(self, "_meteor_controller", None)
            if meteor_controller is not None:
                meteor_controller.shutdown()
            if self._tropical_cyclone_controller is not None:
                self._tropical_cyclone_controller.shutdown()
            if self._jpl_small_body_controller is not None:
                self._jpl_small_body_controller.shutdown()
            if self._terrain_horizon_controller is not None:
                self._terrain_horizon_controller.shutdown()
            if self._water_overlay_controller is not None:
                self._water_overlay_controller.shutdown()
            precipitation_controller = getattr(
                self, "_precipitation_controller", None
            )
            if precipitation_controller is not None:
                precipitation_controller.shutdown()
            road_controller = getattr(self, "_road_night_lights_controller", None)
            if road_controller is not None:
                road_controller.shutdown()
            if self._urban_outline_controller is not None:
                self._urban_outline_controller.shutdown()
            shutdown_gui_worker_pool(wait=True)
            if (
                hasattr(self, "_scheduler_tick_timer")
                and self._scheduler_tick_timer.isActive()
            ):
                self._scheduler_tick_timer.stop()
            if self._interaction_idle_timer.isActive():
                self._interaction_idle_timer.stop()
            if self._viewport_interaction_idle_timer.isActive():
                self._viewport_interaction_idle_timer.stop()
        finally:
            self._hide_shutdown_message()

    def _arm_scheduler_tick_timer(self) -> None:
        if self._is_shutting_down:
            return
        self._scheduler_tick_timer.start(_calendar_second_delay_ms())

    def _on_calendar_second_timeout(self) -> None:
        if self._is_shutting_down:
            return
        try:
            self._on_scheduler_tick()
            self._on_asterism_check_tick()
        finally:
            self._arm_scheduler_tick_timer()

    def _on_asterism_check_tick(self) -> None:
        if self._is_shutting_down:
            return
        if not self.show_asterisms:
            return
        if self.state.mouse_pos is None or self.state.celestial_data is None:
            return

        render_viewer = self._viewer_data_for_render()
        geometry = render_geometry.get_screen_geometry(
            self.client_width(),
            self.client_height(),
            render_viewer.view_center[0],
            edge_fov_deg=render_viewer.edge_fov_deg,
            content_fov_deg=render_viewer.content_fov_deg,
        )
        highlighted = render_stars.find_highlighted_object(
            self.state.celestial_data,
            render_viewer,
            self.state.mouse_pos,
            geometry,
        )
        if highlighted is None:
            return
        obj, _ = highlighted
        if not isinstance(obj, dict):
            return
        source_id = str(obj.get("source_id", "")).strip()
        if not source_id:
            return
        if source_id not in ASTERISM_KEYS_BY_SOURCE_ID:
            return
        self.request_client_update()

    def _predicted_cloud_satellite(self) -> str:
        lat, lon = self.viewer_data.location
        return pick_satellite(lat, lon, ("AUTO",))

    def _geo_satellite_toggle_supported(self) -> bool:
        if not bool(self._geo_satellite_location_resolved):
            return True
        lat, lon = self.viewer_data.location
        return is_within_europe_band(float(lat), float(lon))

    def _satellite_layer_enabled(self) -> bool:
        return self._satellite_toggle_supported and self.satellite_opacity > 0.0

    def _stop_satellite_timers(self) -> None:
        self.state.satellite_next_refresh_utc = None

    def _target_time_utc(self) -> datetime:
        delta_t = self.delta_t
        return target_time_utc_from_delta(delta_t)

    def _satellite_validity_remaining_ms(self) -> int | None:
        refreshed_at_utc = self.satellite_state.refreshed_at_utc
        if refreshed_at_utc is None:
            return None
        age_seconds = abs(
            (datetime.now(timezone.utc) - refreshed_at_utc).total_seconds()
        )
        remaining_seconds = SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS - age_seconds
        return max(0, int(round(remaining_seconds * 1000.0)))

    def _schedule_next_satellite_refresh(self, delay_ms: int | None = None) -> None:
        if not self._satellite_layer_enabled() or self._is_shutting_down:
            return
        interval_ms = delay_ms
        if interval_ms is None:
            interval_ms = self._satellite_validity_remaining_ms()
        if interval_ms is None:
            interval_ms = SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS * 1000
        self.state.satellite_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            milliseconds=max(0, interval_ms)
        )

    def _schedule_satellite_retry_after_failure(self) -> None:
        if not self._satellite_layer_enabled() or self._is_shutting_down:
            return
        self.state.satellite_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=SATELLITE_FAILURE_RETRY_SECONDS
        )

    def _enable_satellite_layer(self, *, reason: str) -> None:
        if not self._satellite_layer_enabled():
            return
        if (
            self.satellite_state.records_by_group
            and self.satellite_state.element_epoch_utc is not None
        ):
            self.reproject_satellite_overlay()
            remaining_ms = self._satellite_validity_remaining_ms()
            if remaining_ms is None:
                self.start_background_satellite_update(reason=reason)
                return
            if remaining_ms <= 0:
                self.start_background_satellite_update(reason=reason)
            else:
                self._schedule_next_satellite_refresh(remaining_ms)
            return
        self.start_background_satellite_update(reason=reason)

    def _aircraft_layer_enabled(self) -> bool:
        return self._aircraft_toggle_supported and self.aircraft_opacity > 0.0

    def _stop_aircraft_timers(self) -> None:
        self.state.aircraft_next_refresh_utc = None

    def _aircraft_cache_age_seconds(self) -> float | None:
        last_success_utc = self.aircraft_state.last_success_utc
        if last_success_utc is None:
            return None
        return max(0.0, (datetime.now(timezone.utc) - last_success_utc).total_seconds())

    def _schedule_next_aircraft_refresh(self, delay_ms: int | None = None) -> None:
        if not self._aircraft_layer_enabled() or self._is_shutting_down:
            return
        interval_ms = (
            AIRCRAFT_REFRESH_INTERVAL_SECONDS * 1000
            if delay_ms is None
            else int(delay_ms)
        )
        self.state.aircraft_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            milliseconds=max(0, interval_ms)
        )

    def _enable_aircraft_layer(self, *, reason: str) -> None:
        if not self._aircraft_layer_enabled():
            return
        if (
            self.aircraft_state.snapshots
            and self.aircraft_state.last_success_utc is not None
        ):
            self.reproject_aircraft_overlay()
            age_seconds = self._aircraft_cache_age_seconds()
            if age_seconds is None:
                self.start_background_aircraft_update(reason=reason)
                return
            remaining_ms = max(
                0,
                int(round((AIRCRAFT_REFRESH_INTERVAL_SECONDS - age_seconds) * 1000.0)),
            )
            if remaining_ms <= 0:
                self.start_background_aircraft_update(reason=reason)
            else:
                self._schedule_next_aircraft_refresh(remaining_ms)
            return
        self.start_background_aircraft_update(reason=reason)

    def closeEvent(self, event: QCloseEvent) -> None:
        g = self.geometry()
        save_window_geometry = self.runtime_options.save_last_window_geometry
        if save_window_geometry is not None:
            save_window_geometry(
                g.x(), g.y(), self.client_width(), self.client_height()
            )
        self._begin_shutdown()
        super().closeEvent(event)


class SkyWindow(SkyWindowCoreMixin, DraggableWindow):
    """Frameless host window using the custom chrome wrapper."""

    FRAMELESS_WINDOW = True

    def _install_window_host(self) -> None:
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.FramelessWindowHint)
        self._frameless_frame = FramelessWindowFrame(self, self._client_widget)
        self.setCentralWidget(self._frameless_frame)
        self.setFocusProxy(self._client_widget)
        self._frameless_frame.setFocusProxy(self._client_widget)
        self.menu_button = self._frameless_frame.menu_button
        self.size_grip = self._frameless_frame.size_grip
        self.add_drag_target(self._client_widget)
        self.add_drag_exclusions(self._frameless_frame.overlay_widgets())


class StandardSkyWindow(SkyWindowCoreMixin, QMainWindow):
    """Standard decorated host window using the platform title bar and menu bar."""

    FRAMELESS_WINDOW = False

    def _install_window_host(self) -> None:
        self.setCentralWidget(self._client_widget)
        self.setFocusProxy(self._client_widget)
        menu_bar = self.menuBar()
        menu_bar.addMenu(self.file_menu)
        menu_bar.addMenu(self.search_menu)
        menu_bar.addMenu(self.display_menu)
        menu_bar.addMenu(self.observer_view_menu)
        menu_bar.addMenu(self.help_menu)
