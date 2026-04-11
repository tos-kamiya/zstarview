# -*- coding: utf-8 -*-
"""
The main window of the ZStarView application.

This module defines the `SkyWindow` class, which is the primary user interface
for the application. It handles rendering the celestial objects, sky background,
clouds, and all user interactions like rotation, zooming, and object highlighting.
"""

import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Tuple, Union

import astropy.time
from PySide6.QtCore import QEvent, QPoint, QRect, Qt, QTimer, Signal
from PySide6.QtGui import (
    QAction,
    QCloseEvent,
    QFont,
    QFontDatabase,
    QGuiApplication,
    QIcon,
    QKeyEvent,
    QKeySequence,
    QMouseEvent,
    QPaintEvent,
    QResizeEvent,
)
from PySide6.QtWidgets import QApplication, QMainWindow, QMenu, QPushButton, QSizeGrip
from PySide6.QtWidgets import QWidget

from ..__about__ import __version__
from ..astro import (
    calculate_visible_stars,
    radec_to_altaz,
)
from ..clouddisc import (
    CloudDisc,
    CloudDiscConfig,
)
from ..clouddisc.providers.select import pick_satellite
from ..overlay_time import overlay_availability_for_delta, target_time_utc_from_delta
from ..satellites import find_satellite_altaz, load_satellite_cache
from ..satellite_constants import (
    SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS,
    SATELLITE_FAILURE_RETRY_SECONDS,
    SATELLITE_ISS_CACHE_KEY,
    SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS,
)
from ..aircraft_constants import AIRCRAFT_REFRESH_INTERVAL_SECONDS
from ..aircraft_constants import AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
from ..paths import (
    APP_ICON_FILE,
    CACHE_PATH,
    CLOUD_UPDATE_INTERVAL,
    GUI_BUTTON_SIZE,
    GUI_MENU_TEXT_COLOR,
    TEXT_FONT_PATH,
    TEXT_FONT_SIZE,
    STATUS_LINE_FONT_SIZE,
    OVERTURE_DERIVED_ROOT_DIR,
    WINDOW_HEIGHT,
    WINDOW_WIDTH,
    CLOUD_MISSING_TINT_RGBA,
)
from ..config import load_last_window_geometry, save_last_window_geometry
from ..location_resolver import project_place_target_to_altaz, search_place_candidates
from ..render import geometry as render_geometry
from ..render import stars as render_stars
from ..render.pipeline import (
    compute_star_render_surface_size,
    compute_star_render_upscale_factor,
)
from ..types import ViewerData
from .draggable_window import DraggableWindow
from .composite import SkyCompositorCache
from .cloud_state import CloudImageState
from .cloud_controller import CloudController
from .satellite_state import SatelliteState
from .satellite_controller import SatelliteController
from .aircraft_state import AircraftState
from .aircraft_controller import AircraftController
from .terrain_state import TerrainHorizonState
from .terrain_controller import TerrainHorizonController
from .famous_star_dialog import NamedStarJumpDialog
from .famous_star_search_dialog import NamedStarSearchDialog
from .place_search_dialog import PlaceSearchDialog
from .famous_star_shortcuts import (
    NamedStarShortcut,
    SearchJumpTarget,
    build_place_search_jump_targets,
)
from ..asterisms import ASTERISM_KEYS_BY_SOURCE_ID
from .sky_worker import SkyDataWorker
from .window_inputs import PreparedWindowCatalogs
from .window_inputs import SkyWindowRuntimeOptions, SkyWindowUserOptions
from .window_render import SkyWindowRenderMixin
from .window_state import SkyWindowState
from .window_updates import SkyWindowUpdatesMixin
from .urban_outline_controller import UrbanOutlineController
from .urban_outline_state import UrbanOutlineState

logger = logging.getLogger(__name__)


WindowGeometryArg = Union[str, Tuple[int, int, int, int]]
DEFAULT_CLOUD_ALT_MIN_DEG = 1.0


def _clamp_window_geometry_to_screen(
    x: int,
    y: int,
    width: int,
    height: int,
    *,
    min_width: int,
    min_height: int,
) -> Tuple[int, int, int, int]:
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


class SkyWindowClientWidget(SkyWindowRenderMixin, QWidget):
    """Client-area widget shared by frameless and decorated host windows."""

    def __init__(self, owner: "SkyWindowCoreMixin") -> None:
        super().__init__(owner)
        self._owner = owner
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    def __getattr__(self, name: str) -> object:
        return getattr(self._owner, name)

    def resizeEvent(self, event: QResizeEvent) -> None:
        self._owner._handle_client_resize(event)
        super().resizeEvent(event)

    def leaveEvent(self, event: QEvent) -> None:
        self._owner._handle_client_leave(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        self._owner._handle_client_mouse_move(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        self._owner._handle_client_key_press(event)


class FramelessWindowFrame(QWidget):
    """Frameless-only window chrome that hosts the client widget and overlay controls."""

    def __init__(
        self, owner: "SkyWindowCoreMixin", client_widget: SkyWindowClientWidget
    ) -> None:
        super().__init__(owner)
        self._owner = owner
        self._client_widget = client_widget
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        self._client_widget.setParent(self)

        self.menu_button = QPushButton("", self)
        self.menu_button.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.menu_button.setStyleSheet(self._owner._menu_button_style_sheet())
        self.menu_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.menu_button.clicked.connect(self._owner.show_menu)

        self.size_grip = QSizeGrip(self)
        self.size_grip.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.size_grip.setStyleSheet(self._owner._size_grip_style_sheet())

        self._layout_chrome()

    def overlay_widgets(self) -> list[QWidget]:
        return [self.menu_button, self.size_grip]

    def _layout_chrome(self) -> None:
        self._client_widget.setGeometry(self.rect())
        button_size = self.menu_button.size()
        self.menu_button.move(self.width() - button_size.width(), 0)
        grip_size = self.size_grip.size()
        self.size_grip.move(
            self.width() - grip_size.width(), self.height() - grip_size.height()
        )
        self.menu_button.raise_()
        self.size_grip.raise_()

    def resizeEvent(self, event: QResizeEvent) -> None:
        self._layout_chrome()
        super().resizeEvent(event)


class SkyWindowCoreMixin(SkyWindowRenderMixin, SkyWindowUpdatesMixin):
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
        """Let the host window paint only its native chrome; the client widget renders the sky view."""
        QMainWindow.paintEvent(self, event)

    def __init__(
        self,
        viewer_data: ViewerData,
        catalogs: PreparedWindowCatalogs,
        user_options: SkyWindowUserOptions = SkyWindowUserOptions(),
        runtime_options: SkyWindowRuntimeOptions = SkyWindowRuntimeOptions(),
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
        self.show_overlay_info: bool = (
            True
            if user_options.show_overlay_info_initial is None
            else bool(user_options.show_overlay_info_initial)
        )
        self._named_stars_by_band = catalogs.named_stars_by_band
        self._named_stars_search_all = catalogs.named_stars_search_all
        self.delta_t = runtime_options.delta_t
        overlay_availability = overlay_availability_for_delta(self.delta_t)
        self.sky_disc_alpha = user_options.sky_disc_alpha
        self._sky_disc_alpha_when_enabled = (
            user_options.sky_disc_alpha if user_options.sky_disc_alpha > 0.0 else 0.3
        )
        requested_satellite_opacity = user_options.satellite_opacity
        self._satellite_toggle_supported = overlay_availability.satellite
        self._satellite_opacity_when_enabled = (
            requested_satellite_opacity if requested_satellite_opacity > 0.0 else 1.0
        )
        self.satellite_opacity = (
            requested_satellite_opacity if self._satellite_toggle_supported else 0.0
        )
        requested_aircraft_opacity = user_options.aircraft_opacity
        self._aircraft_toggle_supported = overlay_availability.aircraft
        self._aircraft_opacity_when_enabled = (
            requested_aircraft_opacity if requested_aircraft_opacity > 0.0 else 1.0
        )
        self.aircraft_opacity = (
            requested_aircraft_opacity if self._aircraft_toggle_supported else 0.0
        )
        self.terrain_horizon_opacity = user_options.terrain_horizon_opacity
        self.urban_outline_opacity = user_options.urban_outline_opacity
        self.ground_tint_opacity = user_options.ground_tint_opacity
        self._terrain_horizon_opacity_when_enabled = (
            user_options.terrain_horizon_opacity
            if user_options.terrain_horizon_opacity > 0.0
            else 0.25
        )
        self._urban_outline_opacity_when_enabled = (
            user_options.urban_outline_opacity
            if user_options.urban_outline_opacity > 0.0
            else 0.2
        )
        self._sky_disc_gui_allowed = bool(user_options.sky_disc_gui_allowed)
        self._cloud_gui_allowed = bool(user_options.cloud_gui_allowed)
        self._satellite_gui_allowed = bool(user_options.satellite_gui_allowed)
        self._aircraft_gui_allowed = bool(user_options.aircraft_gui_allowed)
        self._terrain_horizon_gui_allowed = bool(
            user_options.terrain_horizon_gui_allowed
        )
        self._urban_outline_gui_allowed = bool(user_options.urban_outline_gui_allowed)
        self.show_urban_outline_layer: bool = self.urban_outline_opacity > 0.0
        self.enlarge_moon = user_options.enlarge_moon
        self.star_base_radius = user_options.star_base_radius
        self.vmag_limit = user_options.vmag_limit
        self.sky_update_interval = runtime_options.sky_update_interval
        self.urban_outline_radius_km = float(runtime_options.urban_outline_radius_km)
        self.urban_outline_skyscraper_radius_km = float(
            runtime_options.urban_outline_skyscraper_radius_km
        )
        self.urban_outline_min_height_m = float(
            runtime_options.urban_outline_min_height_m
        )
        self.urban_outline_feature_type = str(
            runtime_options.urban_outline_feature_type
        )
        self.urban_outline_skyscraper_only = bool(
            runtime_options.urban_outline_skyscraper_only
        )
        self.visual_preset = user_options.visual_preset
        self.star_visibility_boost = user_options.star_visibility_boost
        self._star_render_expected_width = runtime_options.star_render_expected_width
        self.content_fov_deg = float(runtime_options.content_fov_deg)
        self._cloud_toggle_supported = overlay_availability.cloud

        # Cloud opacity is disabled if we are looking at a time-shifted view,
        # as we can only fetch current cloud data.
        requested_cloud_alpha = user_options.cloud_disc_alpha
        self._cloud_alpha_when_enabled = (
            requested_cloud_alpha if requested_cloud_alpha > 0.0 else 0.2
        )
        self.cloud_disc_alpha: float = requested_cloud_alpha
        if not self._cloud_toggle_supported:
            self.cloud_disc_alpha = 0.0

        # --- Viewer and Window Setup ---
        self.viewer_data = viewer_data
        self.state = SkyWindowState(
            render_view_center=tuple(self.viewer_data.view_center),
            urban_outlines=None,
            satellite_overlay_points=None,
            aircraft_overlay_points=None,
        )
        self._enabled_satellite_groups: tuple[str, ...] = (SATELLITE_ISS_CACHE_KEY,)
        self._frame_cache_key: object | None = None
        self._frame_cache_image = None
        self._present_frame_cache_key: object | None = None
        self._present_frame_cache_image = None
        self.setWindowTitle(self.viewer_data.city_name)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setFocus(Qt.FocusReason.ActiveWindowFocusReason)
        self.setWindowIcon(QIcon(APP_ICON_FILE))
        initial_x = 100
        initial_y = 100
        initial_width = WINDOW_WIDTH
        initial_height = WINDOW_HEIGHT
        min_width = 400
        min_height = 400
        self.setMinimumSize(min_width, min_height)
        self._window_frame_mode = runtime_options.window_frame_mode
        self._frameless_window = bool(self.FRAMELESS_WINDOW)
        requested_geometry: Optional[Tuple[int, int, int, int]] = None
        if runtime_options.window_geometry_arg == "restore":
            requested_geometry = load_last_window_geometry()
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
        self.setGeometry(initial_x, initial_y, initial_width, initial_height)
        self._client_widget = SkyWindowClientWidget(self)
        self._frameless_frame: Optional[FramelessWindowFrame] = None
        self.menu_button: Optional[QPushButton] = None
        self.size_grip: Optional[QSizeGrip] = None
        self._action_enlarge_moon: Optional[QAction] = None
        self._action_toggle_clouds: Optional[QAction] = None
        self._action_toggle_satellites: Optional[QAction] = None
        self._action_toggle_aircraft: Optional[QAction] = None
        self._action_toggle_terrain_horizon: Optional[QAction] = None
        self._action_toggle_urban_outline: Optional[QAction] = None
        self._action_toggle_dso: Optional[QAction] = None
        self._action_toggle_asterisms: Optional[QAction] = None
        self._action_toggle_guidelines: Optional[QAction] = None
        self._action_toggle_overlay_info: Optional[QAction] = None
        self._action_toggle_sky_disc: Optional[QAction] = None
        self._action_raise_view: Optional[QAction] = None
        self._action_lower_view: Optional[QAction] = None
        self._build_window_menu()
        self._install_window_host()
        self._client_widget.setFocus(Qt.FocusReason.ActiveWindowFocusReason)

        # --- Fonts ---
        text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
        text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
        self.text_font = QFont(text_font_family, TEXT_FONT_SIZE)
        self.status_line_font = QFont(text_font_family, STATUS_LINE_FONT_SIZE)

        # --- Data Update Timers and State ---
        self._is_shutting_down: bool = False
        self._setup_update_infrastructure()

        # --- Cloud Data State and Cache ---
        self.cloud_state = CloudImageState()
        self.satellite_state = SatelliteState()
        self.aircraft_state = AircraftState()
        self.terrain_horizon_state = TerrainHorizonState()
        self.urban_outline_state = UrbanOutlineState()
        self._cloud_controller: Optional[CloudController] = None
        self._satellite_controller: Optional[SatelliteController] = None
        self._aircraft_controller: Optional[AircraftController] = None
        self._terrain_horizon_controller: Optional[TerrainHorizonController] = None
        self._urban_outline_controller: Optional[UrbanOutlineController] = None
        self._cloud_update_timer = QTimer(self)
        self._cloud_update_timer.setInterval(CLOUD_UPDATE_INTERVAL * 1000)
        self._cloud_update_timer.timeout.connect(
            lambda: self.start_background_cloud_update(reason="timer")
        )
        self._satellite_update_timer = QTimer(self)
        self._satellite_update_timer.setSingleShot(True)
        self._satellite_update_timer.setInterval(
            SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS * 1000
        )
        self._satellite_update_timer.timeout.connect(self._on_satellite_refresh_timer)
        self._aircraft_update_timer = QTimer(self)
        self._aircraft_update_timer.setSingleShot(True)
        self._aircraft_update_timer.setInterval(
            AIRCRAFT_REFRESH_INTERVAL_SECONDS * 1000
        )
        self._aircraft_update_timer.timeout.connect(self._on_aircraft_refresh_timer)
        self._overlay_projection_timer = QTimer(self)
        self._overlay_projection_timer.setInterval(
            min(
                SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS,
                AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS,
            )
            * 1000
        )
        self._overlay_projection_timer.timeout.connect(
            self._on_overlay_projection_timer
        )

        # --- CloudDisc Service Initialization ---
        self._clouddisc: Optional[CloudDisc] = None
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
            self._cloud_controller.cloud_ready.connect(self._on_cloud_ready)
            self._cloud_controller.cloud_failed.connect(self._on_cloud_failed)
        self._satellite_controller = SatelliteController(parent=self)
        self._satellite_controller.satellite_started.connect(self._on_satellite_started)
        self._satellite_controller.satellite_ready.connect(self._on_satellite_ready)
        self._satellite_controller.satellite_failed.connect(self._on_satellite_failed)
        self._aircraft_controller = AircraftController(parent=self)
        self._aircraft_controller.aircraft_started.connect(self._on_aircraft_started)
        self._aircraft_controller.aircraft_ready.connect(self._on_aircraft_ready)
        self._aircraft_controller.aircraft_failed.connect(self._on_aircraft_failed)
        if self._action_toggle_clouds is not None:
            self._action_toggle_clouds.setEnabled(
                self._cloud_toggle_supported
                and self._clouddisc is not None
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
        self._urban_outline_controller = UrbanOutlineController(
            derived_root_dir=Path(OVERTURE_DERIVED_ROOT_DIR),
            min_building_height_m=self.urban_outline_min_height_m,
            radius_km=self.urban_outline_radius_km,
            skyscraper_outer_radius_km=self.urban_outline_skyscraper_radius_km,
            feature_type=self.urban_outline_feature_type,
            skyscraper_only=self.urban_outline_skyscraper_only,
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
        self.start_background_sky_data_update(is_initial_load=True)
        if self._clouddisc and self.cloud_disc_alpha > 0.0:
            self.start_background_cloud_update(reason="initial")
        if self.terrain_horizon_opacity > 0.0:
            self.start_background_terrain_horizon_update(reason="initial")
        if self.urban_outline_opacity > 0.0:
            self.start_background_urban_outline_update(reason="initial")
        if self._satellite_layer_enabled():
            self._enable_satellite_layer(reason="initial")
        if self._aircraft_layer_enabled():
            self._enable_aircraft_layer(reason="initial")

    def _setup_update_infrastructure(self) -> None:
        """Initialize timers, worker, and signal wiring for background updates."""
        self._sky_worker = SkyDataWorker(self)
        self._sky_worker.data_ready.connect(self._on_sky_data_calculated)
        self.cloud_repaint_requested.connect(self.request_client_update)

        app = QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(self._begin_shutdown)

        self._sky_data_update_timer = QTimer(self)
        self._sky_data_update_timer.timeout.connect(
            self.start_background_sky_data_update
        )

        self._asterism_check_timer = QTimer(self)
        self._asterism_check_timer.setInterval(1000)
        self._asterism_check_timer.timeout.connect(self._on_asterism_check_tick)
        self._asterism_check_timer.start()

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

    def _install_window_host(self) -> None:
        """Install the host-specific window chrome around the shared client widget."""
        raise NotImplementedError

    def _build_window_menu(self) -> None:
        """Build window actions and the popup menu shared by chrome implementations."""
        self.menu = QMenu("Menu", self)
        self.file_menu = QMenu("File", self)
        self.search_menu = QMenu("Search", self)
        self.observer_view_menu = QMenu("View Direction", self)
        self.display_menu = QMenu("Layers", self)
        if not self._frameless_window:
            self.menu.addMenu(self.file_menu)
        self.menu.addMenu(self.search_menu)
        self.menu.addMenu(self.display_menu)
        self.menu.addMenu(self.observer_view_menu)

        rotate_left = self.observer_view_menu.addAction(
            f"Rotate Left (-{self.state.rotation_step:.0f}°)"
        )
        rotate_left.triggered.connect(
            lambda: self._rotate_view(
                d_az=-self.state.rotation_step, interactive_viewport=True
            )
        )
        rotate_right = self.observer_view_menu.addAction(
            f"Rotate Right (+{self.state.rotation_step:.0f}°)"
        )
        rotate_right.triggered.connect(
            lambda: self._rotate_view(
                d_az=+self.state.rotation_step, interactive_viewport=True
            )
        )
        raise_view = self.observer_view_menu.addAction(
            f"Raise View (+{self.state.rotation_step:.0f}° alt)"
        )
        raise_view.triggered.connect(
            lambda: self._rotate_view(
                d_alt=+self.state.rotation_step, interactive_viewport=True
            )
        )
        self._action_raise_view = raise_view
        lower_view = self.observer_view_menu.addAction(
            f"Lower View (-{self.state.rotation_step:.0f}° alt)"
        )
        lower_view.triggered.connect(
            lambda: self._rotate_view(
                d_alt=-self.state.rotation_step, interactive_viewport=True
            )
        )
        self._action_lower_view = lower_view

        self.observer_view_menu.addSeparator()
        jump_named_star = QAction("Jump to Named Star...", self)
        jump_named_star.setShortcut(QKeySequence("Ctrl+J"))
        jump_named_star.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        jump_named_star.triggered.connect(self._open_named_star_jump_dialog)
        self.search_menu.addAction(jump_named_star)
        self.addAction(jump_named_star)
        search_named_star = QAction("Search Stars and Asterisms...", self)
        search_named_star.setShortcut(QKeySequence("Ctrl+F"))
        search_named_star.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        search_named_star.triggered.connect(self._open_named_star_search_dialog)
        self.search_menu.addAction(search_named_star)
        self.addAction(search_named_star)
        search_place = QAction("Search Places...", self)
        search_place.triggered.connect(self._open_place_search_dialog)
        self.search_menu.addAction(search_place)

        self.display_menu.addSeparator()
        toggle_enlarge_moon_action = QAction("Enlarge Moon", self)
        toggle_enlarge_moon_action.setCheckable(True)
        toggle_enlarge_moon_action.setChecked(self.enlarge_moon)
        toggle_enlarge_moon_action.setShortcut(QKeySequence(Qt.Key.Key_M))
        toggle_enlarge_moon_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_enlarge_moon_action.triggered.connect(self.toggle_enlarge_moon)
        self.display_menu.addAction(toggle_enlarge_moon_action)
        self.addAction(toggle_enlarge_moon_action)
        self._action_enlarge_moon = toggle_enlarge_moon_action
        toggle_dso_action = QAction("DSO", self)
        toggle_dso_action.setCheckable(True)
        toggle_dso_action.setChecked(self.show_dso)
        toggle_dso_action.setEnabled(self.dso_catalog_np is not None)
        toggle_dso_action.setShortcut(QKeySequence(Qt.Key.Key_D))
        toggle_dso_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_dso_action.triggered.connect(self.toggle_dso)
        self.display_menu.addAction(toggle_dso_action)
        self.addAction(toggle_dso_action)
        self._action_toggle_dso = toggle_dso_action
        toggle_asterisms_action = QAction("Asterisms", self)
        toggle_asterisms_action.setCheckable(True)
        toggle_asterisms_action.setChecked(self.show_asterisms)
        toggle_asterisms_action.setShortcut(QKeySequence(Qt.Key.Key_A))
        toggle_asterisms_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_asterisms_action.triggered.connect(self.toggle_asterisms)
        self.display_menu.addAction(toggle_asterisms_action)
        self.addAction(toggle_asterisms_action)
        self._action_toggle_asterisms = toggle_asterisms_action
        toggle_guidelines_action = QAction("Guidelines", self)
        toggle_guidelines_action.setCheckable(True)
        toggle_guidelines_action.setChecked(self.show_guidelines)
        toggle_guidelines_action.setShortcut(QKeySequence(Qt.Key.Key_G))
        toggle_guidelines_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_guidelines_action.triggered.connect(self.toggle_guidelines)
        self.display_menu.addAction(toggle_guidelines_action)
        self.addAction(toggle_guidelines_action)
        self._action_toggle_guidelines = toggle_guidelines_action
        toggle_overlay_info_action = QAction("Observation Info", self)
        toggle_overlay_info_action.setCheckable(True)
        toggle_overlay_info_action.setChecked(self.show_overlay_info)
        toggle_overlay_info_action.triggered.connect(self.toggle_overlay_info)
        self.display_menu.addAction(toggle_overlay_info_action)
        self.addAction(toggle_overlay_info_action)
        self._action_toggle_overlay_info = toggle_overlay_info_action

        self.display_menu.addSeparator()
        toggle_sky_disc_action = QAction("Sky Color Disc", self)
        toggle_sky_disc_action.setCheckable(True)
        toggle_sky_disc_action.setChecked(self.sky_disc_alpha > 0.0)
        toggle_sky_disc_action.setShortcut(QKeySequence(Qt.Key.Key_S))
        toggle_sky_disc_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_sky_disc_action.triggered.connect(self.toggle_sky_disc)
        self.display_menu.addAction(toggle_sky_disc_action)
        self.addAction(toggle_sky_disc_action)
        self._action_toggle_sky_disc = toggle_sky_disc_action
        toggle_clouds_action = QAction("Clouds", self)
        toggle_clouds_action.setCheckable(True)
        toggle_clouds_action.setChecked(self.cloud_disc_alpha > 0.0)
        toggle_clouds_action.setShortcut(QKeySequence(Qt.Key.Key_C))
        toggle_clouds_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_clouds_action.triggered.connect(self.toggle_clouds)
        self.display_menu.addAction(toggle_clouds_action)
        self.addAction(toggle_clouds_action)
        self._action_toggle_clouds = toggle_clouds_action
        toggle_satellites_action = QAction("Satellites", self)
        toggle_satellites_action.setCheckable(True)
        toggle_satellites_action.setChecked(self.satellite_opacity > 0.0)
        toggle_satellites_action.setShortcut(QKeySequence(Qt.Key.Key_I))
        toggle_satellites_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_satellites_action.triggered.connect(self.toggle_satellites)
        self.display_menu.addAction(toggle_satellites_action)
        self.addAction(toggle_satellites_action)
        self._action_toggle_satellites = toggle_satellites_action
        toggle_aircraft_action = QAction("Aircraft", self)
        toggle_aircraft_action.setCheckable(True)
        toggle_aircraft_action.setChecked(self.aircraft_opacity > 0.0)
        toggle_aircraft_action.setShortcut(QKeySequence(Qt.Key.Key_P))
        toggle_aircraft_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_aircraft_action.triggered.connect(self.toggle_aircraft)
        self.display_menu.addAction(toggle_aircraft_action)
        self.addAction(toggle_aircraft_action)
        self._action_toggle_aircraft = toggle_aircraft_action
        toggle_terrain_action = QAction("Terrain Horizon", self)
        toggle_terrain_action.setCheckable(True)
        toggle_terrain_action.setChecked(self.terrain_horizon_opacity > 0.0)
        toggle_terrain_action.setShortcut(QKeySequence(Qt.Key.Key_T))
        toggle_terrain_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_terrain_action.triggered.connect(self.toggle_terrain_horizon)
        self.display_menu.addAction(toggle_terrain_action)
        self.addAction(toggle_terrain_action)
        self._action_toggle_terrain_horizon = toggle_terrain_action
        toggle_urban_outline_action = QAction("Urban Outline", self)
        toggle_urban_outline_action.setCheckable(True)
        toggle_urban_outline_action.setChecked(self.urban_outline_opacity > 0.0)
        toggle_urban_outline_action.setShortcut(QKeySequence(Qt.Key.Key_U))
        toggle_urban_outline_action.setShortcutContext(
            Qt.ShortcutContext.WindowShortcut
        )
        toggle_urban_outline_action.triggered.connect(self.toggle_urban_outline)
        self.display_menu.addAction(toggle_urban_outline_action)
        self.addAction(toggle_urban_outline_action)
        self._action_toggle_urban_outline = toggle_urban_outline_action

        fullscreen_action = QAction("Fullscreen", self)
        fullscreen_action.setShortcut(QKeySequence(Qt.Key.Key_F11))
        fullscreen_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        fullscreen_action.triggered.connect(self.toggle_fullscreen)
        self.file_menu.addAction(fullscreen_action)
        self.addAction(fullscreen_action)

        self.file_menu.addSeparator()
        exit_action = self.file_menu.addAction("Exit")
        exit_action.setShortcut(QKeySequence(Qt.Key.Key_Q))
        exit_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        exit_action.triggered.connect(QApplication.quit)
        self.addAction(exit_action)

        self.display_menu.addSeparator()
        vmag_limit_action = self.display_menu.addAction(self._vmag_limit_menu_text())
        vmag_limit_action.setEnabled(False)
        version_action = self.file_menu.addAction(f"Version {__version__}")
        version_action.setEnabled(False)

        if self._frameless_window:
            self.menu.addSeparator()
            self.menu.addAction(fullscreen_action)
            self.menu.addAction(exit_action)

    def _attach_client_menu_button(self, parent: QWidget) -> None:
        """Attach the legacy popup-menu button directly on the client area."""
        self.menu_button = QPushButton("", parent)
        self.menu_button.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.menu_button.setStyleSheet(self._menu_button_style_sheet())
        self.menu_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.menu_button.clicked.connect(self.show_menu)
        self.menu_button.raise_()

    def _vmag_limit_menu_text(self) -> str:
        return f"Vmag limit {self.vmag_limit:.1f}"

    def _menu_button_style_sheet(self) -> str:
        text = "#%02x%02x%02x" % GUI_MENU_TEXT_COLOR
        return (
            "QPushButton {"
            " border: none;"
            " font-size: 16px;"
            " background: transparent;"
            f" color: {text};"
            "}"
            "QPushButton:hover {"
            " color: white;"
            " background-color: rgba(255, 255, 255, 0.10);"
            "}"
            "QPushButton:pressed {"
            " color: white;"
            " background-color: rgba(255, 255, 255, 0.16);"
            "}"
            "QPushButton:menu-indicator { image: none; }"
        )

    def _size_grip_style_sheet(self) -> str:
        return "QSizeGrip { border: none; background: transparent;}"

    def _raise_overlay_widgets(self) -> None:
        if self.menu_button is not None:
            self.menu_button.raise_()
        if self.size_grip is not None:
            self.size_grip.raise_()

    def _sync_view_altitude_actions(self) -> None:
        alt, _ = self.viewer_data.view_center
        if self._action_raise_view is not None:
            self._action_raise_view.setEnabled(float(alt) < 90.0)
        if self._action_lower_view is not None:
            self._action_lower_view.setEnabled(float(alt) > 0.0)

    def client_width(self) -> int:
        if getattr(self, "_client_widget", None) is not None:
            return self._client_widget.width()
        return (
            self.centralWidget().width()
            if self.centralWidget() is not None
            else super().width()
        )

    def client_height(self) -> int:
        if getattr(self, "_client_widget", None) is not None:
            return self._client_widget.height()
        return (
            self.centralWidget().height()
            if self.centralWidget() is not None
            else super().height()
        )

    def client_size(self):
        if getattr(self, "_client_widget", None) is not None:
            return self._client_widget.size()
        return (
            self.centralWidget().size()
            if self.centralWidget() is not None
            else super().size()
        )

    def client_rect(self):
        if getattr(self, "_client_widget", None) is not None:
            return self._client_widget.rect()
        return (
            self.centralWidget().rect()
            if self.centralWidget() is not None
            else super().rect()
        )

    def request_client_update(self) -> None:
        if getattr(self, "_client_widget", None) is not None:
            self._client_widget.update()
            return
        central = self.centralWidget()
        if central is not None:
            central.update()
            return
        super().update()

    def _handle_client_resize(self, event: QResizeEvent) -> None:
        self._begin_viewport_interaction_mode(preserve_cloud_buffers=True)
        self._disc_generation = int(getattr(self, "_disc_generation", 0)) + 1
        if self._frameless_frame is None and self.menu_button is not None:
            button_size = self.menu_button.size()
            self.menu_button.move(self.client_width() - button_size.width(), 0)
            self._raise_overlay_widgets()

        # Invalidate the composition cache since the size has changed
        self._compositor.invalidate()
        self.request_sky_data_update()
        self.request_client_update()
        self.start_background_cloud_update(reason="resize")

    def _discard_stale_disc_images(self) -> None:
        discarded = False
        if self.state.sky_disc_image is not None:
            self.state.sky_disc_image = None
            discarded = True
        if self.cloud_state.image is not None:
            self.cloud_state.image = None
            discarded = True
        if self.cloud_state.missing_mask is not None:
            self.cloud_state.missing_mask = None
            discarded = True
        if self.cloud_state.cloud_amount_field is not None:
            self.cloud_state.cloud_amount_field = None
            discarded = True
        if discarded:
            self.cloud_state.render_key = None
            self.cloud_state.request_id = None
            self.cloud_state.missing_mask_key = None
            self._compositor.invalidate()

    def _begin_interaction_mode(self) -> None:
        self.state.interaction_mode = True
        self._interaction_idle_timer.start()

    def _end_interaction_mode(self) -> None:
        if not self.state.interaction_mode:
            return
        self.state.interaction_mode = False
        self.request_sky_data_update()
        self.start_background_cloud_update(reason="view-change-idle")
        self.start_background_terrain_horizon_update(reason="view-change-idle")

    def _begin_viewport_interaction_mode(
        self, preserve_cloud_buffers: bool = False
    ) -> None:
        self.state.viewport_interaction_mode = True
        cloud_state = getattr(self, "cloud_state", None)
        cleared_cloud = False
        if cloud_state is not None:
            if not preserve_cloud_buffers:
                if getattr(cloud_state, "image", None) is not None:
                    cloud_state.image = None
                    cleared_cloud = True
                if getattr(cloud_state, "missing_mask", None) is not None:
                    cloud_state.missing_mask = None
                    cleared_cloud = True
                if getattr(cloud_state, "cloud_amount_field", None) is not None:
                    cloud_state.cloud_amount_field = None
                    cleared_cloud = True
                cloud_state.render_key = None
                cloud_state.request_id = None
                cloud_state.missing_mask_key = None
        if cleared_cloud:
            self._compositor.invalidate()
        cloud_controller = getattr(self, "_cloud_controller", None)
        if cloud_controller is not None:
            invalidate = getattr(
                cloud_controller, "invalidate_pending_render_results", None
            )
            if callable(invalidate):
                invalidate()
        self._viewport_interaction_idle_timer.start()

    def _update_viewport_interaction_stars(self) -> None:
        if self.state.celestial_data is None:
            self.state.viewport_interaction_stars = None
            return
        stars, _location = calculate_visible_stars(
            self.star_catalog_np,
            self.viewer_data.location[0],
            self.viewer_data.location[1],
            self.viewer_data.observer_height_m,
            self._current_time_obj(),
            self.state.render_view_center,
            max_vmag=4.0,
            subset_indices=self.star_catalog_lod6_indices,
        )
        self.state.viewport_interaction_stars = stars

    def _end_viewport_interaction_mode(self) -> None:
        if not self.state.viewport_interaction_mode:
            return
        self.state.viewport_interaction_mode = False
        self.state.viewport_interaction_stars = None
        self.request_sky_data_update()
        self.start_background_cloud_update(reason="view-change-idle")
        self.start_background_terrain_horizon_update(reason="view-change-idle")
        self.request_client_update()

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
        dialog = NamedStarSearchDialog(self._named_stars_search_all, self)
        if dialog.exec() == 0:
            return
        target = dialog.selected_target()
        if target is None:
            return
        self._jump_to_search_target(target)

    def _open_place_search_dialog(self) -> None:
        dialog = PlaceSearchDialog(self._search_place_jump_targets, self)
        if dialog.exec() == 0:
            return
        target = dialog.selected_target()
        if target is None:
            return
        self._jump_to_search_target(target)

    def _jump_to_search_target(self, target: SearchJumpTarget) -> None:
        target_kind = getattr(target, "kind", "star")
        if target_kind == "satellite":
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
            projection = project_place_target_to_altaz(
                observer_latitude_deg=self.viewer_data.location[0],
                observer_longitude_deg=self.viewer_data.location[1],
                observer_height_m=getattr(self.viewer_data, "observer_height_m", 1.7),
                target_latitude_deg=target.latitude_deg,
                target_longitude_deg=target.longitude_deg,
            )
            target_alt = float(projection.alt_deg)
            target_az = float(projection.az_deg) % 360.0
        else:
            observer_height_m = getattr(self.viewer_data, "observer_height_m", 1.7)
            alt, az = radec_to_altaz(
                target.ra_hours,
                target.dec_deg,
                self.viewer_data.location[0],
                self.viewer_data.location[1],
                observer_height_m,
                self._current_time_obj(),
            )
            target_alt = float(alt)
            target_az = float(az) % 360.0
        new_alt = max(0.0, min(90.0, target_alt))
        new_az = target_az
        self.viewer_data.view_center = (new_alt, new_az)
        self._sync_view_altitude_actions()

        self.state.jump_highlight_name = target.label
        self.state.jump_highlight_altaz = (target_alt, target_az)
        self.state.jump_highlight_until_ms = (time.monotonic() * 1000.0) + 3000.0

        self._begin_interaction_mode()
        self.request_sky_data_update()
        self.request_client_update()

    def _search_place_jump_targets(self, query: str) -> list[SearchJumpTarget]:
        candidates = search_place_candidates(query)
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

    def _find_satellite_jump_altaz(self, object_key: str) -> tuple[float, float] | None:
        records_by_group = dict(
            getattr(self.satellite_state, "records_by_group", None) or {}
        )
        enabled_groups = tuple(
            getattr(self, "_enabled_satellite_groups", (SATELLITE_ISS_CACHE_KEY,))
        )
        if not records_by_group:
            records_by_group = self._load_cached_satellite_records(enabled_groups)
        if not records_by_group:
            return None
        observer_height_m = getattr(self.viewer_data, "observer_height_m", 1.7)
        altaz = find_satellite_altaz(
            records_by_group,
            object_key=object_key,
            observer_lat=self.viewer_data.location[0],
            observer_lon=self.viewer_data.location[1],
            observer_height_m=observer_height_m,
            time_obj=self._current_time_obj(),
        )
        if altaz is not None:
            return altaz
        merged_records = dict(records_by_group)
        merged_records.update(self._load_cached_satellite_records(enabled_groups))
        if merged_records == records_by_group:
            return None
        return find_satellite_altaz(
            merged_records,
            object_key=object_key,
            observer_lat=self.viewer_data.location[0],
            observer_lon=self.viewer_data.location[1],
            observer_height_m=observer_height_m,
            time_obj=self._current_time_obj(),
        )

    def _load_cached_satellite_records(
        self, enabled_groups: tuple[str, ...]
    ) -> dict[str, list[dict[str, object]]]:
        records_by_group: dict[str, list[dict[str, object]]] = {}
        for group_key in enabled_groups:
            cached = load_satellite_cache(group_key)
            if cached is None:
                continue
            records_by_group[str(group_key)] = list(cached.records)
        return records_by_group

    def _begin_shutdown(self) -> None:
        """Stop scheduling new background work while the app is closing."""
        if self._is_shutting_down:
            return
        self._is_shutting_down = True
        self._sky_worker.shutdown()
        if self._cloud_controller is not None:
            self._cloud_controller.shutdown()
        if self._satellite_controller is not None:
            self._satellite_controller.shutdown()
        if self._aircraft_controller is not None:
            self._aircraft_controller.shutdown()
        if self._terrain_horizon_controller is not None:
            self._terrain_horizon_controller.shutdown()
        if self._urban_outline_controller is not None:
            self._urban_outline_controller.shutdown()
        if self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.stop()
        if self._asterism_check_timer.isActive():
            self._asterism_check_timer.stop()
        if self._cloud_update_timer.isActive():
            self._cloud_update_timer.stop()
        if self._satellite_update_timer.isActive():
            self._satellite_update_timer.stop()
        if self._overlay_projection_timer.isActive():
            self._overlay_projection_timer.stop()
        if self._aircraft_update_timer.isActive():
            self._aircraft_update_timer.stop()
        if self._interaction_idle_timer.isActive():
            self._interaction_idle_timer.stop()
        if self._viewport_interaction_idle_timer.isActive():
            self._viewport_interaction_idle_timer.stop()

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

    def _satellite_layer_enabled(self) -> bool:
        return self._satellite_toggle_supported and self.satellite_opacity > 0.0

    def _stop_satellite_timers(self) -> None:
        if self._satellite_update_timer.isActive():
            self._satellite_update_timer.stop()
        self._sync_overlay_projection_timer()

    def _target_time_utc(self) -> datetime:
        return target_time_utc_from_delta(self.delta_t)

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
        self._satellite_update_timer.start(max(0, interval_ms))

    def _schedule_satellite_retry_after_failure(self) -> None:
        if not self._satellite_layer_enabled() or self._is_shutting_down:
            return
        self._satellite_update_timer.start(SATELLITE_FAILURE_RETRY_SECONDS * 1000)

    def _on_satellite_refresh_timer(self) -> None:
        if not self._satellite_layer_enabled():
            return
        started = self.start_background_satellite_update(reason="timer")
        if not started:
            self._schedule_next_satellite_refresh()

    def _on_overlay_projection_timer(self) -> None:
        if self._satellite_layer_enabled():
            self.refresh_projected_satellite_overlay()
        if self._aircraft_layer_enabled():
            self.refresh_projected_aircraft_overlay()

    def _sync_overlay_projection_timer(self) -> None:
        if self._is_shutting_down:
            if self._overlay_projection_timer.isActive():
                self._overlay_projection_timer.stop()
            return
        any_enabled = self._satellite_layer_enabled() or self._aircraft_layer_enabled()
        if any_enabled:
            if not self._overlay_projection_timer.isActive():
                self._overlay_projection_timer.start()
        elif self._overlay_projection_timer.isActive():
            self._overlay_projection_timer.stop()

    def _enable_satellite_layer(self, *, reason: str) -> None:
        if not self._satellite_layer_enabled():
            return
        self._sync_overlay_projection_timer()
        if (
            self.satellite_state.records_by_group
            and self.satellite_state.element_epoch_utc is not None
        ):
            self.refresh_projected_satellite_overlay()
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
        if self._aircraft_update_timer.isActive():
            self._aircraft_update_timer.stop()
        self._sync_overlay_projection_timer()

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
        self._aircraft_update_timer.start(max(0, interval_ms))

    def _on_aircraft_refresh_timer(self) -> None:
        if not self._aircraft_layer_enabled():
            return
        started = self.start_background_aircraft_update(reason="timer")
        if not started:
            self._schedule_next_aircraft_refresh()

    def _enable_aircraft_layer(self, *, reason: str) -> None:
        if not self._aircraft_layer_enabled():
            return
        self._sync_overlay_projection_timer()
        if (
            self.aircraft_state.snapshots
            and self.aircraft_state.last_success_utc is not None
        ):
            self.refresh_projected_aircraft_overlay()
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

    def toggle_enlarge_moon(self) -> None:
        self.enlarge_moon = not self.enlarge_moon
        if (
            self._action_enlarge_moon is not None
            and self._action_enlarge_moon.isChecked() != self.enlarge_moon
        ):
            self._action_enlarge_moon.setChecked(self.enlarge_moon)
        self.request_client_update()

    def toggle_clouds(self) -> None:
        if (
            not self._cloud_toggle_supported
            or self._clouddisc is None
            or not self._cloud_gui_allowed
        ):
            if self._action_toggle_clouds is not None:
                self._action_toggle_clouds.setChecked(False)
            return

        enable_clouds = self.cloud_disc_alpha <= 0.0
        self.cloud_disc_alpha = self._cloud_alpha_when_enabled if enable_clouds else 0.0
        if (
            self._action_toggle_clouds is not None
            and self._action_toggle_clouds.isChecked() != enable_clouds
        ):
            self._action_toggle_clouds.setChecked(enable_clouds)

        if enable_clouds:
            self.start_background_cloud_update(reason="toggle-on")
            if (
                self._sky_data_update_timer.isActive()
                and not self._cloud_update_timer.isActive()
            ):
                self._cloud_update_timer.start()
        elif self._cloud_update_timer.isActive():
            self._cloud_update_timer.stop()

        self.request_client_update()

    def toggle_satellites(self) -> None:
        if not self._satellite_toggle_supported or not self._satellite_gui_allowed:
            if self._action_toggle_satellites is not None:
                self._action_toggle_satellites.setChecked(self.satellite_opacity > 0.0)
            return

        enable_satellites = self.satellite_opacity <= 0.0
        self.satellite_opacity = (
            self._satellite_opacity_when_enabled if enable_satellites else 0.0
        )
        if (
            self._action_toggle_satellites is not None
            and self._action_toggle_satellites.isChecked() != enable_satellites
        ):
            self._action_toggle_satellites.setChecked(enable_satellites)

        if enable_satellites:
            self._enable_satellite_layer(reason="toggle-on")
        else:
            self._stop_satellite_timers()

        self.request_client_update()

    def toggle_aircraft(self) -> None:
        if not self._aircraft_toggle_supported or not self._aircraft_gui_allowed:
            if self._action_toggle_aircraft is not None:
                self._action_toggle_aircraft.setChecked(self.aircraft_opacity > 0.0)
            return

        enable_aircraft = self.aircraft_opacity <= 0.0
        self.aircraft_opacity = (
            self._aircraft_opacity_when_enabled if enable_aircraft else 0.0
        )
        if (
            self._action_toggle_aircraft is not None
            and self._action_toggle_aircraft.isChecked() != enable_aircraft
        ):
            self._action_toggle_aircraft.setChecked(enable_aircraft)

        if enable_aircraft:
            self._enable_aircraft_layer(reason="toggle-on")
        else:
            self._stop_aircraft_timers()

        self.request_client_update()

    def toggle_dso(self) -> None:
        if self.dso_catalog_np is None:
            self.show_dso = False
            if self._action_toggle_dso is not None:
                self._action_toggle_dso.setChecked(False)
            return
        self.show_dso = not self.show_dso
        if (
            self._action_toggle_dso is not None
            and self._action_toggle_dso.isChecked() != self.show_dso
        ):
            self._action_toggle_dso.setChecked(self.show_dso)
        self.request_client_update()

    def toggle_asterisms(self) -> None:
        self.show_asterisms = not self.show_asterisms
        if (
            self._action_toggle_asterisms is not None
            and self._action_toggle_asterisms.isChecked() != self.show_asterisms
        ):
            self._action_toggle_asterisms.setChecked(self.show_asterisms)
        self.request_client_update()

    def toggle_guidelines(self) -> None:
        self.show_guidelines = not self.show_guidelines
        if (
            self._action_toggle_guidelines is not None
            and self._action_toggle_guidelines.isChecked() != self.show_guidelines
        ):
            self._action_toggle_guidelines.setChecked(self.show_guidelines)
        self.request_client_update()

    def toggle_overlay_info(self) -> None:
        self.show_overlay_info = not self.show_overlay_info
        if (
            self._action_toggle_overlay_info is not None
            and self._action_toggle_overlay_info.isChecked() != self.show_overlay_info
        ):
            self._action_toggle_overlay_info.setChecked(self.show_overlay_info)
        self.request_client_update()

    def toggle_sky_disc(self) -> None:
        if not self._sky_disc_gui_allowed:
            if self._action_toggle_sky_disc is not None:
                self._action_toggle_sky_disc.setChecked(self.sky_disc_alpha > 0.0)
            return
        enable_sky_disc_gradient = self.sky_disc_alpha <= 0.0
        self.sky_disc_alpha = (
            self._sky_disc_alpha_when_enabled if enable_sky_disc_gradient else 0.0
        )
        if (
            self._action_toggle_sky_disc is not None
            and self._action_toggle_sky_disc.isChecked() != enable_sky_disc_gradient
        ):
            self._action_toggle_sky_disc.setChecked(enable_sky_disc_gradient)
        self._compositor.invalidate()
        self.request_sky_data_update()
        self.request_client_update()

    def toggle_terrain_horizon(self) -> None:
        if not self._terrain_horizon_gui_allowed:
            if self._action_toggle_terrain_horizon is not None:
                self._action_toggle_terrain_horizon.setChecked(
                    self.terrain_horizon_opacity > 0.0
                )
            return

        enable_terrain = self.terrain_horizon_opacity <= 0.0
        self.terrain_horizon_opacity = (
            self._terrain_horizon_opacity_when_enabled if enable_terrain else 0.0
        )
        if (
            self._action_toggle_terrain_horizon is not None
            and self._action_toggle_terrain_horizon.isChecked() != enable_terrain
        ):
            self._action_toggle_terrain_horizon.setChecked(enable_terrain)
        self._compositor.invalidate()
        if enable_terrain:
            self.start_background_terrain_horizon_update(reason="toggle-on")
        self.request_client_update()

    def toggle_urban_outline(self) -> None:
        if not self._urban_outline_gui_allowed:
            if self._action_toggle_urban_outline is not None:
                self._action_toggle_urban_outline.setChecked(
                    self.urban_outline_opacity > 0.0
                )
            return
        enable_urban_outline = self.urban_outline_opacity <= 0.0
        self.urban_outline_opacity = (
            self._urban_outline_opacity_when_enabled if enable_urban_outline else 0.0
        )
        self.show_urban_outline_layer = self.urban_outline_opacity > 0.0
        if (
            self._action_toggle_urban_outline is not None
            and self._action_toggle_urban_outline.isChecked() != enable_urban_outline
        ):
            self._action_toggle_urban_outline.setChecked(enable_urban_outline)
        if enable_urban_outline:
            self.start_background_urban_outline_update(reason="toggle-on")
        self.request_client_update()

    def toggle_fullscreen(self) -> None:
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def _handle_client_leave(self, event: QEvent) -> None:
        self.state.mouse_pos = None
        self.request_client_update()
        event.accept()

    def closeEvent(self, event: QCloseEvent) -> None:
        g = self.geometry()
        save_last_window_geometry(g.x(), g.y(), g.width(), g.height())
        self._begin_shutdown()
        super().closeEvent(event)

    def _handle_client_mouse_move(self, event: QMouseEvent) -> None:
        self.state.mouse_pos = event.pos()
        self.request_client_update()
        event.accept()

    def _rotate_view(
        self,
        d_alt: float = 0.0,
        d_az: float = 0.0,
        *,
        interactive_viewport: bool = False,
    ) -> None:
        if interactive_viewport:
            self._begin_viewport_interaction_mode()
        else:
            self._begin_interaction_mode()
        alt, az = self.viewer_data.view_center
        new_alt = max(0.0, min(90.0, alt + d_alt))
        new_az = (az + d_az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)
        self.state.render_view_center = (new_alt, new_az)
        self._sync_view_altitude_actions()
        if interactive_viewport:
            self._update_viewport_interaction_stars()
            self.request_client_update()
            return
        self.request_sky_data_update()

    def _handle_client_key_press(self, event: QKeyEvent) -> None:
        key = event.key()

        # --- View Control ---
        if key == Qt.Key.Key_Left:
            self._rotate_view(d_az=-self.state.rotation_step, interactive_viewport=True)
            event.accept()
        elif key == Qt.Key.Key_Right:
            self._rotate_view(d_az=self.state.rotation_step, interactive_viewport=True)
            event.accept()
        elif key == Qt.Key.Key_Up:
            self._rotate_view(d_alt=self.state.rotation_step, interactive_viewport=True)
            event.accept()
        elif key == Qt.Key.Key_Down:
            self._rotate_view(
                d_alt=-self.state.rotation_step, interactive_viewport=True
            )
            event.accept()

        # --- Toggles ---
        elif key == Qt.Key.Key_M:
            self.toggle_enlarge_moon()
            event.accept()
        elif key == Qt.Key.Key_C:
            self.toggle_clouds()
            event.accept()
        elif key == Qt.Key.Key_I:
            self.toggle_aircraft()
            event.accept()
        elif key == Qt.Key.Key_T:
            self.toggle_terrain_horizon()
            event.accept()
        elif key == Qt.Key.Key_U:
            self.toggle_urban_outline()
            event.accept()
        elif key == Qt.Key.Key_D:
            self.toggle_dso()
            event.accept()
        elif key == Qt.Key.Key_A:
            self.toggle_asterisms()
            event.accept()
        elif key == Qt.Key.Key_G:
            self.toggle_guidelines()
            event.accept()
        elif key == Qt.Key.Key_S:
            self.toggle_sky_disc()
            event.accept()

        # --- Window Control ---
        elif key == Qt.Key.Key_F11:
            self.toggle_fullscreen()
            event.accept()
        elif key == Qt.Key.Key_Escape:
            if self.isFullScreen():
                self.showNormal()
            event.accept()
        elif key == Qt.Key.Key_Q:
            QApplication.quit()
            event.accept()
        else:
            super().keyPressEvent(event)


class FramelessSkyWindow(SkyWindowCoreMixin, DraggableWindow):
    """Frameless host window using the custom chrome wrapper."""

    FRAMELESS_WINDOW = True

    def _install_window_host(self) -> None:
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.FramelessWindowHint)
        self._frameless_frame = FramelessWindowFrame(self, self._client_widget)
        self.setCentralWidget(self._frameless_frame)
        self.menu_button = self._frameless_frame.menu_button
        self.size_grip = self._frameless_frame.size_grip
        self.add_drag_target(self._client_widget)
        self.add_drag_exclusions(self._frameless_frame.overlay_widgets())


class StandardSkyWindow(SkyWindowCoreMixin, QMainWindow):
    """Standard decorated host window using the platform title bar and menu bar."""

    FRAMELESS_WINDOW = False

    def _install_window_host(self) -> None:
        self.setCentralWidget(self._client_widget)
        menu_bar = self.menuBar()
        menu_bar.addMenu(self.file_menu)
        menu_bar.addMenu(self.search_menu)
        menu_bar.addMenu(self.display_menu)
        menu_bar.addMenu(self.observer_view_menu)


SkyWindow = FramelessSkyWindow
