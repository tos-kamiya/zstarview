# -*- coding: utf-8 -*-
"""
The main window of the ZStarView application.

This module defines the `SkyWindow` class, which is the primary user interface
for the application. It handles rendering the celestial objects, sky background,
clouds, and all user interactions like rotation, zooming, and object highlighting.
"""
import logging
import math
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
    QResizeEvent,
)
from PySide6.QtWidgets import QApplication, QMenu, QPushButton, QSizeGrip

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
from ..paths import (
    APP_ICON_FILE,
    APP_DISPLAY_NAME,
    CACHE_PATH,
    CLOUD_UPDATE_INTERVAL,
    GUI_BUTTON_SIZE,
    GUI_MENU_TEXT_COLOR,
    TEXT_FONT_PATH,
    TEXT_FONT_SIZE,
    STATUS_LINE_FONT_SIZE,
    WINDOW_HEIGHT,
    WINDOW_WIDTH,
    CLOUD_MISSING_TINT_RGBA,
)
from ..config import load_last_window_geometry, save_last_window_geometry
from ..render import draw as render_draw
from ..types import ViewerData
from .draggable_window import DraggableWindow
from .composite import SkyCompositorCache
from .cloud_state import CloudImageState
from .cloud_controller import CloudController
from .terrain_state import TerrainHorizonState
from .terrain_controller import TerrainHorizonController
from .famous_star_dialog import NamedStarJumpDialog
from .famous_star_search_dialog import NamedStarSearchDialog
from .famous_star_shortcuts import NamedStarShortcut, SearchJumpTarget
from ..asterisms import ASTERISM_KEYS_BY_SOURCE_ID
from ..urban_skyline_profiles import resolve_urban_skyline_profile_for_city_name
from .sky_worker import SkyDataWorker
from .window_inputs import PreparedWindowCatalogs
from .window_inputs import SkyWindowRuntimeOptions, SkyWindowUserOptions
from .window_render import SkyWindowRenderMixin
from .window_state import SkyWindowState
from .window_updates import SkyWindowUpdatesMixin

logger = logging.getLogger(__name__)


WindowGeometryArg = Union[str, Tuple[int, int, int, int]]


def compute_star_render_surface_size(
    width_px: int,
    height_px: int,
    disc_width_px: int,
    expected_width_px: int,
) -> tuple[int, int]:
    """Return internal star-render surface size.

    - No downsample when `disc_width_px <= expected_width_px`.
    - Above threshold, effective rendered disc width follows:
      `expected_width_px * sqrt(disc_width_px / expected_width_px)`.
    """
    w = max(1, int(width_px))
    h = max(1, int(height_px))
    disc_w = max(1, int(disc_width_px))
    base = max(1, int(expected_width_px))
    if disc_w <= base:
        return (w, h)
    rendered_disc_w = float(base) * math.sqrt(float(disc_w) / float(base))
    scale = rendered_disc_w / float(disc_w)
    return (
        max(1, int(round(w * scale))),
        max(1, int(round(h * scale))),
    )


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
        available_rect = primary.availableGeometry() if primary is not None else screens[0].availableGeometry()

    width = min(width, max(1, available_rect.width()))
    height = min(height, max(1, available_rect.height()))

    max_x = available_rect.right() - width + 1
    max_y = available_rect.bottom() - height + 1
    x = min(max(int(x), available_rect.left()), max_x)
    y = min(max(int(y), available_rect.top()), max_y)
    return x, y, width, height


class SkyWindow(SkyWindowRenderMixin, SkyWindowUpdatesMixin, DraggableWindow):
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
        self.star_catalog_full_np = catalogs.star_catalog_full_np
        self.star_catalog_lod6_np = catalogs.star_catalog_lod6_np
        self.dso_catalog_np = catalogs.dso_catalog_np
        self.show_dso: bool = self.dso_catalog_np is not None
        if user_options.show_dso_initial is not None:
            self.show_dso = bool(user_options.show_dso_initial) and self.dso_catalog_np is not None
        self.show_asterisms: bool = (
            True if user_options.show_asterisms_initial is None else bool(user_options.show_asterisms_initial)
        )
        self._named_stars_by_band = catalogs.named_stars_by_band
        self._named_stars_search_all = catalogs.named_stars_search_all
        self.delta_t = runtime_options.delta_t
        self.sky_disc_alpha = user_options.sky_disc_alpha
        self._sky_disc_alpha_when_enabled = user_options.sky_disc_alpha if user_options.sky_disc_alpha > 0.0 else 0.3
        self.terrain_horizon_opacity = user_options.terrain_horizon_opacity
        self.ground_tint_opacity = user_options.ground_tint_opacity
        self._terrain_horizon_opacity_when_enabled = (
            user_options.terrain_horizon_opacity if user_options.terrain_horizon_opacity > 0.0 else 0.25
        )
        self._terrain_horizon_gui_allowed = bool(user_options.terrain_horizon_gui_allowed)
        self.enlarge_moon = user_options.enlarge_moon
        self.star_base_radius = user_options.star_base_radius
        self.vmag_limit = user_options.vmag_limit
        self.sky_update_interval = runtime_options.sky_update_interval
        self.visual_preset = user_options.visual_preset
        self.star_visibility_boost = user_options.star_visibility_boost
        self._star_render_expected_width = runtime_options.star_render_expected_width
        self._cloud_toggle_supported = runtime_options.delta_t.total_seconds() == 0.0

        # Cloud opacity is disabled if we are looking at a time-shifted view,
        # as we can only fetch current cloud data.
        requested_cloud_alpha = user_options.cloud_disc_alpha
        self._cloud_alpha_when_enabled = requested_cloud_alpha if requested_cloud_alpha > 0.0 else 0.2
        self.cloud_disc_alpha: float = requested_cloud_alpha
        if not self._cloud_toggle_supported:
            self.cloud_disc_alpha = 0.0

        # --- Viewer and Window Setup ---
        self.viewer_data = viewer_data
        self.state = SkyWindowState(
            render_view_center=tuple(self.viewer_data.view_center),
            urban_skyline_profiles=resolve_urban_skyline_profile_for_city_name(
                self.viewer_data.city_name
            ),
        )
        self.setWindowTitle(f"{APP_DISPLAY_NAME} - {self.viewer_data.city_name}")
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
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.FramelessWindowHint)
        requested_geometry: Optional[Tuple[int, int, int, int]] = None
        if runtime_options.window_geometry_arg == "restore":
            requested_geometry = load_last_window_geometry()
        elif isinstance(runtime_options.window_geometry_arg, tuple):
            requested_geometry = runtime_options.window_geometry_arg
        if requested_geometry is not None:
            initial_x, initial_y, initial_width, initial_height = requested_geometry
        initial_x, initial_y, initial_width, initial_height = _clamp_window_geometry_to_screen(
            initial_x,
            initial_y,
            initial_width,
            initial_height,
            min_width=min_width,
            min_height=min_height,
        )
        self.setGeometry(initial_x, initial_y, initial_width, initial_height)
        self.setMouseTracking(True)
        # --- UI Widgets ---
        self.size_grip = QSizeGrip(self)
        self.size_grip.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.size_grip.raise_()
        self._action_enlarge_moon: Optional[QAction] = None
        self._action_toggle_clouds: Optional[QAction] = None
        self._action_toggle_terrain_horizon: Optional[QAction] = None
        self._action_toggle_dso: Optional[QAction] = None
        self._action_toggle_asterisms: Optional[QAction] = None
        self._action_toggle_sky_disc: Optional[QAction] = None
        self._action_raise_view: Optional[QAction] = None
        self._action_lower_view: Optional[QAction] = None
        self._add_hamburger_menu()
        self.add_drag_exclusions([self.menu_button, self.size_grip])

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
        self.terrain_horizon_state = TerrainHorizonState()
        self._cloud_controller: Optional[CloudController] = None
        self._terrain_horizon_controller: Optional[TerrainHorizonController] = None
        self._cloud_update_timer = QTimer(self)
        self._cloud_update_timer.setInterval(CLOUD_UPDATE_INTERVAL * 1000)
        self._cloud_update_timer.timeout.connect(lambda: self.start_background_cloud_update(reason="timer"))

        # --- CloudDisc Service Initialization ---
        self._clouddisc: Optional[CloudDisc] = None
        clouddisc_config = CloudDiscConfig(
            cache_dir=CACHE_PATH,
            sat_priority=("AUTO",),
            bt_warm_k=310.0,
            bt_cold_k=190.0,
            alt_min_deg=0.0,
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
        if self._action_toggle_clouds is not None:
            self._action_toggle_clouds.setEnabled(self._cloud_toggle_supported and self._clouddisc is not None)

        terrain_cache_dir = Path(CACHE_PATH) / "copernicus-dem"
        self._terrain_horizon_controller = TerrainHorizonController(
            cache_dir=terrain_cache_dir,
            parent=self,
        )
        self._terrain_horizon_controller.terrain_started.connect(self._on_terrain_horizon_started)
        self._terrain_horizon_controller.terrain_ready.connect(self._on_terrain_horizon_ready)
        self._terrain_horizon_controller.terrain_failed.connect(self._on_terrain_horizon_failed)
        if self._action_toggle_terrain_horizon is not None:
            self._action_toggle_terrain_horizon.setEnabled(self._terrain_horizon_gui_allowed)

        # --- Composition Cache (moved to dedicated class) ---
        target_stripes, width_factor = runtime_options.cloud_stripe_style
        missing_tint_alpha = int(round(255.0 * runtime_options.cloud_missing_tint_opacity))
        missing_tint_rgba = (
            int(CLOUD_MISSING_TINT_RGBA[0]),
            int(CLOUD_MISSING_TINT_RGBA[1]),
            int(CLOUD_MISSING_TINT_RGBA[2]),
            missing_tint_alpha,
        )
        self._compositor = SkyCompositorCache(
            cloud_target_stripes=int(target_stripes),
            cloud_stripe_width_factor=float(width_factor),
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

    def _setup_update_infrastructure(self) -> None:
        """Initialize timers, worker, and signal wiring for background updates."""
        self._sky_worker = SkyDataWorker(self)
        self._sky_worker.data_ready.connect(self._on_sky_data_calculated)
        self.cloud_repaint_requested.connect(self.update)

        app = QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(self._begin_shutdown)

        self._sky_data_update_timer = QTimer(self)
        self._sky_data_update_timer.timeout.connect(self.start_background_sky_data_update)

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
        self._viewport_interaction_idle_timer.setInterval(self.state.viewport_interaction_idle_ms)
        self._viewport_interaction_idle_timer.timeout.connect(self._end_viewport_interaction_mode)

    def _add_hamburger_menu(self) -> None:
        """Adds a hamburger menu button and its corresponding actions."""
        self.menu_button = QPushButton("☰", self)
        self.menu_button.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.menu_button.setStyleSheet(
            "QPushButton { border: none; font-size: 18px; background-color: transparent; color: " + "#%02x%02x%02x" % GUI_MENU_TEXT_COLOR + "; }"
            "QPushButton:hover { color: white; }"
            "QPushButton:menu-indicator { image: none; }"
        )
        self.menu_button.clicked.connect(self.show_menu)

        self.menu = QMenu(self)
        rotate_left = self.menu.addAction(f"Rotate Left (-{self.state.rotation_step:.0f}°)")
        rotate_left.triggered.connect(lambda: self._rotate_view(d_az=-self.state.rotation_step))
        rotate_right = self.menu.addAction(f"Rotate Right (+{self.state.rotation_step:.0f}°)")
        rotate_right.triggered.connect(lambda: self._rotate_view(d_az=+self.state.rotation_step))
        raise_view = self.menu.addAction(f"Raise View (+{self.state.rotation_step:.0f}° alt)")
        raise_view.triggered.connect(lambda: self._rotate_view(d_alt=+self.state.rotation_step))
        self._action_raise_view = raise_view
        lower_view = self.menu.addAction(f"Lower View (-{self.state.rotation_step:.0f}° alt)")
        lower_view.triggered.connect(lambda: self._rotate_view(d_alt=-self.state.rotation_step))
        self._action_lower_view = lower_view

        self.menu.addSeparator()
        jump_named_star = self.menu.addAction("Jump to Named Star...")
        jump_named_star.setShortcut(QKeySequence("Ctrl+J"))
        jump_named_star.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        jump_named_star.triggered.connect(self._open_named_star_jump_dialog)
        self.addAction(jump_named_star)
        search_named_star = self.menu.addAction("Search Stars and Asterisms...")
        search_named_star.setShortcut(QKeySequence("Ctrl+F"))
        search_named_star.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        search_named_star.triggered.connect(self._open_named_star_search_dialog)
        self.addAction(search_named_star)

        self.menu.addSeparator()
        toggle_enlarge_moon_action = QAction("Enlarge Moon", self)
        toggle_enlarge_moon_action.setCheckable(True)
        toggle_enlarge_moon_action.setChecked(self.enlarge_moon)
        toggle_enlarge_moon_action.setShortcut(QKeySequence(Qt.Key.Key_M))
        toggle_enlarge_moon_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_enlarge_moon_action.triggered.connect(self.toggle_enlarge_moon)
        self.menu.addAction(toggle_enlarge_moon_action)
        self.addAction(toggle_enlarge_moon_action)
        self._action_enlarge_moon = toggle_enlarge_moon_action
        toggle_dso_action = QAction("DSO", self)
        toggle_dso_action.setCheckable(True)
        toggle_dso_action.setChecked(self.show_dso)
        toggle_dso_action.setEnabled(self.dso_catalog_np is not None)
        toggle_dso_action.setShortcut(QKeySequence(Qt.Key.Key_D))
        toggle_dso_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_dso_action.triggered.connect(self.toggle_dso)
        self.menu.addAction(toggle_dso_action)
        self.addAction(toggle_dso_action)
        self._action_toggle_dso = toggle_dso_action
        toggle_asterisms_action = QAction("Asterisms", self)
        toggle_asterisms_action.setCheckable(True)
        toggle_asterisms_action.setChecked(self.show_asterisms)
        toggle_asterisms_action.setShortcut(QKeySequence(Qt.Key.Key_A))
        toggle_asterisms_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_asterisms_action.triggered.connect(self.toggle_asterisms)
        self.menu.addAction(toggle_asterisms_action)
        self.addAction(toggle_asterisms_action)
        self._action_toggle_asterisms = toggle_asterisms_action

        self.menu.addSeparator()
        toggle_sky_disc_action = QAction("Sky Color Disc", self)
        toggle_sky_disc_action.setCheckable(True)
        toggle_sky_disc_action.setChecked(self.sky_disc_alpha > 0.0)
        toggle_sky_disc_action.setShortcut(QKeySequence(Qt.Key.Key_S))
        toggle_sky_disc_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_sky_disc_action.triggered.connect(self.toggle_sky_disc)
        self.menu.addAction(toggle_sky_disc_action)
        self.addAction(toggle_sky_disc_action)
        self._action_toggle_sky_disc = toggle_sky_disc_action
        toggle_clouds_action = QAction("Clouds", self)
        toggle_clouds_action.setCheckable(True)
        toggle_clouds_action.setChecked(self.cloud_disc_alpha > 0.0)
        toggle_clouds_action.setShortcut(QKeySequence(Qt.Key.Key_C))
        toggle_clouds_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_clouds_action.triggered.connect(self.toggle_clouds)
        self.menu.addAction(toggle_clouds_action)
        self.addAction(toggle_clouds_action)
        self._action_toggle_clouds = toggle_clouds_action
        toggle_terrain_action = QAction("Terrain Horizon", self)
        toggle_terrain_action.setCheckable(True)
        toggle_terrain_action.setChecked(self.terrain_horizon_opacity > 0.0)
        toggle_terrain_action.setShortcut(QKeySequence(Qt.Key.Key_T))
        toggle_terrain_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        toggle_terrain_action.triggered.connect(self.toggle_terrain_horizon)
        self.menu.addAction(toggle_terrain_action)
        self.addAction(toggle_terrain_action)
        self._action_toggle_terrain_horizon = toggle_terrain_action

        self.menu.addSeparator()
        fullscreen_action = self.menu.addAction("Fullscreen")
        fullscreen_action.setShortcut(QKeySequence(Qt.Key.Key_F11))
        fullscreen_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        fullscreen_action.triggered.connect(self.toggle_fullscreen)
        self.addAction(fullscreen_action)

        self.menu.addSeparator()
        exit_action = self.menu.addAction("Exit")
        exit_action.setShortcut(QKeySequence(Qt.Key.Key_Q))
        exit_action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
        exit_action.triggered.connect(QApplication.quit)
        self.addAction(exit_action)

        self.menu.addSeparator()
        version_action = self.menu.addAction(f"Version {__version__}")
        version_action.setEnabled(False)

    def _sync_view_altitude_actions(self) -> None:
        alt, _ = self.viewer_data.view_center
        if self._action_raise_view is not None:
            self._action_raise_view.setEnabled(float(alt) < 90.0)
        if self._action_lower_view is not None:
            self._action_lower_view.setEnabled(float(alt) > 0.0)

    def resizeEvent(self, event: QResizeEvent) -> None:
        self._begin_viewport_interaction_mode()
        grip_size = self.size_grip.size()
        self.size_grip.move(self.width() - grip_size.width(), self.height() - grip_size.height())

        button_size = self.menu_button.size()
        self.menu_button.move(self.width() - button_size.width() - 8, 8)

        # Invalidate the composition cache since the size has changed
        self._compositor.invalidate()

        super().resizeEvent(event)

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

    def _begin_viewport_interaction_mode(self) -> None:
        self.state.viewport_interaction_mode = True
        self._viewport_interaction_idle_timer.start()

    def _update_viewport_interaction_stars(self) -> None:
        if self.state.celestial_data is None:
            self.state.viewport_interaction_stars = None
            return
        stars, _location = calculate_visible_stars(
            self.star_catalog_lod6_np,
            self.viewer_data.location[0],
            self.viewer_data.location[1],
            self.viewer_data.observer_height_m,
            self._current_time_obj(),
            self.state.render_view_center,
            max_vmag=4.0,
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
        self.update()

    def show_menu(self) -> None:
        self._sync_view_altitude_actions()
        menu_pos = self.menu_button.mapToGlobal(QPoint(0, self.menu_button.height()))
        self.menu.exec(menu_pos)

    def _current_time_obj(self) -> astropy.time.Time:
        now = datetime.now(timezone.utc) + self.delta_t
        return astropy.time.Time(now)

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

    def _jump_to_search_target(self, target: SearchJumpTarget) -> None:
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
        self.update()

    def _jump_to_named_star(self, star: NamedStarShortcut) -> None:
        self._jump_to_search_target(
            SearchJumpTarget(
                label=star.name,
                ra_hours=star.ra_hours,
                dec_deg=star.dec_deg,
                kind="star",
                subtitle=f"Vmag {star.vmag:.2f}",
                sort_key=(star.vmag, star.name.casefold()),
            )
        )

    def _begin_shutdown(self) -> None:
        """Stop scheduling new background work while the app is closing."""
        if self._is_shutting_down:
            return
        self._is_shutting_down = True
        self._sky_worker.shutdown()
        if self._cloud_controller is not None:
            self._cloud_controller.shutdown()
        if self._terrain_horizon_controller is not None:
            self._terrain_horizon_controller.shutdown()
        if self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.stop()
        if self._asterism_check_timer.isActive():
            self._asterism_check_timer.stop()
        if self._cloud_update_timer.isActive():
            self._cloud_update_timer.stop()
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
        geometry = render_draw.get_screen_geometry(self.width(), self.height(), render_viewer.view_center[0])
        highlighted = render_draw.find_highlighted_object(
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
        self.update()

    def _predicted_cloud_satellite(self) -> str:
        lat, lon = self.viewer_data.location
        return pick_satellite(lat, lon, ("AUTO",))

    def toggle_enlarge_moon(self) -> None:
        self.enlarge_moon = not self.enlarge_moon
        if self._action_enlarge_moon is not None and self._action_enlarge_moon.isChecked() != self.enlarge_moon:
            self._action_enlarge_moon.setChecked(self.enlarge_moon)
        self.update()  # Redraw with the new setting

    def toggle_clouds(self) -> None:
        if not self._cloud_toggle_supported or self._clouddisc is None:
            if self._action_toggle_clouds is not None:
                self._action_toggle_clouds.setChecked(False)
            return

        enable_clouds = self.cloud_disc_alpha <= 0.0
        self.cloud_disc_alpha = self._cloud_alpha_when_enabled if enable_clouds else 0.0
        if self._action_toggle_clouds is not None and self._action_toggle_clouds.isChecked() != enable_clouds:
            self._action_toggle_clouds.setChecked(enable_clouds)

        if enable_clouds:
            self.start_background_cloud_update(reason="toggle-on")
            if self._sky_data_update_timer.isActive() and not self._cloud_update_timer.isActive():
                self._cloud_update_timer.start()
        elif self._cloud_update_timer.isActive():
            self._cloud_update_timer.stop()

        self.update()

    def toggle_dso(self) -> None:
        if self.dso_catalog_np is None:
            self.show_dso = False
            if self._action_toggle_dso is not None:
                self._action_toggle_dso.setChecked(False)
            return
        self.show_dso = not self.show_dso
        if self._action_toggle_dso is not None and self._action_toggle_dso.isChecked() != self.show_dso:
            self._action_toggle_dso.setChecked(self.show_dso)
        self.update()

    def toggle_asterisms(self) -> None:
        self.show_asterisms = not self.show_asterisms
        if self._action_toggle_asterisms is not None and self._action_toggle_asterisms.isChecked() != self.show_asterisms:
            self._action_toggle_asterisms.setChecked(self.show_asterisms)
        self.update()

    def toggle_sky_disc(self) -> None:
        enable_sky_disc_gradient = self.sky_disc_alpha <= 0.0
        self.sky_disc_alpha = self._sky_disc_alpha_when_enabled if enable_sky_disc_gradient else 0.0
        if (
            self._action_toggle_sky_disc is not None
            and self._action_toggle_sky_disc.isChecked() != enable_sky_disc_gradient
        ):
            self._action_toggle_sky_disc.setChecked(enable_sky_disc_gradient)
        self._compositor.invalidate()
        self.request_sky_data_update()
        self.update()

    def toggle_terrain_horizon(self) -> None:
        if not self._terrain_horizon_gui_allowed:
            if self._action_toggle_terrain_horizon is not None:
                self._action_toggle_terrain_horizon.setChecked(self.terrain_horizon_opacity > 0.0)
            return

        enable_terrain = self.terrain_horizon_opacity <= 0.0
        self.terrain_horizon_opacity = self._terrain_horizon_opacity_when_enabled if enable_terrain else 0.0
        if self._action_toggle_terrain_horizon is not None and self._action_toggle_terrain_horizon.isChecked() != enable_terrain:
            self._action_toggle_terrain_horizon.setChecked(enable_terrain)
        self._compositor.invalidate()
        if enable_terrain:
            self.start_background_terrain_horizon_update(reason="toggle-on")
        self.update()

    def toggle_fullscreen(self) -> None:
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def leaveEvent(self, event: QEvent) -> None:
        self.state.mouse_pos = None
        self.update()
        event.accept()

    def closeEvent(self, event: QCloseEvent) -> None:
        g = self.geometry()
        save_last_window_geometry(g.x(), g.y(), g.width(), g.height())
        self._begin_shutdown()
        super().closeEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        self.state.mouse_pos = event.pos()
        self.update()  # Trigger a repaint to show hover effects
        # We accept the event to prevent it from propagating further.
        # This is why the manual drag in DraggableWindow does not work.
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
            self.update()
            return
        self.request_sky_data_update()

    def keyPressEvent(self, event: QKeyEvent) -> None:
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
            self._rotate_view(d_alt=-self.state.rotation_step, interactive_viewport=True)
            event.accept()

        # --- Toggles ---
        elif key == Qt.Key.Key_M:
            self.toggle_enlarge_moon()
            event.accept()
        elif key == Qt.Key.Key_C:
            self.toggle_clouds()
            event.accept()
        elif key == Qt.Key.Key_T:
            self.toggle_terrain_horizon()
            event.accept()
        elif key == Qt.Key.Key_D:
            self.toggle_dso()
            event.accept()
        elif key == Qt.Key.Key_A:
            self.toggle_asterisms()
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
