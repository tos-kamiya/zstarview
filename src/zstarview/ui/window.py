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
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, Optional, Tuple, Union

import astropy.time
import polars as pl
from PySide6.QtCore import QEvent, QPoint, QPointF, QRect, QRectF, Qt, QTimer, Signal
from PySide6.QtGui import (
    QAction,
    QCloseEvent,
    QFont,
    QFontDatabase,
    QGuiApplication,
    QIcon,
    QImage,
    QKeyEvent,
    QMouseEvent,
    QPainter,
    QPaintEvent,
    QResizeEvent,
)
from PySide6.QtWidgets import QApplication, QMenu, QPushButton, QSizeGrip

from ..__about__ import __version__
from ..astro import (
    DeepSkyCatalogArrays,
    StarCatalogArrays,
    prepare_deep_sky_catalog_arrays,
    prepare_star_catalog_arrays,
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
from ..types import CelestialData, ViewerData
from .draggable_window import DraggableWindow
from .composite import SkyCompositorCache
from .cloud_state import CloudImageState
from .cloud_controller import CloudController
from .famous_star_dialog import NamedStarJumpDialog
from .famous_star_search_dialog import NamedStarSearchDialog
from .famous_star_shortcuts import (
    NamedStarShortcut,
    SearchJumpTarget,
    build_search_jump_targets,
    build_named_star_shortcuts,
)
from ..asterisms import ASTERISM_KEYS_BY_SOURCE_ID
from .sky_worker import SkyDataWorker
from .window_render import SkyWindowRenderMixin
from .window_state import SkyWindowState
from .window_updates import SkyWindowUpdatesMixin

logger = logging.getLogger(__name__)


DEBUG_ECLIPSES = True
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
        city_name: str,
        city_data: Tuple[float, float, str],
        view_center: Tuple[float, float],
        star_catalog: pl.DataFrame,
        dso_catalog: Optional[pl.DataFrame] = None,
        delta_t: timedelta = timedelta(days=0, hours=0),
        sky_disc_alpha: float = 0.3,
        cloud_disc_alpha: float = 0.6,
        enlarge_moon: bool = False,
        star_base_radius: float = 4.0,
        vmag_limit: float = 6.0,
        sky_update_interval: int = 60,  # sec
        visual_preset: str = "night",
        star_visibility_boost: float = 1.0,
        vmag_brightness_scale: float = -0.39,
        cloud_stripe_style: Tuple[int, float] = (50, 0.2),
        cloud_missing_tint_opacity: float = float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0,
        star_render_expected_width: int = 600,
        window_geometry_arg: Optional[WindowGeometryArg] = None,
        show_dso_initial: Optional[bool] = None,
        show_asterisms_initial: Optional[bool] = None,
    ) -> None:
        """
        Initializes the SkyWindow.

        Args:
            city_name: The name of the observer's city.
            city_data: A tuple of (latitude, longitude, timezone_name).
            view_center: The initial view center (altitude, azimuth) in degrees.
            star_catalog: A Polars DataFrame containing the star catalog.
            delta_t: Time offset from the current time.
            sky_disc_alpha: Opacity of the daytime sky color overlay.
            cloud_disc_alpha: Opacity of the cloud layer.
            enlarge_moon: Whether to draw the Moon larger than its true scale.
            star_base_radius: Base size for 2nd-magnitude stars.
            vmag_limit: The faintest star magnitude to display.
            visual_preset: UI visual preset name.
            star_visibility_boost: Multiplier for star visibility on bright presets.
            vmag_brightness_scale: Slope applied to magnitude for color intensity (negative number).
        """
        super().__init__()
        self.star_catalog = star_catalog
        self.vmag_brightness_scale = vmag_brightness_scale
        self.star_catalog_full_np: StarCatalogArrays = prepare_star_catalog_arrays(
            star_catalog, vmag_brightness_scale=self.vmag_brightness_scale
        )
        self.star_catalog_lod6_np: StarCatalogArrays = prepare_star_catalog_arrays(
            star_catalog,
            max_vmag=6.0,
            vmag_brightness_scale=self.vmag_brightness_scale,
        )
        self.dso_catalog_np: Optional[DeepSkyCatalogArrays]
        if dso_catalog is None:
            self.dso_catalog_np = None
        else:
            self.dso_catalog_np = prepare_deep_sky_catalog_arrays(dso_catalog)
        self.show_dso: bool = self.dso_catalog_np is not None
        if show_dso_initial is not None:
            self.show_dso = bool(show_dso_initial) and self.dso_catalog_np is not None
        self.show_asterisms: bool = True if show_asterisms_initial is None else bool(show_asterisms_initial)
        self._named_stars_by_band = build_named_star_shortcuts(star_catalog, max_vmag=2.0)
        self._named_stars_search_all = build_search_jump_targets(star_catalog)
        self.delta_t = delta_t
        self.sky_disc_alpha = sky_disc_alpha
        self.enlarge_moon = enlarge_moon
        self.star_base_radius = star_base_radius
        self.vmag_limit = vmag_limit
        self.sky_update_interval = sky_update_interval
        self.visual_preset = visual_preset
        self.star_visibility_boost = star_visibility_boost
        self._star_render_expected_width = max(1, int(star_render_expected_width))
        self._cloud_toggle_supported = delta_t.total_seconds() == 0.0

        # Cloud opacity is disabled if we are looking at a time-shifted view,
        # as we can only fetch current cloud data.
        requested_cloud_alpha = max(0.0, min(1.0, cloud_disc_alpha))
        self._cloud_alpha_when_enabled = requested_cloud_alpha if requested_cloud_alpha > 0.0 else 0.2
        self.cloud_disc_alpha: float = requested_cloud_alpha
        if not self._cloud_toggle_supported:
            self.cloud_disc_alpha = 0.0

        # --- Viewer and Window Setup ---
        lat, lon, tz_name = city_data
        self.viewer_data = ViewerData(
            location=(lat, lon),
            timezone_name=tz_name,
            city_name=city_name,
            view_center=view_center,
        )
        self.state = SkyWindowState(
            render_view_center=tuple(self.viewer_data.view_center),
        )
        self._debug_eclipses = DEBUG_ECLIPSES
        self.setWindowTitle(f"{APP_DISPLAY_NAME} - {self.viewer_data.city_name.title()}")
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
        if window_geometry_arg == "restore":
            requested_geometry = load_last_window_geometry()
        elif isinstance(window_geometry_arg, tuple):
            requested_geometry = window_geometry_arg
        if requested_geometry is not None:
            initial_x, initial_y, initial_width, initial_height = requested_geometry
        initial_x, initial_y, initial_width, initial_height = self._clamp_window_geometry_to_screen(
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
        self._action_toggle_dso: Optional[QAction] = None
        self._action_toggle_asterisms: Optional[QAction] = None
        self._add_hamburger_menu()
        self.add_drag_exclusions([self.menu_button, self.size_grip])

        # --- Fonts ---
        text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
        text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
        self.text_font = QFont(text_font_family, TEXT_FONT_SIZE)
        self.status_line_font = QFont(text_font_family, STATUS_LINE_FONT_SIZE)

        # --- Data Update Timers and State ---
        self._sky_worker = SkyDataWorker(self)
        self._sky_worker.data_ready.connect(self._on_sky_data_calculated)
        self.cloud_repaint_requested.connect(self.update)
        self._is_shutting_down: bool = False
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

        # --- Cloud Data State and Cache ---
        self.cloud_state = CloudImageState()
        self._cloud_controller: Optional[CloudController] = None
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

        # --- Composition Cache (moved to dedicated class) ---
        target_stripes, width_factor = cloud_stripe_style
        missing_tint_alpha = int(round(255.0 * max(0.0, min(1.0, float(cloud_missing_tint_opacity)))))
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
        )

        # Cloud error banner is kept inside CloudImageState

        # --- Initial Data Load ---
        self.start_background_sky_data_update(is_initial_load=True)
        if self._clouddisc and self.cloud_disc_alpha > 0.0:
            self.start_background_cloud_update(reason="initial")

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
        jump_named_star = self.menu.addAction("Jump to Named Star...")
        jump_named_star.triggered.connect(self._open_named_star_jump_dialog)
        search_named_star = self.menu.addAction("Search Named Stars...")
        search_named_star.triggered.connect(self._open_named_star_search_dialog)

        toggle_enlarge_moon_action = QAction("Enlarge Moon (M)", self)
        toggle_enlarge_moon_action.setCheckable(True)
        toggle_enlarge_moon_action.setChecked(self.enlarge_moon)
        toggle_enlarge_moon_action.triggered.connect(self.toggle_enlarge_moon)
        self.menu.addAction(toggle_enlarge_moon_action)
        self._action_enlarge_moon = toggle_enlarge_moon_action
        toggle_clouds_action = QAction("Clouds", self)
        toggle_clouds_action.setCheckable(True)
        toggle_clouds_action.setChecked(self.cloud_disc_alpha > 0.0)
        toggle_clouds_action.triggered.connect(self.toggle_clouds)
        self.menu.addAction(toggle_clouds_action)
        self._action_toggle_clouds = toggle_clouds_action
        toggle_dso_action = QAction("DSO", self)
        toggle_dso_action.setCheckable(True)
        toggle_dso_action.setChecked(self.show_dso)
        toggle_dso_action.setEnabled(self.dso_catalog_np is not None)
        toggle_dso_action.triggered.connect(self.toggle_dso)
        self.menu.addAction(toggle_dso_action)
        self._action_toggle_dso = toggle_dso_action
        toggle_asterisms_action = QAction("Asterisms", self)
        toggle_asterisms_action.setCheckable(True)
        toggle_asterisms_action.setChecked(self.show_asterisms)
        toggle_asterisms_action.triggered.connect(self.toggle_asterisms)
        self.menu.addAction(toggle_asterisms_action)
        self._action_toggle_asterisms = toggle_asterisms_action

        self.menu.addSeparator()
        fullscreen_action = self.menu.addAction("Fullscreen (F11)")
        fullscreen_action.triggered.connect(self.toggle_fullscreen)

        self.menu.addSeparator()
        exit_action = self.menu.addAction("Exit (Q)")
        exit_action.triggered.connect(QApplication.quit)

        self.menu.addSeparator()
        version_action = self.menu.addAction(f"Version {__version__}")
        version_action.setEnabled(False)

    def resizeEvent(self, event: QResizeEvent) -> None:
        self._begin_interaction_mode()
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

    def show_menu(self) -> None:
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
        alt, az = radec_to_altaz(
            target.ra_hours,
            target.dec_deg,
            self.viewer_data.location[0],
            self.viewer_data.location[1],
            self._current_time_obj(),
        )
        new_alt = max(0.0, min(90.0, float(alt)))
        new_az = float(az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)

        self.state.jump_highlight_name = target.label
        self.state.jump_highlight_altaz = (new_alt, new_az)
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
        if self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.stop()
        if self._asterism_check_timer.isActive():
            self._asterism_check_timer.stop()
        if self._cloud_update_timer.isActive():
            self._cloud_update_timer.stop()

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

    def _safe_request_cloud_repaint(self) -> None:
        return SkyWindowUpdatesMixin._safe_request_cloud_repaint(self)

    def _predicted_cloud_satellite(self) -> str:
        lat, lon = self.viewer_data.location
        return pick_satellite(lat, lon, ("AUTO",))

    def _cloud_status_line(self) -> str:
        return SkyWindowUpdatesMixin._cloud_status_line(self)

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

    @staticmethod
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

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        self.state.mouse_pos = event.pos()
        self.update()  # Trigger a repaint to show hover effects
        # We accept the event to prevent it from propagating further.
        # This is why the manual drag in DraggableWindow does not work.
        event.accept()

    def set_sky_data(self, data: CelestialData, *, trigger_update: bool = True) -> None:
        self.state.celestial_data = data
        if trigger_update:
            self.update()

    def _viewer_data_for_render(self) -> ViewerData:
        return SkyWindowRenderMixin._viewer_data_for_render(self)

    def _active_jump_highlight_object(
        self,
        geometry,
    ):
        return SkyWindowRenderMixin._active_jump_highlight_object(self, geometry)

    def paintEvent(self, event: QPaintEvent) -> None:
        return SkyWindowRenderMixin.paintEvent(self, event)

    def _clear_background_layer(self, painter: QPainter) -> None:
        return SkyWindowRenderMixin._clear_background_layer(self, painter)

    def _draw_background_layer(self, painter: QPainter, geometry: render_draw.ScreenGeometry) -> None:
        return SkyWindowRenderMixin._draw_background_layer(self, painter, geometry)

    def _draw_sky_cloud_layers(self, painter: QPainter, geometry: render_draw.ScreenGeometry) -> None:
        return SkyWindowRenderMixin._draw_sky_cloud_layers(self, painter, geometry)

    def _draw_terrain_layers(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
        highlighted_object: Any | None,
        highlighted_dso: Any | None,
        label_reservations: list[QRectF],
        label_candidates: list[dict[str, Any]],
    ) -> None:
        return SkyWindowRenderMixin._draw_terrain_layers(
            self,
            painter,
            geometry,
            celestial_data,
            render_viewer,
            highlighted_object,
            highlighted_dso,
            label_reservations,
            label_candidates,
        )

    def _draw_star_layer(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
    ) -> None:
        return SkyWindowRenderMixin._draw_star_layer(
            self,
            painter,
            geometry,
            celestial_data,
            render_viewer,
        )

    def _draw_planet_layer(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
        enlarge_moon: bool,
        label_candidates: list[dict[str, Any]],
    ) -> None:
        return SkyWindowRenderMixin._draw_planet_layer(
            self,
            painter,
            geometry,
            celestial_data,
            render_viewer,
            enlarge_moon,
            label_candidates,
        )

    def _draw_overlay_layer(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        celestial_data: CelestialData,
        render_viewer: ViewerData,
        highlighted_object: Any | None,
        highlighted_dso: Any | None,
        enlarge_moon: bool,
        label_reservations: list[QRectF],
        label_candidates: list[dict[str, Any]],
    ) -> None:
        return SkyWindowRenderMixin._draw_overlay_layer(
            self,
            painter,
            geometry,
            celestial_data,
            render_viewer,
            highlighted_object,
            highlighted_dso,
            enlarge_moon,
            label_reservations,
            label_candidates,
        )

    def _draw_label_layer(self, painter: QPainter, label_candidates: list[dict[str, Any]]) -> None:
        return SkyWindowRenderMixin._draw_label_layer(self, painter, label_candidates)

    def _draw_status_line(self, painter: QPainter) -> None:
        return SkyWindowRenderMixin._draw_status_line(self, painter)

    # Composited drawing handled by SkyCompositorCache

    def _on_sky_data_calculated(self, payload: Dict) -> None:
        return SkyWindowUpdatesMixin._on_sky_data_calculated(self, payload)

    def request_sky_data_update(
        self,
        star_vmag_limit: Optional[float] = None,
    ) -> None:
        return SkyWindowUpdatesMixin.request_sky_data_update(self, star_vmag_limit)

    def start_background_sky_data_update(
        self,
        is_initial_load: bool = False,
        star_vmag_limit: Optional[float] = None,
    ) -> bool:
        return SkyWindowUpdatesMixin.start_background_sky_data_update(
            self,
            is_initial_load=is_initial_load,
            star_vmag_limit=star_vmag_limit,
        )

    def start_background_cloud_update(self, reason: str = "manual") -> None:
        return SkyWindowUpdatesMixin.start_background_cloud_update(self, reason)

    def _on_cloud_started(self, payload: Dict) -> None:
        return SkyWindowUpdatesMixin._on_cloud_started(self, payload)

    def _on_cloud_ready(self, payload: Dict) -> None:
        return SkyWindowUpdatesMixin._on_cloud_ready(self, payload)

    def _on_cloud_failed(self, payload: Dict) -> None:
        return SkyWindowUpdatesMixin._on_cloud_failed(self, payload)

    def _rotate_view(self, d_alt: float = 0.0, d_az: float = 0.0) -> None:
        self._begin_interaction_mode()
        alt, az = self.viewer_data.view_center
        new_alt = max(0.0, min(90.0, alt + d_alt))
        new_az = (az + d_az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)

        # Recalculate sky data for the current full quality on each key input.
        self.request_sky_data_update()

    def keyPressEvent(self, event: QKeyEvent) -> None:
        key = event.key()

        # --- View Control ---
        if key == Qt.Key.Key_Left:
            self._rotate_view(d_az=-self.state.rotation_step)
            event.accept()
        elif key == Qt.Key.Key_Right:
            self._rotate_view(d_az=self.state.rotation_step)
            event.accept()
        elif key == Qt.Key.Key_Up:
            self._rotate_view(d_alt=self.state.rotation_step)
            event.accept()
        elif key == Qt.Key.Key_Down:
            self._rotate_view(d_alt=-self.state.rotation_step)
            event.accept()

        # --- Toggles ---
        elif key == Qt.Key.Key_M:
            self.toggle_enlarge_moon()
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
