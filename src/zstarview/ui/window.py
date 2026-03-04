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
from PySide6.QtCore import QEvent, QPoint, QPointF, QRect, Qt, QTimer, Signal
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
    build_named_star_shortcuts,
    flatten_named_star_shortcuts,
)
from .sky_worker import SkyDataWorker

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


class SkyWindow(DraggableWindow):
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
        star_render_expected_width: int = 600,
        window_geometry_arg: Optional[WindowGeometryArg] = None,
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
        self._rotation_step: float = 5.0  # degrees
        self._interaction_idle_ms: int = 300
        self._interaction_mode: bool = False

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
        self._named_stars_by_band = build_named_star_shortcuts(star_catalog, max_vmag=2.0)
        self._named_stars_search_all = flatten_named_star_shortcuts(
            build_named_star_shortcuts(star_catalog, max_vmag=None)
        )
        self._jump_highlight_name: Optional[str] = None
        self._jump_highlight_altaz: Optional[Tuple[float, float]] = None
        self._jump_highlight_until_ms: float = 0.0
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
        self.mouse_pos: Optional[QPoint] = None

        # --- UI Widgets ---
        self.size_grip = QSizeGrip(self)
        self.size_grip.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.size_grip.raise_()
        self._action_enlarge_moon: Optional[QAction] = None
        self._action_toggle_clouds: Optional[QAction] = None
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
        self._interaction_idle_timer = QTimer(self)
        self._interaction_idle_timer.setSingleShot(True)
        self._interaction_idle_timer.setInterval(self._interaction_idle_ms)
        self._interaction_idle_timer.timeout.connect(self._end_interaction_mode)
        self._sky_update_pending: bool = False
        self._pending_star_vmag_limit: Optional[float] = None
        self._cloud_repaint_deferred: bool = False
        self._last_star_render_stats: Optional[tuple[int, int, int, int]] = None
        self.celestial_data: Optional[CelestialData] = None
        self._sky_disc_base_size: int = 1024
        self._sky_disc_image: Optional[QImage] = None

        # --- Cloud Data State and Cache ---
        self._cloud_base_size: int = 256
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
        self._compositor = SkyCompositorCache(
            cloud_target_stripes=int(target_stripes),
            cloud_stripe_width_factor=float(width_factor),
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
        rotate_left = self.menu.addAction(f"Rotate Left (-{self._rotation_step:.0f}°)")
        rotate_left.triggered.connect(lambda: self._rotate_view(d_az=-self._rotation_step))
        rotate_right = self.menu.addAction(f"Rotate Right (+{self._rotation_step:.0f}°)")
        rotate_right.triggered.connect(lambda: self._rotate_view(d_az=+self._rotation_step))
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
        self._interaction_mode = True
        self._interaction_idle_timer.start()

    def _end_interaction_mode(self) -> None:
        if not self._interaction_mode:
            return
        self._interaction_mode = False
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
        star = dialog.selected_star()
        if star is None:
            return
        self._jump_to_named_star(star)

    def _jump_to_named_star(self, star: NamedStarShortcut) -> None:
        lat, lon = self.viewer_data.location
        alt, az = radec_to_altaz(
            star.ra_hours,
            star.dec_deg,
            lat,
            lon,
            self._current_time_obj(),
        )
        new_alt = max(0.0, min(90.0, float(alt)))
        new_az = float(az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)

        # Mark the selected result for 3s using the same overlay highlight style.
        self._jump_highlight_name = star.name
        self._jump_highlight_altaz = (new_alt, new_az)
        self._jump_highlight_until_ms = (time.monotonic() * 1000.0) + 3000.0

        self._begin_interaction_mode()
        self.request_sky_data_update()
        self.update()

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
        if self._cloud_update_timer.isActive():
            self._cloud_update_timer.stop()

    def _safe_request_cloud_repaint(self) -> None:
        """Best-effort repaint request; ignores teardown-time signal errors."""
        if self._is_shutting_down:
            return
        try:
            self.cloud_repaint_requested.emit()
        except RuntimeError:
            logger.debug("Skip cloud repaint emit during shutdown.")

    def _predicted_cloud_satellite(self) -> str:
        lat, lon = self.viewer_data.location
        return pick_satellite(lat, lon, ("AUTO",))

    def _cloud_status_line(self) -> str:
        sat = self.cloud_state.current_satellite or self._predicted_cloud_satellite()
        if self.cloud_state.banner_text:
            detail = self.cloud_state.banner_text.removeprefix("Clouds:").strip()
            return f"Clouds [{sat}]: {detail}"
        meta = self.cloud_state.meta
        if meta is not None:
            try:
                t = meta.time_utc.strftime("%H:%MZ")
                return f"Clouds [{meta.satellite}]: {meta.product} {t}"
            except Exception:
                pass
        return f"Clouds [{sat}]: idle"

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

    def toggle_fullscreen(self) -> None:
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def leaveEvent(self, event: QEvent) -> None:
        self.mouse_pos = None
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
        self.mouse_pos = event.pos()
        self.update()  # Trigger a repaint to show hover effects
        # We accept the event to prevent it from propagating further.
        # This is why the manual drag in DraggableWindow does not work.
        event.accept()

    def set_sky_data(self, data: CelestialData) -> None:
        self.celestial_data = data
        self.update()

    def _active_jump_highlight_object(
        self,
        geometry,
    ):
        if not self._jump_highlight_name:
            return None
        if time.monotonic() * 1000.0 > self._jump_highlight_until_ms:
            self._jump_highlight_name = None
            self._jump_highlight_altaz = None
            self._jump_highlight_until_ms = 0.0
            return None
        if not self.celestial_data:
            return None

        target_name = self._jump_highlight_name
        stars = self.celestial_data.stars
        best_idx = None
        best_vmag = float("inf")
        for idx, raw_name in enumerate(stars["name"]):
            name = str(raw_name).strip()
            if name != target_name:
                continue
            vmag = float(stars["vmag"][idx])
            if vmag < best_vmag:
                best_vmag = vmag
                best_idx = idx

        if best_idx is not None:
            alt = float(stars["alt"][best_idx])
            az = float(stars["az"][best_idx])
        elif self._jump_highlight_altaz is not None:
            alt, az = self._jump_highlight_altaz
        else:
            return None

        nx, ny = render_draw.altaz_to_normalized_xy(alt, az, self.viewer_data.view_center)
        px, py = render_draw.normalized_to_screen_xy(nx, ny, geometry)
        return ({"name": target_name}, QPointF(px, py))

    def paintEvent(self, event: QPaintEvent) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)

        # If data is not yet loaded, show a loading message.
        if not self.celestial_data:
            loading_color, _ = render_draw.get_text_style(self.visual_preset)
            painter.setPen(loading_color)
            painter.setFont(self.text_font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Loading celestial data...")
            return

        alt = self.viewer_data.view_center[0]
        geometry = render_draw.get_screen_geometry(self.width(), self.height(), alt)

        highlighted_object = None
        highlighted_dso = None
        if self.mouse_pos is not None:
            highlighted_object = render_draw.find_highlighted_object(self.celestial_data, self.viewer_data, self.mouse_pos, geometry)
            highlighted_dso = render_draw.find_highlighted_dso(self.celestial_data, self.viewer_data, self.mouse_pos, geometry)
        jump_highlight = self._active_jump_highlight_object(geometry)
        if jump_highlight is not None:
            highlighted_object = jump_highlight

        self._clear_background_layer(painter)
        self._draw_background_layer(painter, geometry)
        self._draw_sky_cloud_layers(painter, geometry)
        self._draw_terrain_layers(painter, geometry, highlighted_dso)
        self._draw_star_layer(painter, geometry)

        enlarge_moon = self.enlarge_moon
        if highlighted_object is not None:
            obj = highlighted_object[0]
            name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
            enlarge_moon = enlarge_moon or name == "moon"

        self._draw_planet_layer(painter, geometry, enlarge_moon)
        self._draw_overlay_layer(painter, geometry, highlighted_object, highlighted_dso, enlarge_moon)
        self._draw_status_line(painter)

    def _clear_background_layer(self, painter: QPainter) -> None:
        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode_Clear)
        painter.fillRect(self.rect(), Qt.transparent)
        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
        painter.restore()

    def _draw_background_layer(self, painter: QPainter, geometry: render_draw.ScreenGeometry) -> None:
        render_draw.draw_radial_background(painter, self.rect(), geometry, preset=self.visual_preset)

    def _draw_sky_cloud_layers(self, painter: QPainter, geometry: render_draw.ScreenGeometry) -> None:
        self._compositor.draw(
            painter,
            geometry,
            self._sky_disc_image,
            self.cloud_state.image,
            cloud_alpha=self.cloud_disc_alpha,
            stripe_density=self.cloud_state.stripe_density,
        )

    def _draw_terrain_layers(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        highlighted_dso: Any | None,
    ) -> None:
        render_draw.draw_deep_sky_shapes(
            painter,
            geometry,
            self.celestial_data,
            self.viewer_data,
            preset=self.visual_preset,
        )
        render_draw.draw_dso_hover_info(
            painter,
            geometry,
            self.viewer_data,
            highlighted_dso,
            self.text_font,
            preset=self.visual_preset,
        )
        render_draw.draw_sky_reference_lines(painter, geometry, self.celestial_data)
        render_draw.draw_direction_labels(
            painter,
            geometry,
            self.viewer_data.view_center,
            self.text_font,
            preset=self.visual_preset,
        )
        render_draw.draw_zenith_marker(painter, geometry, self.viewer_data.view_center)

    def _draw_star_layer(self, painter: QPainter, geometry: render_draw.ScreenGeometry) -> None:
        win_w = self.width()
        win_h = self.height()
        low_w, low_h = compute_star_render_surface_size(
            win_w,
            win_h,
            geometry.radius * 2,
            self._star_render_expected_width,
        )
        stats = (win_w, win_h, low_w, low_h)
        if stats != self._last_star_render_stats:
            logger.info(
                "Star render resolution: window=%dx%d draw=%dx%d",
                win_w,
                win_h,
                low_w,
                low_h,
            )
            self._last_star_render_stats = stats

        if low_w == win_w and low_h == win_h:
            render_draw.draw_stars(
                painter,
                geometry,
                self.celestial_data,
                self.viewer_data,
                self.star_base_radius,
                visibility_boost=self.star_visibility_boost,
                draw_vmag_limit=self.vmag_limit,
                viewport_size=(win_w, win_h),
            )
            return

        low_img = QImage(low_w, low_h, QImage.Format.Format_ARGB32_Premultiplied)
        low_img.fill(Qt.GlobalColor.transparent)
        low_painter = QPainter(low_img)
        low_painter.setRenderHint(QPainter.Antialiasing, False)
        low_painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
        sx = low_w / max(1.0, float(win_w))
        sy = low_h / max(1.0, float(win_h))
        low_geometry = render_draw.ScreenGeometry(
            center=(
                int(round(geometry.center[0] * sx)),
                int(round(geometry.center[1] * sy)),
            ),
            radius=max(1, int(round(geometry.radius * min(sx, sy)))),
        )
        render_draw.draw_stars(
            low_painter,
            low_geometry,
            self.celestial_data,
            self.viewer_data,
            self.star_base_radius,
            visibility_boost=self.star_visibility_boost,
            draw_vmag_limit=self.vmag_limit,
            viewport_size=(low_w, low_h),
        )
        low_painter.end()

        painter.save()
        painter.setRenderHint(QPainter.SmoothPixmapTransform, False)
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
        painter.drawImage(self.rect(), low_img)
        painter.restore()

    def _draw_planet_layer(self, painter: QPainter, geometry: render_draw.ScreenGeometry, enlarge_moon: bool) -> None:
        render_draw.draw_solar_system_bodies(
            painter,
            geometry,
            self.celestial_data,
            self.viewer_data,
            enlarge_moon,
            text_font=self.text_font,
            preset=self.visual_preset,
        )

    def _draw_overlay_layer(
        self,
        painter: QPainter,
        geometry: render_draw.ScreenGeometry,
        highlighted_object: Any | None,
        highlighted_dso: Any | None,
        enlarge_moon: bool,
    ) -> None:
        render_draw.draw_overlay_info(
            painter,
            geometry,
            self.celestial_data,
            self.viewer_data,
            self.vmag_limit,
            enlarge_moon,
            highlighted_dso,
            highlighted_object,
            self.text_font,
            preset=self.visual_preset,
        )

    def _draw_status_line(self, painter: QPainter) -> None:
        if self.cloud_disc_alpha > 0.0 or self.cloud_state.banner_text:
            render_draw.draw_status_line_text(
                painter=painter,
                message=self._cloud_status_line(),
                status_line_font=self.status_line_font,
                viewport_rect=self.rect(),
                preset=self.visual_preset,
            )

    # Composited drawing handled by SkyCompositorCache

    def _on_sky_data_calculated(self, payload: Dict) -> None:
        """
        Slot to handle the arrival of newly calculated sky data from the worker thread.

        Args:
            payload: A dictionary containing the 'celestial' data and 'sky_disc' image.
        """
        self.set_sky_data(payload["celestial"])
        self._sky_disc_image = payload["sky_disc"]

        # Invalidate the composition cache
        self._compositor.invalidate()

        # On the very first load, start the periodic update timers.
        if not self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.start(self.sky_update_interval * 1000)
            if self._clouddisc and self.cloud_disc_alpha > 0.0 and not self._cloud_update_timer.isActive():
                self._cloud_update_timer.start()
            self.initial_data_loaded.emit()

        # If a request came in while worker was busy, run one catch-up update
        # using the latest view center/time once current payload is applied.
        if self._sky_update_pending and not self._is_shutting_down:
            self.request_sky_data_update(
                self._pending_star_vmag_limit,
            )

        if self._cloud_repaint_deferred and not self._interaction_mode:
            self._cloud_repaint_deferred = False
            self._safe_request_cloud_repaint()

    def request_sky_data_update(
        self,
        star_vmag_limit: Optional[float] = None,
    ) -> None:
        """Requests a sky data update if one is not already in progress."""
        if self.start_background_sky_data_update(
            star_vmag_limit=star_vmag_limit,
        ):
            self._sky_update_pending = False
            self._pending_star_vmag_limit = None
            return
        self._sky_update_pending = True
        self._pending_star_vmag_limit = star_vmag_limit
        logger.debug("Sky data update deferred; worker is busy.")

    def start_background_sky_data_update(
        self,
        is_initial_load: bool = False,
        star_vmag_limit: Optional[float] = None,
    ) -> bool:
        lat, lon = self.viewer_data.location
        use_lod6_catalog = star_vmag_limit is not None and float(star_vmag_limit) <= 6.0
        star_catalog = self.star_catalog_lod6_np if use_lod6_catalog else self.star_catalog_full_np
        worker_star_vmag_limit = None if use_lod6_catalog else star_vmag_limit
        started = self._sky_worker.update(
            lat=lat,
            lon=lon,
            view_center=self.viewer_data.view_center,
            star_catalog=star_catalog,
            dso_catalog=self.dso_catalog_np,
            star_vmag_limit=worker_star_vmag_limit,
            delta_t=self.delta_t,
            sky_disc_alpha=self.sky_disc_alpha,
            sky_disc_base_size=self._sky_disc_base_size,
            debug_eclipses=DEBUG_ECLIPSES,
        )
        if started:
            if is_initial_load:
                logger.info("Calculating initial sky data...")
            else:
                logger.info("Updating sky data...")
        return started

    def start_background_cloud_update(self, reason: str = "manual") -> None:
        if self._is_shutting_down:
            return
        if not (self._cloud_controller and self.cloud_disc_alpha > 0.0):
            return
        lat, lon = self.viewer_data.location
        alt, az = self.viewer_data.view_center
        self._cloud_controller.update(
            lat=lat,
            lon=lon,
            alt=alt,
            az=az,
            radius_px=self._cloud_base_size,
            reason=reason,
        )

    def _on_cloud_started(self, payload: Dict) -> None:
        sat = str(payload.get("satellite", "")).strip()
        banner = str(payload.get("banner", "")).strip()
        if sat:
            self.cloud_state.current_satellite = sat
        if banner:
            self.cloud_state.set_error_banner(banner)

    def _on_cloud_ready(self, payload: Dict) -> None:
        self.cloud_state.set_result(
            payload["image"],
            payload.get("meta"),
            az=float(payload["az"]),
            time_utc=payload["time_utc"],
            stripe_density=payload.get("stripe_density"),
        )
        self._compositor.invalidate()
        if self._interaction_mode:
            self._cloud_repaint_deferred = True
            return
        self._safe_request_cloud_repaint()

    def _on_cloud_failed(self, payload: Dict) -> None:
        banner = str(payload.get("banner", "")).strip()
        self.cloud_state.image = None
        self.cloud_state.stripe_density = None
        self._compositor.invalidate()
        if banner:
            self.cloud_state.set_error_banner(banner)
        if self._interaction_mode:
            self._cloud_repaint_deferred = True
            return
        self._safe_request_cloud_repaint()

    def _rotate_view(self, d_alt: float = 0.0, d_az: float = 0.0) -> None:
        self._begin_interaction_mode()
        alt, az = self.viewer_data.view_center
        new_alt = max(0.0, min(90.0, alt + d_alt))
        new_az = (az + d_az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)

        # Recalculate sky data for the current full quality on each key input.
        self.request_sky_data_update()

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if not event:
            super().keyPressEvent(event)
            return

        key = event.key()

        # --- View Control ---
        if key == Qt.Key.Key_Left:
            self._rotate_view(d_az=-self._rotation_step)
            event.accept()
        elif key == Qt.Key.Key_Right:
            self._rotate_view(d_az=self._rotation_step)
            event.accept()
        elif key == Qt.Key.Key_Up:
            self._rotate_view(d_alt=self._rotation_step)
            event.accept()
        elif key == Qt.Key.Key_Down:
            self._rotate_view(d_alt=-self._rotation_step)
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
