# -*- coding: utf-8 -*-
"""
The main window of the ZStarView application.

This module defines the `SkyWindow` class, which is the primary user interface
for the application. It handles rendering the celestial objects, sky background,
clouds, and all user interactions like rotation, zooming, and object highlighting.
"""
import logging
import threading
from datetime import datetime, timedelta, timezone
from typing import Dict, Optional, Tuple

import polars as pl
from PySide6.QtCore import QEvent, QPoint, Qt, QTimer, Signal
from PySide6.QtGui import (
    QAction,
    QFont,
    QFontDatabase,
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
from ..clouddisc import (
    CloudDisc,
    CloudDiscConfig,
    CloudDiscError,
    DataNotFoundError,
    DownloadError,
    TimeoutError,
    RenderError,
    VisibilityError,
    cleanup_satellite_cache,
)
from ..paths import (
    APP_ICON_FILE,
    APP_DISPLAY_NAME,
    CACHE_PATH,
    CLOUD_SHELL_KM,
    CLOUD_UPDATE_INTERVAL,
    EMOJI_FONT_PATH,
    EMOJI_FONT_SIZE,
    GUI_BUTTON_SIZE,
    GUI_MENU_TEXT_COLOR,
    TEXT_FONT_PATH,
    TEXT_FONT_SIZE,
    STATUS_LINE_FONT_SIZE,
    WINDOW_HEIGHT,
    WINDOW_WIDTH,
)
from ..render import draw as render_draw
from ..types import CelestialData, ViewerData
from ..utils.qt import pil_to_qimage
from .draggable_window import DraggableWindow
from .composite import SkyCompositorCache
from .cloud_state import CloudImageState
from .sky_worker import SkyDataWorker

logger = logging.getLogger(__name__)


DEBUG_ECLIPSES = True


# compositing helpers and cache moved to ui/composite.py


class SkyWindow(DraggableWindow):
    """
    The main application window, displaying the sky view.

    This class orchestrates the entire application. It manages fetching and
    processing of celestial and cloud data, renders all visual elements, and
    handles user input for interaction.
    """

    # Signal emitted when initial sky data is fully loaded and window can be shown.
    initial_data_loaded = Signal()

    # TODO(CloudController): Extract cloud fetching to a QObject controller
    #   (e.g., ui/cloud_controller.py) with:
    #   - Signals: cloud_ready(QImage|None, dict|None), status(str|None)
    #   - Methods: update(lat, lon, view_center, radius_px, params)
    #   - Responsibility: call CloudDisc.render_now, manage running/pending state,
    #     cleanup cadence, and set status messages; no UI painting.
    #   Wiring here:
    #   - Replace start_background_cloud_update/_update_cloud_in_background with
    #     controller.update(...); on signals: update cloud_state + compositor.invalidate()
    #   - Optionally move cleanup timer here; for now Window retains QTimer

    def __init__(
        self,
        city_name: str,
        city_data: Tuple[float, float, str],
        view_center: Tuple[float, float],
        star_catalog: pl.DataFrame,
        delta_t: timedelta = timedelta(days=0, hours=0),
        sky_disc_alpha: float = 0.3,
        cloud_disc_alpha: float = 0.6,
        enlarge_moon: bool = False,
        star_base_radius: float = 8.0,
        vmag_limit: float = 6.0,
        sky_update_interval: int = 3 * 60  # sec
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
            star_base_radius: Base radius for drawing stars.
            vmag_limit: The faintest star magnitude to display.
        """
        super().__init__()
        self._rotation_step: float = 5.0  # degrees

        self.star_catalog = star_catalog
        self.delta_t = delta_t
        self.sky_disc_alpha = sky_disc_alpha
        self.enlarge_moon = enlarge_moon
        self.star_base_radius = star_base_radius
        self.vmag_limit = vmag_limit
        self.sky_update_interval = sky_update_interval

        # Cloud opacity is disabled if we are looking at a time-shifted view,
        # as we can only fetch current cloud data.
        self.cloud_disc_alpha: float = max(0.0, min(1.0, cloud_disc_alpha))
        if delta_t.total_seconds() != 0.0:
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
        self.setMinimumSize(400, 400)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.FramelessWindowHint)
        self.setGeometry(100, 100, WINDOW_WIDTH, WINDOW_HEIGHT)
        self.setMouseTracking(True)
        self.mouse_pos: Optional[QPoint] = None

        # --- UI Widgets ---
        self.size_grip = QSizeGrip(self)
        self.size_grip.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.size_grip.raise_()
        self._action_enlarge_moon: Optional[QAction] = None
        self._add_hamburger_menu()
        self.add_drag_exclusions([self.menu_button, self.size_grip])

        # --- Fonts ---
        emoji_font_id = QFontDatabase.addApplicationFont(EMOJI_FONT_PATH)
        emoji_font_family = QFontDatabase.applicationFontFamilies(emoji_font_id)[0]
        self.emoji_font = QFont(emoji_font_family, EMOJI_FONT_SIZE)
        text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
        text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
        self.text_font = QFont(text_font_family, TEXT_FONT_SIZE)
        self.status_line_font = QFont(text_font_family, STATUS_LINE_FONT_SIZE)

        # --- Data Update Timers and State ---
        self._sky_worker = SkyDataWorker(self)
        self._sky_worker.data_ready.connect(self._on_sky_data_calculated)
        self._sky_data_update_timer = QTimer(self)
        self._sky_data_update_timer.timeout.connect(self.start_background_sky_data_update)
        self.celestial_data: Optional[CelestialData] = None
        self._sky_disc_base_size: int = 1024
        self._sky_disc_image: Optional[QImage] = None

        # --- Cloud Data State and Cache ---
        self._cloud_base_size: int = 256
        self.cloud_state = CloudImageState()
        self._is_cloud_update_running: bool = False
        self._cloud_update_pending: bool = False
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
            alt_min_deg=-1.0,
            search_back_minutes=120,
        )
        try:
            self._clouddisc = CloudDisc(clouddisc_config)
        except Exception as e:
            logger.warning(f"CloudDisc init failed: {e}")

        # --- Composition Cache (moved to dedicated class) ---
        self._compositor = SkyCompositorCache()

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

        toggle_enlarge_moon_action = QAction("Enlarge Moon (M)", self)
        toggle_enlarge_moon_action.setCheckable(True)
        toggle_enlarge_moon_action.setChecked(self.enlarge_moon)
        toggle_enlarge_moon_action.triggered.connect(self.toggle_enlarge_moon)
        self.menu.addAction(toggle_enlarge_moon_action)
        self._action_enlarge_moon = toggle_enlarge_moon_action

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
        grip_size = self.size_grip.size()
        self.size_grip.move(self.width() - grip_size.width(), self.height() - grip_size.height())

        button_size = self.menu_button.size()
        self.menu_button.move(self.width() - button_size.width() - 8, 8)

        # Invalidate the composition cache since the size has changed
        self._compositor.invalidate()

        super().resizeEvent(event)

    def show_menu(self) -> None:
        menu_pos = self.menu_button.mapToGlobal(QPoint(0, self.menu_button.height()))
        self.menu.exec(menu_pos)

    def toggle_enlarge_moon(self) -> None:
        self.enlarge_moon = not self.enlarge_moon
        if self._action_enlarge_moon is not None and self._action_enlarge_moon.isChecked() != self.enlarge_moon:
            self._action_enlarge_moon.setChecked(self.enlarge_moon)
        self.update()  # Redraw with the new setting

    def toggle_fullscreen(self) -> None:
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def leaveEvent(self, event: QEvent) -> None:
        self.mouse_pos = None
        self.update()
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        self.mouse_pos = event.pos()
        self.update()  # Trigger a repaint to show hover effects
        # We accept the event to prevent it from propagating further.
        # This is why the manual drag in DraggableWindow does not work.
        event.accept()

    def set_star_base_radius(self, star_base_radius: float) -> None:
        """Sets the base radius for drawing stars."""
        self.star_base_radius = star_base_radius

    def set_enlarge_moon(self, enlarge_moon: bool) -> None:
        """Sets the enlarge moon flag."""
        self.enlarge_moon = enlarge_moon

    def set_sky_data(self, data: CelestialData) -> None:
        self.celestial_data = data
        self.update()

    def paintEvent(self, event: QPaintEvent) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)

        # If data is not yet loaded, show a loading message.
        if not self.celestial_data:
            painter.setPen(Qt.GlobalColor.white)
            painter.setFont(self.text_font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Loading celestial data...")
            return

        # --- Main Drawing Sequence ---
        alt = self.viewer_data.view_center[0]
        geometry = render_draw.get_screen_geometry(self.width(), self.height(), alt)

        # Determine if any object is under the mouse cursor
        highlighted_object = None
        if self.mouse_pos is not None:
            highlighted_object = render_draw.find_highlighted_object(self.celestial_data, self.viewer_data, self.mouse_pos, geometry)

        # 1. Clear the background to be fully transparent
        painter.setCompositionMode(QPainter.CompositionMode_Clear)
        painter.fillRect(self.rect(), Qt.transparent)
        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)

        # 2. Draw the radial gradient background
        render_draw.draw_radial_background(painter, self.rect(), geometry)

        # 3. Draw the composited sky and cloud disc
        self._compositor.draw(
            painter,
            geometry,
            self._sky_disc_image,
            self.cloud_state.image,
            cloud_alpha=self.cloud_disc_alpha,
        )

        # 4. Draw reference lines (horizon, equator, etc.)
        render_draw.draw_sky_reference_lines(painter, geometry, self.celestial_data)
        render_draw.draw_direction_labels(painter, geometry, self.viewer_data.view_center, self.text_font)
        render_draw.draw_zenith_marker(painter, geometry, self.viewer_data.view_center)

        # 5. Draw celestial objects
        render_draw.draw_stars(painter, geometry, self.celestial_data, self.viewer_data, self.star_base_radius)

        # Enlarge moon if the global flag is set or if it's being hovered over.
        enlarge_moon = self.enlarge_moon
        if highlighted_object is not None:
            obj = highlighted_object[0]
            name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
            enlarge_moon = enlarge_moon or name == "moon"
        render_draw.draw_planets(painter, geometry, self.celestial_data, self.viewer_data, enlarge_moon, self.emoji_font)

        # 6. Draw overlay information text
        render_draw.draw_overlay_info(
            painter,
            self.celestial_data,
            self.viewer_data,
            self.vmag_limit,
            enlarge_moon,
            highlighted_object,
            self.text_font,
        )

        # 7. Draw persistent cloud error message (bottom-left), if any
        if self.cloud_state.banner_text:
            render_draw.draw_status_line_text(
                painter=painter,
                message=self.cloud_state.banner_text,
                status_line_font=self.status_line_font,
                viewport_rect=self.rect(),
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

    def request_sky_data_update(self) -> None:
        """Requests a sky data update if one is not already in progress."""
        if not self.start_background_sky_data_update():
            logger.warning("Sky data update request ignored; an update is already running.")

    def start_background_sky_data_update(self, is_initial_load: bool = False) -> bool:
        lat, lon = self.viewer_data.location
        started = self._sky_worker.update(
            lat=lat,
            lon=lon,
            view_center=self.viewer_data.view_center,
            star_catalog=self.star_catalog,
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
        if not (self._clouddisc and self.cloud_disc_alpha > 0.0):
            return

        # Run cache cleanup in a background thread periodically.
        if self.cloud_state.tick_cleanup():
            cleanup_thread = threading.Thread(target=self._cleanup_cache)
            cleanup_thread.daemon = True
            cleanup_thread.start()

        if self._is_cloud_update_running:
            # If an update is already running, flag that another one is pending.
            # This prevents dropping requests during rapid view changes.
            self._cloud_update_pending = True
            return

        self.cloud_state.set_error_banner("Clouds: downloading…")  # e.g., "Clouds: Downloading..."
        self.update()

        self._is_cloud_update_running = True
        t = threading.Thread(target=self._update_cloud_in_background, args=(reason,))
        t.daemon = True
        t.start()

    def _update_cloud_in_background(self, reason: str) -> None:
        try:
            lat, lon = self.viewer_data.location
            alt, az = self.viewer_data.view_center

            if reason == "initial":
                logger.info("Fetching initial cloud data...")
            else:
                logger.info("Updating cloud data...")

            try:
                pil_img, meta = self._clouddisc.render_now(
                    lat=lat,
                    lon=lon,
                    alt=float(alt),
                    az=float(az),
                    radius_px=self._cloud_base_size,
                    edge_fov_deg=90,
                    mask_fov_deg=93,
                    cloud_shell_km=CLOUD_SHELL_KM,
                )
                qimg = pil_to_qimage(pil_img)
                self.cloud_state.set_result(
                    qimg, meta, az=az, time_utc=datetime.now(timezone.utc)
                )

                # Invalidate composition cache to force a redraw with new clouds
                self._compositor.invalidate()
            except VisibilityError as e:
                logger.error("Invalid params for cloud-disc image generation: %s", e)
            except DownloadError as e:
                logger.warning("Network/S3 download error: %s", e)
                self.cloud_state.set_error_banner("Clouds: Network/S3 download error")
                self.update()
            except TimeoutError as e:
                logger.warning("Clouds fetch timed out: %s", e)
                self.cloud_state.set_error_banner("Clouds: Clouds fetch timed out")
                self.update()
            except DataNotFoundError as e:
                logger.info("Satellite data not found in search window: %s", e)
                self.cloud_state.set_error_banner("Clouds: Satellite data not found in search window")
            except RenderError as e:
                logger.error("Failed to decode/render satellite data: %s", e)
                self.cloud_state.set_error_banner("Clouds: Failed to decode/render satellite data")
            except CloudDiscError as e:
                logger.error("Unexpected clouddisc error: %s", e)
                self.cloud_state.set_error_banner("Clouds: Unexpected clouddisc error")

            
            # Trigger a repaint on the main thread
            self.update()
        except Exception as e:
            logger.error(f"Cloud update failed: {e}", exc_info=True)
        finally:
            self._is_cloud_update_running = False
            # If another update was requested while this one was running, start it now.
            if self._cloud_update_pending:
                self._cloud_update_pending = False
                self.start_background_cloud_update(reason="pending")

    def _cleanup_cache(self) -> None:
        if not self._clouddisc:
            return
        try:
            logger.info("Running satellite cache cleanup...")
            cleanup_satellite_cache(self._clouddisc.cfg.cache_root())
            logger.info("Done: Satellite cache cleanup.")
        except Exception as e:
            logger.error(f"Error during cache cleanup: {e}", exc_info=True)

    def _rotate_view(self, d_alt: float = 0.0, d_az: float = 0.0) -> None:
        alt, az = self.viewer_data.view_center
        new_alt = max(0.0, min(90.0, alt + d_alt))
        new_az = (az + d_az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)

        # Request updates for both sky and clouds since the view has changed.
        self.request_sky_data_update()
        self.start_background_cloud_update(reason="az-change")

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
