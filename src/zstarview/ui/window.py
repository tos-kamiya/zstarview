from datetime import datetime, timedelta, timezone
import threading
import sys
from typing import Optional, Tuple

from PySide6.QtCore import Qt, QPoint, QTimer
from PySide6.QtGui import QFont, QFontDatabase, QIcon, QPainter, QImage
from PySide6.QtGui import QAction, QKeyEvent, QMouseEvent, QPaintEvent, QResizeEvent
from PySide6.QtWidgets import QMainWindow, QSizeGrip, QApplication, QPushButton, QMenu

import astropy
import polars as pl

import zstarview.render.draw_sky_disc

from ..__about__ import __version__

from ..paths import EMOJI_FONT_PATH, EMOJI_FONT_SIZE, TEXT_FONT_PATH, TEXT_FONT_SIZE
from ..paths import APP_ICON_FILE, GUI_BUTTON_SIZE, GUI_MENU_TEXT_COLOR, WINDOW_HEIGHT, WINDOW_WIDTH
from ..types import SkyData, ViewerData
from ..astro import (
    calculate_visible_stars,
    calculate_planets,
    calculate_celestial_equator_points,
    calculate_ecliptic_points,
    calculate_horizon_points,
)
from ..render import draw as render_draw


class SkyWindow(QMainWindow):
    # Using Qt signal objects requires attribute creation at runtime; avoid type hints here
    from PySide6.QtCore import Signal

    sky_data_calculated = Signal(object)
    initial_data_loaded = Signal()

    def __init__(
        self,
        city_name: str,
        city_data: Tuple[float, float, str],
        view_center: Tuple[float, float],
        star_catalog: pl.DataFrame,
        delta_t: timedelta = timedelta(days=0, hours=0),
        enlarge_moon: bool = False,
        star_base_radius: float = 8.0,
        vmag_limit: float = 6.0,
    ):
        super().__init__()
        self._rotation_step: float = 5.0  # degrees

        self.star_catalog = star_catalog
        self.delta_t = delta_t
        self.enlarge_moon = enlarge_moon
        self.star_base_radius = star_base_radius
        self.vmag_limit = vmag_limit

        lat, lon, tz_name = city_data
        self.viewer_data = ViewerData(
            location=(lat, lon),
            timezone_name=tz_name,
            city_name=city_name,
            view_center=view_center,
        )
        self.setWindowTitle(f"Zenith Star View - {self.viewer_data.city_name.title()}")

        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setFocus(Qt.FocusReason.ActiveWindowFocusReason)

        self.setWindowIcon(QIcon(APP_ICON_FILE))

        self.setMinimumSize(400, 400)
        self.setAttribute(Qt.WA_TranslucentBackground)

        # Preserve existing flags and add Frameless; improves drag behavior on Qt6
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.FramelessWindowHint)
        self.setGeometry(100, 100, WINDOW_WIDTH, WINDOW_HEIGHT)

        # Size grip
        self.size_grip = QSizeGrip(self)
        self.size_grip.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.size_grip.raise_()

        self.setMouseTracking(True)
        self.mouse_pos: Optional[QPoint] = None
        self._drag_active: bool = False
        self._drag_pos: QPoint = QPoint(0, 0)

        # Menu
        self._action_enlarge_moon = None
        self._add_hamburger_menu()

        # Fonts for drawing
        emoji_font_id = QFontDatabase.addApplicationFont(EMOJI_FONT_PATH)
        emoji_font_family = QFontDatabase.applicationFontFamilies(emoji_font_id)[0]
        self.emoji_font = QFont(emoji_font_family, EMOJI_FONT_SIZE)
        text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
        text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
        self.text_font = QFont(text_font_family, TEXT_FONT_SIZE)

        self.sky_data_calculated.connect(self._on_sky_data_calculated)

        self.sky_data: Optional[SkyData] = None
        self._is_sky_data_calculation_running: bool = False
        self._sky_data_update_timer = QTimer(self)
        self._sky_data_update_timer.timeout.connect(self.start_background_update)

        self._sky_disc_cache_image: Optional[QImage] = None
        self._sky_disc_cache_key: Optional[tuple] = None

        self._resize_ongoing_timer = QTimer(self)
        self._resize_ongoing_timer.setSingleShot(True)
        self._resize_ongoing_timer.setInterval(1000)  # 1000ms
        self._resize_ongoing_timer.timeout.connect(self.update)

        self.start_background_update(is_initial_load=True)

    def _add_hamburger_menu(self):
        self.menu_button = QPushButton("☰", self)
        self.menu_button.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.menu_button.setStyleSheet(
            "QPushButton { border: none; font-size: 18px; background-color: transparent; color: "
            + GUI_MENU_TEXT_COLOR
            + "; }"
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

    def resizeEvent(self, event: QResizeEvent):
        grip_size = self.size_grip.size()
        self.size_grip.move(self.width() - grip_size.width(), self.height() - grip_size.height())

        button_size = self.menu_button.size()
        self.menu_button.move(self.width() - button_size.width() - 8, 8)

        # invalidate sky-disc image
        self._sky_disc_cache_image = None
        self._sky_disc_cache_key = None

        self._resize_ongoing_timer.start()

        super().resizeEvent(event)

    def show_menu(self):
        menu_pos = self.menu_button.mapToGlobal(QPoint(0, self.menu_button.height()))
        self.menu.exec(menu_pos)

    def toggle_enlarge_moon(self):
        self.enlarge_moon = not self.enlarge_moon
        if self._action_enlarge_moon is not None and self._action_enlarge_moon.isChecked() != self.enlarge_moon:
            self._action_enlarge_moon.setChecked(self.enlarge_moon)
        self.update()  # Redraw

    def toggle_fullscreen(self):
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.LeftButton:
            # Ignore drags initiated on interactive child widgets
            child = self.childAt(event.pos())
            if child in (self.menu_button, self.size_grip):
                super().mousePressEvent(event)
                return

            # On Wayland (and some platforms), top-level window move must be delegated
            try:
                wh = self.windowHandle()
                if wh is not None:
                    try:
                        if bool(wh.startSystemMove()):
                            event.accept()
                            return
                    except Exception:
                        pass
            except Exception:
                pass

            self._drag_active = True
            try:
                global_pos = event.globalPosition().toPoint()  # Qt6 API
            except AttributeError:
                global_pos = event.globalPos()  # Fallback for Qt5-style
            self._drag_pos = global_pos - self.frameGeometry().topLeft()
            event.accept()

    def leaveEvent(self, event):
        self.mouse_pos = None
        self.update()
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent):
        if getattr(self, "_drag_active", False) and event.buttons() & Qt.LeftButton:
            try:
                global_pos = event.globalPosition().toPoint()
            except AttributeError:
                global_pos = event.globalPos()
            target_pos = global_pos - self._drag_pos
            self.move(target_pos)
            event.accept()
        else:
            self.mouse_pos = event.pos()
            self.update()
            event.accept()

    def mouseReleaseEvent(self, event: QMouseEvent):
        self._drag_active = False
        event.accept()

    def set_star_base_radius(self, star_base_radius: float):
        self.star_base_radius = star_base_radius

    def set_enlarge_moon(self, enlarge_moon: bool):
        self.enlarge_moon = enlarge_moon

    def set_sky_data(self, data: SkyData):
        self.sky_data = data
        self.update()

    def paintEvent(self, event: QPaintEvent):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        if not self.sky_data:
            painter.setPen(Qt.GlobalColor.white)
            painter.setFont(self.text_font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Loading sky data...")
            return

        alt = self.viewer_data.view_center[0]
        geometry = render_draw.get_screen_geometry(self.width(), self.height(), alt)

        highlighted_object = None
        if self.mouse_pos is not None:
            highlighted_object = render_draw.find_highlighted_object(
                self.sky_data, self.viewer_data, self.mouse_pos, geometry
            )

        painter.setCompositionMode(QPainter.CompositionMode_Clear)
        painter.fillRect(self.rect(), Qt.transparent)

        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)

        render_draw.draw_radial_background(painter, self.rect(), geometry)

        if not self._resize_ongoing_timer.isActive():  # skip sky-disc image generation during window resizing
            self._handle_sky_disc_drawing(painter, geometry)

        render_draw.draw_sky_reference_lines(painter, geometry, self.sky_data)
        render_draw.draw_direction_labels(painter, geometry, self.viewer_data.view_center, self.text_font)
        render_draw.draw_zenith_marker(painter, geometry, self.viewer_data.view_center)

        render_draw.draw_stars(painter, geometry, self.sky_data, self.viewer_data, self.star_base_radius)

        enlarge_moon = self.enlarge_moon
        if highlighted_object is not None:
            obj = highlighted_object[0]
            name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
            enlarge_moon = enlarge_moon or name == "moon"
        render_draw.draw_planets(painter, geometry, self.sky_data, self.viewer_data, enlarge_moon, self.emoji_font)

        render_draw.draw_overlay_info(
            painter, self.sky_data, self.viewer_data, self.vmag_limit, enlarge_moon, highlighted_object, self.text_font
        )

    def _handle_sky_disc_drawing(self, painter, geometry):
        sun_altaz = None
        for body in self.sky_data.planets:
            if body.name == "sun":
                sun_altaz = (body.alt, body.az)
        if sun_altaz:
            current_cache_key = (
                geometry.center,
                geometry.radius,
                self.viewer_data.view_center,
                sun_altaz,
            )
            if (
                self._sky_disc_cache_image is None
                or self._sky_disc_cache_key != current_cache_key
            ):
                print("Updating sky disc...")
                self._sky_disc_cache_image = zstarview.render.draw_sky_disc.draw_sky_color_disc(
                    geometry, self.viewer_data.view_center, sun_altaz
                )
                self._sky_disc_cache_key = current_cache_key

            # Draw the cached image
            painter.drawImage(
                int(geometry.center[0] - geometry.radius),
                int(geometry.center[1] - geometry.radius),
                self._sky_disc_cache_image,
            )

    def _on_sky_data_calculated(self, sky_data: SkyData):
        self.set_sky_data(sky_data)
        if not self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.start(5 * 60 * 1000)
            self.initial_data_loaded.emit()
        self._is_sky_data_calculation_running = False

    def request_sky_update(self):
        """Request a sky update, but only if no update is currently running."""
        f = self.start_background_update()
        if not f:
            print("Warning: sky-data updating canceled.")

    def update_sky_data_in_background(self):
        try:
            now = datetime.now(timezone.utc) + self.delta_t
            time_obj = astropy.time.Time(now)
            lat, lon = self.viewer_data.location
            # Pass the latest view_center to the calculation function
            stars, loc = calculate_visible_stars(self.star_catalog, lat, lon, time_obj, self.viewer_data.view_center)
            planets = calculate_planets(lat, lon, time_obj, self.viewer_data.view_center)
            celestial_equator_points = calculate_celestial_equator_points(loc, time_obj, self.viewer_data.view_center)
            ecliptic_points = calculate_ecliptic_points(loc, time_obj, self.viewer_data.view_center)
            horizon_points = calculate_horizon_points(loc, time_obj, self.viewer_data.view_center)
            sky_data = SkyData(
                time=time_obj,
                planets=planets,
                stars=stars,
                celestial_equator_points=celestial_equator_points,
                ecliptic_points=ecliptic_points,
                horizon_points=horizon_points,
            )
            self.sky_data_calculated.emit(sky_data)
        except Exception as e:
            print(f"Error in background update thread: {e}", file=sys.stderr)
            import traceback

            traceback.print_exc()

    def start_background_update(self, is_initial_load: bool = False) -> bool:
        if self._is_sky_data_calculation_running:
            return False
        self._is_sky_data_calculation_running = True

        if is_initial_load:
            print("Calculating initial sky data...")
        else:
            print("Updating sky data...")
        thread = threading.Thread(target=self.update_sky_data_in_background)
        thread.daemon = True
        thread.start()
        return True

    def _rotate_view(self, d_alt: float = 0.0, d_az: float = 0.0):
        alt, az = self.viewer_data.view_center
        alt = max(0.0, min(90.0, alt + d_alt))
        az = (az + d_az) % 360.0
        self.viewer_data.view_center = (alt, az)
        self.request_sky_update()

    def keyPressEvent(self, event: QKeyEvent):
        if not event:
            super().keyPressEvent(event)
            return

        key = event.key()

        # View control
        if key == Qt.Key.Key_Left:
            self._rotate_view(d_az=-self._rotation_step)
            event.accept()
        elif key == Qt.Key.Key_Right:
            self._rotate_view(d_az=self._rotation_step)
            event.accept()

        # Moon size
        if key == Qt.Key.Key_M:
            self.toggle_enlarge_moon()
            event.accept()

        # Window control
        if key == Qt.Key.Key_F11:
            self.toggle_fullscreen()
        elif key == Qt.Key.Key_Escape:
            if self.isFullScreen():
                self.showNormal()
        elif key == Qt.Key.Key_Q:
            QApplication.quit()
        else:
            super().keyPressEvent(event)
