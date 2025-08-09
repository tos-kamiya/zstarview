from datetime import datetime, timedelta, timezone
import threading
import sys
from typing import Any, Dict, List, Optional, Tuple

from PyQt5.QtCore import Qt, QPoint, QTimer
from PyQt5.QtGui import (
    QFont,
    QFontDatabase,
    QIcon,
    QPainter,
    QResizeEvent,
    QPaintEvent,
    QMouseEvent,
    QKeyEvent,
)
from PyQt5.QtWidgets import QMainWindow, QSizeGrip, QApplication

import astropy

from ..paths import EMOJI_FONT_PATH, TEXT_FONT_SIZE, APP_ICON_FILE
from ..types import SkyData, ViewerData
from ..astro import (
    calculate_visible_stars,
    calculate_visible_planets,
    calculate_celestial_equator_points,
    calculate_ecliptic_points,
    calculate_horizon_points,
)
from ..render import draw as render_draw
from ..render import draw_o as render_draw_o


class SkyWindow(QMainWindow):
    # Using Qt signal objects requires attribute creation at runtime; avoid type hints here
    from PyQt5.QtCore import pyqtSignal

    data_updated = pyqtSignal(object)
    initial_data_loaded = pyqtSignal()

    def __init__(
        self,
        city_name: str,
        city_data: Tuple[float, float, str],
        star_catalog: List[Dict[str, Any]],
        delta_t: timedelta,
        enlarge_moon: bool,
        star_base_radius: float,
        view_center: Tuple[float, float],
    ):
        super().__init__()
        self.setWindowIcon(QIcon(APP_ICON_FILE))

        lat, lon, tz_name = city_data
        self.viewer_data = ViewerData(
            location=(lat, lon),
            timezone_name=tz_name,
            city_name=city_name,
            view_center=view_center,
        )

        self.star_catalog = star_catalog
        self.delta_t = delta_t
        self.enlarge_moon = enlarge_moon
        self.star_base_radius = star_base_radius

        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setWindowFlags(Qt.FramelessWindowHint)
        self.setWindowTitle(f"Zenith Star View - {self.viewer_data.city_name.title()}")
        self.setGeometry(100, 100, 800, 800)

        # Size grip
        self.size_grip = QSizeGrip(self)
        self.size_grip.setFixedSize(24, 24)
        self.size_grip.raise_()

        # Fonts for drawing
        font_id = QFontDatabase.addApplicationFont(EMOJI_FONT_PATH)
        font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
        self.emoji_font = QFont(font_family, TEXT_FONT_SIZE + 4)
        self.text_font = QFont("Arial", TEXT_FONT_SIZE)

        self.setMouseTracking(True)
        self.setMinimumSize(400, 400)
        self.sky_data: Optional[SkyData] = None
        self.mouse_pos: Optional[QPoint] = None

        self.data_updated.connect(self.on_data_updated)
        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.start_background_update)

        self.start_background_update(is_initial_load=True)

    def resizeEvent(self, event: QResizeEvent):
        grip_size = self.size_grip.size()
        self.size_grip.move(self.width() - grip_size.width(), self.height() - grip_size.height())
        super().resizeEvent(event)

    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.LeftButton:
            self._drag_active = True
            self._drag_pos = event.globalPos() - self.frameGeometry().topLeft()
            event.accept()

    def leaveEvent(self, event):
        self.mouse_pos = None
        self.update()
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent):
        if getattr(self, "_drag_active", False) and event.buttons() & Qt.LeftButton:
            self.move(event.globalPos() - self._drag_pos)
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
            painter.setFont(QFont("Arial", 16))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Loading sky data...")
            return

        alt = self.viewer_data.view_center[0]
        geometry = render_draw.get_screen_geometry(self.width(), self.height(), alt)

        painter.setCompositionMode(QPainter.CompositionMode_Clear)
        painter.fillRect(self.rect(), Qt.transparent)

        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)

        render_draw.draw_radial_background(painter, self.rect(), geometry)

        render_draw.draw_sky_reference_lines(painter, geometry, self.sky_data)
        render_draw.draw_direction_labels(painter, geometry, self.viewer_data.view_center, self.text_font)

        render_draw_o.draw_stars_fully_vectorized(painter, geometry, self.sky_data, self.viewer_data, self.star_base_radius)
        # render_draw.draw_stars(painter, geometry, self.sky_data, self.viewer_data, self.star_base_radius)

        render_draw.draw_planets(painter, geometry, self.sky_data, self.viewer_data, self.enlarge_moon, self.emoji_font)

        highlighted_object = None
        if self.mouse_pos is not None:
            highlighted_object = render_draw.find_highlighted_object(self.sky_data, self.viewer_data, self.mouse_pos, geometry)
        render_draw.draw_overlay_info(painter, self.sky_data, self.viewer_data, highlighted_object, self.text_font)

    def on_data_updated(self, sky_data: SkyData):
        self.set_sky_data(sky_data)
        if not self.update_timer.isActive():
            self.update_timer.start(5 * 60 * 1000)
            self.initial_data_loaded.emit()

    def update_sky_data_in_background(self):
        try:
            now = datetime.now(timezone.utc) + self.delta_t
            time_obj = astropy.time.Time(now)
            lat, lon = self.viewer_data.location
            stars, loc = calculate_visible_stars(self.star_catalog, lat, lon, time_obj, self.viewer_data.view_center)
            planets = calculate_visible_planets(lat, lon, time_obj, self.viewer_data.view_center)
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
            self.data_updated.emit(sky_data)
        except Exception as e:
            print(f"Error in background update thread: {e}", file=sys.stderr)
            import traceback

            traceback.print_exc()

    def start_background_update(self, is_initial_load: bool = False):
        if is_initial_load:
            print("Calculating initial sky data..")
        else:
            print("Updating sky data...")
        thread = threading.Thread(target=self.update_sky_data_in_background)
        thread.daemon = True
        thread.start()

    def keyPressEvent(self, event: QKeyEvent):
        if event and event.key() == Qt.Key.Key_F11:
            if self.isFullScreen():
                self.showNormal()
            else:
                self.showFullScreen()
        elif event and event.key() == Qt.Key.Key_Escape:
            if self.isFullScreen():
                self.showNormal()
        elif event and event.key() == Qt.Key.Key_Q:
            QApplication.quit()
        else:
            super().keyPressEvent(event)
