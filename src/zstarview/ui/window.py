from datetime import datetime, timedelta, timezone
import threading
import sys
from typing import Optional, Tuple

from PySide6.QtCore import Qt, QPoint, QTimer, QRect, Signal
from PySide6.QtGui import QFont, QFontDatabase, QIcon, QImage, QPainter, QPainterPath
from PySide6.QtGui import QAction, QKeyEvent, QMouseEvent, QPaintEvent, QResizeEvent
from PySide6.QtWidgets import QMainWindow, QSizeGrip, QApplication, QPushButton, QMenu

import astropy
import polars as pl

from ..clouddisc import CloudDisc, CloudDiscConfig

from ..__about__ import __version__
from ..astro import (
    calculate_visible_stars,
    calculate_planets,
    calculate_celestial_equator_points,
    calculate_ecliptic_points,
    calculate_horizon_points,
    eclipse_factor_from_info,
)
from .draggable_window import DraggableWindow
from ..paths import CACHE_PATH
from ..paths import EMOJI_FONT_PATH, EMOJI_FONT_SIZE, TEXT_FONT_PATH, TEXT_FONT_SIZE
from ..paths import APP_ICON_FILE, GUI_BUTTON_SIZE, GUI_MENU_TEXT_COLOR, WINDOW_HEIGHT, WINDOW_HEIGHT, WINDOW_WIDTH
from ..paths import SKY_UPDATE_INTERVAL, CLOUD_UPDATE_INTERVAL
from ..render import draw as render_draw
from ..render import draw_sky_disc
from ..types import CelestialData, ViewerData
from ..utils.qt import pil_to_qimage

DEBUG_ECLIPSES = True


class SkyWindow(DraggableWindow):
    sky_data_calculated = Signal(object)
    initial_data_loaded = Signal()

    def __init__(
        self,
        city_name: str,
        city_data: Tuple[float, float, str],
        view_center: Tuple[float, float],
        star_catalog: pl.DataFrame,
        delta_t: timedelta = timedelta(days=0, hours=0),
        sky_disc_alpha: float = 0.3,
        enlarge_moon: bool = False,
        star_base_radius: float = 8.0,
        vmag_limit: float = 6.0,
        clouds_alpha: float = 1.0,
    ):
        """Initializes the SkyWindow."""
        super().__init__()
        self._rotation_step: float = 5.0  # degrees

        self.star_catalog = star_catalog
        self.delta_t = delta_t
        self.sky_disc_alpha = sky_disc_alpha
        self.enlarge_moon = enlarge_moon
        self.star_base_radius = star_base_radius
        self.vmag_limit = vmag_limit

        # --- added: clouds alpha
        self.clouds_alpha: float = max(0.0, min(1.0, clouds_alpha))
        if delta_t.total_seconds() != 0.0:
            self.clouds_alpha = 0.0  # can not obtain time-shifted cloud data, but just the current one

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

        self.setMouseTracking(True)
        self.mouse_pos: Optional[QPoint] = None

        # Size grip
        self.size_grip = QSizeGrip(self)
        self.size_grip.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.size_grip.raise_()

        # Menu
        self._action_enlarge_moon = None
        self._add_hamburger_menu()

        self.add_drag_exclusions([self.menu_button, self.size_grip])

        # Fonts for drawing
        emoji_font_id = QFontDatabase.addApplicationFont(EMOJI_FONT_PATH)
        emoji_font_family = QFontDatabase.applicationFontFamilies(emoji_font_id)[0]
        self.emoji_font = QFont(emoji_font_family, EMOJI_FONT_SIZE)
        text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
        text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
        self.text_font = QFont(text_font_family, TEXT_FONT_SIZE)

        self.sky_data_calculated.connect(self._on_sky_data_calculated)
        self._is_sky_data_calculation_running: bool = False
        self._sky_data_update_timer = QTimer(self)
        self._sky_data_update_timer.timeout.connect(self.start_background_sky_data_update)

        self.celestial_data: Optional[CelestialData] = None

        self._sky_disc_base_size: int = 1024
        self._sky_disc_image: Optional[QImage] = None

        # --- added: CloudDisc & timers/flags/cache
        self._cloud_base_size: int = 128
        self._cloud_img: Optional[QImage] = None
        self._cloud_meta: Optional[dict] = None
        self._is_cloud_update_running: bool = False
        self._cloud_update_pending: bool = False
        self._last_cloud_az: Optional[float] = None
        self._last_cloud_time_utc: Optional[datetime] = None

        self._cloud_update_timer = QTimer(self)
        self._cloud_update_timer.setInterval(CLOUD_UPDATE_INTERVAL * 1000)
        self._cloud_update_timer.timeout.connect(lambda: self.start_background_cloud_update(reason="timer"))

        self._clouddisc = None
        clouddisc_config = CloudDiscConfig(
            cache_dir=CACHE_PATH,
            sat_priority=("AUTO",),
            bt_warm_k=310.0,
            bt_cold_k=190.0,
            gamma=1.6,
            alt_min_deg=-2.0,
            search_back_minutes=120,
            edge_antialias=True,
        )
        try:
            self._clouddisc = CloudDisc(clouddisc_config)
        except Exception as e:
            print(f"[warn] CloudDisc init failed: {e}", file=sys.stderr)

        # 初期ロード開始
        self.start_background_sky_data_update(is_initial_load=True)
        # 雲も初回生成（非同期）
        if self._clouddisc and self.clouds_alpha > 0.0:
            self.start_background_cloud_update(reason="initial")

    def _add_hamburger_menu(self):
        """Adds a hamburger menu to the window."""
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

    def resizeEvent(self, event: QResizeEvent):
        """Handles the window resize event."""
        grip_size = self.size_grip.size()
        self.size_grip.move(self.width() - grip_size.width(), self.height() - grip_size.height())

        button_size = self.menu_button.size()
        self.menu_button.move(self.width() - button_size.width() - 8, 8)

        super().resizeEvent(event)

    def show_menu(self):
        """Shows the hamburger menu."""
        menu_pos = self.menu_button.mapToGlobal(QPoint(0, self.menu_button.height()))
        self.menu.exec(menu_pos)

    def toggle_enlarge_moon(self):
        """Toggles the moon enlargement."""
        self.enlarge_moon = not self.enlarge_moon
        if self._action_enlarge_moon is not None and self._action_enlarge_moon.isChecked() != self.enlarge_moon:
            self._action_enlarge_moon.setChecked(self.enlarge_moon)
        self.update()  # Redraw

    def toggle_fullscreen(self):
        """Toggles fullscreen mode."""
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def leaveEvent(self, event):
        """Handles the mouse leave event."""
        self.mouse_pos = None
        self.update()
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent):
        """Handles the mouse move event."""
        self.mouse_pos = event.pos()
        self.update()
        event.accept()

    def set_star_base_radius(self, star_base_radius: float):
        """Sets the base radius for stars."""
        self.star_base_radius = star_base_radius

    def set_enlarge_moon(self, enlarge_moon: bool):
        """Sets the enlarge moon flag."""
        self.enlarge_moon = enlarge_moon

    def set_sky_data(self, data: CelestialData):
        """Sets the celestial data."""
        self.celestial_data = data
        self.update()

    def paintEvent(self, event: QPaintEvent):
        """Handles the paint event."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)

        if not self.celestial_data:
            painter.setPen(Qt.GlobalColor.white)
            painter.setFont(self.text_font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Loading celestial data...")
            return

        alt = self.viewer_data.view_center[0]
        geometry = render_draw.get_screen_geometry(self.width(), self.height(), alt)

        highlighted_object = None
        if self.mouse_pos is not None:
            highlighted_object = render_draw.find_highlighted_object(self.celestial_data, self.viewer_data, self.mouse_pos, geometry)

        painter.setCompositionMode(QPainter.CompositionMode_Clear)
        painter.fillRect(self.rect(), Qt.transparent)

        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)

        render_draw.draw_radial_background(painter, self.rect(), geometry)

        # 背景の空色ディスク
        self._draw_sky_disc_scaled(painter, geometry)

        # --- added: clouds disc (over the sky color, under stars)
        self._draw_clouds_scaled(painter, geometry)

        render_draw.draw_sky_reference_lines(painter, geometry, self.celestial_data)
        render_draw.draw_direction_labels(painter, geometry, self.viewer_data.view_center, self.text_font)
        render_draw.draw_zenith_marker(painter, geometry, self.viewer_data.view_center)

        render_draw.draw_stars(painter, geometry, self.celestial_data, self.viewer_data, self.star_base_radius)

        enlarge_moon = self.enlarge_moon
        if highlighted_object is not None:
            obj = highlighted_object[0]
            name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
            enlarge_moon = enlarge_moon or name == "moon"
        render_draw.draw_planets(painter, geometry, self.celestial_data, self.viewer_data, enlarge_moon, self.emoji_font)

        render_draw.draw_overlay_info(
            painter,
            self.celestial_data,
            self.viewer_data,
            self.vmag_limit,
            enlarge_moon,
            highlighted_object,
            self.text_font,
        )

    def _draw_sky_disc_scaled(self, painter: QPainter, geometry):
        """Draws the scaled sky disc."""
        img = self._sky_disc_image
        if img is None or self.sky_disc_alpha <= 0.0:
            return

        x = int(geometry.center[0] - geometry.radius)
        y = int(geometry.center[1] - geometry.radius)
        w = h = int(geometry.radius * 2)

        painter.drawImage(QRect(x, y, w, h), img)

    # --- added: draw clouds
    def _draw_clouds_scaled(self, painter: QPainter, geometry):
        """Draws the scaled cloud disc with alpha."""
        if self._cloud_img is None or self.clouds_alpha <= 0.0:
            return

        x = int(geometry.center[0] - geometry.radius)
        y = int(geometry.center[1] - geometry.radius)
        w = h = int(geometry.radius * 2)

        painter.save()
        painter.setOpacity(self.clouds_alpha)
        path = QPainterPath()
        path.addEllipse(QRect(x, y, w, h))
        painter.setClipPath(path)
        painter.drawImage(QRect(x, y, w, h), self._cloud_img)
        painter.restore()

    def _on_sky_data_calculated(self, payload):
        """Handles the sky data calculated signal."""
        self.set_sky_data(payload["celestial"])
        self._sky_disc_image = payload["sky_disc"]

        # Emit the "initial_data_loaded" signal only once:
        if not self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.start(SKY_UPDATE_INTERVAL * 1000)
            # --- start clouds periodic timer once window is up
            if self._clouddisc and self.clouds_alpha > 0.0 and not self._cloud_update_timer.isActive():
                self._cloud_update_timer.start()
            self.initial_data_loaded.emit()

        self._is_sky_data_calculation_running = False

    def request_sky_data_update(self):
        """Request a celestial data/sky disc image update, but only if no update is currently running."""
        f = self.start_background_sky_data_update()
        if not f:
            print("Warning: celestial data/sky disc image updating canceled.")

    def update_sky_data_in_background(self):
        """Updates the sky data in a background thread."""
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
            celestial_data = CelestialData(
                time=time_obj,
                planets=planets,
                stars=stars,
                celestial_equator_points=celestial_equator_points,
                ecliptic_points=ecliptic_points,
                horizon_points=horizon_points,
            )

            if DEBUG_ECLIPSES:
                for body in planets:
                    if body.name == "sun" and body.solar_eclipse_info.is_eclipse:
                        print("Info: Solar eclipse detected")
                    if body.name == "moon" and body.lunar_eclipse_info.is_eclipse:
                        print("Info: Lunar eclipse detected")

            sun_altaz = None
            solar_eclipse_info = None
            for body in planets:
                if body.name == "sun":
                    sun_altaz = (body.alt, body.az)
                    solar_eclipse_info = body.solar_eclipse_info
                    break

            sky_disc_img = None
            if self.sky_disc_alpha > 0.0 and sun_altaz is not None:
                base = self._sky_disc_base_size
                fixed_geom = render_draw.get_screen_geometry(base, base, base // 2)

                ef = eclipse_factor_from_info(solar_eclipse_info)
                sky_disc_img = draw_sky_disc.draw_sky_color_disc(
                    fixed_geom,
                    self.viewer_data.view_center,
                    sun_altaz,
                    alpha=self.sky_disc_alpha,
                    eclipse_factor=ef,
                )

            payload = {"celestial": celestial_data, "sky_disc": sky_disc_img}
            self.sky_data_calculated.emit(payload)
        except Exception as e:
            print(f"Error in background update thread: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()

    def start_background_sky_data_update(self, is_initial_load: bool = False) -> bool:
        """Starts a background thread to update the sky data."""
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

    def start_background_cloud_update(self, reason: str = "manual"):
        """Kick a background cloud generation if conditions allow."""
        if not (self._clouddisc and self.clouds_alpha > 0.0):
            return

        if self._is_cloud_update_running:
            # 連打時の取りこぼしを避けるためペンディングを立てる
            self._cloud_update_pending = True
            return

        self._is_cloud_update_running = True
        t = threading.Thread(target=self._update_cloud_in_background, args=(reason,))
        t.daemon = True
        t.start()

    def _update_cloud_in_background(self, reason: str):
        """Generate cloud disc via CloudDisc in a worker thread."""
        try:
            lat, lon = self.viewer_data.location
            alt, az = self.viewer_data.view_center

            print(f"[clouds] updating (reason={reason}, alt={alt:.1f}, az={az:.1f})")

            # CloudDisc は内部で時刻を選んでくれる render_now を利用
            pil_img, meta = self._clouddisc.render_now(
                lat=lat, lon=lon, alt=float(alt), az=float(az), radius_px=self._cloud_base_size
            )

            qimg = pil_to_qimage(pil_img)

            # 反映
            self._cloud_img = qimg
            self._cloud_meta = meta
            self._last_cloud_az = az
            self._last_cloud_time_utc = datetime.now(timezone.utc)

            # 画面更新
            self.update()
        except Exception as e:
            print(f"[clouds] update failed: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
        finally:
            self._is_cloud_update_running = False
            # ペンディングが立っていたらもう一度キック
            if self._cloud_update_pending:
                self._cloud_update_pending = False
                self.start_background_cloud_update(reason="pending")

    def _rotate_view(self, d_alt: float = 0.0, d_az: float = 0.0):
        """Rotates the view."""
        alt, az = self.viewer_data.view_center
        new_alt = max(0.0, min(90.0, alt + d_alt))
        new_az = (az + d_az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)

        self.request_sky_data_update()

        self.start_background_cloud_update(reason="az-change")

    def keyPressEvent(self, event: QKeyEvent):
        """Handles key press events."""
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
