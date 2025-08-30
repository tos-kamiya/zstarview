import logging
logger = logging.getLogger(__name__)

from datetime import datetime, timedelta, timezone
import threading
import sys
from typing import Optional, Tuple

from PySide6.QtCore import Qt, QPoint, QTimer, QRect, Signal
from PySide6.QtGui import QFont, QFontDatabase, QIcon, QImage, QPainter, QPainterPath
from PySide6.QtGui import QAction, QKeyEvent, QMouseEvent, QPaintEvent, QResizeEvent
from PySide6.QtWidgets import QMainWindow, QSizeGrip, QApplication, QPushButton, QMenu

import astropy
import numpy as np
from PIL import Image
import polars as pl

from ..clouddisc import CloudDisc, CloudDiscConfig
from ..clouddisc import VisibilityError, CloudDiscError, DataNotFoundError, DownloadError, RenderError

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
from ..utils.qt import pil_to_qimage, qimage_to_pil

DEBUG_ECLIPSES = True


def tint_cloud_gray_with_skyR_qimage(
    sky_qimg: QImage,
    cloud_gray_qimg: QImage,
    factor: float = 0.3,
    warm_rgb=(255, 156, 82),
    keep_alpha: bool = False,
) -> QImage:
    """
    グレー雲画像(L=雲の濃さ)に、空画像のR成分で暖色をブレンドして着色する。
    出力はRGBA (Aはデフォルト255)。Plus合成を想定。

    Args:
        sky_qimg: 背景の空(QImage)。Rチャンネルを参照する。
        cloud_gray_qimg: 雲のグレー画像(QImage.Format_Grayscale8相当)。
                         画素値0..255が雲の強さ(白=強)。
        factor: 空の赤さ→暖色の寄り具合の係数。0..1くらいで調整。
        warm_rgb: 暖色のRGB。白(255,255,255)とこの色の間を補間する。
        keep_alpha: cloud_gray_qimgにαがあるならそれをAに使う（なければ255）。

    Returns:
        QImage (RGBA Premultiplied推奨)
    """
    # --- QImage -> PIL
    sky_pil   = qimage_to_pil(sky_qimg).convert("RGB")
    cloud_pil = qimage_to_pil(cloud_gray_qimg)
    Wc, Hc = cloud_pil.size

    # 空を雲サイズへ
    sky_resized = sky_pil.resize((Wc, Hc), Image.BILINEAR)
    sky_arr = np.array(sky_resized, dtype=np.uint8)  # (H,W,3)

    # 雲の強さ（L）。必要ならαも取得
    if cloud_pil.mode == "LA":
        L_arr = np.array(cloud_pil.getchannel("L"), dtype=np.uint8)
        A_arr = np.array(cloud_pil.getchannel("A"), dtype=np.uint8)
    else:
        L_arr = np.array(cloud_pil, dtype=np.uint8)
        A_arr = None

    # ---- 色づけロジック ----
    # 空のRを0..1に正規化し、factorでスケール → k: 白↔暖色の混合率
    sky_R = sky_arr[..., 0].astype(np.float32) / 255.0
    k = np.clip(factor * sky_R, 0.0, 1.0)[..., None]  # (H,W,1)

    white = np.array([255, 255, 255], dtype=np.float32)
    warm  = np.array(warm_rgb,       dtype=np.float32)
    base_color = white * (1.0 - k) + warm * k  # (H,W,3) 0..255

    # 雲の明るさL(0..255)でスケール（L=0で黒、L=255でbase_colorそのまま）
    Lf = (L_arr.astype(np.float32) / 255.0)[..., None]  # (H,W,1)
    rgb = (base_color * Lf).clip(0, 255).astype(np.uint8)

    # アルファ：既定は255（合成時はsetOpacityやCompositionで調整）
    if keep_alpha and A_arr is not None:
        A = A_arr
    else:
        A = np.full((Hc, Wc), 255, dtype=np.uint8)

    out = np.dstack([rgb, A])  # (H,W,4)
    tinted_pil = Image.fromarray(out, "RGBA")
    return pil_to_qimage(tinted_pil, premultiplied=True)

from typing import Optional
from PySide6.QtGui import QImage, QPainter, QBrush, QPen, QColor, QPixmap, QPainterPath
from PySide6.QtCore import Qt, QRect, QPoint
import math

# 互換：alphaChannel未使用の実装で上書き
def cloud_with_hatched_alpha_no_alphaChannel(
    cloud_img: QImage,
    disc_rect: QRect,
    tile_px: int = 12,
    line_px: int = 6,
    angle_deg: float = 45.0,
    strength: int = 255,
    # thickness_scale: float = 1.0,
    phase_px: float = 0.0,
) -> QImage:
    out = cloud_img.convertToFormat(QImage.Format_ARGB32_Premultiplied)
    w, h = disc_rect.width(), disc_rect.height()

    # ストライプ（透明RGBA）
    hatch = QImage(w, h, QImage.Format_ARGB32_Premultiplied)
    hatch.fill(Qt.transparent)

    # 円盤クリップ付きのペイント面
    painter = QPainter(hatch)
    painter.setRenderHint(QPainter.Antialiasing, True)
    path = QPainterPath(); path.addEllipse(0, 0, w, h)
    painter.setClipPath(path)

    # 連続直線ストライプを描く
    L = int(math.hypot(w, h))
    painter.save()
    painter.translate(w/2, h/2)
    painter.rotate(angle_deg)
    painter.setPen(Qt.NoPen)
    color = QColor(0, 0, 0, max(0, min(255, strength)))
    pitch = max(1, int(tile_px))
    band  = max(1, line_px)
    # band  = max(1, int(round(line_px * thickness_scale)))
    y = -L + (phase_px % pitch)
    while y <= L:
        painter.fillRect(-L, int(y - band/2), 2*L, band, color)
        y += pitch
    painter.restore()
    painter.end()

    # “抜き”適用：DestinationOut（= 雲のαからhatch分を減算）
    p = QPainter(out)
    p.setCompositionMode(QPainter.CompositionMode_DestinationOut)
    p.drawImage(disc_rect.topLeft(), hatch)
    p.end()

    return out


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
        cloud_disc_alpha: float = 0.6,
        enlarge_moon: bool = False,
        star_base_radius: float = 8.0,
        vmag_limit: float = 6.0,
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
        self.cloud_disc_alpha: float = max(0.0, min(1.0, cloud_disc_alpha))
        if delta_t.total_seconds() != 0.0:
            self.cloud_disc_alpha = 0.0  # can not obtain time-shifted cloud data, but just the current one

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
        self._cloud_base_size: int = 512
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
            alt_min_deg=-1.0,
            search_back_minutes=120,
        )
        try:
            self._clouddisc = CloudDisc(clouddisc_config)
        except Exception as e:
            logger.warning(f"CloudDisc init failed: {e}")

        # 初期ロード開始
        self.start_background_sky_data_update(is_initial_load=True)
        # 雲も初回生成（非同期）
        if self._clouddisc and self.cloud_disc_alpha > 0.0:
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
        if self._cloud_img is None or self.cloud_disc_alpha <= 0.0:
            return

        img = self._cloud_img
        # if self._sky_disc_image is not None:
        #     img = tint_cloud_gray_with_skyR_qimage(self._sky_disc_image, self._cloud_img, factor=0.5)

        x = int(geometry.center[0] - geometry.radius)
        y = int(geometry.center[1] - geometry.radius)
        w = h = int(geometry.radius * 2)

        disc = QRect(0, 0, img.width(), img.height())

        # 元の img はそのまま、新しい画像 new_img を作る
        new_img = cloud_with_hatched_alpha_no_alphaChannel(
            img, disc,
            tile_px=8,               # ストライプ間隔（ピッチ）
            line_px=4,                # 基本太さ
            angle_deg=45,
            strength=255,
            phase_px=0.0,             # 必要なら微調整
        )

        path = QPainterPath()
        path.addEllipse(QRect(x, y, w, h))
        painter.save()
        painter.setClipPath(path)
        painter.setOpacity(self.cloud_disc_alpha)
        painter.setCompositionMode(QPainter.CompositionMode_Plus)
        painter.drawImage(QRect(x, y, w, h), new_img)
        painter.restore()

    def _on_sky_data_calculated(self, payload):
        """Handles the sky data calculated signal."""
        self.set_sky_data(payload["celestial"])
        self._sky_disc_image = payload["sky_disc"]

        # Emit the "initial_data_loaded" signal only once:
        if not self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.start(SKY_UPDATE_INTERVAL * 1000)
            # --- start clouds periodic timer once window is up
            if self._clouddisc and self.cloud_disc_alpha > 0.0 and not self._cloud_update_timer.isActive():
                self._cloud_update_timer.start()
            self.initial_data_loaded.emit()

        self._is_sky_data_calculation_running = False

    def request_sky_data_update(self):
        """Request a celestial data/sky disc image update, but only if no update is currently running."""
        f = self.start_background_sky_data_update()
        if not f:
            logger.warning("Warning: celestial data/sky disc image updating canceled.")

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
                        logger.debug("Solar eclipse detected")
                    if body.name == "moon" and body.lunar_eclipse_info.is_eclipse:
                        logger.debug("Lunar eclipse detected")

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
                    mask_fov_deg=93,
                )

            payload = {"celestial": celestial_data, "sky_disc": sky_disc_img}
            self.sky_data_calculated.emit(payload)
        except Exception as e:
            logger.error(f"Error in background update thread: {e}")
            import traceback
            traceback.print_exc()

    def start_background_sky_data_update(self, is_initial_load: bool = False) -> bool:
        """Starts a background thread to update the sky data."""
        if self._is_sky_data_calculation_running:
            return False
        self._is_sky_data_calculation_running = True

        if is_initial_load:
            logger.info("Calculating initial sky data...")
        else:
            logger.info("Updating sky data...")
        thread = threading.Thread(target=self.update_sky_data_in_background)
        thread.daemon = True
        thread.start()
        return True

    def start_background_cloud_update(self, reason: str = "manual"):
        """Kick a background cloud generation if conditions allow."""
        if not (self._clouddisc and self.cloud_disc_alpha > 0.0):
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

            if reason == "initial":
                logger.info("Fetching initial could data...")
            else:
                logger.info("Updating cloud data...")

            try:
                pil_img, meta = self._clouddisc.render_now(
                    lat=lat, lon=lon, alt=float(alt), az=float(az),
                    radius_px=self._cloud_base_size, edge_fov_deg=90, mask_fov_deg=93,
                )
                qimg = pil_to_qimage(pil_img)
                self._cloud_img = qimg
                self._cloud_meta = meta
            except VisibilityError as e:
                logger.error("Invalid params for cloud-disc image generation: %s", e)
            except DownloadError as e:
                logger.warning("Network/S3 download error (transient=%s): %s", e.transient, e)
            except DataNotFoundError as e:
                logger.info("Satellite data not found in search window: %s", e)
            except RenderError as e:
                logger.error("Failed to decode/render satellite data: %s", e)
            except CloudDiscError as e:  # 保険
                logger.error("Unexpected clouddisc error: %s", e)

            self._last_cloud_az = az
            self._last_cloud_time_utc = datetime.now(timezone.utc)

            # 画面更新
            self.update()
        except Exception as e:
            logger.error(f"Warn: clouds update failed: {e}")
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
