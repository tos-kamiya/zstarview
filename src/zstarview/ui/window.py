import logging
logger = logging.getLogger(__name__)

from collections import OrderedDict
from datetime import datetime, timedelta, timezone
import math
import threading
import sys
from typing import Any, Dict, Hashable, Optional, Tuple


from PySide6.QtCore import Qt, QPoint, QTimer, QRect, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase, QIcon, QImage, QPainter, QPainterPath
from PySide6.QtGui import QAction, QKeyEvent, QMouseEvent, QPaintEvent, QResizeEvent
from PySide6.QtWidgets import QMainWindow, QSizeGrip, QApplication, QPushButton, QMenu

import astropy
import numpy as np
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
from ..paths import HatchConfig, CLOUD_HATCH_DEFAULT
from ..render import draw as render_draw
from ..render import draw_sky_disc
from ..types import CelestialData, ViewerData
from ..utils.qt import pil_to_qimage, qimage_to_np_rgba, np_rgba_to_qimage

DEBUG_ECLIPSES = True


class _LRU:
    """Tiny LRU cache for QImage/ndarray blobs keyed by hashable tuples."""
    def __init__(self, max_items: int = 16):
        self.max = max_items
        self._d: "OrderedDict[Hashable, Any]" = OrderedDict()

    def get(self, key: Hashable):
        v = self._d.get(key)
        if v is not None:
            self._d.move_to_end(key)
        return v

    def set(self, key: Hashable, val: Any):
        self._d[key] = val
        self._d.move_to_end(key)
        if len(self._d) > self.max:
            self._d.popitem(last=False)


def cloud_with_hatched_alpha(
    cloud_img: QImage,
    disc_rect: QRect,
    hatch_cfg: HatchConfig,
    hatch_cache: Optional["_LRU"] = None,
) -> QImage:
    """
    雲画像にハッチ状のアルファ抜きを適用する。
    回転AAによる幅ブレを避けるため、NumPyで最終サイズの二値（＋端1px勾配）マスクを生成する。
    """
    assert cloud_img.width() == disc_rect.width() and cloud_img.height() == disc_rect.height(), \
        "Hatch must be generated at the final display size before composition."

    # 出力をPremultipliedに
    out = cloud_img if cloud_img.format() == QImage.Format_ARGB32_Premultiplied \
        else cloud_img.convertToFormat(QImage.Format_ARGB32_Premultiplied)

    w, h = disc_rect.width(), disc_rect.height()

    (tile_px, line_px, ang_key, strength_i, pha_key, edge) = hatch_cfg.as_key()
    pitch = max(1.0, float(tile_px))
    band  = float(max(1, int(line_px)))
    ang   = float(ang_key)
    pha   = float(pha_key)
    strength_i = int(strength_i)

    hatch_key = ("hatch_np", w, h, pitch, band, ang_key, strength_i, pha_key, edge)

    cache = hatch_cache
    hatch = cache.get(hatch_key) if cache is not None else None
    if hatch is None:
        # --- NumPyでハッチαマスクを生成（0..255）
        # 画像中心原点（ピクセル中心合わせ）
        cy = (h - 1) * 0.5
        cx = (w - 1) * 0.5

        # 回転座標 u による縞: (u + phase) % pitch < band
        theta = np.deg2rad(ang)
        ct, st = -math.cos(theta), math.sin(theta)

        # グリッド（float32で十分）
        y = (np.arange(h, dtype=np.float32) - cy)[:, None]
        x = (np.arange(w, dtype=np.float32) - cx)[None, :]

        # 目的の回転方向に合わせて基準座標を作る
        # 画面座標(x,y) → パターン座標(u,v)
        # u =  x*ct + y*st   （u軸に沿って平行縞）
        u = x * ct + y * st
        u = u + pha  # 位相

        # 0..pitch の範囲に折り返し
        umod = np.mod(u, pitch)

        # 中央部（硬い二値）
        inside = umod < band

        # ★ 境界を“わずかに”ソフトに（見た目を立たせるが、幅はブレさせない）
        if edge > 0.0:
            # 0 付近と band 付近で 0..1 の線形ランプ
            near_low  = np.clip((edge - umod) / max(1e-6, edge), 0.0, 1.0)
            near_high = np.clip((edge - (band - umod)) / max(1e-6, edge), 0.0, 1.0)
            edge_ramp = np.maximum(near_low, near_high)  # 端から内側に1px分だけ勾配

            # αは「中央は最大、端は勾配」で作る
            alpha = np.where(inside, 1.0, 0.0)
            # 端の勾配を上書き（中央の1px近辺のみ有効）
            edge_zone = (umod < edge) | (umod > (band - edge))
            alpha = np.where(edge_zone & inside, np.maximum(alpha, edge_ramp), alpha)
        else:
            alpha = np.where(inside, 1.0, 0.0)

        alpha_u8 = (alpha * int(np.clip(strength_i, 0, 255))).astype(np.uint8)

        # 円盤外は完全透明に（プレマルチ前提でRGB=0, A=0）
        # ここで円マスクを作ってαに適用
        rr = min(cx, cy)
        r2 = (x - 0.0)**2 + (y - 0.0)**2
        disc_mask = (r2 <= (rr + 0.25)**2)
        alpha_u8 = np.where(disc_mask, alpha_u8, 0)

        # RGBAにしてQImage化（RGB=0でOK。DestinationOutでαだけ使う）
        hatch_rgba = np.zeros((h, w, 4), dtype=np.uint8)
        hatch_rgba[..., 3] = alpha_u8
        hatch = np_rgba_to_qimage(hatch_rgba)  # Premultipliedに変換される実装ならそれでOK

        if cache is not None:
            cache.set(hatch_key, hatch)

    # --- αカットアウト
    p = QPainter(out)
    p.setCompositionMode(QPainter.CompositionMode_DestinationOut)
    p.drawImage(disc_rect.topLeft(), hatch)  # サイズ一致前提
    p.end()

    return out


def _compose_cloud_over_sky_no_global_cache(
    sky_img: QImage,
    cloud_img_rgba: QImage,
    dest_rect: QRect,
    cloud_opacity: float = 1.0,
    gray_mix: float = 1.0,
) -> QImage:
    """グローバルキャッシュを使わない合成（描画サイズ変化/元画像更新時のみ呼ばれる前提）。"""
    w, h = dest_rect.width(), dest_rect.height()

    # スケール（呼び出し側で合わせている想定だが、保険でサイズ確認）
    if sky_img.width() != w or sky_img.height() != h:
        sky_img = sky_img.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)
    if cloud_img_rgba.width() != w or cloud_img_rgba.height() != h:
        cloud_img_rgba = cloud_img_rgba.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

    sky_np   = qimage_to_np_rgba(sky_img).astype(np.uint8)     # (h,w,4)
    cloud_np = qimage_to_np_rgba(cloud_img_rgba).astype(np.uint8)

    # グレースケール（整数近似）
    r = sky_np[..., 0].astype(np.uint16)
    g = sky_np[..., 1].astype(np.uint16)
    b = sky_np[..., 2].astype(np.uint16)
    gray_u8 = ((77 * r + 150 * g + 29 * b) >> 8).astype(np.uint8)  # (h,w)

    a = (cloud_np[..., 3].astype(np.float32) / 255.0) * float(np.clip(gray_mix, 0.0, 1.0))
    a8 = (np.clip(a, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint16)  # (h,w)
    inv_a8 = (255 - a8).astype(np.uint16)

    sky_rgb_u16  = sky_np[..., :3].astype(np.uint16)
    gray_rgb_u16 = np.repeat(gray_u8[:, :, None], 3, axis=2).astype(np.uint16)

    base_u16 = (inv_a8[:, :, None] * sky_rgb_u16 + a8[:, :, None] * gray_rgb_u16) // 255

    cop = float(np.clip(cloud_opacity, 0.0, 1.0))
    if cop > 0.0:
        add_u16 = (cloud_np[..., :3].astype(np.uint16) * int(round(cop * 255))) // 255
        out_u16 = base_u16 + add_u16
        np.minimum(out_u16, 255, out=out_u16)
    else:
        out_u16 = base_u16

    out = np.empty((h, w, 4), dtype=np.uint8)
    out[..., :3] = out_u16.astype(np.uint8)
    out[..., 3]  = 255

    cy = (h - 1) * 0.5
    cx = (w - 1) * 0.5
    rr = min(cx, cy)  # 画像サイズぴったりの円
    # ピクセル中心合わせ（+0.5/-0.5）でジャギ減
    y = np.arange(h, dtype=np.float32)[:, None]
    x = np.arange(w, dtype=np.float32)[None, :]
    r2 = (x - cx)**2 + (y - cy)**2
    mask = r2 <= (rr + 0.25)**2   # ほんのわずか内側に寄せてリング抑制
    # 外側はアルファ=0、RGBもプレマルチ的に0へ
    out[..., 3][~mask] = 0
    out[..., :3][~mask] = 0

    return np_rgba_to_qimage(out)


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

        self._composited_img: Optional[QImage] = None
        self._composite_key: Optional[Tuple] = None  # ("comp", sky_ck, cloud_ck, w, h, cloud_alpha, hatch_params)
        self._hatch_cache = _LRU(max_items=8)

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

        self._composite_key = None
        self._composited_img = None

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

        self._draw_sky_and_clouds_scaled(painter, geometry)

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

    def _draw_sky_and_clouds_scaled(self, painter: QPainter, geometry):
        """空色ディスクと雲ディスクを一括で合成し、1回の drawImage で描画。"""
        # どちらも無効なら描かない
        if (self._sky_disc_image is None or self.sky_disc_alpha <= 0.0) and \
        (self._cloud_img is None or self.cloud_disc_alpha <= 0.0):
            return

        x = int(geometry.center[0] - geometry.radius)
        y = int(geometry.center[1] - geometry.radius)
        w = h = int(geometry.radius * 2)

        sky_ck   = int(self._sky_disc_image.cacheKey()) if self._sky_disc_image else 0
        cloud_ck = int(self._cloud_img.cacheKey()) if self._cloud_img else 0

        hatch_cfg = CLOUD_HATCH_DEFAULT
        comp_key = ("comp", sky_ck, cloud_ck, w, h, float(self.cloud_disc_alpha), hatch_cfg.as_key())

        # キーが同じなら再計算しない
        if self._composite_key != comp_key or self._composited_img is None:
            # 1) 必要サイズへスケーリング
            def _scaled(qimg: Optional[QImage]) -> Optional[QImage]:
                if qimg is None:
                    return None
                if qimg.width() == w and qimg.height() == h:
                    return qimg
                return qimg.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

            sky_s   = _scaled(self._sky_disc_image)
            cloud_s = _scaled(self._cloud_img)

            # 2) 雲のハッチ（ある場合のみ）
            if cloud_s is not None and self.cloud_disc_alpha > 0.0:
                disc = QRect(0, 0, w, h)
                cloud_s = cloud_with_hatched_alpha(cloud_s, disc, hatch_cfg, self._hatch_cache)

            # 3) 合成
            if cloud_s is None or self.cloud_disc_alpha <= 0.0:
                # 雲なし → 空だけ
                composited = sky_s if sky_s is not None else QImage(w, h, QImage.Format_ARGB32_Premultiplied)
                if composited is None or composited.isNull():
                    # 念のため透過で埋める
                    composited = QImage(w, h, QImage.Format_ARGB32_Premultiplied)
                    composited.fill(Qt.transparent)
            else:
                # 雲あり → グレーブレンド＋雲のRGB加算（従来のロジック）
                # compose_cloud_over_sky_scaled を「キャッシュ非依存」版に差し替えたものを使う
                composited = _compose_cloud_over_sky_no_global_cache(
                    sky_img = sky_s if sky_s is not None else QImage(w, h, QImage.Format_ARGB32_Premultiplied),
                    cloud_img_rgba = cloud_s,
                    dest_rect = QRect(0, 0, w, h),
                    cloud_opacity = float(self.cloud_disc_alpha * 0.5),
                    gray_mix = 1.0,
                )

            self._composited_img = composited
            self._composite_key = comp_key

        # 4) 描画
        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        painter.drawImage(QRect(x, y, w, h), self._composited_img)  # クリップしない（合成画像に円形アルファが入っている）
        painter.restore()

    def _on_sky_data_calculated(self, payload):
        """Handles the sky data calculated signal."""
        self.set_sky_data(payload["celestial"])
        self._sky_disc_image = payload["sky_disc"]

        self._composite_key = None
        self._composited_img = None

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

                self._composite_key = None
                self._composited_img = None
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
