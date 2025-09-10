# -*- coding: utf-8 -*-
"""
The main window of the ZStarView application.

This module defines the `SkyWindow` class, which is the primary user interface
for the application. It handles rendering the celestial objects, sky background,
clouds, and all user interactions like rotation, zooming, and object highlighting.
"""
import logging
import math
import sys
import threading
from datetime import datetime, timedelta, timezone
from typing import Any, Optional, Tuple

import astropy
import numpy as np
import polars as pl
from PySide6.QtCore import QEvent, QPoint, QRect, Qt, QTimer, Signal
from PySide6.QtGui import (
    QAction,
    QColor,
    QFont,
    QFontDatabase,
    QIcon,
    QImage,
    QKeyEvent,
    QMouseEvent,
    QPainter,
    QPaintEvent,
    QPixmap,
    QResizeEvent,
)
from PySide6.QtWidgets import QApplication, QMainWindow, QMenu, QPushButton, QSizeGrip

from ..__about__ import __version__
from ..astro import (
    calculate_celestial_equator_points,
    calculate_ecliptic_points,
    calculate_horizon_points,
    calculate_planets,
    calculate_visible_stars,
    eclipse_factor_from_info,
)
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
    CACHE_PATH,
    CLOUD_HATCH_DEFAULT,
    CLOUD_SHELL_KM,
    CLOUD_UPDATE_INTERVAL,
    EMOJI_FONT_PATH,
    EMOJI_FONT_SIZE,
    GUI_BUTTON_SIZE,
    GUI_MENU_TEXT_COLOR,
    HatchConfig,
    SKY_UPDATE_INTERVAL,
    TEXT_FONT_PATH,
    TEXT_FONT_SIZE,
    STATUS_LINE_FONT_SIZE,
    WINDOW_HEIGHT,
    WINDOW_WIDTH,
)
from ..render import draw as render_draw
from ..render import draw_sky_disc
from ..types import CelestialData, ScreenGeometry, ViewerData
from ..utils.qt import np_rgba_to_qimage, pil_to_qimage, qimage_to_np_rgba
from .draggable_window import DraggableWindow

logger = logging.getLogger(__name__)


DEBUG_ECLIPSES = True


def make_hatch_tile_qimage(W: int, H: int, line_px: int, strength: int) -> QImage:
    """
    Generates a hatch tile image for masking.

    Creates a QImage with a diagonal line pattern (hatch). This tile can be
    used to create a stylized, less realistic look for UI elements like clouds.
    The format is ARGB32_Premultiplied for direct use with QPainter.

    Args:
        W: The width of the tile in pixels.
        H: The height of the tile in pixels.
        line_px: The thickness of the hatch lines in pixels.
        strength: The alpha value (0-255) of the hatch lines.

    Returns:
        A QImage containing the generated hatch pattern.
    """
    norm = math.sqrt(W * W + H * H)
    P = W * H
    band_u = max(1, int(round(line_px * norm)))

    # Create a coordinate grid
    xs = np.arange(W, dtype=np.int32)[None, :]
    ys = np.arange(H, dtype=np.int32)[:, None]

    # Calculate the distance from the diagonal line for each pixel
    u = H * xs - W * ys
    u_mod = np.mod(u, P)
    dist = np.minimum(u_mod, P - u_mod)
    mask = dist <= (band_u / 2)

    # Create the image buffer
    arr = np.zeros((H, W, 4), dtype=np.uint8)
    arr[..., 0:3] = 0  # Black color
    arr[..., 3] = 0  # Transparent background
    arr[..., 3][mask] = np.uint8(np.clip(strength, 0, 255))  # Apply alpha to lines

    # Convert numpy array to QImage
    buf = arr.tobytes()
    qimg = QImage(buf, W, H, QImage.Format_ARGB32_Premultiplied)
    return qimg.copy()


def cloud_with_hatched_alpha(cloud_img: QImage, hatch_cfg: HatchConfig) -> QImage:
    """
    Applies a hatch mask to a cloud image's alpha channel.

    This function takes a cloud image and uses a generated hatch pattern to
    modify its transparency. The hatch pattern is drawn using the
    `CompositionMode_DestinationOut` mode, effectively "erasing" parts of the
    cloud image's alpha channel to create a stylized, hatched appearance.

    Args:
        cloud_img: The source cloud image (must be convertible to ARGB32_Premultiplied).
        hatch_cfg: Configuration for the hatch pattern.

    Returns:
        A new QImage with the hatch pattern applied to its alpha channel.
    """
    # Ensure the image is in the correct format for painter operations
    out = cloud_img if cloud_img.format() == QImage.Format_ARGB32_Premultiplied else cloud_img.convertToFormat(QImage.Format_ARGB32_Premultiplied)

    hatch_tile = make_hatch_tile_qimage(
        hatch_cfg.tile_w_px,
        hatch_cfg.tile_h_px,
        hatch_cfg.line_px,
        hatch_cfg.strength,
    )

    # Paint the hatch tile over the image to modify its alpha channel
    p = QPainter(out)
    p.setCompositionMode(QPainter.CompositionMode_DestinationOut)
    p.drawTiledPixmap(out.rect(), QPixmap.fromImage(hatch_tile))
    p.end()

    return out


def _compose_cloud_over_sky_no_global_cache(
    sky_img: QImage,
    cloud_img_rgba: QImage,
    dest_rect: QRect,
    cloud_opacity: float = 1.0,
    gray_mix: float = 1.0,
) -> QImage:
    """
    Composes cloud and sky images without using a global cache.

    This function is intended for use when the draw size changes or when the
    source images are updated. It performs a complex blend:
    1. The sky is partially desaturated to grayscale based on `gray_mix`.
    2. The cloud's color is added on top, controlled by `cloud_opacity`.
    3. The final image is clipped to a circle.

    Args:
        sky_img: The background sky image.
        cloud_img_rgba: The foreground cloud image with an alpha channel.
        dest_rect: The target rectangle for the composition.
        cloud_opacity: The opacity of the cloud color overlay (0.0 to 1.0).
        gray_mix: The factor for desaturating the sky behind the clouds (0.0 to 1.0).

    Returns:
        The composed QImage.
    """
    w, h = dest_rect.width(), dest_rect.height()

    # Ensure images are scaled to the destination size (as a safeguard)
    if sky_img.width() != w or sky_img.height() != h:
        sky_img = sky_img.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)
    if cloud_img_rgba.width() != w or cloud_img_rgba.height() != h:
        cloud_img_rgba = cloud_img_rgba.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

    sky_np = qimage_to_np_rgba(sky_img).astype(np.uint8)
    cloud_np = qimage_to_np_rgba(cloud_img_rgba).astype(np.uint8)

    # --- Grayscale conversion of the sky (integer approximation) ---
    r = sky_np[..., 0].astype(np.uint16)
    g = sky_np[..., 1].astype(np.uint16)
    b = sky_np[..., 2].astype(np.uint16)
    gray_u8 = ((77 * r + 150 * g + 29 * b) >> 8).astype(np.uint8)

    # --- Blend sky with its grayscale version based on cloud alpha ---
    a = (cloud_np[..., 3].astype(np.float32) / 255.0) * float(np.clip(gray_mix, 0.0, 1.0))
    a8 = (np.clip(a, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint16)
    inv_a8 = 255 - a8

    sky_rgb_u16 = sky_np[..., :3].astype(np.uint16)
    gray_rgb_u16 = np.repeat(gray_u8[:, :, None], 3, axis=2).astype(np.uint16)

    # This creates the base layer: sky that becomes grayer where clouds are denser.
    base_u16 = (inv_a8[:, :, None] * sky_rgb_u16 + a8[:, :, None] * gray_rgb_u16) // 255

    # --- Add cloud color on top ---
    cop = float(np.clip(cloud_opacity, 0.0, 1.0))
    if cop > 0.0:
        add_u16 = (cloud_np[..., :3].astype(np.uint16) * int(round(cop * 255))) // 255
        out_u16 = base_u16 + add_u16
        np.minimum(out_u16, 255, out=out_u16)  # Clamp to 255
    else:
        out_u16 = base_u16

    # --- Final assembly and circular mask ---
    out = np.empty((h, w, 4), dtype=np.uint8)
    out[..., :3] = out_u16.astype(np.uint8)
    out[..., 3] = 255  # Start with a fully opaque alpha channel

    # Create a circular mask to make the final disc shape
    cy, cx = (h - 1) * 0.5, (w - 1) * 0.5
    rr = min(cx, cy)
    y, x = np.arange(h, dtype=np.float32)[:, None], np.arange(w, dtype=np.float32)[None, :]
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    mask = r2 <= (rr + 0.25) ** 2  # Anti-aliasing trick to reduce jagged edges

    # Apply the mask: pixels outside the circle become transparent
    out[..., 3][~mask] = 0
    out[..., :3][~mask] = 0  # Premultiply alpha

    return np_rgba_to_qimage(out)


class SkyWindow(DraggableWindow):
    """
    The main application window, displaying the sky view.

    This class orchestrates the entire application. It manages fetching and
    processing of celestial and cloud data, renders all visual elements, and
    handles user input for interaction.
    """

    # Signals for cross-thread communication
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
        self.setWindowTitle(f"Zenith Star View - {self.viewer_data.city_name.title()}")
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
        self.sky_data_calculated.connect(self._on_sky_data_calculated)
        self._is_sky_data_calculation_running: bool = False
        self._sky_data_update_timer = QTimer(self)
        self._sky_data_update_timer.timeout.connect(self.start_background_sky_data_update)
        self.celestial_data: Optional[CelestialData] = None
        self._sky_disc_base_size: int = 1024
        self._sky_disc_image: Optional[QImage] = None

        # --- Cloud Data State and Cache ---
        self._cloud_base_size: int = 256
        self._cloud_img: Optional[QImage] = None
        self._cloud_meta: Optional[dict] = None
        self._is_cloud_update_running: bool = False
        self._cloud_update_pending: bool = False
        self._last_cloud_az: Optional[float] = None
        self._last_cloud_time_utc: Optional[datetime] = None
        self._cloud_update_timer = QTimer(self)
        self._cloud_update_timer.setInterval(CLOUD_UPDATE_INTERVAL * 1000)
        self._cloud_update_timer.timeout.connect(lambda: self.start_background_cloud_update(reason="timer"))

        # --- Cache Cleanup Counter ---
        self._cleanup_counter: int = 0
        self._cleanup_interval: int = 10  # Run cleanup every 10 cloud updates

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

        # --- Composition Cache ---
        self._composited_img: Optional[QImage] = None
        self._composite_key: Optional[Tuple] = None

        # --- Cloud error banner (cleared when next cloud update starts) ---
        self._cloud_banner_text: Optional[str] = None

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
        """
        Handles the window resize event. Repositions UI elements like the
        size grip and menu button.

        Args:
            event: The QResizeEvent.
        """
        grip_size = self.size_grip.size()
        self.size_grip.move(self.width() - grip_size.width(), self.height() - grip_size.height())

        button_size = self.menu_button.size()
        self.menu_button.move(self.width() - button_size.width() - 8, 8)

        # Invalidate the composition cache since the size has changed
        self._composite_key = None
        self._composited_img = None

        super().resizeEvent(event)

    def show_menu(self) -> None:
        """Shows the hamburger menu at the button's position."""
        menu_pos = self.menu_button.mapToGlobal(QPoint(0, self.menu_button.height()))
        self.menu.exec(menu_pos)

    def toggle_enlarge_moon(self) -> None:
        """Toggles the moon enlargement feature and triggers a redraw."""
        self.enlarge_moon = not self.enlarge_moon
        if self._action_enlarge_moon is not None and self._action_enlarge_moon.isChecked() != self.enlarge_moon:
            self._action_enlarge_moon.setChecked(self.enlarge_moon)
        self.update()  # Redraw with the new setting

    def toggle_fullscreen(self) -> None:
        """Toggles the window's fullscreen state."""
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def leaveEvent(self, event: QEvent) -> None:
        """
        Handles the mouse leaving the window area.

        Args:
            event: The QEvent.
        """
        self.mouse_pos = None
        self.update()
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        """
        Handles the mouse move event to track the cursor for highlighting.

        Note: This implementation overrides the base class's dragging
        functionality. Dragging is handled by the native system move or is
        unavailable if that fails.

        Args:
            event: The QMouseEvent.
        """
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
        """
        Sets the celestial data and triggers a repaint.

        Args:
            data: The new CelestialData to display.
        """
        self.celestial_data = data
        self.update()

    def paintEvent(self, event: QPaintEvent) -> None:
        """
        The main paint event handler for the window.

        This method is responsible for drawing all visual elements in the correct
        order: background, sky/clouds, reference lines, celestial objects, and
        overlay text.

        Args:
            event: The QPaintEvent.
        """
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
        self._draw_sky_and_clouds_scaled(painter, geometry)

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
        if self._cloud_banner_text:
            render_draw.draw_status_line_text(
                painter=painter,
                message=self._cloud_banner_text,
                status_line_font=self.status_line_font,
                viewport_rect=self.rect(),
            )

    def _draw_sky_and_clouds_scaled(self, painter: QPainter, geometry: ScreenGeometry) -> None:
        """
        Composes and draws the sky and cloud discs in a single operation.

        This method handles scaling, caching, and compositing of the sky and
        cloud images to optimize rendering performance. The final composited
        image is drawn once.

        Args:
            painter: The QPainter to draw with.
            geometry: The screen geometry for the current view.
        """
        # If both layers are disabled, there's nothing to draw.
        if (self._sky_disc_image is None or self.sky_disc_alpha <= 0.0) and (self._cloud_img is None or self.cloud_disc_alpha <= 0.0):
            return

        # Define the destination rectangle for the disc
        x = int(geometry.center[0] - geometry.radius)
        y = int(geometry.center[1] - geometry.radius)
        w = h = int(geometry.radius * 2)

        # --- Caching Logic ---
        # Create a key to identify the current state. If the key matches the
        # cached key, we can reuse the previously composited image.
        sky_ck = int(self._sky_disc_image.cacheKey()) if self._sky_disc_image else 0
        cloud_ck = int(self._cloud_img.cacheKey()) if self._cloud_img else 0
        comp_key = ("comp", sky_ck, cloud_ck, w, h, float(self.cloud_disc_alpha))

        if self._composite_key != comp_key or self._composited_img is None:
            # --- 1. Scale source images to the required size ---
            def _scaled(qimg: Optional[QImage]) -> Optional[QImage]:
                if qimg is None:
                    return None
                if qimg.width() == w and qimg.height() == h:
                    return qimg
                return qimg.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

            sky_s = _scaled(self._sky_disc_image)
            cloud_s = _scaled(self._cloud_img)

            # --- 2. Apply hatch effect to clouds if enabled ---
            if cloud_s is not None and self.cloud_disc_alpha > 0.0:
                cloud_s = cloud_with_hatched_alpha(cloud_s, CLOUD_HATCH_DEFAULT)

            # --- 3. Composite the layers ---
            if cloud_s is None or self.cloud_disc_alpha <= 0.0:
                # No clouds, just use the sky image.
                composited = sky_s if sky_s is not None else QImage(w, h, QImage.Format_ARGB32_Premultiplied)
                if composited is None or composited.isNull():
                    composited = QImage(w, h, QImage.Format_ARGB32_Premultiplied)
                    composited.fill(Qt.transparent)
            else:
                # Clouds are present, perform the full composition.
                composited = _compose_cloud_over_sky_no_global_cache(
                    sky_img=sky_s if sky_s is not None else QImage(w, h, QImage.Format_ARGB32_Premultiplied),
                    cloud_img_rgba=cloud_s,
                    dest_rect=QRect(0, 0, w, h),
                    cloud_opacity=self.cloud_disc_alpha * 0.7,
                    gray_mix=1.0,
                )

            # Cache the result
            self._composited_img = composited
            self._composite_key = comp_key

        # --- 4. Draw the final composited image ---
        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        painter.drawImage(QRect(x, y, w, h), self._composited_img)
        painter.restore()

    def _on_sky_data_calculated(self, payload: dict) -> None:
        """
        Slot to handle the arrival of newly calculated sky data from the worker thread.

        Args:
            payload: A dictionary containing the 'celestial' data and 'sky_disc' image.
        """
        self.set_sky_data(payload["celestial"])
        self._sky_disc_image = payload["sky_disc"]

        # Invalidate the composition cache
        self._composite_key = None
        self._composited_img = None

        # On the very first load, start the periodic update timers.
        if not self._sky_data_update_timer.isActive():
            self._sky_data_update_timer.start(SKY_UPDATE_INTERVAL * 1000)
            if self._clouddisc and self.cloud_disc_alpha > 0.0 and not self._cloud_update_timer.isActive():
                self._cloud_update_timer.start()
            self.initial_data_loaded.emit()

        self._is_sky_data_calculation_running = False

    def request_sky_data_update(self) -> None:
        """Requests a sky data update if one is not already in progress."""
        if not self.start_background_sky_data_update():
            logger.warning("Sky data update request ignored; an update is already running.")

    def update_sky_data_in_background(self) -> None:
        """
        Performs the calculation of celestial data in a background thread.

        This function calculates the positions of stars, planets, and other
        celestial phenomena. It also generates the sky color disc. The results
        are emitted via a signal to the main thread.
        """
        try:
            now = datetime.now(timezone.utc) + self.delta_t
            time_obj = astropy.time.Time(now)
            lat, lon = self.viewer_data.location

            # Calculate all celestial object positions based on the current time and view
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

            # Generate the sky color disc based on the Sun's position
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

            # Emit the results to the main thread
            payload = {"celestial": celestial_data, "sky_disc": sky_disc_img}
            self.sky_data_calculated.emit(payload)
        except Exception as e:
            logger.error(f"Error in background sky update thread: {e}", exc_info=True)

    def start_background_sky_data_update(self, is_initial_load: bool = False) -> bool:
        """
        Starts a background thread to update the sky data.

        Args:
            is_initial_load: True if this is the first data load.

        Returns:
            True if the update was started, False if an update was already running.
        """
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

    def start_background_cloud_update(self, reason: str = "manual") -> None:
        """
        Kicks off a background cloud image generation if conditions allow.

        Args:
            reason: A string indicating why the update was triggered (e.g., "initial", "timer").
        """
        if not (self._clouddisc and self.cloud_disc_alpha > 0.0):
            return

        # Run cache cleanup in a background thread periodically.
        if self._cleanup_counter % self._cleanup_interval == 0:
            cleanup_thread = threading.Thread(target=self._cleanup_cache)
            cleanup_thread.daemon = True
            cleanup_thread.start()
        self._cleanup_counter += 1

        if self._is_cloud_update_running:
            # If an update is already running, flag that another one is pending.
            # This prevents dropping requests during rapid view changes.
            self._cloud_update_pending = True
            return

        self._cloud_banner_text = "Clouds: downloading…"  # 例: "雲データ: ダウンロード中…"
        self.update()

        self._is_cloud_update_running = True
        t = threading.Thread(target=self._update_cloud_in_background, args=(reason,))
        t.daemon = True
        t.start()

    def _update_cloud_in_background(self, reason: str) -> None:
        """
        Fetches and renders the cloud disc in a worker thread.

        This method uses the CloudDisc service to get the latest satellite
        imagery and render it for the current view.

        Args:
            reason: The reason for the update.
        """
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
                self._cloud_img = qimg
                self._cloud_meta = meta

                # Invalidate composition cache to force a redraw with new clouds
                self._composite_key = None
                self._composited_img = None

                self._cloud_banner_text = None
            except VisibilityError as e:
                logger.error("Invalid params for cloud-disc image generation: %s", e)
            except DownloadError as e:
                logger.warning("Network/S3 download error: %s", e)
                self._cloud_banner_text = "Clouds: Network/S3 download error"
                self.update()
            except TimeoutError as e:
                logger.warning("Clouds fetch timed out: %s", e)
                self._cloud_banner_text = "Clouds: Clouds fetch timed out"
                self.update()
            except DataNotFoundError as e:
                logger.info("Satellite data not found in search window: %s", e)
                self._cloud_banner_text = "Clouds: Satellite data not found in search window"
            except RenderError as e:
                logger.error("Failed to decode/render satellite data: %s", e)
                self._cloud_banner_text = "Clouds: Failed to decode/render satellite data"
            except CloudDiscError as e:
                logger.error("Unexpected clouddisc error: %s", e)
                self._cloud_banner_text = "Clouds: Unexpected clouddisc error"

            self._last_cloud_az = az
            self._last_cloud_time_utc = datetime.now(timezone.utc)

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
        """
        Rotates the view by a given delta in altitude and azimuth.

        Args:
            d_alt: The change in altitude (degrees).
            d_az: The change in azimuth (degrees).
        """
        alt, az = self.viewer_data.view_center
        new_alt = max(0.0, min(90.0, alt + d_alt))
        new_az = (az + d_az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)

        # Request updates for both sky and clouds since the view has changed.
        self.request_sky_data_update()
        self.start_background_cloud_update(reason="az-change")

    def keyPressEvent(self, event: QKeyEvent) -> None:
        """
        Handles key press events for controlling the application.

        Args:
            event: The QKeyEvent.
        """
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
