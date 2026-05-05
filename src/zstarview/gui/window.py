# -*- coding: utf-8 -*-
"""
The main window of the ZStarView application.

This module defines the `SkyWindow` class, which is the primary user interface
for the application. It handles rendering the celestial objects, sky background,
clouds, and all user interactions like rotation, zooming, and object highlighting.
"""

import logging
import os
import sys
import time
from dataclasses import replace
from datetime import datetime, timezone
from datetime import timedelta
from pathlib import Path
from typing import Callable, Optional, Tuple, Union

import astropy.time
from PySide6.QtCore import QEvent, QPoint, QPointF, QRect, QSize, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import (
    QAction,
    QCloseEvent,
    QColor,
    QDesktopServices,
    QFont,
    QFontDatabase,
    QGuiApplication,
    QIcon,
    QKeyEvent,
    QKeySequence,
    QMouseEvent,
    QPaintEvent,
    QPainter,
    QPen,
    QResizeEvent,
)
from PySide6.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QMenu,
    QPushButton,
)
from PySide6.QtWidgets import QWidget

from ..__about__ import __version__
from ..astro import (
    calculate_visible_stars,
    radec_to_altaz as _radec_to_altaz,
)
from ..clouddisc import (
    CloudDisc,
    CloudDiscConfig,
)
from ..clouddisc.providers.select import pick_satellite
from ..cli.args import SKY_OPACITY_DEFAULT
from ..overlay_time import overlay_availability_for_delta, target_time_utc_from_delta
from ..satellites import (
    fetch_horizons_lookup,
    find_satellite_altaz,
    load_satellite_cache,
    satellite_cache_scope_key,
)
from ..satellite_constants import (
    SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS,
    SATELLITE_FAILURE_RETRY_SECONDS,
    SATELLITE_HORIZONS_CACHE_KEY,
    SATELLITE_ISS_CACHE_KEY,
    SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS,
)
from ..aircraft_constants import AIRCRAFT_REFRESH_INTERVAL_SECONDS
from ..aircraft_constants import AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS
from ..paths import (
    APP_ICON_FILE,
    CACHE_PATH,
    CLOUD_UPDATE_INTERVAL,
    GUI_BUTTON_SIZE,
    TEXT_FONT_PATH,
    TEXT_FONT_SIZE,
    STATUS_LINE_FONT_SIZE,
    OVERTURE_DERIVED_ROOT_DIR,
    WINDOW_HEIGHT,
    WINDOW_WIDTH,
    CLOUD_MISSING_TINT_RGBA,
    OBSERVER_MIN_ALT_DEG,
    THEME_STYLES_BY_PRESET,
)
from ..config import load_last_window_geometry, save_last_window_geometry
from ..location_resolver import (
    project_place_target_to_altaz as _project_place_target_to_altaz,
    search_place_candidates,
)
from ..render import geometry as render_geometry
from ..render import stars as render_stars
from ..render.pipeline import (
    compute_star_render_surface_size,
    compute_star_render_upscale_factor,
)
from ..types import ViewerData
from .draggable_window import DraggableWindow
from .composite import SkyCompositorCache
from .cloud_state import CloudImageState
from .cloud_controller import CloudController
from .satellite_state import SatelliteState
from .satellite_controller import SatelliteController
from .aircraft_state import AircraftState
from .aircraft_controller import AircraftController
from .jpl_small_body_controller import JplSmallBodyController
from .terrain_state import TerrainHorizonState
from .terrain_controller import TerrainHorizonController
from .famous_star_dialog import NamedStarJumpDialog
from .famous_star_search_dialog import NamedStarSearchDialog
from .place_search_dialog import PlaceSearchDialog
from .famous_star_shortcuts import (
    NamedStarShortcut,
    build_place_search_jump_targets,
)
from ..search.jpl import search_jpl_targets
from ..search.jpl import project_jpl_target_altaz_from_state_vector
from ..search.jpl import resolve_jpl_target_state_vector
from ..search.satellites import search_satellite_targets
from ..search.models import SearchJumpTarget
from ..asterisms import ASTERISM_KEYS_BY_SOURCE_ID
from .sky_worker import SkyDataWorker
from .window_inputs import PreparedWindowCatalogs
from .window_inputs import SkyWindowRuntimeOptions, SkyWindowUserOptions
from .window_render import SkyWindowRenderMixin
from .window_state import SkyWindowState
from .window_updates import SkyWindowUpdatesMixin
from .urban_outline_controller import UrbanOutlineController
from .urban_outline_state import UrbanOutlineState

logger = logging.getLogger(__name__)


def _resize_event_size(event: QResizeEvent, attr: str) -> tuple[int, int]:
    value = getattr(event, attr, None)
    if callable(value):
        try:
            size = value()
            return int(size.width()), int(size.height())
        except Exception:
            pass
    return (-1, -1)


def _widget_rect_tuple(widget: QWidget) -> tuple[int, int, int, int]:
    geometry = getattr(widget, "geometry", None)
    if callable(geometry):
        try:
            rect = geometry()
            return rect.getRect()
        except Exception:
            pass
    return (-1, -1, -1, -1)


def _widget_visible(widget: QWidget) -> bool:
    visible = getattr(widget, "isVisible", None)
    if callable(visible):
        try:
            return bool(visible())
        except Exception:
            pass
    return False


def _widget_parent_name(widget: QWidget) -> str:
    parent_widget = getattr(widget, "parentWidget", None)
    if callable(parent_widget):
        try:
            parent = parent_widget()
            if parent is not None:
                return type(parent).__name__
        except Exception:
            pass
    return "None"


def _chrome_debug_print(message: str) -> None:
    if _chrome_debug_enabled():
        print(message, file=sys.stderr, flush=True)


def _chrome_event_name(event_type: QEvent.Type) -> str:
    name = getattr(event_type, "name", None)
    return str(name) if name is not None else str(int(event_type))


def _draw_resize_grip_marker(painter: QPainter, rect: QRect) -> None:
    if rect.width() <= 2 or rect.height() <= 2:
        return
    color = QColor(220, 220, 220, 225)
    pen = QPen(color, 2.0)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setPen(pen)
    right = float(rect.right()) - 6.0
    bottom = float(rect.bottom()) - 6.0
    painter.drawLine(
        QPointF(right - 16.0, bottom),
        QPointF(right, bottom - 16.0),
    )
    painter.restore()

_JPL_DEBUG_ENV = "ZSTARVIEW_DEBUG_JPL_SEARCH"
_CHROME_DEBUG_ENV = "ZSTARVIEW_DEBUG_WINDOW_CHROME"


def _jpl_debug_enabled() -> bool:
    raw = os.getenv(_JPL_DEBUG_ENV, "").strip().casefold()
    return raw in {"1", "true", "yes", "on"}


def _jpl_debug_print(message: str) -> None:
    if not _jpl_debug_enabled():
        return
    print(f"[jpl-debug] {message}", file=sys.stderr, flush=True)


def _chrome_debug_enabled() -> bool:
    raw = os.getenv(_CHROME_DEBUG_ENV, "").strip().casefold()
    return raw in {"1", "true", "yes", "on"}

GITHUB_CODE_DATA_LICENSES_AND_CREDITS_URL = (
    "https://github.com/tos-kamiya/zstarview#code-data-licenses-and-credits"
)


def open_code_data_licenses_and_credits() -> None:
    QDesktopServices.openUrl(QUrl(GITHUB_CODE_DATA_LICENSES_AND_CREDITS_URL))


def radec_to_altaz(*args, **kwargs):
    return _radec_to_altaz(*args, **kwargs)


def project_place_target_to_altaz(*args, **kwargs):
    return _project_place_target_to_altaz(*args, **kwargs)


WindowGeometryArg = Union[str, Tuple[int, int, int, int]]
DEFAULT_CLOUD_ALT_MIN_DEG = 1.0


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
        available_rect = (
            primary.availableGeometry()
            if primary is not None
            else screens[0].availableGeometry()
        )

    width = min(width, max(1, available_rect.width()))
    height = min(height, max(1, available_rect.height()))

    max_x = available_rect.right() - width + 1
    max_y = available_rect.bottom() - height + 1
    x = min(max(int(x), available_rect.left()), max_x)
    y = min(max(int(y), available_rect.top()), max_y)
    return x, y, width, height


class SkyWindowClientWidget(SkyWindowRenderMixin, QWidget):
    """Client-area widget shared by frameless and decorated host windows."""

    def __init__(self, owner: "SkyWindowCoreMixin") -> None:
        super().__init__(owner)
        self._owner = owner
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    def __getattr__(self, name: str) -> object:
        return getattr(self._owner, name)

    def resizeEvent(self, event: QResizeEvent) -> None:
        self._owner._handle_client_resize(event)
        super().resizeEvent(event)

    def leaveEvent(self, event: QEvent) -> None:
        self._owner._handle_client_leave(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        self._owner._handle_client_mouse_move(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        self._owner._handle_client_key_press(event)

    def keyReleaseEvent(self, event: QKeyEvent) -> None:
        self._owner._handle_client_key_release(event)


class ShutdownMessageOverlay(QLabel):
    """Centered shutdown message shown while background workers are stopping."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__("Shutting down (closing sub-processes)... please wait.", parent)
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setWordWrap(True)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self.setStyleSheet(
            "QLabel {"
            " background-color: rgba(0, 0, 0, 150);"
            " color: rgba(255, 255, 255, 225);"
            " border: 3px solid rgba(255, 255, 255, 55);"
            " padding: 10px 18px;"
            " font-size: 16px;"
            "}"
        )


class ResizeGripWidget(QWidget):
    """Bottom-right resize affordance with explicit paint and drag handling."""

    def __init__(self, owner: "SkyWindowCoreMixin", parent: QWidget) -> None:
        super().__init__(parent)
        self._owner = owner
        self._drag_active = False
        self._drag_start_global_pos: QPoint | None = None
        self._drag_start_size: tuple[int, int] | None = None
        self.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.setCursor(Qt.CursorShape.SizeFDiagCursor)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, False)
        self.setAutoFillBackground(False)

    def _start_system_resize(self, event: QMouseEvent) -> bool:
        window_handle = self.windowHandle()
        if window_handle is None:
            return False
        start_system_resize = getattr(window_handle, "startSystemResize", None)
        if not callable(start_system_resize):
            return False
        try:
            return bool(
                start_system_resize(
                    Qt.Edge.BottomEdge | Qt.Edge.RightEdge,
                )
            )
        except Exception:
            return False

    def _begin_manual_resize(self, event: QMouseEvent) -> None:
        self._drag_active = True
        self._drag_start_global_pos = event.globalPosition().toPoint()
        window = self.window()
        self._drag_start_size = (int(window.width()), int(window.height()))
        event.accept()

    def _update_manual_resize(self, event: QMouseEvent) -> bool:
        if not self._drag_active:
            return False
        if not (event.buttons() & Qt.MouseButton.LeftButton):
            return False
        if self._drag_start_global_pos is None or self._drag_start_size is None:
            return False
        window = self.window()
        delta = event.globalPosition().toPoint() - self._drag_start_global_pos
        min_size = window.minimumSize()
        new_width = max(int(min_size.width()), self._drag_start_size[0] + int(delta.x()))
        new_height = max(int(min_size.height()), self._drag_start_size[1] + int(delta.y()))
        window.resize(new_width, new_height)
        event.accept()
        return True

    def paintEvent(self, event: QPaintEvent) -> None:
        super().paintEvent(event)
        painter = QPainter(self)
        try:
            _draw_resize_grip_marker(painter, self.rect())
        finally:
            painter.end()

    def mousePressEvent(self, event: QMouseEvent) -> None:
        if event.button() != Qt.MouseButton.LeftButton:
            super().mousePressEvent(event)
            return
        if self._start_system_resize(event):
            event.accept()
            return
        self._begin_manual_resize(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        if self._update_manual_resize(event):
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton and self._drag_active:
            self._drag_active = False
            self._drag_start_global_pos = None
            self._drag_start_size = None
            event.accept()
            return
        super().mouseReleaseEvent(event)


class FramelessWindowFrame(QWidget):
    """Frameless-only window chrome that hosts the client widget and overlay controls."""

    def __init__(
        self, owner: "SkyWindowCoreMixin", client_widget: SkyWindowClientWidget
    ) -> None:
        super().__init__(owner)
        self._owner = owner
        self._client_widget = client_widget
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        self._client_widget.setParent(self)

        self.menu_button = QPushButton("", self)
        self.menu_button.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.menu_button.setStyleSheet(self._owner._menu_button_style_sheet(self._owner.theme))
        self.menu_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.menu_button.clicked.connect(self._owner.show_menu)
        self.menu_button.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        self.size_grip = ResizeGripWidget(self._owner, self)
        self.size_grip.setStyleSheet(self._owner._size_grip_style_sheet())
        self.size_grip.installEventFilter(self)
        self._client_widget.installEventFilter(self)

        self._layout_chrome()

    def overlay_widgets(self) -> list[QWidget]:
        return [self.menu_button, self.size_grip]

    def _log_chrome_state(self, prefix: str) -> None:
        if not _chrome_debug_enabled():
            return
        menu_geo = self.menu_button.geometry()
        grip_geo = self.size_grip.geometry()
        _chrome_debug_print(
            "%s: frame=%sx%s client=%sx%s menu=%s visible=%s grip=%s visible=%s grip_parent=%s"
            % (
                prefix,
                self.width(),
                self.height(),
                self._client_widget.width(),
                self._client_widget.height(),
                menu_geo.getRect(),
                self.menu_button.isVisible(),
                grip_geo.getRect(),
                self.size_grip.isVisible(),
                type(self.size_grip.parentWidget()).__name__,
            )
        )

    def eventFilter(self, watched: object, event: QEvent) -> bool:
        if _chrome_debug_enabled() and isinstance(watched, QWidget):
            if watched in {self.size_grip, self._client_widget}:
                event_type = event.type()
                if event_type in {
                    QEvent.Type.Resize,
                    QEvent.Type.Move,
                    QEvent.Type.Show,
                    QEvent.Type.Hide,
                    QEvent.Type.Paint,
                    QEvent.Type.MouseButtonPress,
                    QEvent.Type.MouseButtonRelease,
                    QEvent.Type.MouseMove,
                }:
                    _chrome_debug_print(
                        "chrome-event: watched=%s type=%s geom=%s visible=%s"
                        % (
                            "size_grip"
                            if watched is self.size_grip
                            else "client_widget",
                            _chrome_event_name(event_type),
                            _widget_rect_tuple(watched),
                            _widget_visible(watched),
                        )
                    )
        return super().eventFilter(watched, event)

    def _layout_chrome(self) -> None:
        self._client_widget.setGeometry(self.rect())
        button_size = self.menu_button.size()
        self.menu_button.move(self.width() - button_size.width(), 0)
        grip_size = self.size_grip.size()
        self.size_grip.move(
            self.width() - grip_size.width(), self.height() - grip_size.height()
        )
        self.menu_button.raise_()
        self.size_grip.raise_()
        self._log_chrome_state("layout-chrome")

    def resizeEvent(self, event: QResizeEvent) -> None:
        if _chrome_debug_enabled():
            old_w, old_h = _resize_event_size(event, "oldSize")
            new_w, new_h = _resize_event_size(event, "size")
            _chrome_debug_print(
                "frameless-frame-resize: old=%sx%s new=%sx%s"
                % (old_w, old_h, new_w, new_h)
            )
        super().resizeEvent(event)
        self._layout_chrome()


class SkyWindowCoreMixin(SkyWindowRenderMixin, SkyWindowUpdatesMixin):
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
    compute_star_render_upscale_factor = staticmethod(
        compute_star_render_upscale_factor
    )

    FRAMELESS_WINDOW = False

    def paintEvent(self, event: QPaintEvent) -> None:
        """Keep the host window from clearing the client area during repaints."""
        event.accept()

    def eventFilter(self, watched: object, event: QEvent) -> bool:
        if isinstance(watched, QApplication):
            if event.type() in (QEvent.Type.KeyPress, QEvent.Type.KeyRelease):
                key_event = event if isinstance(event, QKeyEvent) else None
                if key_event is not None and key_event.key() in {
                    Qt.Key.Key_Left,
                    Qt.Key.Key_Right,
                    Qt.Key.Key_Up,
                    Qt.Key.Key_Down,
                }:
                    if QApplication.activePopupWidget() is None and self.isActiveWindow():
                        if event.type() == QEvent.Type.KeyPress:
                            self._handle_client_key_press(key_event)
                        else:
                            self._handle_client_key_release(key_event)
                        return True
        return super().eventFilter(watched, event)

    def showEvent(self, event) -> None:
        super().showEvent(event)
        if self._frameless_window or self._client_geometry_sync_done:
            return
        self._client_geometry_sync_done = True
        target_client_width, target_client_height = self._target_client_size
        self._resize_client_area(target_client_width, target_client_height)

    def __init__(
        self,
        viewer_data: ViewerData,
        catalogs: PreparedWindowCatalogs,
        user_options: SkyWindowUserOptions,
        runtime_options: SkyWindowRuntimeOptions,
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
        self.star_catalog_np = catalogs.star_catalog_np
        self.star_catalog_lod6_indices = catalogs.star_catalog_lod6_indices
        self.star_catalog_meta = catalogs.star_catalog_meta
        self.dso_catalog_np = catalogs.dso_catalog_np
        self.show_dso: bool = self.dso_catalog_np is not None
        if user_options.show_dso_initial is not None:
            self.show_dso = (
                bool(user_options.show_dso_initial) and self.dso_catalog_np is not None
            )
        self.show_asterisms: bool = (
            True
            if user_options.show_asterisms_initial is None
            else bool(user_options.show_asterisms_initial)
        )
        self.show_guidelines: bool = (
            True
            if user_options.show_guidelines_initial is None
            else bool(user_options.show_guidelines_initial)
        )
        # Determine observation info startup mode. The user option may be a mode string
        # ('auto','top','bottom','off') or None.
        raw_observation_mode = user_options.observation_info_mode
        if isinstance(raw_observation_mode, str):
            observation_mode = raw_observation_mode
        elif raw_observation_mode is None:
            observation_mode = "auto"
        else:
            observation_mode = "auto" if bool(raw_observation_mode) else "off"
        self.observation_info_mode = observation_mode
        # Pinned when explicitly top or bottom
        self.observation_info_pinned = observation_mode in ("top", "bottom")
        # Visible unless 'off'
        self.show_observation_info: bool = observation_mode != "off"
        # Initialize overlay_info_bottom_left state: True means bottom-left, False means top-left
        # When pinned, set the HUD state immediately so rendering respects the pin.
        if observation_mode == "bottom":
            try:
                self.state.overlay_info_bottom_left = True
            except Exception:
                # state may not yet be fully initialized; ignore and rely on the render path
                pass
        elif observation_mode == "top":
            try:
                self.state.overlay_info_bottom_left = False
            except Exception:
                pass

        self._named_stars_by_band = catalogs.named_stars_by_band
        self._named_stars_search_all = catalogs.named_stars_search_all
        self.delta_t = runtime_options.delta_t
        overlay_availability = overlay_availability_for_delta(self.delta_t)
        self.sky_disc_alpha = user_options.sky_disc_alpha
        self._sky_disc_alpha_when_enabled = (
            user_options.sky_disc_alpha
            if user_options.sky_disc_alpha > 0.0
            else SKY_OPACITY_DEFAULT
        )
        requested_satellite_opacity = user_options.satellite_opacity
        self._satellite_toggle_supported = overlay_availability.satellite
        self._satellite_opacity_when_enabled = (
            requested_satellite_opacity if requested_satellite_opacity > 0.0 else 1.0
        )
        self.satellite_opacity = (
            requested_satellite_opacity if self._satellite_toggle_supported else 0.0
        )
        requested_aircraft_opacity = user_options.aircraft_opacity
        self._aircraft_toggle_supported = overlay_availability.aircraft
        self._aircraft_opacity_when_enabled = (
            requested_aircraft_opacity if requested_aircraft_opacity > 0.0 else 1.0
        )
        self.aircraft_opacity = (
            requested_aircraft_opacity if self._aircraft_toggle_supported else 0.0
        )
        self.terrain_horizon_opacity = user_options.terrain_horizon_opacity
        self.earth_guide_opacity = user_options.earth_guide_opacity
        self.urban_outline_opacity = user_options.urban_outline_opacity
        self.ground_tint_opacity = user_options.ground_tint_opacity
        self._terrain_horizon_opacity_when_enabled = (
            user_options.terrain_horizon_opacity
            if user_options.terrain_horizon_opacity > 0.0
            else 0.25
        )
        self._earth_guide_opacity_when_enabled = (
            user_options.earth_guide_opacity
            if user_options.earth_guide_opacity > 0.0
            else 0.028
        )
        self._urban_outline_opacity_when_enabled = (
            user_options.urban_outline_opacity
            if user_options.urban_outline_opacity > 0.0
            else 0.2
        )
        self._sky_disc_gui_allowed = bool(user_options.sky_disc_gui_allowed)
        self._cloud_gui_allowed = bool(user_options.cloud_gui_allowed)
        self._satellite_gui_allowed = bool(user_options.satellite_gui_allowed)
        self._aircraft_gui_allowed = bool(user_options.aircraft_gui_allowed)
        self._terrain_horizon_gui_allowed = bool(
            user_options.terrain_horizon_gui_allowed
        )
        self._earth_guide_gui_allowed = bool(user_options.earth_guide_gui_allowed)
        self._urban_outline_gui_allowed = bool(user_options.urban_outline_gui_allowed)
        self.show_urban_outline_layer: bool = self.urban_outline_opacity > 0.0
        self.enlarge_moon = user_options.enlarge_moon
        self.bright_bodies_mode = user_options.bright_bodies_mode
        self.star_base_radius = user_options.star_base_radius
        self.vmag_limit = user_options.vmag_limit
        self.sky_update_interval = runtime_options.sky_update_interval
        self.urban_outline_radius_km = float(runtime_options.urban_outline_radius_km)
        self.urban_outline_skyscraper_radius_km = float(
            runtime_options.urban_outline_skyscraper_radius_km
        )
        self.urban_outline_min_height_m = float(
            runtime_options.urban_outline_min_height_m
        )
        self.urban_outline_feature_type = str(
            runtime_options.urban_outline_feature_type
        )
        self.urban_outline_skyscraper_only = bool(
            runtime_options.urban_outline_skyscraper_only
        )
        self.visual_preset = user_options.visual_preset
        self.theme = THEME_STYLES_BY_PRESET.get(
            self.visual_preset,
            THEME_STYLES_BY_PRESET["night"],
        )
        self.star_visibility_boost = user_options.star_visibility_boost
        self.asterism_visibility_boost = user_options.asterism_visibility_boost
        self.earth_guide_visibility_boost = user_options.earth_guide_visibility_boost
        self._star_render_expected_width = runtime_options.star_render_expected_width
        self.content_fov_deg = float(runtime_options.content_fov_deg)
        self._cloud_toggle_supported = overlay_availability.cloud

        # Cloud opacity is disabled if we are looking at a time-shifted view,
        # as we can only fetch current cloud data.
        requested_cloud_alpha = user_options.cloud_disc_alpha
        self._cloud_alpha_when_enabled = (
            requested_cloud_alpha if requested_cloud_alpha > 0.0 else 0.2
        )
        self.cloud_disc_alpha: float = requested_cloud_alpha
        if not self._cloud_toggle_supported:
            self.cloud_disc_alpha = 0.0

        # --- Viewer and Window Setup ---
        self.viewer_data = viewer_data
        self._search_view_center_base = tuple(self.viewer_data.view_center)
        self._search_view_center_alt_specified = False
        self._search_view_center_az_specified = False
        self.state = SkyWindowState(
            render_view_center=tuple(self.viewer_data.view_center),
            urban_outlines=None,
            satellite_overlay_points=None,
            aircraft_overlay_points=None,
        )
        self._viewport_rotation_keys_down: set[int] = set()
        # Ensure overlay_info_bottom_left reflects the startup mode now that
        # the mutable state object exists. True==bottom-left, False==top-left.
        try:
            if self.observation_info_mode == "bottom":
                self.state.overlay_info_bottom_left = True
            elif self.observation_info_mode == "top":
                self.state.overlay_info_bottom_left = False
        except Exception:
            pass

        self._enabled_satellite_groups: tuple[str, ...] = (
            SATELLITE_ISS_CACHE_KEY,
            SATELLITE_HORIZONS_CACHE_KEY,
        )
        self._disc_generation = 0
        self._frame_cache_key: object | None = None
        self._frame_cache_image = None
        self._present_frame_cache_key: object | None = None
        self._present_frame_cache_image = None
        self.setWindowTitle(self.viewer_data.city_name)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setFocus(Qt.FocusReason.ActiveWindowFocusReason)
        self.setWindowIcon(QIcon(APP_ICON_FILE))
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        initial_x = 100
        initial_y = 100
        initial_width = WINDOW_WIDTH
        initial_height = WINDOW_HEIGHT
        min_width = 400
        min_height = 400
        self.setMinimumSize(min_width, min_height)
        self._window_frame_mode = runtime_options.window_frame_mode
        self._frameless_window = bool(self.FRAMELESS_WINDOW)
        requested_geometry: Optional[Tuple[int, int, int, int]] = None
        if runtime_options.window_geometry_arg == "restore":
            requested_geometry = load_last_window_geometry()
        elif isinstance(runtime_options.window_geometry_arg, tuple):
            requested_geometry = runtime_options.window_geometry_arg
        if requested_geometry is not None:
            initial_x, initial_y, initial_width, initial_height = requested_geometry
        initial_x, initial_y, initial_width, initial_height = (
            _clamp_window_geometry_to_screen(
                initial_x,
                initial_y,
                initial_width,
                initial_height,
                min_width=min_width,
                min_height=min_height,
            )
        )
        self._target_client_size = (int(initial_width), int(initial_height))
        self._client_geometry_sync_done = False
        self.setGeometry(initial_x, initial_y, initial_width, initial_height)
        self._client_widget = SkyWindowClientWidget(self)
        self._shutdown_overlay: Optional[ShutdownMessageOverlay] = None
        self._frameless_frame: Optional[FramelessWindowFrame] = None
        self.menu_button: Optional[QPushButton] = None
        self.size_grip: Optional[QWidget] = None
        self._action_enlarge_moon: Optional[QAction] = None
        self._action_toggle_clouds: Optional[QAction] = None
        self._action_toggle_satellites: Optional[QAction] = None
        self._action_toggle_aircraft: Optional[QAction] = None
        self._action_toggle_terrain_horizon: Optional[QAction] = None
        self._action_toggle_earth_guide: Optional[QAction] = None
        self._action_toggle_urban_outline: Optional[QAction] = None
        self._action_toggle_dso: Optional[QAction] = None
        self._action_toggle_asterisms: Optional[QAction] = None
        self._action_toggle_guidelines: Optional[QAction] = None
        self._action_toggle_observation_info: Optional[QAction] = None
        self._action_toggle_sky_disc: Optional[QAction] = None
        self._action_raise_view: Optional[QAction] = None
        self._action_lower_view: Optional[QAction] = None
        self._build_window_menu()
        self._install_window_host()
        app = QApplication.instance()
        if app is not None:
            app.installEventFilter(self)
        self._client_widget.setFocus(Qt.FocusReason.ActiveWindowFocusReason)

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
        self.satellite_state = SatelliteState()
        self.aircraft_state = AircraftState()
        self.terrain_horizon_state = TerrainHorizonState()
        self.urban_outline_state = UrbanOutlineState()
        self._cloud_controller: Optional[CloudController] = None
        self._satellite_controller: Optional[SatelliteController] = None
        self._aircraft_controller: Optional[AircraftController] = None
        self._jpl_small_body_controller: Optional[JplSmallBodyController] = None
        self._terrain_horizon_controller: Optional[TerrainHorizonController] = None
        self._urban_outline_controller: Optional[UrbanOutlineController] = None
        self._cloud_update_timer = QTimer(self)
        self._cloud_update_timer.setInterval(CLOUD_UPDATE_INTERVAL * 1000)
        self._cloud_update_timer.timeout.connect(
            lambda: self.start_background_cloud_update(reason="timer")
        )
        self._satellite_update_timer = QTimer(self)
        self._satellite_update_timer.setSingleShot(True)
        self._satellite_update_timer.setInterval(
            SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS * 1000
        )
        self._satellite_update_timer.timeout.connect(self._on_satellite_refresh_timer)
        self._aircraft_update_timer = QTimer(self)
        self._aircraft_update_timer.setSingleShot(True)
        self._aircraft_update_timer.setInterval(
            AIRCRAFT_REFRESH_INTERVAL_SECONDS * 1000
        )
        self._aircraft_update_timer.timeout.connect(self._on_aircraft_refresh_timer)
        self._overlay_projection_timer = QTimer(self)
        self._overlay_projection_timer.setInterval(
            min(
                SATELLITE_POSITION_REFRESH_INTERVAL_SECONDS,
                AIRCRAFT_PREDICTION_REFRESH_INTERVAL_SECONDS,
            )
            * 1000
        )
        self._overlay_projection_timer.timeout.connect(
            self._on_overlay_projection_timer
        )

        # --- CloudDisc Service Initialization ---
        self._clouddisc: Optional[CloudDisc] = None
        clouddisc_config = CloudDiscConfig(
            cache_dir=CACHE_PATH,
            sat_priority=("AUTO",),
            bt_warm_k=310.0,
            bt_cold_k=190.0,
            alt_min_deg=DEFAULT_CLOUD_ALT_MIN_DEG,
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
        self._satellite_controller = SatelliteController(parent=self)
        self._satellite_controller.satellite_started.connect(self._on_satellite_started)
        self._satellite_controller.satellite_ready.connect(self._on_satellite_ready)
        self._satellite_controller.satellite_failed.connect(self._on_satellite_failed)
        self._aircraft_controller = AircraftController(parent=self)
        self._aircraft_controller.aircraft_started.connect(self._on_aircraft_started)
        self._aircraft_controller.aircraft_ready.connect(self._on_aircraft_ready)
        self._aircraft_controller.aircraft_failed.connect(self._on_aircraft_failed)
        self._jpl_small_body_controller = JplSmallBodyController(parent=self)
        self._jpl_small_body_controller.jpl_started.connect(self._on_jpl_started)
        self._jpl_small_body_controller.jpl_ready.connect(self._on_jpl_ready)
        self._jpl_small_body_controller.jpl_failed.connect(self._on_jpl_failed)
        self._persistent_search_update_timer = QTimer(self)
        self._persistent_search_update_timer.setSingleShot(True)
        self._persistent_search_update_timer.timeout.connect(
            self._on_persistent_search_refresh_timer
        )
        if self._action_toggle_clouds is not None:
            self._action_toggle_clouds.setEnabled(
                self._cloud_toggle_supported
                and self._clouddisc is not None
                and self._cloud_gui_allowed
            )
        if self._action_toggle_satellites is not None:
            self._action_toggle_satellites.setEnabled(
                self._satellite_toggle_supported and self._satellite_gui_allowed
            )
        if self._action_toggle_aircraft is not None:
            self._action_toggle_aircraft.setEnabled(
                self._aircraft_toggle_supported and self._aircraft_gui_allowed
            )

        terrain_cache_dir = Path(CACHE_PATH) / "copernicus-dem"
        self._terrain_horizon_controller = TerrainHorizonController(
            cache_dir=terrain_cache_dir,
            parent=self,
        )
        self._terrain_horizon_controller.terrain_started.connect(
            self._on_terrain_horizon_started
        )
        self._terrain_horizon_controller.terrain_ready.connect(
            self._on_terrain_horizon_ready
        )
        self._terrain_horizon_controller.terrain_failed.connect(
            self._on_terrain_horizon_failed
        )
        self._urban_outline_controller = UrbanOutlineController(
            derived_root_dir=Path(OVERTURE_DERIVED_ROOT_DIR),
            min_building_height_m=self.urban_outline_min_height_m,
            radius_km=self.urban_outline_radius_km,
            skyscraper_outer_radius_km=self.urban_outline_skyscraper_radius_km,
            feature_type=self.urban_outline_feature_type,
            skyscraper_only=self.urban_outline_skyscraper_only,
            parent=self,
        )
        self._urban_outline_controller.urban_started.connect(
            self._on_urban_outline_started
        )
        self._urban_outline_controller.urban_ready.connect(self._on_urban_outline_ready)
        self._urban_outline_controller.urban_failed.connect(
            self._on_urban_outline_failed
        )
        if self._action_toggle_terrain_horizon is not None:
            self._action_toggle_terrain_horizon.setEnabled(
                self._terrain_horizon_gui_allowed
            )
        if self._action_toggle_earth_guide is not None:
            self._action_toggle_earth_guide.setEnabled(self._earth_guide_gui_allowed)
        if self._action_toggle_sky_disc is not None:
            self._action_toggle_sky_disc.setEnabled(self._sky_disc_gui_allowed)
        if self._action_toggle_urban_outline is not None:
            self._action_toggle_urban_outline.setEnabled(
                self._urban_outline_gui_allowed
            )

        # --- Composition Cache (moved to dedicated class) ---
        target_stripes, width_factor = runtime_options.cloud_stripe_style
        missing_tint_alpha = int(
            round(255.0 * runtime_options.cloud_missing_tint_opacity)
        )
        missing_tint_rgba = (
            int(CLOUD_MISSING_TINT_RGBA[0]),
            int(CLOUD_MISSING_TINT_RGBA[1]),
            int(CLOUD_MISSING_TINT_RGBA[2]),
            missing_tint_alpha,
        )
        self._compositor = SkyCompositorCache(
            cloud_target_stripes=int(target_stripes),
            cloud_stripe_width_factor=float(width_factor),
            cloud_stripe_mode=runtime_options.cloud_stripe_mode,
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
        if self.urban_outline_opacity > 0.0:
            self.start_background_urban_outline_update(reason="initial")
        if self._satellite_layer_enabled():
            self._enable_satellite_layer(reason="initial")
        if self._aircraft_layer_enabled():
            self._enable_aircraft_layer(reason="initial")

    def _setup_update_infrastructure(self) -> None:
        """Initialize timers, worker, and signal wiring for background updates."""
        self._sky_worker = SkyDataWorker(self)
        self._sky_worker.data_ready.connect(self._on_sky_data_calculated)
        self.cloud_repaint_requested.connect(self.request_client_update)

        app = QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(self._begin_shutdown)

        self._sky_data_update_timer = QTimer(self)
        self._sky_data_update_timer.timeout.connect(
            lambda: self.start_background_sky_data_update(reason="timer")
        )

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
        self._viewport_interaction_idle_timer.setInterval(
            self.state.viewport_interaction_idle_ms
        )
        self._viewport_interaction_idle_timer.timeout.connect(
            self._end_viewport_interaction_mode
        )

    def _resize_client_area(self, target_client_width: int, target_client_height: int) -> None:
        """Resize the host so the client widget reaches the requested size."""
        if self.isFullScreen() or self.isMaximized():
            self.showNormal()
        self._target_client_size = (int(target_client_width), int(target_client_height))
        current_client_width = int(self.client_width())
        current_client_height = int(self.client_height())
        delta_width = int(target_client_width) - current_client_width
        delta_height = int(target_client_height) - current_client_height
        if delta_width == 0 and delta_height == 0:
            return
        self.resize(self.width() + delta_width, self.height() + delta_height)

    def _fit_client_area_to_screen(self) -> None:
        """Resize the client area to the largest aspect-preserving size on screen."""
        current_client_width = max(1, int(self.client_width()))
        current_client_height = max(1, int(self.client_height()))
        if self.isFullScreen() or self.isMaximized():
            self.showNormal()

        screen = None
        window_handle = self.windowHandle()
        if window_handle is not None:
            screen = window_handle.screen()
        if screen is None:
            screen = self.screen()
        if screen is None:
            screen = QGuiApplication.primaryScreen()
        if screen is None:
            return

        available_geometry = screen.availableGeometry()
        if not available_geometry.isValid():
            return

        frame_width = max(0, int(self.width()) - int(self.client_width()))
        frame_height = max(0, int(self.height()) - int(self.client_height()))
        max_window_width = max(1, int(available_geometry.width()) - frame_width)
        max_window_height = max(1, int(available_geometry.height()) - frame_height)

        target_window_size = (
            QSize(current_client_width, current_client_height)
            .scaled(max_window_width, max_window_height, Qt.AspectRatioMode.KeepAspectRatio)
        )
        target_client_width = max(1, int(target_window_size.width()))
        target_client_height = max(1, int(target_window_size.height()))
        if (
            target_client_width == current_client_width
            and target_client_height == current_client_height
        ):
            return
        self._resize_client_area(target_client_width, target_client_height)

    def square_client_area(self) -> None:
        """Resize the client area so width and height match."""
        side = max(1, int(self.client_width()))
        self._resize_client_area(side, side)

    def _install_window_host(self) -> None:
        """Install the host-specific window chrome around the shared client widget."""
        raise NotImplementedError

    def _build_window_menu(self) -> None:
        """Build window actions and the popup menu shared by chrome implementations."""
        self.menu = QMenu("Menu", self)
        self.file_menu = QMenu("File", self)
        self.search_menu = QMenu("Search", self)
        self.observer_view_menu = QMenu("View Direction", self)
        self.display_menu = QMenu("Layers", self)
        self.help_menu = QMenu("Help", self)
        if not self._frameless_window:
            self.menu.addMenu(self.file_menu)
        self.menu.addMenu(self.search_menu)
        self.menu.addMenu(self.display_menu)
        self.menu.addMenu(self.observer_view_menu)

        self._add_menu_action(
            self.observer_view_menu,
            f"Rotate Left (-{self.state.rotation_step:.0f}°)",
            triggered=lambda: self._rotate_view(
                d_az=-self.state.rotation_step, interactive_viewport=True
            ),
        )
        self._add_menu_action(
            self.observer_view_menu,
            f"Rotate Right (+{self.state.rotation_step:.0f}°)",
            triggered=lambda: self._rotate_view(
                d_az=+self.state.rotation_step, interactive_viewport=True
            ),
        )
        self._action_raise_view = self._add_menu_action(
            self.observer_view_menu,
            f"Raise View (+{self.state.rotation_step:.0f}° alt)",
            triggered=lambda: self._rotate_view(
                d_alt=+self.state.rotation_step, interactive_viewport=True
            ),
        )
        self._action_lower_view = self._add_menu_action(
            self.observer_view_menu,
            f"Lower View (-{self.state.rotation_step:.0f}° alt)",
            triggered=lambda: self._rotate_view(
                d_alt=-self.state.rotation_step, interactive_viewport=True
            ),
        )

        self.observer_view_menu.addSeparator()
        self._add_menu_action(
            self.search_menu,
            "Jump to Named Star...",
            shortcut=QKeySequence("Ctrl+J"),
            triggered=self._open_named_star_jump_dialog,
        )
        self._add_menu_action(
            self.search_menu,
            "Search Objects...",
            shortcut=QKeySequence("Ctrl+F"),
            triggered=self._open_named_star_search_dialog,
        )
        self._add_menu_action(
            self.search_menu,
            "Search Places...",
            triggered=self._open_place_search_dialog,
        )

        self.display_menu.addSeparator()
        self._action_enlarge_moon = self._add_checkable_menu_action(
            self.display_menu,
            "Enlarge Moon",
            checked=self.enlarge_moon,
            shortcut=QKeySequence(Qt.Key.Key_M),
            triggered=self.toggle_enlarge_moon,
        )
        self._action_toggle_dso = self._add_checkable_menu_action(
            self.display_menu,
            "DSO",
            checked=self.show_dso,
            enabled=self.dso_catalog_np is not None,
            shortcut=QKeySequence(Qt.Key.Key_D),
            triggered=self.toggle_dso,
        )
        self._action_toggle_asterisms = self._add_checkable_menu_action(
            self.display_menu,
            "Asterisms",
            checked=self.show_asterisms,
            shortcut=QKeySequence(Qt.Key.Key_A),
            triggered=self.toggle_asterisms,
        )
        self._action_toggle_guidelines = self._add_checkable_menu_action(
            self.display_menu,
            "Guidelines",
            checked=self.show_guidelines,
            shortcut=QKeySequence(Qt.Key.Key_G),
            triggered=self.toggle_guidelines,
        )
        self._action_toggle_observation_info = self._add_checkable_menu_action(
            self.display_menu,
            "Observation Info",
            checked=self.show_observation_info,
            enabled=self.observation_info_mode != "off",
            triggered=self.toggle_observation_info,
        )

        self.display_menu.addSeparator()
        self._action_toggle_sky_disc = self._add_checkable_menu_action(
            self.display_menu,
            "Sky Color Disc",
            checked=self.sky_disc_alpha > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_S),
            triggered=self.toggle_sky_disc,
        )
        self._action_toggle_clouds = self._add_checkable_menu_action(
            self.display_menu,
            "Clouds",
            checked=self.cloud_disc_alpha > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_C),
            triggered=self.toggle_clouds,
        )
        self._action_toggle_satellites = self._add_checkable_menu_action(
            self.display_menu,
            "Satellites",
            checked=self.satellite_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_I),
            triggered=self.toggle_satellites,
        )
        self._action_toggle_aircraft = self._add_checkable_menu_action(
            self.display_menu,
            "Aircraft",
            checked=self.aircraft_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_P),
            triggered=self.toggle_aircraft,
        )
        self._action_toggle_terrain_horizon = self._add_checkable_menu_action(
            self.display_menu,
            "Terrain Horizon",
            checked=self.terrain_horizon_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_T),
            triggered=self.toggle_terrain_horizon,
        )
        self._action_toggle_earth_guide = self._add_checkable_menu_action(
            self.display_menu,
            "Earth Guide",
            checked=self.earth_guide_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_E),
            triggered=self.toggle_earth_guide,
        )
        self._action_toggle_urban_outline = self._add_checkable_menu_action(
            self.display_menu,
            "Urban Outline",
            checked=self.urban_outline_opacity > 0.0,
            shortcut=QKeySequence(Qt.Key.Key_U),
            triggered=self.toggle_urban_outline,
        )

        square_client_area_action = self._add_menu_action(
            self.file_menu,
            "Square Client Area",
            triggered=self.square_client_area,
        )
        fit_to_screen_action = self._add_menu_action(
            self.file_menu,
            "Fit to Screen",
            triggered=self._fit_client_area_to_screen,
        )
        fullscreen_action = self._add_menu_action(
            self.file_menu,
            "Fullscreen",
            shortcut=QKeySequence(Qt.Key.Key_F11),
            triggered=self.toggle_fullscreen,
        )

        self.file_menu.addSeparator()
        self.file_menu.addMenu(self.help_menu)
        self.file_menu.addSeparator()
        request_quit = getattr(self, "_request_application_quit", None)
        exit_action = self._add_menu_action(
            self.file_menu,
            "Exit",
            shortcut=QKeySequence(Qt.Key.Key_Q),
            triggered=request_quit if callable(request_quit) else None,
        )

        self.display_menu.addSeparator()
        vmag_limit_action = self._add_menu_action(
            self.display_menu,
            self._vmag_limit_menu_text(),
        )
        vmag_limit_action.setEnabled(False)
        version_action = self._add_menu_action(
            self.help_menu,
            f"Version {__version__}",
        )
        version_action.setEnabled(False)
        self._add_menu_action(
            self.help_menu,
            "Code, Data Licenses, and Credits ->",
            triggered=open_code_data_licenses_and_credits,
        )

        if self._frameless_window:
            self.menu.addSeparator()
            self.menu.addAction(square_client_area_action)
            self.menu.addAction(fit_to_screen_action)
            self.menu.addAction(fullscreen_action)
            self.menu.addSeparator()
            self.menu.addMenu(self.help_menu)
            self.menu.addSeparator()
            self.menu.addAction(exit_action)

    def _create_menu_action(self, menu: QMenu, text: str) -> QAction:
        """Create a menu action and add it to the target menu."""
        action = QAction(text, self)
        menu.addAction(action)
        return action

    def _add_menu_action(
        self,
        menu: QMenu,
        text: str,
        *,
        shortcut: QKeySequence | str | int | None = None,
        enabled: bool = True,
        triggered: Callable[..., object] | None = None,
    ) -> QAction:
        """Create a non-checkable menu action."""
        action = self._create_menu_action(menu, text)
        if shortcut is not None:
            action.setShortcut(QKeySequence(shortcut))
            action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
            self.addAction(action)
        action.setEnabled(enabled)
        if triggered is not None:
            action.triggered.connect(triggered)
        return action

    def _add_checkable_menu_action(
        self,
        menu: QMenu,
        text: str,
        *,
        checked: bool = False,
        shortcut: QKeySequence | str | int | None = None,
        enabled: bool = True,
        triggered: Callable[..., object] | None = None,
    ) -> QAction:
        """Create a checkable menu action."""
        action = self._create_menu_action(menu, text)
        action.setCheckable(True)
        action.setChecked(checked)
        if shortcut is not None:
            action.setShortcut(QKeySequence(shortcut))
            action.setShortcutContext(Qt.ShortcutContext.WindowShortcut)
            self.addAction(action)
        action.setEnabled(enabled)
        if triggered is not None:
            action.triggered.connect(triggered)
        return action

    def _ensure_shutdown_overlay(self) -> ShutdownMessageOverlay:
        if self._shutdown_overlay is None:
            self._shutdown_overlay = ShutdownMessageOverlay(self)
        self._shutdown_overlay.setGeometry(self._shutdown_message_geometry())
        return self._shutdown_overlay

    def _shutdown_message_geometry(self) -> QRect:
        overlay = self._shutdown_overlay
        client_width = max(1, int(self.client_width()))
        width = max(1, int(round(client_width * 0.95)))
        width = min(width, client_width)
        if overlay is None:
            height = 64
        else:
            height = max(1, int(overlay.sizeHint().height()))
        x = max(0, (client_width - width) // 2)
        y = max(0, (int(self.height()) - height) // 2)
        return QRect(x, y, width, height)

    def _show_shutdown_message(self) -> None:
        overlay = self._ensure_shutdown_overlay()
        overlay.raise_()
        overlay.show()
        overlay.setGeometry(self._shutdown_message_geometry())
        QApplication.processEvents()

    def _hide_shutdown_message(self) -> None:
        if self._shutdown_overlay is not None:
            self._shutdown_overlay.hide()

    def _request_application_quit(self) -> None:
        self._show_shutdown_message()
        QApplication.quit()

    def _add_help_menu_actions(self, menu: QMenu) -> None:
        version_action = self._add_menu_action(
            menu,
            f"Version {__version__}",
        )
        version_action.setEnabled(False)
        self._add_menu_action(
            menu,
            "Code, Data Licenses, and Credits ->",
            triggered=open_code_data_licenses_and_credits,
        )

    def _attach_help_menu_to_file_menu(self) -> None:
        self.file_menu.addSeparator()
        self.file_menu.addMenu(self.help_menu)
        self.file_menu.addSeparator()

    def _attach_help_menu_to_frameless_menu(
        self,
        square_client_area_action: QAction,
        fullscreen_action: QAction,
        exit_action: QAction,
    ) -> None:
        self.menu.addSeparator()
        self.menu.addAction(square_client_area_action)
        self.menu.addAction(fullscreen_action)
        self.menu.addSeparator()
        self.menu.addMenu(self.help_menu)
        self.menu.addSeparator()
        self.menu.addAction(exit_action)

    def _attach_client_menu_button(self, parent: QWidget) -> None:
        """Attach the legacy popup-menu button directly on the client area."""
        self.menu_button = QPushButton("", parent)
        self.menu_button.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
        self.menu_button.setStyleSheet(self._menu_button_style_sheet(self.theme))
        self.menu_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.menu_button.clicked.connect(self.show_menu)
        self.menu_button.raise_()

    def _vmag_limit_menu_text(self) -> str:
        return f"Vmag limit {self.vmag_limit:.1f}"

    def _menu_button_style_sheet(self, theme) -> str:
        chrome = theme.window_chrome
        text = "#%02x%02x%02x" % chrome.menu_button_text_rgb
        return (
            "QPushButton {"
            " border: none;"
            " font-size: 16px;"
            " background: transparent;"
            f" color: {text};"
            "}"
            "QPushButton:hover {"
            f" color: rgb({chrome.menu_hover_text_rgb[0]}, {chrome.menu_hover_text_rgb[1]}, {chrome.menu_hover_text_rgb[2]});"
            f" background-color: rgba({chrome.menu_hover_bg_rgba[0]}, {chrome.menu_hover_bg_rgba[1]}, {chrome.menu_hover_bg_rgba[2]}, {chrome.menu_hover_bg_rgba[3]});"
            "}"
            "QPushButton:pressed {"
            f" color: rgb({chrome.menu_pressed_text_rgb[0]}, {chrome.menu_pressed_text_rgb[1]}, {chrome.menu_pressed_text_rgb[2]});"
            f" background-color: rgba({chrome.menu_pressed_bg_rgba[0]}, {chrome.menu_pressed_bg_rgba[1]}, {chrome.menu_pressed_bg_rgba[2]}, {chrome.menu_pressed_bg_rgba[3]});"
            "}"
            "QPushButton:menu-indicator { image: none; }"
        )

    def _size_grip_style_sheet(self) -> str:
        return "QWidget { border: none; background: transparent;}"

    def _raise_overlay_widgets(self) -> None:
        if self.menu_button is not None:
            self.menu_button.raise_()
        if self.size_grip is not None:
            self.size_grip.raise_()

    def _sync_view_altitude_actions(self) -> None:
        alt, _ = self.viewer_data.view_center
        if self._action_raise_view is not None:
            self._action_raise_view.setEnabled(float(alt) < 90.0)
        if self._action_lower_view is not None:
            self._action_lower_view.setEnabled(float(alt) > OBSERVER_MIN_ALT_DEG)

    def client_width(self) -> int:
        if self._client_widget is not None:
            return self._client_widget.width()
        return (
            self.centralWidget().width()
            if self.centralWidget() is not None
            else super().width()
        )

    def client_height(self) -> int:
        if self._client_widget is not None:
            return self._client_widget.height()
        return (
            self.centralWidget().height()
            if self.centralWidget() is not None
            else super().height()
        )

    def client_size(self):
        if self._client_widget is not None:
            return self._client_widget.size()
        return (
            self.centralWidget().size()
            if self.centralWidget() is not None
            else super().size()
        )

    def client_rect(self):
        if self._client_widget is not None:
            return self._client_widget.rect()
        return (
            self.centralWidget().rect()
            if self.centralWidget() is not None
            else super().rect()
        )

    def request_client_update(self) -> None:
        if self._client_widget is not None:
            self._client_widget.update()
            return
        central = self.centralWidget()
        if central is not None:
            central.update()
            return
        super().update()

    def _handle_client_resize(self, event: QResizeEvent) -> None:
        old_w, old_h = _resize_event_size(event, "oldSize")
        new_w, new_h = _resize_event_size(event, "size")
        if _chrome_debug_enabled():
            _chrome_debug_print(
                "client-resize: old=%sx%s new=%sx%s frameless=%s"
                % (old_w, old_h, new_w, new_h, self._frameless_frame is not None)
            )
        self._begin_viewport_interaction_mode()
        self._disc_generation = int(self._disc_generation) + 1
        if self._frameless_frame is None and self.menu_button is not None:
            button_size = self.menu_button.size()
            self.menu_button.move(self.client_width() - button_size.width(), 0)

        # Invalidate the composition cache since the size has changed
        self._compositor.invalidate()
        self.request_sky_data_update()
        self.request_client_update()
        self.start_background_cloud_update(reason="resize")
        self._raise_overlay_widgets()
        if self.size_grip is not None:
            if _chrome_debug_enabled():
                _chrome_debug_print(
                    "client-resize-post: grip=%s visible=%s parent=%s"
                    % (
                        _widget_rect_tuple(self.size_grip),
                        _widget_visible(self.size_grip),
                        _widget_parent_name(self.size_grip),
                    )
                )

    def _discard_stale_disc_images(self) -> None:
        discarded = False
        if self.state.sky_disc_image is not None:
            self.state.sky_disc_image = None
            discarded = True
        if self.cloud_state.image is not None:
            self.cloud_state.image = None
            discarded = True
        if self.cloud_state.missing_mask is not None:
            self.cloud_state.missing_mask = None
            discarded = True
        if self.cloud_state.cloud_amount_field is not None:
            self.cloud_state.cloud_amount_field = None
            discarded = True
        if discarded:
            self.cloud_state.render_key = None
            self.cloud_state.request_id = None
            self.cloud_state.missing_mask_key = None
            self._compositor.invalidate()

    def _begin_interaction_mode(self) -> None:
        self.state.interaction_mode = True
        self._interaction_idle_timer.start()

    def _end_interaction_mode(self) -> None:
        if not self.state.interaction_mode:
            return
        self.state.interaction_mode = False
        self.request_sky_data_update(reason="interaction-idle")
        self.start_background_cloud_update(reason="view-change-idle")
        self.start_background_terrain_horizon_update(reason="view-change-idle")

    def _begin_viewport_interaction_mode(
        self,
        preserve_cloud_buffers: bool = False,
        *,
        start_idle_timer: bool = True,
    ) -> None:
        self.state.viewport_interaction_mode = True
        self.state.viewport_interaction_release_pending = False
        cloud_state = self.cloud_state
        cloud_controller = self._cloud_controller
        preserve_cloud_buffers = bool(preserve_cloud_buffers)
        cleared_cloud = False
        if cloud_state is not None:
            if not preserve_cloud_buffers:
                if cloud_state.image is not None:
                    cloud_state.image = None
                    cleared_cloud = True
                if cloud_state.missing_mask is not None:
                    cloud_state.missing_mask = None
                    cleared_cloud = True
                if cloud_state.cloud_amount_field is not None:
                    cloud_state.cloud_amount_field = None
                    cleared_cloud = True
                cloud_state.render_key = None
                cloud_state.request_id = None
                cloud_state.missing_mask_key = None
        if cleared_cloud:
            self._compositor.invalidate()
        if cloud_controller is not None:
            invalidate = cloud_controller.invalidate_pending_render_results
            if callable(invalidate):
                invalidate()
        if start_idle_timer:
            self._viewport_interaction_idle_timer.start()

    def _update_viewport_interaction_stars(self) -> None:
        if self.state.celestial_data is None:
            self.state.viewport_interaction_stars = None
            return
        stars, _location = calculate_visible_stars(
            self.star_catalog_np,
            self.viewer_data.location[0],
            self.viewer_data.location[1],
            self.viewer_data.observer_height_m,
            self._current_time_obj(),
            self.state.render_view_center,
            max_vmag=4.0,
            subset_indices=self.star_catalog_lod6_indices,
        )
        self.state.viewport_interaction_stars = stars

    def _end_viewport_interaction_mode(
        self,
        reason: str = "viewport-interaction-idle",
    ) -> None:
        if not self.state.viewport_interaction_mode:
            return
        self.request_sky_data_update(reason=reason)
        if reason.endswith("release"):
            self.state.viewport_interaction_release_pending = True
            return
        self.state.viewport_interaction_mode = False
        self.state.viewport_interaction_stars = None
        refresh_reason = (
            "view-change-release"
            if reason.endswith("release")
            else "view-change-idle"
        )
        self.start_background_cloud_update(reason=refresh_reason)
        self.start_background_terrain_horizon_update(reason=refresh_reason)
        self.request_client_update()

    def show_menu(self) -> None:
        if self.menu_button is None:
            return
        self._sync_view_altitude_actions()
        menu_pos = self.menu_button.mapToGlobal(QPoint(0, self.menu_button.height()))
        self.menu.exec(menu_pos)

    def _current_time_obj(self) -> astropy.time.Time:
        return astropy.time.Time(self._target_time_utc())

    def _open_named_star_jump_dialog(self) -> None:
        dialog = NamedStarJumpDialog(self._named_stars_by_band, self)
        if dialog.exec() == 0:
            return
        star = dialog.selected_star()
        if star is None:
            return
        self._jump_to_named_star(star)

    def _open_named_star_search_dialog(self) -> None:
        dialog = NamedStarSearchDialog(
            self._named_stars_search_all,
            self,
            cli_view_center_alt_specified=bool(
                getattr(self, "_search_view_center_alt_specified", False)
            ),
            cli_view_center_az_specified=bool(
                getattr(self, "_search_view_center_az_specified", False)
            ),
            satellite_search_callback=self._search_satellite_targets,
            jpl_search_callback=self._search_jpl_targets,
        )
        if dialog.exec() == 0:
            return
        if dialog.should_clear_persistent_marker():
            clear_persistent_search = getattr(self, "_clear_persistent_search", None)
            if callable(clear_persistent_search):
                clear_persistent_search()
            return
        target = dialog.selected_target()
        if target is None:
            return
        self._jump_to_search_target(target)

    def _open_place_search_dialog(self) -> None:
        dialog = PlaceSearchDialog(
            self._search_place_jump_targets,
            self,
            cli_view_center_alt_specified=bool(
                getattr(self, "_search_view_center_alt_specified", False)
            ),
            cli_view_center_az_specified=bool(
                getattr(self, "_search_view_center_az_specified", False)
            ),
        )
        if dialog.exec() == 0:
            return
        target = dialog.selected_target()
        if target is None:
            return
        self._jump_to_search_target(target)

    def _search_jpl_targets(self, query: str) -> list[SearchJumpTarget]:
        target_time_utc = self._target_time_utc()
        delta_t = getattr(self, "delta_t", None)
        _jpl_debug_print(
            "search "
            f"query={query!r} "
            f"delta_t={delta_t!r} "
            f"target_time_utc={target_time_utc.astimezone(timezone.utc).isoformat()}"
        )
        return search_jpl_targets(
            query,
            target_time_utc=target_time_utc,
            lookup_fetch=fetch_horizons_lookup,
        )

    def _search_satellite_targets(self, query: str) -> list[SearchJumpTarget]:
        return search_satellite_targets(query, target_time_utc=self._target_time_utc())

    def _find_satellite_jump_altaz(self, object_key: str) -> tuple[float, float] | None:
        records = self.satellite_state.records_by_group or None
        if not records:
            records = self._load_cached_satellite_records(tuple(self._enabled_satellite_groups))
        if not records:
            return None
        return find_satellite_altaz(
            records,
            object_key=object_key,
            observer_lat=float(self.viewer_data.location[0]),
            observer_lon=float(self.viewer_data.location[1]),
            observer_height_m=float(self.viewer_data.observer_height_m),
            time_obj=self._current_time_obj(),
        )

    def _build_horizons_command(self, target: dict[str, object], *, group: str) -> str:
        spkid = str(target.get("spkid", "")).strip()
        if spkid:
            if group == "sb":
                return f"DES={spkid};"
            return spkid
        pdes = str(target.get("pdes", "")).strip()
        if pdes:
            if group == "sb":
                return f"DES={pdes};"
            return pdes
        name = str(target.get("name", "")).strip()
        if name:
            if group == "sb":
                return name if name.endswith(";") else f"{name};"
            return name
        return ""

    def _extract_horizons_altaz(self, rows: list[list[str]]) -> tuple[float, float] | None:
        for row in rows:
            numeric_values: list[float] = []
            for value in row:
                try:
                    numeric_values.append(float(str(value).strip()))
                except (TypeError, ValueError):
                    continue
            if len(numeric_values) >= 2:
                # Horizons observer CSV reports azimuth first and elevation second.
                return numeric_values[1], numeric_values[0]
        return None

    def _jpl_small_body_persistent_target(self) -> SearchJumpTarget | None:
        target = self.state.persistent_search_target
        if target is None:
            return None
        if not bool(getattr(target, "persistent_keep_marker", False)):
            return None
        return target

    def _clear_persistent_search(self) -> None:
        self.state.persistent_search_target = None
        self.state.persistent_search_reference_time_utc = None
        self.state.persistent_search_next_refresh_utc = None
        self.state.persistent_search_last_refresh_utc = None
        self.state.persistent_search_last_error = None
        if hasattr(self, "_persistent_search_update_timer") and self._persistent_search_update_timer.isActive():
            self._persistent_search_update_timer.stop()

    def _schedule_persistent_search_refresh(self) -> None:
        if self._is_shutting_down:
            return
        target = self._jpl_small_body_persistent_target()
        if target is None:
            self._clear_persistent_search()
            return
        next_refresh_utc = self.state.persistent_search_next_refresh_utc
        if next_refresh_utc is None:
            reference_time_utc = target.target_time_utc or self._target_time_utc()
            self.state.persistent_search_reference_time_utc = reference_time_utc
            next_refresh_utc = reference_time_utc + timedelta(hours=1)
            self.state.persistent_search_next_refresh_utc = next_refresh_utc
        delay_ms = int(
            max(
                0.0,
                round(
                    (next_refresh_utc - datetime.now(timezone.utc)).total_seconds()
                    * 1000.0
                ),
            )
        )
        self._persistent_search_update_timer.start(delay_ms)

    def _on_persistent_search_refresh_timer(self) -> None:
        if self._is_shutting_down:
            return
        target = self._jpl_small_body_persistent_target()
        if target is None:
            self._clear_persistent_search()
            return
        query_time_utc = self.state.persistent_search_next_refresh_utc
        if query_time_utc is None:
            query_time_utc = target.target_time_utc or self._target_time_utc()
        if self._jpl_small_body_controller is None:
            return
        started = self._jpl_small_body_controller.update(
            observer_lat=float(self.viewer_data.location[0]),
            observer_lon=float(self.viewer_data.location[1]),
            observer_height_m=float(self.viewer_data.observer_height_m),
            target=target,
            target_time_utc=query_time_utc,
            reason="timer",
        )
        if not started:
            return

    def _on_jpl_started(self, payload: object) -> None:
        banner = ""
        if isinstance(payload, dict):
            banner = str(payload.get("banner", "")).strip()
        if banner:
            self.state.persistent_search_last_error = None
        self.request_client_update()

    def _log_persistent_search_target_update(
        self,
        *,
        action: str,
        target: SearchJumpTarget,
        target_time_utc: datetime,
        alt_deg: float,
        az_deg: float,
    ) -> None:
        logger.info(
            "JPL persistent target %s: label=%s kind=%s group=%s target_time_utc=%s alt=%.1f az=%.1f command=%s",
            action,
            str(getattr(target, "label", "")).strip() or "<unnamed>",
            str(getattr(target, "kind", "")).strip() or "<unknown>",
            str(getattr(target, "jpl_group", "")).strip() or "<none>",
            target_time_utc.astimezone(timezone.utc).isoformat(),
            float(alt_deg),
            float(az_deg) % 360.0,
            str(getattr(target, "command", "")).strip() or "<missing>",
        )

    def _on_jpl_ready(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        target = payload.get("target")
        if not isinstance(target, SearchJumpTarget):
            return
        current_target = self.state.persistent_search_target
        if current_target is None:
            return
        if current_target.command != target.command or current_target.label != target.label:
            return
        target_time_utc = payload.get("target_time_utc")
        if not isinstance(target_time_utc, datetime):
            target_time_utc = current_target.target_time_utc or self._target_time_utc()
        horizons_epoch_utc = payload.get("horizons_epoch_utc")
        horizons_position_km = payload.get("horizons_position_km")
        horizons_velocity_km_s = payload.get("horizons_velocity_km_s")
        viewer_data = getattr(self, "viewer_data", None)
        projected_altaz = None
        if not (
            viewer_data is not None
            and isinstance(horizons_epoch_utc, datetime)
            and isinstance(horizons_position_km, (list, tuple))
            and isinstance(horizons_velocity_km_s, (list, tuple))
        ):
            _jpl_debug_print(
                "ready-skip "
                f"label={target.label} command={target.command} "
                f"target_time_utc={target_time_utc.astimezone(timezone.utc).isoformat()} "
                f"epoch={horizons_epoch_utc!r} pos={horizons_position_km!r} vel={horizons_velocity_km_s!r}"
            )
            return
        vector_target = replace(
            current_target,
            horizons_epoch_utc=horizons_epoch_utc,
            horizons_position_km=tuple(horizons_position_km),
            horizons_velocity_km_s=tuple(horizons_velocity_km_s),
        )
        projected_altaz = project_jpl_target_altaz_from_state_vector(
            vector_target,
            observer_lat=float(viewer_data.location[0]),
            observer_lon=float(viewer_data.location[1]),
            observer_height_m=float(viewer_data.observer_height_m),
            time_obj=self._current_time_obj(),
        )
        if projected_altaz is None:
            _jpl_debug_print(
                "ready-project-none "
                f"label={target.label} command={target.command} "
                f"target_time_utc={target_time_utc.astimezone(timezone.utc).isoformat()} "
                f"epoch={horizons_epoch_utc.astimezone(timezone.utc).isoformat()} "
                f"pos={tuple(horizons_position_km)!r} vel={tuple(horizons_velocity_km_s)!r}"
            )
            return
        alt_deg, az_deg = projected_altaz
        _jpl_debug_print(
            "ready "
            f"label={target.label} command={target.command} "
            f"target_time_utc={target_time_utc.astimezone(timezone.utc).isoformat()} "
            f"epoch={horizons_epoch_utc.astimezone(timezone.utc).isoformat()} "
            f"pos={tuple(horizons_position_km)!r} vel={tuple(horizons_velocity_km_s)!r} "
            f"projected_alt={float(alt_deg):.3f} projected_az={float(az_deg) % 360.0:.3f}"
        )
        updated_target = replace(
            vector_target,
            alt_deg=alt_deg,
            az_deg=az_deg,
            target_time_utc=target_time_utc,
        )
        self._log_persistent_search_target_update(
            action="refreshed",
            target=updated_target,
            target_time_utc=target_time_utc,
            alt_deg=alt_deg,
            az_deg=az_deg,
        )
        self.state.persistent_search_target = updated_target
        self.state.persistent_search_reference_time_utc = target_time_utc
        self.state.persistent_search_last_refresh_utc = payload.get("refreshed_at_utc")
        self.state.persistent_search_last_error = None
        self.state.persistent_search_next_refresh_utc = target_time_utc + timedelta(hours=1)
        self.request_client_update()
        self._schedule_persistent_search_refresh()

    def _on_jpl_failed(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return
        target = payload.get("target")
        if not isinstance(target, SearchJumpTarget):
            return
        current_target = self.state.persistent_search_target
        if current_target is None:
            return
        if current_target.command != target.command or current_target.label != target.label:
            return
        refreshed_at_utc = payload.get("refreshed_at_utc")
        if not isinstance(refreshed_at_utc, datetime):
            refreshed_at_utc = datetime.now(timezone.utc)
        error_text = str(payload.get("error", "")).strip()
        banner_text = str(payload.get("banner", "")).strip()
        last_error = error_text or banner_text
        if last_error.casefold() == "none":
            last_error = ""
        self.state.persistent_search_last_error = last_error
        self.state.persistent_search_last_refresh_utc = refreshed_at_utc
        self.state.persistent_search_next_refresh_utc = refreshed_at_utc + timedelta(hours=1)
        self.request_client_update()
        self._schedule_persistent_search_refresh()

    def _jump_to_search_target(self, target: SearchJumpTarget) -> None:
        target_kind = target.kind
        target_time_utc_fn = getattr(self, "_target_time_utc", None)
        if callable(target_time_utc_fn):
            current_time = target_time_utc_fn()
        else:
            current_time = datetime.now(timezone.utc)
        current_time_obj_fn = getattr(self, "_current_time_obj", None)
        if callable(current_time_obj_fn):
            current_time_obj = current_time_obj_fn()
        else:
            current_time_obj = astropy.time.Time(current_time)
        if target_kind == "satellite":
            if target.alt_deg is not None and target.az_deg is not None:
                target_altaz = (float(target.alt_deg), float(target.az_deg) % 360.0)
            else:
                target_altaz = self._find_satellite_jump_altaz(target.object_key or target.label)
            if target_altaz is None:
                self.satellite_state.set_banner(
                    f"Satellites: {target.label} not available"
                )
                self.request_client_update()
                return
            target_alt = float(target_altaz[0])
            target_az = float(target_altaz[1]) % 360.0
        elif target_kind == "place":
            if target.latitude_deg is None or target.longitude_deg is None:
                return
            projection = project_place_target_to_altaz(
                observer_latitude_deg=float(self.viewer_data.location[0]),
                observer_longitude_deg=float(self.viewer_data.location[1]),
                observer_height_m=float(self.viewer_data.observer_height_m),
                target_latitude_deg=float(target.latitude_deg),
                target_longitude_deg=float(target.longitude_deg),
            )
            target_alt = float(projection.alt_deg)
            target_az = float(projection.az_deg) % 360.0
        elif target.kind in ("jpl_small_body", "jpl_body"):
            state_vector_target = target
            if (
                target.horizons_epoch_utc is None
                or target.horizons_position_km is None
                or target.horizons_velocity_km_s is None
            ):
                state_vector = resolve_jpl_target_state_vector(
                    target,
                    target_time_utc=target.target_time_utc or current_time,
                )
                if state_vector is None:
                    _jpl_debug_print(
                        "jump-resolve-none "
                        f"label={target.label} command={target.command} "
                        f"target_time_utc={(target.target_time_utc or current_time).astimezone(timezone.utc).isoformat()}"
                    )
                    return
                horizons_epoch_utc, horizons_position_km, horizons_velocity_km_s = state_vector
                state_vector_target = replace(
                    target,
                    horizons_epoch_utc=horizons_epoch_utc,
                    horizons_position_km=horizons_position_km,
                    horizons_velocity_km_s=horizons_velocity_km_s,
                )
            target_altaz = project_jpl_target_altaz_from_state_vector(
                state_vector_target,
                observer_lat=float(self.viewer_data.location[0]),
                observer_lon=float(self.viewer_data.location[1]),
                observer_height_m=float(self.viewer_data.observer_height_m),
                time_obj=current_time_obj,
            )
            if target_altaz is None:
                _jpl_debug_print(
                    "jump-project-none "
                    f"label={target.label} command={target.command} "
                    f"target_time_utc={(target.target_time_utc or current_time).astimezone(timezone.utc).isoformat()} "
                    f"epoch={state_vector_target.horizons_epoch_utc.astimezone(timezone.utc).isoformat() if state_vector_target.horizons_epoch_utc else '<none>'} "
                    f"pos={state_vector_target.horizons_position_km!r} vel={state_vector_target.horizons_velocity_km_s!r}"
                )
                return
            target_alt = float(target_altaz[0])
            target_az = float(target_altaz[1]) % 360.0
            _jpl_debug_print(
                "jump "
                f"label={target.label} command={target.command} "
                f"target_time_utc={(target.target_time_utc or current_time).astimezone(timezone.utc).isoformat()} "
                f"epoch={state_vector_target.horizons_epoch_utc.astimezone(timezone.utc).isoformat() if state_vector_target.horizons_epoch_utc else '<none>'} "
                f"pos={state_vector_target.horizons_position_km!r} vel={state_vector_target.horizons_velocity_km_s!r} "
                f"projected_alt={target_alt:.3f} projected_az={target_az:.3f}"
            )
        else:
            target_alt, target_az = radec_to_altaz(
                target.ra_hours,
                target.dec_deg,
                float(self.viewer_data.location[0]),
                float(self.viewer_data.location[1]),
                float(self.viewer_data.observer_height_m),
                current_time,
            )
            target_alt = float(target_alt)
            target_az = float(target_az) % 360.0
        # When jumping to a target, keep the jump highlight altitude as-reported
        # (may be negative), but clamp the actual view center to the horizon (0°)
        # so the view doesn't go below the horizon automatically.
        base_alt, base_az = tuple(
            getattr(self, "_search_view_center_base", self.viewer_data.view_center)
        )
        preserve_cli_view_center = getattr(
            target, "preserve_cli_view_center", None
        )
        fixed_alt = bool(getattr(self, "_search_view_center_alt_specified", False))
        fixed_az = bool(getattr(self, "_search_view_center_az_specified", False))
        if preserve_cli_view_center is False:
            fixed_alt = False
            fixed_az = False
        new_alt = float(base_alt) if fixed_alt else max(
            OBSERVER_MIN_ALT_DEG, min(90.0, target_alt)
        )
        new_az = float(base_az) % 360.0 if fixed_az else target_az
        self.viewer_data.view_center = (new_alt, new_az)
        self._sync_view_altitude_actions()

        self.state.jump_highlight_name = target.label
        self.state.jump_highlight_altaz = (target_alt, target_az)
        self.state.jump_highlight_until_ms = (time.monotonic() * 1000.0) + 3000.0
        if bool(getattr(target, "persistent_keep_marker", False)):
            reference_time_utc = target.target_time_utc or current_time
            horizons_epoch_utc = state_vector_target.horizons_epoch_utc
            horizons_position_km = state_vector_target.horizons_position_km
            horizons_velocity_km_s = state_vector_target.horizons_velocity_km_s
            updated_target = replace(
                state_vector_target,
                alt_deg=target_alt,
                az_deg=target_az,
                persistent_keep_marker=True,
                horizons_epoch_utc=horizons_epoch_utc,
                horizons_position_km=horizons_position_km,
                horizons_velocity_km_s=horizons_velocity_km_s,
            )
            self._log_persistent_search_target_update(
                action="set",
                target=updated_target,
                target_time_utc=reference_time_utc,
                alt_deg=target_alt,
                az_deg=target_az,
            )
            self.state.persistent_search_target = updated_target
            self.state.persistent_search_reference_time_utc = reference_time_utc
            self.state.persistent_search_last_error = None
            self.state.persistent_search_last_refresh_utc = target.target_time_utc
            if target.kind in {"jpl_small_body", "jpl_body"}:
                self.state.persistent_search_next_refresh_utc = (
                    reference_time_utc + timedelta(hours=1)
                )
                self._schedule_persistent_search_refresh()
            else:
                self.state.persistent_search_next_refresh_utc = None
        else:
            clear_persistent_search = getattr(self, "_clear_persistent_search", None)
            if callable(clear_persistent_search):
                clear_persistent_search()

        self._begin_interaction_mode()
        self.request_sky_data_update()
        self.request_client_update()

    def _search_place_jump_targets(self, query: str) -> list[SearchJumpTarget]:
        candidates = search_place_candidates(query)
        return build_place_search_jump_targets(candidates)

    def _jump_to_named_star(self, star: NamedStarShortcut) -> None:
        self._jump_to_search_target(
            SearchJumpTarget(
                label=star.name,
                ra_hours=star.ra_hours,
                dec_deg=star.dec_deg,
                kind=star.kind,
                subtitle=star.subtitle or f"Vmag {star.vmag:.2f}",
                sort_key=(star.vmag, star.name.casefold()),
                object_key=star.name if star.kind == "satellite" else "",
            )
        )

    def _load_cached_satellite_records(
        self, enabled_groups: tuple[str, ...]
    ) -> dict[str, list[dict[str, object]]]:
        records_by_group: dict[str, list[dict[str, object]]] = {}
        for group_key in enabled_groups:
            cache_scope = None
            if group_key == SATELLITE_HORIZONS_CACHE_KEY:
                cache_scope = satellite_cache_scope_key(
                    observer_lat=float(self.viewer_data.location[0]),
                    observer_lon=float(self.viewer_data.location[1]),
                    observer_height_m=float(self.viewer_data.observer_height_m),
                )
            if cache_scope is None:
                cached = load_satellite_cache(group_key)
            else:
                cached = load_satellite_cache(group_key, cache_scope_key=cache_scope)
            if cached is None:
                continue
            records_by_group[str(group_key)] = list(cached.records)
        return records_by_group

    def _begin_shutdown(self) -> None:
        """Stop scheduling new background work while the app is closing."""
        if self._is_shutting_down:
            return
        self._is_shutting_down = True
        self._show_shutdown_message()
        try:
            self._sky_worker.shutdown()
            if self._cloud_controller is not None:
                self._cloud_controller.shutdown()
            if self._satellite_controller is not None:
                self._satellite_controller.shutdown()
            if self._aircraft_controller is not None:
                self._aircraft_controller.shutdown()
            if self._jpl_small_body_controller is not None:
                self._jpl_small_body_controller.shutdown()
            if self._terrain_horizon_controller is not None:
                self._terrain_horizon_controller.shutdown()
            if self._urban_outline_controller is not None:
                self._urban_outline_controller.shutdown()
            if self._sky_data_update_timer.isActive():
                self._sky_data_update_timer.stop()
            if self._asterism_check_timer.isActive():
                self._asterism_check_timer.stop()
            if self._cloud_update_timer.isActive():
                self._cloud_update_timer.stop()
            if self._satellite_update_timer.isActive():
                self._satellite_update_timer.stop()
            if self._overlay_projection_timer.isActive():
                self._overlay_projection_timer.stop()
            if self._aircraft_update_timer.isActive():
                self._aircraft_update_timer.stop()
            if hasattr(self, "_persistent_search_update_timer") and self._persistent_search_update_timer.isActive():
                self._persistent_search_update_timer.stop()
            if self._interaction_idle_timer.isActive():
                self._interaction_idle_timer.stop()
            if self._viewport_interaction_idle_timer.isActive():
                self._viewport_interaction_idle_timer.stop()
        finally:
            self._hide_shutdown_message()

    def _on_asterism_check_tick(self) -> None:
        if self._is_shutting_down:
            return
        if not self.show_asterisms:
            return
        if self.state.mouse_pos is None or self.state.celestial_data is None:
            return

        render_viewer = self._viewer_data_for_render()
        geometry = render_geometry.get_screen_geometry(
            self.client_width(),
            self.client_height(),
            render_viewer.view_center[0],
        )
        highlighted = render_stars.find_highlighted_object(
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
        self.request_client_update()

    def _predicted_cloud_satellite(self) -> str:
        lat, lon = self.viewer_data.location
        return pick_satellite(lat, lon, ("AUTO",))

    def _satellite_layer_enabled(self) -> bool:
        return self._satellite_toggle_supported and self.satellite_opacity > 0.0

    def _stop_satellite_timers(self) -> None:
        if self._satellite_update_timer.isActive():
            self._satellite_update_timer.stop()
        self._sync_overlay_projection_timer()

    def _target_time_utc(self) -> datetime:
        delta_t = getattr(self, "delta_t", None)
        if delta_t is None:
            delta_t = timedelta(0)
        target_time_utc = target_time_utc_from_delta(delta_t)
        _jpl_debug_print(
            "target-time "
            f"delta_t={delta_t!r} "
            f"target_time_utc={target_time_utc.astimezone(timezone.utc).isoformat()}"
        )
        return target_time_utc

    def _satellite_validity_remaining_ms(self) -> int | None:
        refreshed_at_utc = self.satellite_state.refreshed_at_utc
        if refreshed_at_utc is None:
            return None
        age_seconds = abs(
            (datetime.now(timezone.utc) - refreshed_at_utc).total_seconds()
        )
        remaining_seconds = SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS - age_seconds
        return max(0, int(round(remaining_seconds * 1000.0)))

    def _schedule_next_satellite_refresh(self, delay_ms: int | None = None) -> None:
        if not self._satellite_layer_enabled() or self._is_shutting_down:
            return
        interval_ms = delay_ms
        if interval_ms is None:
            interval_ms = self._satellite_validity_remaining_ms()
        if interval_ms is None:
            interval_ms = SATELLITE_ELEMENT_REFRESH_INTERVAL_SECONDS * 1000
        self._satellite_update_timer.start(max(0, interval_ms))

    def _schedule_satellite_retry_after_failure(self) -> None:
        if not self._satellite_layer_enabled() or self._is_shutting_down:
            return
        self._satellite_update_timer.start(SATELLITE_FAILURE_RETRY_SECONDS * 1000)

    def _on_satellite_refresh_timer(self) -> None:
        if not self._satellite_layer_enabled():
            return
        started = self.start_background_satellite_update(reason="timer")
        if not started:
            self._schedule_next_satellite_refresh()

    def _on_overlay_projection_timer(self) -> None:
        refresh_persistent_search = getattr(
            self,
            "refresh_projected_persistent_search_target",
            None,
        )
        if callable(refresh_persistent_search):
            refresh_persistent_search()
        if self._satellite_layer_enabled():
            self.refresh_projected_satellite_overlay()
        if self._aircraft_layer_enabled():
            self.refresh_projected_aircraft_overlay()

    def _sync_overlay_projection_timer(self) -> None:
        if self._is_shutting_down:
            if self._overlay_projection_timer.isActive():
                self._overlay_projection_timer.stop()
            return
        any_enabled = self._satellite_layer_enabled() or self._aircraft_layer_enabled()
        if any_enabled:
            if not self._overlay_projection_timer.isActive():
                self._overlay_projection_timer.start()
        elif self._overlay_projection_timer.isActive():
            self._overlay_projection_timer.stop()

    def _enable_satellite_layer(self, *, reason: str) -> None:
        if not self._satellite_layer_enabled():
            return
        self._sync_overlay_projection_timer()
        if (
            self.satellite_state.records_by_group
            and self.satellite_state.element_epoch_utc is not None
        ):
            self.refresh_projected_satellite_overlay()
            remaining_ms = self._satellite_validity_remaining_ms()
            if remaining_ms is None:
                self.start_background_satellite_update(reason=reason)
                return
            if remaining_ms <= 0:
                self.start_background_satellite_update(reason=reason)
            else:
                self._schedule_next_satellite_refresh(remaining_ms)
            return
        self.start_background_satellite_update(reason=reason)

    def _aircraft_layer_enabled(self) -> bool:
        return self._aircraft_toggle_supported and self.aircraft_opacity > 0.0

    def _stop_aircraft_timers(self) -> None:
        if self._aircraft_update_timer.isActive():
            self._aircraft_update_timer.stop()
        self._sync_overlay_projection_timer()

    def _aircraft_cache_age_seconds(self) -> float | None:
        last_success_utc = self.aircraft_state.last_success_utc
        if last_success_utc is None:
            return None
        return max(0.0, (datetime.now(timezone.utc) - last_success_utc).total_seconds())

    def _schedule_next_aircraft_refresh(self, delay_ms: int | None = None) -> None:
        if not self._aircraft_layer_enabled() or self._is_shutting_down:
            return
        interval_ms = (
            AIRCRAFT_REFRESH_INTERVAL_SECONDS * 1000
            if delay_ms is None
            else int(delay_ms)
        )
        self._aircraft_update_timer.start(max(0, interval_ms))

    def _on_aircraft_refresh_timer(self) -> None:
        if not self._aircraft_layer_enabled():
            return
        started = self.start_background_aircraft_update(reason="timer")
        if not started:
            self._schedule_next_aircraft_refresh()

    def _enable_aircraft_layer(self, *, reason: str) -> None:
        if not self._aircraft_layer_enabled():
            return
        self._sync_overlay_projection_timer()
        if (
            self.aircraft_state.snapshots
            and self.aircraft_state.last_success_utc is not None
        ):
            self.refresh_projected_aircraft_overlay()
            age_seconds = self._aircraft_cache_age_seconds()
            if age_seconds is None:
                self.start_background_aircraft_update(reason=reason)
                return
            remaining_ms = max(
                0,
                int(round((AIRCRAFT_REFRESH_INTERVAL_SECONDS - age_seconds) * 1000.0)),
            )
            if remaining_ms <= 0:
                self.start_background_aircraft_update(reason=reason)
            else:
                self._schedule_next_aircraft_refresh(remaining_ms)
            return
        self.start_background_aircraft_update(reason=reason)

    def toggle_enlarge_moon(self) -> None:
        self.enlarge_moon = not self.enlarge_moon
        if (
            self._action_enlarge_moon is not None
            and self._action_enlarge_moon.isChecked() != self.enlarge_moon
        ):
            self._action_enlarge_moon.setChecked(self.enlarge_moon)
        self.request_client_update()

    def toggle_clouds(self) -> None:
        if (
            not self._cloud_toggle_supported
            or self._clouddisc is None
            or not self._cloud_gui_allowed
        ):
            if self._action_toggle_clouds is not None:
                self._action_toggle_clouds.setChecked(False)
            return

        enable_clouds = self.cloud_disc_alpha <= 0.0
        self.cloud_disc_alpha = self._cloud_alpha_when_enabled if enable_clouds else 0.0
        if (
            self._action_toggle_clouds is not None
            and self._action_toggle_clouds.isChecked() != enable_clouds
        ):
            self._action_toggle_clouds.setChecked(enable_clouds)

        if enable_clouds:
            self.start_background_cloud_update(reason="toggle-on")
            if (
                self._sky_data_update_timer.isActive()
                and not self._cloud_update_timer.isActive()
            ):
                self._cloud_update_timer.start()
        elif self._cloud_update_timer.isActive():
            self._cloud_update_timer.stop()

        self.request_client_update()

    def toggle_satellites(self) -> None:
        if not self._satellite_toggle_supported or not self._satellite_gui_allowed:
            if self._action_toggle_satellites is not None:
                self._action_toggle_satellites.setChecked(self.satellite_opacity > 0.0)
            return

        enable_satellites = self.satellite_opacity <= 0.0
        self.satellite_opacity = (
            self._satellite_opacity_when_enabled if enable_satellites else 0.0
        )
        if (
            self._action_toggle_satellites is not None
            and self._action_toggle_satellites.isChecked() != enable_satellites
        ):
            self._action_toggle_satellites.setChecked(enable_satellites)

        if enable_satellites:
            self._enable_satellite_layer(reason="toggle-on")
        else:
            self._stop_satellite_timers()

        self.request_client_update()

    def toggle_aircraft(self) -> None:
        if not self._aircraft_toggle_supported or not self._aircraft_gui_allowed:
            if self._action_toggle_aircraft is not None:
                self._action_toggle_aircraft.setChecked(self.aircraft_opacity > 0.0)
            return

        enable_aircraft = self.aircraft_opacity <= 0.0
        self.aircraft_opacity = (
            self._aircraft_opacity_when_enabled if enable_aircraft else 0.0
        )
        if (
            self._action_toggle_aircraft is not None
            and self._action_toggle_aircraft.isChecked() != enable_aircraft
        ):
            self._action_toggle_aircraft.setChecked(enable_aircraft)

        if enable_aircraft:
            self._enable_aircraft_layer(reason="toggle-on")
        else:
            self._stop_aircraft_timers()

        self.request_client_update()

    def toggle_dso(self) -> None:
        if self.dso_catalog_np is None:
            self.show_dso = False
            if self._action_toggle_dso is not None:
                self._action_toggle_dso.setChecked(False)
            return
        self.show_dso = not self.show_dso
        if (
            self._action_toggle_dso is not None
            and self._action_toggle_dso.isChecked() != self.show_dso
        ):
            self._action_toggle_dso.setChecked(self.show_dso)
        self.request_client_update()

    def toggle_asterisms(self) -> None:
        self.show_asterisms = not self.show_asterisms
        if (
            self._action_toggle_asterisms is not None
            and self._action_toggle_asterisms.isChecked() != self.show_asterisms
        ):
            self._action_toggle_asterisms.setChecked(self.show_asterisms)
        self.request_client_update()

    def toggle_guidelines(self) -> None:
        self.show_guidelines = not self.show_guidelines
        if (
            self._action_toggle_guidelines is not None
            and self._action_toggle_guidelines.isChecked() != self.show_guidelines
        ):
            self._action_toggle_guidelines.setChecked(self.show_guidelines)
        self.request_client_update()

    def toggle_observation_info(self) -> None:
        # If CLI explicitly disabled the overlay, the menu toggle should be inert.
        if self.observation_info_mode == "off":
            self.show_observation_info = False
            if (
                self._action_toggle_observation_info is not None
                and self._action_toggle_observation_info.isChecked() != self.show_observation_info
            ):
                self._action_toggle_observation_info.setChecked(self.show_observation_info)
            return

        # Toggle visibility. When re-enabling and the overlay is pinned, ensure
        # the pinned position is enforced.
        self.show_observation_info = not self.show_observation_info

        if self.show_observation_info and self.observation_info_pinned:
            mode = self.observation_info_mode
            try:
                if mode == "bottom":
                    self.state.overlay_info_bottom_left = True
                elif mode == "top":
                    self.state.overlay_info_bottom_left = False
            except Exception:
                pass

        if (
            self._action_toggle_observation_info is not None
            and self._action_toggle_observation_info.isChecked()
            != self.show_observation_info
        ):
            self._action_toggle_observation_info.setChecked(self.show_observation_info)
        self.request_client_update()

    def toggle_sky_disc(self) -> None:
        if not self._sky_disc_gui_allowed:
            if self._action_toggle_sky_disc is not None:
                self._action_toggle_sky_disc.setChecked(self.sky_disc_alpha > 0.0)
            return
        enable_sky_disc_gradient = self.sky_disc_alpha <= 0.0
        self.sky_disc_alpha = (
            self._sky_disc_alpha_when_enabled if enable_sky_disc_gradient else 0.0
        )
        if (
            self._action_toggle_sky_disc is not None
            and self._action_toggle_sky_disc.isChecked() != enable_sky_disc_gradient
        ):
            self._action_toggle_sky_disc.setChecked(enable_sky_disc_gradient)
        self._compositor.invalidate()
        self.request_sky_data_update()
        self.request_client_update()

    def toggle_terrain_horizon(self) -> None:
        if not self._terrain_horizon_gui_allowed:
            if self._action_toggle_terrain_horizon is not None:
                self._action_toggle_terrain_horizon.setChecked(
                    self.terrain_horizon_opacity > 0.0
                )
            return

        enable_terrain = self.terrain_horizon_opacity <= 0.0
        self.terrain_horizon_opacity = (
            self._terrain_horizon_opacity_when_enabled if enable_terrain else 0.0
        )
        if (
            self._action_toggle_terrain_horizon is not None
            and self._action_toggle_terrain_horizon.isChecked() != enable_terrain
        ):
            self._action_toggle_terrain_horizon.setChecked(enable_terrain)
        self._compositor.invalidate()
        if enable_terrain:
            self.start_background_terrain_horizon_update(reason="toggle-on")
        self.request_client_update()

    def toggle_earth_guide(self) -> None:
        if not self._earth_guide_gui_allowed:
            if self._action_toggle_earth_guide is not None:
                self._action_toggle_earth_guide.setChecked(self.earth_guide_opacity > 0.0)
            return

        enable_earth_guide = self.earth_guide_opacity <= 0.0
        self.earth_guide_opacity = (
            self._earth_guide_opacity_when_enabled if enable_earth_guide else 0.0
        )
        if (
            self._action_toggle_earth_guide is not None
            and self._action_toggle_earth_guide.isChecked() != enable_earth_guide
        ):
            self._action_toggle_earth_guide.setChecked(enable_earth_guide)
        self._compositor.invalidate()
        self.request_client_update()

    def toggle_urban_outline(self) -> None:
        if not self._urban_outline_gui_allowed:
            if self._action_toggle_urban_outline is not None:
                self._action_toggle_urban_outline.setChecked(
                    self.urban_outline_opacity > 0.0
                )
            return
        enable_urban_outline = self.urban_outline_opacity <= 0.0
        self.urban_outline_opacity = (
            self._urban_outline_opacity_when_enabled if enable_urban_outline else 0.0
        )
        self.show_urban_outline_layer = self.urban_outline_opacity > 0.0
        if (
            self._action_toggle_urban_outline is not None
            and self._action_toggle_urban_outline.isChecked() != enable_urban_outline
        ):
            self._action_toggle_urban_outline.setChecked(enable_urban_outline)
        if enable_urban_outline:
            self.start_background_urban_outline_update(reason="toggle-on")
        self.request_client_update()

    def toggle_fullscreen(self) -> None:
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def _handle_client_leave(self, event: QEvent) -> None:
        self.state.mouse_pos = None
        self.request_client_update()
        event.accept()

    def closeEvent(self, event: QCloseEvent) -> None:
        g = self.geometry()
        save_last_window_geometry(g.x(), g.y(), self.client_width(), self.client_height())
        self._begin_shutdown()
        super().closeEvent(event)

    def _handle_client_mouse_move(self, event: QMouseEvent) -> None:
        self.state.mouse_pos = event.pos()
        self.request_client_update()
        event.accept()

    def _rotate_view(
        self,
        d_alt: float = 0.0,
        d_az: float = 0.0,
        *,
        interactive_viewport: bool = False,
        start_viewport_idle_timer: bool = True,
    ) -> None:
        if interactive_viewport:
            if not bool(getattr(self.state, "viewport_interaction_mode", False)):
                self._begin_viewport_interaction_mode(
                    start_idle_timer=start_viewport_idle_timer
                )
        else:
            self._begin_interaction_mode()
        alt, az = self.viewer_data.view_center
        new_alt = max(OBSERVER_MIN_ALT_DEG, min(90.0, alt + d_alt))
        new_az = (az + d_az) % 360.0
        self.viewer_data.view_center = (new_alt, new_az)
        self.state.render_view_center = (new_alt, new_az)
        self._sync_view_altitude_actions()
        if interactive_viewport:
            self._update_viewport_interaction_stars()
            self.request_client_update()
            return
        self.request_sky_data_update()

    def _viewport_rotation_keys(self) -> set[int]:
        keys = getattr(self, "_viewport_rotation_keys_down", None)
        if keys is None:
            keys = set()
            self._viewport_rotation_keys_down = keys
        return keys

    def _handle_client_key_press(self, event: QKeyEvent) -> None:
        key = event.key()

        # --- View Control ---
        if key == Qt.Key.Key_Left:
            if not event.isAutoRepeat():
                self._viewport_rotation_keys().add(key)
            self._rotate_view(
                d_az=-self.state.rotation_step,
                interactive_viewport=True,
                start_viewport_idle_timer=False,
            )
            event.accept()
        elif key == Qt.Key.Key_Right:
            if not event.isAutoRepeat():
                self._viewport_rotation_keys().add(key)
            self._rotate_view(
                d_az=self.state.rotation_step,
                interactive_viewport=True,
                start_viewport_idle_timer=False,
            )
            event.accept()
        elif key == Qt.Key.Key_Up:
            if not event.isAutoRepeat():
                self._viewport_rotation_keys().add(key)
            self._rotate_view(
                d_alt=self.state.rotation_step,
                interactive_viewport=True,
                start_viewport_idle_timer=False,
            )
            event.accept()
        elif key == Qt.Key.Key_Down:
            if not event.isAutoRepeat():
                self._viewport_rotation_keys().add(key)
            self._rotate_view(
                d_alt=-self.state.rotation_step,
                interactive_viewport=True,
                start_viewport_idle_timer=False,
            )
            event.accept()

        # --- Toggles ---
        elif key == Qt.Key.Key_M:
            self.toggle_enlarge_moon()
            event.accept()
        elif key == Qt.Key.Key_C:
            self.toggle_clouds()
            event.accept()
        elif key == Qt.Key.Key_I:
            self.toggle_aircraft()
            event.accept()
        elif key == Qt.Key.Key_T:
            self.toggle_terrain_horizon()
            event.accept()
        elif key == Qt.Key.Key_E:
            self.toggle_earth_guide()
            event.accept()
        elif key == Qt.Key.Key_U:
            self.toggle_urban_outline()
            event.accept()
        elif key == Qt.Key.Key_D:
            self.toggle_dso()
            event.accept()
        elif key == Qt.Key.Key_A:
            self.toggle_asterisms()
            event.accept()
        elif key == Qt.Key.Key_G:
            self.toggle_guidelines()
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
            self._request_application_quit()
            event.accept()
        else:
            super().keyPressEvent(event)

    def _handle_client_key_release(self, event: QKeyEvent) -> None:
        key = event.key()
        if key not in {
            Qt.Key.Key_Left,
            Qt.Key.Key_Right,
            Qt.Key.Key_Up,
            Qt.Key.Key_Down,
        }:
            return
        if event.isAutoRepeat():
            event.accept()
            return
        keys_down = self._viewport_rotation_keys()
        keys_down.discard(key)
        if not keys_down:
            self._end_viewport_interaction_mode(
                reason="viewport-interaction-release"
            )
        event.accept()


class FramelessSkyWindow(SkyWindowCoreMixin, DraggableWindow):
    """Frameless host window using the custom chrome wrapper."""

    FRAMELESS_WINDOW = True

    def _install_window_host(self) -> None:
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.FramelessWindowHint)
        self._frameless_frame = FramelessWindowFrame(self, self._client_widget)
        self.setCentralWidget(self._frameless_frame)
        self.setFocusProxy(self._client_widget)
        self._frameless_frame.setFocusProxy(self._client_widget)
        self.menu_button = self._frameless_frame.menu_button
        self.size_grip = self._frameless_frame.size_grip
        self.add_drag_target(self._client_widget)
        self.add_drag_exclusions(self._frameless_frame.overlay_widgets())


class StandardSkyWindow(SkyWindowCoreMixin, QMainWindow):
    """Standard decorated host window using the platform title bar and menu bar."""

    FRAMELESS_WINDOW = False

    def _install_window_host(self) -> None:
        self.setCentralWidget(self._client_widget)
        self.setFocusProxy(self._client_widget)
        menu_bar = self.menuBar()
        menu_bar.addMenu(self.file_menu)
        menu_bar.addMenu(self.search_menu)
        menu_bar.addMenu(self.display_menu)
        menu_bar.addMenu(self.observer_view_menu)
        menu_bar.addMenu(self.help_menu)


SkyWindow = FramelessSkyWindow
