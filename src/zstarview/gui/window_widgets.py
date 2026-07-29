from __future__ import annotations

import logging
from html import escape
from typing import Any

from PySide6.QtCore import QEvent, QPoint, QPointF, QRect, Qt, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QKeyEvent,
    QMouseEvent,
    QPainter,
    QPaintEvent,
    QPen,
    QResizeEvent,
)
from PySide6.QtWidgets import QLabel, QTextEdit, QWidget

from ..paths import GUI_BUTTON_SIZE
from .window_render import SkyWindowRenderMixin


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


class SkyWindowClientWidget(SkyWindowRenderMixin, QWidget):
    """Client-area widget shared by frameless and decorated host windows."""

    def __init__(self, owner: Any) -> None:
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
        self._owner._layout_startup_log_overlay()
        super().resizeEvent(event)

    def leaveEvent(self, event: QEvent) -> None:
        if self._owner._startup_input_blocked():
            event.accept()
            return
        self._owner._handle_client_leave(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        if self._owner._startup_input_blocked():
            event.accept()
            return
        self._owner._handle_client_mouse_move(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if self._owner._startup_input_blocked():
            event.accept()
            return
        self._owner._handle_client_key_press(event)

    def keyReleaseEvent(self, event: QKeyEvent) -> None:
        if self._owner._startup_input_blocked():
            event.accept()
            return
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


class StartupLogOverlay(QTextEdit):
    """Temporary log view shown in the main window during startup."""

    line_received = Signal(str, int)

    def __init__(
        self,
        parent: QWidget,
        *,
        text_rgb: tuple[int, int, int] | None = None,
        background_rgba: tuple[int, int, int, int] | None = None,
    ) -> None:
        super().__init__(parent)
        if text_rgb is None:
            text_rgb = (245, 245, 245)
        if background_rgba is None:
            background_rgba = (0, 0, 0, 180)
        red, green, blue = (int(text_rgb[0]), int(text_rgb[1]), int(text_rgb[2]))
        bg_red, bg_green, bg_blue, bg_alpha = (
            int(background_rgba[0]),
            int(background_rgba[1]),
            int(background_rgba[2]),
            int(background_rgba[3]),
        )
        border_rgb = (
            (0, 0, 0)
            if ((77 * bg_red + 150 * bg_green + 29 * bg_blue) >> 8) >= 128
            else (255, 255, 255)
        )
        self._info_color = f"#{red:02x}{green:02x}{blue:02x}"
        self.setReadOnly(True)
        self.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)
        self.setUndoRedoEnabled(False)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self.setStyleSheet(
            "QTextEdit {"
            f" background-color: rgba({bg_red}, {bg_green}, {bg_blue}, {bg_alpha});"
            f" color: rgba({red}, {green}, {blue}, 230);"
            f" border: 2px solid rgba({border_rgb[0]}, {border_rgb[1]}, {border_rgb[2]}, 60);"
            " padding: 10px 12px;"
            " font-family: monospace;"
            " font-size: 10px;"
            "}"
        )
        fixed_font = QFont("Monospace")
        fixed_font.setStyleHint(QFont.StyleHint.Monospace)
        self.setFont(fixed_font)
        self.document().setMaximumBlockCount(200)
        self.line_received.connect(self._append_line)

    def append_line(self, line: str, levelno: int) -> None:
        self.line_received.emit(line, int(levelno))

    def _line_color(self, levelno: int) -> str:
        if levelno >= logging.ERROR:
            return "#ff8080"
        if levelno >= logging.WARNING:
            return "#ffcf80"
        return self._info_color

    def _append_line(self, line: str, levelno: int) -> None:
        color = self._line_color(levelno)
        self.append(f"<span style=\"color: {color}\">{escape(line)}</span>")
        scrollbar = self.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())


class ResizeGripWidget(QWidget):
    """Resize affordance with explicit paint and drag handling."""

    def __init__(self, owner: Any, parent: QWidget) -> None:
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


class MenuButtonWidget(QWidget):
    """Menu button with a custom drawn hamburger icon."""

    clicked = Signal()

    def __init__(self, chrome: Any, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._chrome = chrome
        self._hovered = False
        self._pressed = False
        self.setAttribute(Qt.WidgetAttribute.WA_Hover, True)
        self.setMouseTracking(True)

    def enterEvent(self, event) -> None:
        self._hovered = True
        self.update()
        super().enterEvent(event)

    def leaveEvent(self, event) -> None:
        self._hovered = False
        self._pressed = False
        self.update()
        super().leaveEvent(event)

    def mousePressEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            self._pressed = True
            self.update()
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            was_pressed = self._pressed
            self._pressed = False
            self.update()
            if was_pressed and self.rect().contains(event.position().toPoint()):
                self.clicked.emit()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def paintEvent(self, event: QPaintEvent) -> None:
        painter = QPainter(self)
        try:
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            chrome = self._chrome
            if self._pressed:
                background_rgba = chrome.menu_pressed_bg_rgba
            elif self._hovered:
                background_rgba = chrome.menu_hover_bg_rgba
            else:
                background_rgba = chrome.menu_fill_rgba
            background = QColor(
                background_rgba[0],
                background_rgba[1],
                background_rgba[2],
                background_rgba[3],
            )
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(background)
            painter.drawRoundedRect(self.rect().adjusted(0, 0, -1, -1), 5.0, 5.0)

            color = QColor(
                chrome.menu_button_text_rgb[0],
                chrome.menu_button_text_rgb[1],
                chrome.menu_button_text_rgb[2],
            )
            if not self.isEnabled():
                color.setAlpha(128)
            pen = QPen(color, 2.0)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            width = float(self.width())
            center_x = width / 2.0
            line_half = max(4.0, min(8.0, width / 4.0))
            center_y = float(self.height()) / 2.0
            offsets = (-5.0, 0.0, 5.0)
            for offset in offsets:
                y = center_y + offset
                painter.drawLine(
                    QPointF(center_x - line_half, y),
                    QPointF(center_x + line_half, y),
                )
        finally:
            painter.end()


class FramelessWindowFrame(QWidget):
    """Frameless-only window chrome that hosts the client widget and overlay controls."""

    def __init__(
        self, owner: Any, client_widget: SkyWindowClientWidget
    ) -> None:
        super().__init__(owner)
        self._owner = owner
        self._client_widget = client_widget
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAutoFillBackground(False)
        self._client_widget.setParent(self)

        self.menu_button = MenuButtonWidget(self._owner.theme.window_chrome, self)
        self.menu_button.setFixedSize(GUI_BUTTON_SIZE, GUI_BUTTON_SIZE)
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

    def eventFilter(self, watched: object, event: QEvent) -> bool:
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

    def resizeEvent(self, event: QResizeEvent) -> None:
        super().resizeEvent(event)
        self._layout_chrome()
