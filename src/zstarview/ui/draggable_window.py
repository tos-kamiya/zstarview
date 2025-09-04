# -*- coding: utf-8 -*-
"""
Provides a draggable, frameless main window base class.
"""
from typing import Any, Set

from PySide6.QtCore import Qt, QPoint
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QMainWindow, QWidget


class DraggableWindow(QMainWindow):
    """
    A frameless main window that can be moved by clicking and dragging.

    This class provides a base for creating custom-styled windows without the
    standard OS window frame. It implements dragging behavior that works across
    different windowing systems.

    Features:
        - Uses `startSystemMove()` when available (e.g., on Wayland) for
          native, smooth window dragging.
        - Falls back to a manual implementation of dragging on other systems.
        - Allows specific child widgets to be excluded from initiating a drag,
          so they can handle their own mouse events (e.g., buttons, sliders).
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """Initializes the DraggableWindow."""
        super().__init__(*args, **kwargs)
        self._drag_active: bool = False
        self._drag_pos: QPoint = QPoint(0, 0)
        self._drag_exclusions: Set[QWidget] = set()

    # ---- public API -------------------------------------------------
    def add_drag_exclusion(self, widget: QWidget) -> None:
        """
        Registers a child widget that should not trigger a window drag.

        When a mouse press occurs on an excluded widget or one of its children,
        the window will not initiate a drag, allowing the widget to process
        the event.

        Args:
            widget: The widget to exclude from drag initiation.
        """
        if widget is not None:
            self._drag_exclusions.add(widget)

    def add_drag_exclusions(self, widgets: list[QWidget]) -> None:
        """
        Registers multiple child widgets to be excluded from dragging.

        Args:
            widgets: A list of QWidgets to exclude.
        """
        for w in widgets:
            self.add_drag_exclusion(w)

    # ---- internals --------------------------------------------------
    def _is_in_exclusions(self, w: QWidget | None) -> bool:
        """
        Checks if a widget or any of its ancestors are in the exclusion set.

        This allows a click on a complex widget (e.g., a button with a label)
        to be correctly ignored for dragging.

        Args:
            w: The widget to check.

        Returns:
            True if the widget or any of its parents are in the exclusion list.
        """
        while w is not None:
            if w in self._drag_exclusions:
                return True
            w = w.parentWidget()
        return False

    # ---- Qt events --------------------------------------------------
    def mousePressEvent(self, event: QMouseEvent) -> None:
        """
        Handles the mouse press event to initiate a window drag.

        If the click is on the main window background (not an excluded widget),
        it attempts to start a system-native move. If that fails or is not
        available, it prepares for manual dragging.

        Args:
            event: The QMouseEvent from the mouse press.
        """
        if event.button() == Qt.LeftButton:
            child = self.childAt(event.pos())
            if self._is_in_exclusions(child):
                # The click was on an excluded widget, so let it handle the event.
                super().mousePressEvent(event)
                return

            # --- Native Drag ---
            # Prefer using the native window system's move functionality if available.
            # This provides a better user experience (e.g., on Wayland).
            wh = self.windowHandle()
            if wh and wh.startSystemMove():
                event.accept()
                return

            # --- Manual Drag Fallback ---
            # If native move is not available, fall back to manual implementation.
            self._drag_active = True
            # Get the position of the click relative to the window's top-left corner.
            try:
                global_pos = event.globalPosition().toPoint()  # Qt6
            except AttributeError:
                global_pos = event.globalPos()  # Qt5-style fallback
            self._drag_pos = global_pos - self.frameGeometry().topLeft()
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        """
        Handles the mouse move event to drag the window.

        If a manual drag is active, this method calculates the new window
        position based on the mouse movement.

        Args:
            event: The QMouseEvent from the mouse move.
        """
        # This check ensures we only move the window when a manual drag has been initiated.
        if self._drag_active and (event.buttons() & Qt.LeftButton):
            try:
                global_pos = event.globalPosition().toPoint()  # Qt6
            except AttributeError:
                global_pos = event.globalPos()  # Qt5-style fallback
            # Move the window to the new global position, adjusted by the initial click offset.
            self.move(global_pos - self._drag_pos)
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        """
        Handles the mouse release event to end the window drag.

        Args:
            event: The QMouseEvent from the mouse release.
        """
        self._drag_active = False
        event.accept()
