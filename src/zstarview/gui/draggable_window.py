# -*- coding: utf-8 -*-
"""Provides a draggable main window base class."""
from typing import Any, Set

from PySide6.QtCore import QEvent, QObject, QPoint, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QMainWindow, QWidget


class DraggableWindow(QMainWindow):
    """
    A main window that can optionally be moved by dragging a registered widget.

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
        self._drag_targets: Set[QWidget] = set()

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

    def add_drag_target(self, widget: QWidget) -> None:
        """Register a widget whose background should initiate window dragging."""
        if widget is None:
            return
        self._drag_targets.add(widget)
        widget.installEventFilter(self)

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

    def _begin_drag(self, child: QWidget | None, event: QMouseEvent) -> bool:
        if event.button() == Qt.LeftButton:
            if self._is_in_exclusions(child):
                return False

            wh = self.windowHandle()
            if wh and wh.startSystemMove():
                event.accept()
                return True

            self._drag_active = True
            try:
                global_pos = event.globalPosition().toPoint()  # Qt6
            except AttributeError:
                global_pos = event.globalPos()  # Qt5-style fallback
            self._drag_pos = global_pos - self.frameGeometry().topLeft()
            event.accept()
            return True
        return False

    def _update_drag(self, event: QMouseEvent) -> bool:
        if self._drag_active and (event.buttons() & Qt.LeftButton):
            try:
                global_pos = event.globalPosition().toPoint()  # Qt6
            except AttributeError:
                global_pos = event.globalPos()  # Qt5-style fallback
            self.move(global_pos - self._drag_pos)
            event.accept()
            return True
        return False

    def _end_drag(self, event: QMouseEvent) -> bool:
        self._drag_active = False
        event.accept()
        return True

    # ---- Qt events --------------------------------------------------
    def eventFilter(self, watched: QObject, event: QEvent) -> bool:
        if isinstance(watched, QWidget) and watched in self._drag_targets:
            if event.type() == QEvent.Type.MouseButtonPress and isinstance(event, QMouseEvent):
                child = watched.childAt(event.position().toPoint())
                if self._begin_drag(child, event):
                    return True
            elif event.type() == QEvent.Type.MouseMove and isinstance(event, QMouseEvent):
                if self._update_drag(event):
                    return True
            elif event.type() == QEvent.Type.MouseButtonRelease and isinstance(event, QMouseEvent):
                if self._drag_active:
                    return self._end_drag(event)
        return super().eventFilter(watched, event)

    def mousePressEvent(self, event: QMouseEvent) -> None:
        child = self.childAt(event.pos())
        if self._begin_drag(child, event):
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        """
        Handles the mouse move event to drag the window.

        If a manual drag is active, this method calculates the new window
        position based on the mouse movement.

        Args:
            event: The QMouseEvent from the mouse move.
        """
        if self._update_drag(event):
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        """
        Handles the mouse release event to end the window drag.

        Args:
            event: The QMouseEvent from the mouse release.
        """
        if self._drag_active:
            self._end_drag(event)
            return
        super().mouseReleaseEvent(event)
