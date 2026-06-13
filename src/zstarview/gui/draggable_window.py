# -*- coding: utf-8 -*-
"""Provides a draggable main window base class."""
from typing import Any, Set

from PySide6.QtCore import QEvent, QObject, QPoint, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QApplication, QMainWindow, QWidget


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
        self._drag_press_pending: bool = False
        self._drag_press_pos: QPoint = QPoint(0, 0)
        self._drag_press_child: QWidget | None = None
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

    def _press_hits_exclusion(self, root: QWidget, pos: QPoint) -> bool:
        for widget in self._drag_exclusions:
            if widget is None or not widget.isVisible():
                continue
            try:
                global_pos = root.mapToGlobal(pos)
                local_pos = widget.mapFromGlobal(global_pos)
            except Exception:
                continue
            if widget.rect().contains(local_pos):
                return True
        return False

    def _set_press_pending_state(self, active: bool) -> None:
        owner = getattr(self, "_on_background_press_state_changed", None)
        if callable(owner):
            owner(bool(active))

    def _begin_drag(
        self,
        child: QWidget | None,
        event: QMouseEvent,
        *,
        root: QWidget | None = None,
    ) -> bool:
        if event.button() == Qt.LeftButton:
            if self._is_in_exclusions(child):
                return False
            try:
                position = event.position()  # Qt6
            except AttributeError:
                position = event.pos()
            local_pos = position.toPoint() if hasattr(position, "toPoint") else position
            root_widget = self if root is None else root
            if self._press_hits_exclusion(root_widget, local_pos):
                return False
            try:
                global_pos = event.globalPosition().toPoint()  # Qt6
            except AttributeError:
                global_pos = event.globalPos()  # Qt5-style fallback
            self._drag_press_pending = True
            self._drag_press_pos = global_pos
            self._drag_press_child = child
            self._drag_pos = global_pos - self.frameGeometry().topLeft()
            self._set_press_pending_state(True)
            event.accept()
            return True
        return False

    def _start_drag(self, event: QMouseEvent) -> bool:
        if not self._drag_press_pending or self._drag_active:
            return False
        self._drag_press_pending = False
        self._drag_press_child = None
        self._set_press_pending_state(False)
        wh = self.windowHandle()
        if wh and wh.startSystemMove():
            self._drag_active = True
        else:
            self._drag_active = True
        event.accept()
        return True

    def _maybe_start_drag(self, event: QMouseEvent) -> bool:
        if not self._drag_press_pending or self._drag_active:
            return False
        try:
            global_pos = event.globalPosition().toPoint()  # Qt6
        except AttributeError:
            global_pos = event.globalPos()  # Qt5-style fallback
        delta = global_pos - self._drag_press_pos
        if int(delta.manhattanLength()) < int(QApplication.startDragDistance()):
            return False
        return self._start_drag(event)

    def _update_drag(self, event: QMouseEvent) -> bool:
        if self._maybe_start_drag(event):
            return True
        if self._drag_active and (event.buttons() & Qt.LeftButton):
            try:
                global_pos = event.globalPosition().toPoint()  # Qt6
            except AttributeError:
                global_pos = event.globalPos()  # Qt5-style fallback
            new_pos = global_pos - self._drag_pos
            if new_pos != self.frameGeometry().topLeft():
                self.move(new_pos)
            event.accept()
            return True
        return False

    def _end_drag(self, event: QMouseEvent) -> bool:
        self._drag_active = False
        self._drag_press_pending = False
        self._drag_press_child = None
        self._set_press_pending_state(False)
        event.accept()
        return True

    def _end_pending_press(self, event: QMouseEvent) -> bool:
        if not self._drag_press_pending:
            return False
        self._drag_press_pending = False
        self._drag_press_child = None
        self._set_press_pending_state(False)
        event.accept()
        return True

    # ---- Qt events --------------------------------------------------
    def eventFilter(self, watched: QObject, event: QEvent) -> bool:
        if isinstance(watched, QWidget) and watched in self._drag_targets:
            if event.type() == QEvent.Type.MouseButtonPress and isinstance(event, QMouseEvent):
                child = watched.childAt(event.position().toPoint())
                if self._begin_drag(child, event, root=watched):
                    return True
            elif event.type() == QEvent.Type.MouseMove and isinstance(event, QMouseEvent):
                if self._update_drag(event):
                    return True
            elif event.type() == QEvent.Type.MouseButtonRelease and isinstance(event, QMouseEvent):
                if self._drag_active:
                    return self._end_drag(event)
                if self._end_pending_press(event):
                    return True
        return super().eventFilter(watched, event)

    def mousePressEvent(self, event: QMouseEvent) -> None:
        child = self.childAt(event.pos())
        if self._begin_drag(child, event, root=self):
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
        if self._end_pending_press(event):
            return
        super().mouseReleaseEvent(event)
