from PySide6.QtCore import Qt, QPoint
from PySide6.QtWidgets import QMainWindow, QWidget


class DraggableWindow(QMainWindow):
    """
    Frameless window base with click-drag move.
    - Uses startSystemMove() when available (Wayland etc.).
    - Falls back to manual dragging otherwise.
    - Supports registering child widgets that disable dragging when clicked.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._drag_active = False
        self._drag_pos = QPoint(0, 0)
        self._drag_exclusions: set[QWidget] = set()

    # ---- public API -------------------------------------------------
    def add_drag_exclusion(self, widget: QWidget) -> None:
        """Exclude a child widget from initiating window drag."""
        if widget is not None:
            self._drag_exclusions.add(widget)

    def add_drag_exclusions(self, widgets) -> None:
        for w in widgets:
            self.add_drag_exclusion(w)

    # ---- internals --------------------------------------------------
    def _is_in_exclusions(self, w: QWidget | None) -> bool:
        """Return True if w or any ancestor is in the exclusion set."""
        while w is not None:
            if w in self._drag_exclusions:
                return True
            w = w.parentWidget()
        return False

    # ---- Qt events --------------------------------------------------
    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            child = self.childAt(event.pos())
            if self._is_in_exclusions(child):
                # Let the child handle the event
                super().mousePressEvent(event)
                return

            # Prefer native system move when available (Wayland, etc.)
            wh = None
            try:
                wh = self.windowHandle()
            except Exception:
                wh = None

            if wh is not None:
                try:
                    if bool(wh.startSystemMove()):
                        event.accept()
                        return
                except Exception:
                    # Fallback to manual drag
                    pass

            # Manual dragging fallback
            self._drag_active = True
            try:
                global_pos = event.globalPosition().toPoint()  # Qt6
            except AttributeError:
                global_pos = event.globalPos()  # Qt5-style fallback
            self._drag_pos = global_pos - self.frameGeometry().topLeft()
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._drag_active and (event.buttons() & Qt.LeftButton):
            try:
                global_pos = event.globalPosition().toPoint()
            except AttributeError:
                global_pos = event.globalPos()
            self.move(global_pos - self._drag_pos)
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        self._drag_active = False
        event.accept()
