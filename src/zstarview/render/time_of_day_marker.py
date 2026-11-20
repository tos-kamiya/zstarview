"""Time-of-day marker drawn in the HUD-free corner of the viewport."""

from __future__ import annotations

from PySide6.QtCore import QPointF, QRect, Qt
from PySide6.QtGui import QColor, QPainter, QPolygonF

TIME_OF_DAY_MARKER_SIZE_PX = 31.0
TIME_OF_DAY_MARKER_MARGIN_PX = 4.0
TIME_OF_DAY_MARKER_OPACITY = 0.7


def draw_time_of_day_marker(
    painter: QPainter,
    viewport_rect: QRect,
    *,
    bottom_left: bool,
    tint_rgba: tuple[int, int, int, int],
) -> None:
    """Draw a triangular time-of-day marker in the HUD-free left corner."""
    left = float(viewport_rect.left()) + TIME_OF_DAY_MARKER_MARGIN_PX
    right = left + TIME_OF_DAY_MARKER_SIZE_PX
    top = float(viewport_rect.top()) + TIME_OF_DAY_MARKER_MARGIN_PX
    bottom = float(viewport_rect.bottom()) - TIME_OF_DAY_MARKER_MARGIN_PX
    if bottom_left:
        points = (
            QPointF(left, bottom),
            QPointF(right, bottom),
            QPointF(left, bottom - TIME_OF_DAY_MARKER_SIZE_PX),
        )
    else:
        points = (
            QPointF(left, top),
            QPointF(right, top),
            QPointF(left, top + TIME_OF_DAY_MARKER_SIZE_PX),
        )
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    color = QColor(*tint_rgba)
    color.setAlpha(round(color.alpha() * TIME_OF_DAY_MARKER_OPACITY))
    painter.setBrush(color)
    painter.drawPolygon(QPolygonF(points))
    painter.restore()
