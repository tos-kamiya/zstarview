from __future__ import annotations

from PySide6.QtCore import QPointF, QRect
from PySide6.QtGui import QColor, QPainter, QPolygonF
from PySide6.QtCore import Qt

from ..paths import ThemeStyle
from .background import atlas_background_tint_rgba

ATLAS_TIME_OF_DAY_MARKER_SIZE_PX = 31.0
ATLAS_TIME_OF_DAY_MARKER_MARGIN_PX = 4.0


def draw_instrument_background(
    painter: QPainter,
    viewport_rect: QRect,
    *,
    theme: ThemeStyle,
) -> None:
    """Draw the flat white background used by Atlas."""
    painter.fillRect(viewport_rect, QColor(*theme.window_background.inner_rgba))


def draw_instrument_time_of_day_marker(
    painter: QPainter,
    viewport_rect: QRect,
    *,
    sun_alt_deg: float | None,
    bottom_left: bool,
    tint_rgba: tuple[int, int, int, int] | None = None,
) -> None:
    """Draw a time-of-day marker in the HUD-free left corner."""
    if tint_rgba is None:
        tint_rgba = atlas_background_tint_rgba(sun_alt_deg)
    if tint_rgba is None:
        return
    left = float(viewport_rect.left()) + ATLAS_TIME_OF_DAY_MARKER_MARGIN_PX
    right = left + ATLAS_TIME_OF_DAY_MARKER_SIZE_PX
    top = float(viewport_rect.top()) + ATLAS_TIME_OF_DAY_MARKER_MARGIN_PX
    bottom = float(viewport_rect.bottom()) - ATLAS_TIME_OF_DAY_MARKER_MARGIN_PX
    if bottom_left:
        points = (
            QPointF(left, bottom),
            QPointF(right, bottom),
            QPointF(left, bottom - ATLAS_TIME_OF_DAY_MARKER_SIZE_PX),
        )
    else:
        points = (
            QPointF(left, top),
            QPointF(right, top),
            QPointF(left, top + ATLAS_TIME_OF_DAY_MARKER_SIZE_PX),
        )
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QColor(*tint_rgba))
    painter.drawPolygon(QPolygonF(points))
    painter.restore()
