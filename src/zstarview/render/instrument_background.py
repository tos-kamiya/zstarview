from __future__ import annotations

from PySide6.QtCore import QRect
from PySide6.QtGui import QColor, QPainter

from ..paths import ThemeStyle


def draw_instrument_background(
    painter: QPainter,
    viewport_rect: QRect,
    *,
    theme: ThemeStyle,
) -> None:
    """Draw the solid background used by Atlas."""
    painter.fillRect(viewport_rect, QColor(*theme.window_background.inner_rgba))
