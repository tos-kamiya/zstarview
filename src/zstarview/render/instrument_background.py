from __future__ import annotations

from PySide6.QtCore import QRect
from PySide6.QtGui import QColor, QImage, QPainter, QRadialGradient

from ..paths import ThemeStyle


_PAPER_TEXTURE_CACHE: dict[tuple[int, int, tuple[int, int, int]], QImage] = {}


def _paper_texture_image(width: int, height: int, base_rgb: tuple[int, int, int]) -> QImage:
    key = (int(width), int(height), tuple(int(value) for value in base_rgb))
    cached = _PAPER_TEXTURE_CACHE.get(key)
    if cached is not None:
        return cached

    image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(QColor(*base_rgb, 255))
    texture_painter = QPainter(image)
    texture_painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

    # A low-contrast radial edge darkening keeps the map-like paper visible
    # without introducing a scenic sky gradient.
    center = (float(width) * 0.5, float(height) * 0.47)
    radius = max(float(width), float(height)) * 0.76
    edge = QRadialGradient(center[0], center[1], radius)
    edge.setColorAt(0.0, QColor(255, 255, 255, 0))
    edge.setColorAt(0.72, QColor(92, 62, 32, 0))
    edge.setColorAt(1.0, QColor(92, 62, 32, 34))
    texture_painter.fillRect(image.rect(), edge)

    # Deterministic, sparse mottling avoids a repeated bitmap pattern and is
    # intentionally subtle at normal viewing size.
    state = (width * 73856093) ^ (height * 19349663)
    for index in range(72):
        state = (state * 1664525 + 1013904223 + index) & 0xFFFFFFFF
        x = float(state % max(1, width))
        state = (state * 1664525 + 1013904223) & 0xFFFFFFFF
        y = float(state % max(1, height))
        state = (state * 1664525 + 1013904223) & 0xFFFFFFFF
        radius_px = 8.0 + float(state % 26)
        alpha = 4 + int(state % 7)
        texture_painter.setBrush(QColor(111, 77, 39, alpha))
        texture_painter.setPen(QColor(0, 0, 0, 0))
        texture_painter.drawEllipse(
            int(x - radius_px),
            int(y - radius_px * 0.65),
            int(radius_px * 2.0),
            int(radius_px * 1.3),
        )
    texture_painter.end()

    if len(_PAPER_TEXTURE_CACHE) >= 8:
        _PAPER_TEXTURE_CACHE.pop(next(iter(_PAPER_TEXTURE_CACHE)))
    _PAPER_TEXTURE_CACHE[key] = image
    return image


def draw_instrument_background(
    painter: QPainter,
    viewport_rect: QRect,
    *,
    theme: ThemeStyle,
) -> None:
    """Draw the paper background used by Atlas."""
    if not theme.paper_texture:
        painter.fillRect(viewport_rect, QColor(*theme.window_background.inner_rgba))
        return
    image = _paper_texture_image(
        max(1, viewport_rect.width()),
        max(1, viewport_rect.height()),
        theme.window_background.base_rgb,
    )
    painter.drawImage(viewport_rect.topLeft(), image)
