from typing import Any, Dict, List, Optional, Tuple

from PySide6.QtCore import QPointF, QRect, QRectF
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPainterPath, QPen

from ..paths import STATUS_LINE_STYLES_BY_PRESET, TEXT_STYLES_BY_PRESET


def _qcolor_from_rgba(color: tuple[int, ...]) -> QColor:
    if len(color) == 3:
        return QColor(*color)
    return QColor(*color[:3], color[3])


def _text_bounds_at_baseline(text: str, font: QFont, baseline_pos: QPointF) -> QRectF:
    fm = QFontMetrics(font)
    bounds = fm.tightBoundingRect(text)
    return QRectF(
        baseline_pos.x() + bounds.x(),
        baseline_pos.y() + bounds.y(),
        bounds.width(),
        bounds.height(),
    )


def _rect_overlap_count(rect: QRectF, others: List[QRectF], pad_px: float = 2.0) -> int:
    if not others:
        return 0
    test = rect.adjusted(-pad_px, -pad_px, pad_px, pad_px)
    return sum(
        1
        for other in others
        if test.intersects(other.adjusted(-pad_px, -pad_px, pad_px, pad_px))
    )


def _clamp_baseline_pos_to_viewport(
    text: str,
    font: QFont,
    baseline_pos: QPointF,
    viewport: QRectF,
    *,
    margin_px: float = 2.0,
) -> QPointF:
    """Shift baseline position so text bounds stay inside viewport."""
    rect = _text_bounds_at_baseline(text, font, baseline_pos)
    dx = 0.0
    dy = 0.0
    left = viewport.left() + margin_px
    right = viewport.right() - margin_px
    top = viewport.top() + margin_px
    bottom = viewport.bottom() - margin_px

    if rect.left() < left:
        dx += left - rect.left()
    elif rect.right() > right:
        dx -= rect.right() - right

    if rect.top() < top:
        dy += top - rect.top()
    elif rect.bottom() > bottom:
        dy -= rect.bottom() - bottom

    return QPointF(baseline_pos.x() + dx, baseline_pos.y() + dy)


def get_text_style(preset: str = "night", *, status_line: bool = False) -> Tuple[QColor, QColor]:
    """Return (text_color, outline_color) for the selected preset and text role."""
    style_table = STATUS_LINE_STYLES_BY_PRESET if status_line else TEXT_STYLES_BY_PRESET
    style = style_table.get(preset, style_table["night"])
    return _qcolor_from_rgba(style.text), _qcolor_from_rgba(style.outline)


def draw_outlined_text(
    painter: QPainter,
    text: str,
    pos: QPointF,
    font: QFont,
    text_color: QColor = QColor(255, 255, 255),
    outline_color: QColor = QColor.fromRgbF(0, 0, 0, 0.3),
    outline_width: float = 3.0,
) -> None:
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)

    path = QPainterPath()
    path.addText(pos, font, text)

    pen = QPen(outline_color, outline_width)
    painter.setPen(pen)
    painter.drawPath(path)
    painter.fillPath(path, text_color)
    painter.restore()


def _draw_label_candidates(
    painter: QPainter,
    candidates: List[Dict[str, Any]],
    text_font: QFont,
) -> None:
    """Draw label candidates as the final label layer with overlap avoidance."""
    if not candidates:
        return
    reservations: List[QRectF] = []
    offsets = (
        (0.0, 0.0),
        (12.0, 0.0),
        (-12.0, 0.0),
        (0.0, -12.0),
        (0.0, 12.0),
        (20.0, -10.0),
        (-20.0, -10.0),
        (24.0, 0.0),
        (-24.0, 0.0),
        (0.0, -24.0),
        (0.0, 24.0),
        (30.0, -14.0),
        (-30.0, -14.0),
        (30.0, 14.0),
        (-30.0, 14.0),
    )
    painter.save()
    painter.setFont(text_font)
    viewport = QRectF(painter.viewport())
    ordered = sorted(candidates, key=lambda c: int(c.get("priority", 999)))
    for cand in ordered:
        text = str(cand.get("text", "")).strip()
        if not text:
            continue
        anchor = cand.get("pos")
        if not isinstance(anchor, QPointF):
            continue
        text_color = cand.get("text_color")
        outline_color = cand.get("outline_color")
        if not isinstance(text_color, QColor) or not isinstance(outline_color, QColor):
            continue
        hide_on_overlap = bool(cand.get("hide_on_overlap", False))
        outline_width = float(cand.get("outline_width", 3.0))
        placed = False
        best_nonfree: Optional[Tuple[int, float, QPointF, QRectF]] = None
        for dx, dy in offsets:
            pos = QPointF(anchor.x() + dx, anchor.y() + dy)
            pos = _clamp_baseline_pos_to_viewport(text, text_font, pos, viewport)
            rect = _text_bounds_at_baseline(text, text_font, pos)
            overlap_count = _rect_overlap_count(rect, reservations)
            if overlap_count > 0:
                if not hide_on_overlap:
                    distance2 = (pos.x() - anchor.x()) ** 2 + (pos.y() - anchor.y()) ** 2
                    score = (overlap_count, distance2, pos, rect)
                    if best_nonfree is None or score[:2] < best_nonfree[:2]:
                        best_nonfree = score
                continue
            draw_outlined_text(
                painter,
                text,
                pos,
                text_font,
                text_color,
                outline_color,
                outline_width=outline_width,
            )
            reservations.append(rect)
            placed = True
            break
        if not placed and not hide_on_overlap and best_nonfree is not None:
            _, _, pos, rect = best_nonfree
            draw_outlined_text(
                painter,
                text,
                pos,
                text_font,
                text_color,
                outline_color,
                outline_width=outline_width,
            )
            reservations.append(rect)
        elif not placed:
            continue
    painter.restore()


def _draw_status_line_text(
    painter: QPainter,
    message: str,
    status_line_font: QFont,
    viewport_rect: QRect,
    *,
    preset: str = "night",
) -> None:
    """Draw a single-line status message at the bottom-left corner."""
    if not message:
        return

    color, outline_color = get_text_style(preset, status_line=True)

    painter.save()
    painter.setFont(status_line_font)
    fm = painter.fontMetrics()
    margin = fm.lineSpacing()
    baseline_y = viewport_rect.bottom() - margin // 4
    x = margin
    draw_outlined_text(
        painter,
        "> " + message,
        QPointF(x, baseline_y),
        status_line_font,
        color,
        outline_color,
        outline_width=3.0,
    )
    painter.restore()
