from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from PySide6.QtCore import QPointF, QRect, QRectF
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPainterPath, QPen

from ..paths import ThemeStyle


@dataclass(frozen=True, slots=True)
class ResolvedTextStyle:
    font: QFont
    text_color: QColor
    outline_color: QColor
    outline_width: float


def _qcolor_from_rgba(color: tuple[int, ...]) -> QColor:
    if len(color) == 3:
        return QColor(*color)
    return QColor(*color[:3], color[3])


def blend_color_toward_white(color: QColor, amount: float = 0.1) -> QColor:
    """Blend a color slightly toward white while preserving alpha."""
    t = max(0.0, min(1.0, float(amount)))
    return QColor(
        int(round(color.red() * (1.0 - t) + 255.0 * t)),
        int(round(color.green() * (1.0 - t) + 255.0 * t)),
        int(round(color.blue() * (1.0 - t) + 255.0 * t)),
        color.alpha(),
    )


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


def _rects_overlap(rect_a: QRectF, rect_b: QRectF, pad_px: float = 2.0) -> bool:
    return rect_a.adjusted(-pad_px, -pad_px, pad_px, pad_px).intersects(
        rect_b.adjusted(-pad_px, -pad_px, pad_px, pad_px)
    )


def _label_candidate_offsets() -> Tuple[Tuple[float, float], ...]:
    """Return an ordered search pattern for label placement."""
    offsets = [
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
    ]
    for radius in (42.0, 54.0, 66.0):
        offsets.extend(
            [
                (radius, 0.0),
                (-radius, 0.0),
                (0.0, -radius),
                (0.0, radius),
                (radius, -radius * 0.5),
                (-radius, -radius * 0.5),
                (radius, radius * 0.5),
                (-radius, radius * 0.5),
            ]
        )
    return tuple(offsets)


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


def _label_candidate_layout(
    text: str,
    font: QFont,
    anchor: QPointF,
    viewport: QRectF,
) -> tuple[QPointF, QRectF]:
    baseline_pos = _clamp_baseline_pos_to_viewport(text, font, anchor, viewport)
    rect = _text_bounds_at_baseline(text, font, baseline_pos)
    return baseline_pos, rect


def _cluster_label_candidate_groups(
    items: List[Dict[str, Any]],
    *,
    pad_px: float = 2.0,
) -> List[List[int]]:
    if len(items) <= 1:
        return [list(range(len(items)))] if items else []

    parent = list(range(len(items)))

    def find(index: int) -> int:
        root = index
        while parent[root] != root:
            root = parent[root]
        while parent[index] != index:
            parent[index], index = root, parent[index]
        return root

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for left_index in range(len(items)):
        left_rect = items[left_index]["base_rect"]
        if not isinstance(left_rect, QRectF):
            continue
        for right_index in range(left_index + 1, len(items)):
            right_rect = items[right_index]["base_rect"]
            if not isinstance(right_rect, QRectF):
                continue
            if _rects_overlap(left_rect, right_rect, pad_px=pad_px):
                union(left_index, right_index)

    grouped_indices: dict[int, List[int]] = {}
    for index in range(len(items)):
        root = find(index)
        grouped_indices.setdefault(root, []).append(index)

    ordered_groups: List[List[int]] = []
    seen_roots: set[int] = set()
    for index in range(len(items)):
        root = find(index)
        if root in seen_roots:
            continue
        seen_roots.add(root)
        ordered_groups.append(grouped_indices[root])
    return ordered_groups


def _order_label_candidate_group(
    items: List[Dict[str, Any]],
    group: List[int],
) -> List[int]:
    if len(group) <= 1:
        return list(group)

    if len(group) == 2:
        return sorted(
            group,
            key=lambda index: (
                float(items[index]["base_rect"].top()),
                int(items[index]["priority"]),
                float(items[index]["base_rect"].left()),
                index,
            ),
        )

    cx = sum(float(items[index]["base_rect"].center().x()) for index in group) / float(len(group))
    cy = sum(float(items[index]["base_rect"].center().y()) for index in group) / float(len(group))
    return sorted(
        group,
        key=lambda index: (
            (
                float(items[index]["base_rect"].center().x()) - cx
            )
            ** 2
            + (
                float(items[index]["base_rect"].center().y()) - cy
            )
            ** 2,
            float(items[index]["base_rect"].top()),
            float(items[index]["base_rect"].left()),
            int(items[index]["priority"]),
            index,
        ),
    )


def get_text_style(
    theme: ThemeStyle,
    *,
    status_line: bool = False,
) -> Tuple[QColor, QColor]:
    """Return (text_color, outline_color) for the selected theme and text role."""
    style = theme.status_text if status_line else theme.text
    return _qcolor_from_rgba(style.foreground_rgb), _qcolor_from_rgba(style.outline_rgba)


def get_text_outline_width(
    theme: ThemeStyle,
    *,
    status_line: bool = False,
) -> float:
    """Return the outline stroke width for the selected theme and text role."""
    style = theme.status_text if status_line else theme.text
    return float(style.outline_width)


def resolve_text_style(
    theme: ThemeStyle,
    font: QFont,
    *,
    status_line: bool = False,
    opacity: float = 1.0,
) -> ResolvedTextStyle:
    """Resolve a theme into a concrete text render style."""
    text_color, outline_color = get_text_style(theme, status_line=status_line)
    alpha_scale = max(0.0, min(1.0, float(opacity)))
    if alpha_scale < 1.0:
        text_color = QColor(text_color)
        outline_color = QColor(outline_color)
        text_color.setAlpha(max(0, min(255, int(round(text_color.alpha() * alpha_scale)))))
        outline_color.setAlpha(max(0, min(255, int(round(outline_color.alpha() * alpha_scale)))))
    return ResolvedTextStyle(
        font=font,
        text_color=text_color,
        outline_color=outline_color,
        outline_width=get_text_outline_width(theme, status_line=status_line),
    )


def resolve_label_text_style(
    theme: ThemeStyle,
    font: QFont,
    *,
    opacity: float = 1.0,
) -> ResolvedTextStyle:
    """Resolve a label style and suppress outlines in bright themes."""
    style = resolve_text_style(theme, font, opacity=opacity)
    if not theme.label_outline_suppressed:
        return style
    outline_color = QColor(style.outline_color)
    outline_color.setAlpha(0)
    return ResolvedTextStyle(
        font=style.font,
        text_color=style.text_color,
        outline_color=outline_color,
        outline_width=0.0,
    )


def recolor_text_style(
    style: ResolvedTextStyle,
    rgb: tuple[int, int, int],
) -> ResolvedTextStyle:
    """Return a copy of a resolved text style with a new text color."""
    text_color = QColor(*rgb)
    text_color.setAlpha(style.text_color.alpha())
    return ResolvedTextStyle(
        font=style.font,
        text_color=text_color,
        outline_color=QColor(style.outline_color),
        outline_width=style.outline_width,
    )


def draw_outlined_text(
    painter: QPainter,
    text: str,
    pos: QPointF,
    font: QFont | None = None,
    text_color: QColor = QColor(255, 255, 255),
    outline_color: QColor = QColor.fromRgbF(0, 0, 0, 0.3),
    outline_width: float = 3.0,
    *,
    style: ResolvedTextStyle | None = None,
) -> None:
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    if style is not None:
        font = style.font
        text_color = style.text_color
        outline_color = style.outline_color
        outline_width = style.outline_width
    if font is None:
        raise ValueError("draw_outlined_text requires either a font or a style")

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
    painter.save()
    painter.setFont(text_font)
    viewport = QRectF(painter.viewport())
    ordered = sorted(candidates, key=lambda c: int(c.get("priority", 999)))
    items: List[Dict[str, Any]] = []
    for order_index, cand in enumerate(ordered):
        text = str(cand.get("text", "")).strip()
        if not text:
            continue
        anchor = cand.get("pos")
        if not isinstance(anchor, QPointF):
            continue
        style = cand.get("style")
        if not isinstance(style, ResolvedTextStyle):
            text_color = cand.get("text_color")
            outline_color = cand.get("outline_color")
            if not isinstance(text_color, QColor) or not isinstance(outline_color, QColor):
                continue
            style = ResolvedTextStyle(
                font=text_font,
                text_color=text_color,
                outline_color=outline_color,
                outline_width=float(cand.get("outline_width", 3.0)),
            )
        base_pos, base_rect = _label_candidate_layout(text, text_font, anchor, viewport)
        items.append(
            {
                "cand": cand,
                "text": text,
                "anchor": anchor,
                "style": style,
                "base_pos": base_pos,
                "base_rect": base_rect,
                "priority": int(cand.get("priority", 999)),
                "hide_on_overlap": bool(cand.get("hide_on_overlap", False)),
                "order_index": order_index,
            }
        )

    reservations: List[QRectF] = []
    offsets = _label_candidate_offsets()
    placement_order: List[Dict[str, Any]] = []
    for group in _cluster_label_candidate_groups(items):
        for index in _order_label_candidate_group(items, group):
            placement_order.append(items[index])

    for item in placement_order:
        text = item["text"]
        anchor = item["anchor"]
        style = item["style"]
        hide_on_overlap = bool(item["hide_on_overlap"])
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
                style=style,
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
                style=style,
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
    theme: ThemeStyle,
) -> None:
    """Draw a single-line status message at the bottom-left corner."""
    if not message:
        return
    text_color, outline_color = get_text_style(theme, status_line=True)
    style = ResolvedTextStyle(
        font=status_line_font,
        text_color=text_color,
        outline_color=outline_color,
        outline_width=get_text_outline_width(theme, status_line=True),
    )

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
        style=style,
    )
    painter.restore()
