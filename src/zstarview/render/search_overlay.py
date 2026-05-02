from __future__ import annotations

from typing import Any, Dict, List, Optional

from PySide6.QtCore import QPointF, QRectF
from PySide6.QtGui import QPainter

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import ThemeStyle
from ..search.models import SearchJumpTarget
from . import guides as render_guides
from . import text as render_text
from .geometry import normalized_to_screen_xy


def _project_altaz_to_normalized_xy(
    alt_deg: float,
    az_deg: float,
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
) -> tuple[float, float]:
    return altaz_to_normalized_xy(
        alt_deg,
        az_deg,
        view_center,
        edge_fov_deg=edge_fov_deg,
    )


def draw_search_target_overlay(
    painter: QPainter,
    geometry,
    target: SearchJumpTarget,
    *,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    text_font,
    draw_marker: bool = True,
    draw_label: bool = True,
    marker_scale: float = 1.0,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    label_reservations: Optional[List[QRectF]] = None,
    theme: ThemeStyle,
) -> None:
    alt = getattr(target, "alt_deg", None)
    az = getattr(target, "az_deg", None)
    if alt is None or az is None:
        return
    if not is_in_fov(float(alt), float(az), view_center, fov_deg=float(content_fov_deg)):
        return

    nx, ny = _project_altaz_to_normalized_xy(
        float(alt),
        float(az),
        view_center,
        edge_fov_deg=edge_fov_deg,
    )
    px, py = normalized_to_screen_xy(nx, ny, geometry)
    pos = QPointF(px, py)
    color, _outline = render_text.get_text_style(theme)

    if draw_marker:
        render_guides.draw_gauge_cross(
            painter,
            color,
            pos,
            scale=0.42 * max(1.0, float(marker_scale)),
            pen_width=2.0,
        )

    if not draw_label:
        return

    label = str(getattr(target, "label", "")).strip()
    if not label:
        return
    label_pos = QPointF(pos.x() + 12.0, pos.y() - 10.0)
    style = render_text.resolve_label_text_style(theme, text_font)
    if label_candidates is not None:
        label_candidates.append(
            {
                "text": label,
                "pos": label_pos,
                "style": style,
                "priority": 15,
            }
        )
        return
    label_pos = render_text._clamp_baseline_pos_to_viewport(
        label,
        text_font,
        label_pos,
        QRectF(painter.viewport()),
    )
    if label_reservations is not None:
        label_reservations.append(
            render_text._text_bounds_at_baseline(label, text_font, label_pos)
        )
    render_text.draw_outlined_text(
        painter,
        label,
        label_pos,
        style=style,
    )
