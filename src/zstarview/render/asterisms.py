from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Tuple

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF

from ..asterisms import ASTERISMS, pick_rotating_asterism
from ..astro import altaz_to_normalized_xy, resolve_star_source_ids
from ..paths import ThemeStyle
from ..types import CelestialData, CelestialObject, ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy
from .guides import _clip_polyline_to_radius, _great_circle_altaz_points, split_by_gaps
from .stars import _content_fov_deg_from_viewer
from .text import (
    _text_bounds_at_baseline,
    draw_outlined_text,
    recolor_text_style,
    resolve_text_style,
)

ASTERISM_BASE_MID_WIDTH = 2.0
ASTERISM_HIGHLIGHT_MID_WIDTH = 2.2
ASTERISM_HIGHLIGHT_CORE_WIDTH = 1.0


def draw_asterisms(
    painter: QPainter,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    celestial_data: CelestialData,
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
    label_reservations: Optional[List[QRectF]] = None,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    *,
    theme: ThemeStyle,
    line_width_scale: float = 1.0,
    base_line_width_scale: float = 1.0,
    base_line_alpha_scale: float = 1.0,
    content_fov_deg: float | None = None,
    draw_base: bool = True,
    draw_highlight: bool = True,
) -> None:
    """Draw dim asterisms always, and brighten the hovered selection with a label."""

    stars = celestial_data.stars
    source_ids = resolve_star_source_ids(stars, celestial_data.star_catalog_meta)
    if source_ids.size == 0:
        return

    star_altaz_by_source: Dict[str, Tuple[float, float]] = {}
    for idx, raw_source in enumerate(source_ids):
        source_id = str(raw_source).strip()
        if not source_id:
            continue
        if source_id in star_altaz_by_source:
            continue
        star_altaz_by_source[source_id] = (float(stars["alt"][idx]), float(stars["az"][idx]))
    if not star_altaz_by_source:
        return

    asterism_style = theme.overlays.asterism
    is_bright_theme = theme.label_outline_suppressed
    asterism_rgb = asterism_style.rgb
    asterism_label_rgb = asterism_style.label_rgb or asterism_rgb
    highlight_mid_color = QColor(*asterism_rgb, 150 if is_bright_theme else 120)

    painter.save()
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)
    clip_radius = effective_fov_deg / max(1.0e-6, float(viewer_data.edge_fov_deg))
    width_scale = max(0.5, float(line_width_scale) * asterism_style.width_scale)
    base_width_scale = max(0.5, float(base_line_width_scale) * asterism_style.width_scale)
    base_alpha_scale = max(0.0, float(base_line_alpha_scale))

    def _make_pen(color: QColor, width: float) -> QPen:
        pen = QPen(color, width)
        pen.setCosmetic(True)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        return pen

    def _draw_segments(
        segments: Iterable[Tuple[str, str]],
        pens: Iterable[QPen],
    ) -> List[QPointF]:
        label_points: List[QPointF] = []
        for source_a, source_b in segments:
            pos_a = star_altaz_by_source.get(source_a)
            pos_b = star_altaz_by_source.get(source_b)
            if pos_a is None or pos_b is None:
                continue
            arc_altaz = _great_circle_altaz_points(pos_a[0], pos_a[1], pos_b[0], pos_b[1])
            arc_points: List[Tuple[float, float]] = []
            for alt_i, az_i in arc_altaz:
                nx_i, ny_i = altaz_to_normalized_xy(
                    alt_i,
                    az_i,
                    viewer_data.view_center,
                    edge_fov_deg=float(viewer_data.edge_fov_deg),
                )
                arc_points.append((nx_i, ny_i))
            for raw_frag in split_by_gaps(arc_points):
                clipped_frags = _clip_polyline_to_radius(raw_frag, clip_radius)
                for frag in clipped_frags:
                    if len(frag) < 2:
                        continue
                    poly = QPolygonF([QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag])
                    for pen in pens:
                        painter.setPen(pen)
                        painter.drawPolyline(poly)
                    label_points.extend(poly)
        return label_points

    def _base_pass() -> QPen:
        base_alpha = (
            asterism_style.line_alpha
            if asterism_style.line_alpha is not None
            else (32 if is_bright_theme else 24)
        )
        core_color = QColor(*asterism_rgb, min(255, int(round(base_alpha * base_alpha_scale))))
        core_width_scale = base_width_scale if base_width_scale > 1.0 else width_scale
        return _make_pen(core_color, ASTERISM_BASE_MID_WIDTH * core_width_scale)

    highlighted_asterism = None
    if draw_highlight and highlighted_object is not None:
        hovered_obj, _ = highlighted_object
        if isinstance(hovered_obj, dict):
            hovered_source_id = str(hovered_obj.get("source_id", "")).strip()
            if hovered_source_id:
                second_slot = int(datetime.now().timestamp()) // 3
                highlighted_asterism = pick_rotating_asterism(hovered_source_id, second_slot)

    label_points: List[QPointF] = []
    if draw_base:
        base_segments: set[Tuple[str, str]] = set()
        for asterism in ASTERISMS:
            for source_a, source_b in asterism.segments():
                if source_a == source_b:
                    continue
                base_segments.add(tuple(sorted((source_a, source_b))))
        _draw_segments(sorted(base_segments), (_base_pass(),))

    if highlighted_asterism is not None:
        highlight_pen = _make_pen(highlight_mid_color, ASTERISM_HIGHLIGHT_MID_WIDTH * width_scale)
        label_points = _draw_segments(
            highlighted_asterism.segments(),
            (highlight_pen,),
        )

    if label_points:
        cx = sum(pt.x() for pt in label_points) / len(label_points)
        cy = sum(pt.y() for pt in label_points) / len(label_points)
        label_pos = QPointF(cx + 8.0, cy - 8.0)
        text_style = recolor_text_style(resolve_text_style(theme, text_font), asterism_label_rgb)
        if label_candidates is not None:
            label_candidates.append(
                {
                    "text": highlighted_asterism.name,
                    "pos": label_pos,
                    "style": text_style,
                    "priority": 30,
                }
            )
        else:
            draw_outlined_text(
                painter,
                highlighted_asterism.name,
                label_pos,
                text_font,
                style=text_style,
            )
            if label_reservations is not None:
                label_reservations.append(_text_bounds_at_baseline(highlighted_asterism.name, text_font, label_pos))

    painter.restore()
