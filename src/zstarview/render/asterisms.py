from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Tuple

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy, resolve_star_source_ids
from ..asterisms import ASTERISMS, pick_rotating_asterism
from ..types import CelestialData, CelestialObject, ScreenGeometry, ViewerData
from .stars import _content_fov_deg_from_viewer
from .geometry import normalized_to_screen_xy
from .guides import _clip_polyline_to_radius, _great_circle_altaz_points, split_by_gaps
from .text import _text_bounds_at_baseline, draw_outlined_text, get_text_outline_width, get_text_style


def draw_asterisms(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
    label_reservations: Optional[List[QRectF]] = None,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    *,
    preset: str = "night",
    line_width_scale: float = 1.0,
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

    if preset in ("white", "day"):
        base_line_color = QColor(26, 114, 214, 25)
        base_outline_color = QColor(190, 220, 250, 12)
        highlight_line_color = QColor(26, 114, 214, 124)
        highlight_outline_color = QColor(190, 220, 250, 52)
    else:
        base_line_color = QColor(82, 142, 214, 21)
        base_outline_color = QColor(24, 48, 86, 9)
        highlight_line_color = QColor(120, 190, 255, 92)
        highlight_outline_color = QColor(32, 76, 130, 44)

    painter.save()
    effective_fov_deg = _content_fov_deg_from_viewer(viewer_data) if content_fov_deg is None else float(content_fov_deg)
    clip_radius = effective_fov_deg / 90.0
    width_scale = max(1.0, float(line_width_scale))

    def _make_pen(color: QColor, width: float) -> QPen:
        pen = QPen(color, width)
        pen.setCosmetic(True)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        return pen

    def _draw_segments(
        segments: Iterable[Tuple[str, str]],
        outline_pen: QPen,
        line_pen: QPen,
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
                nx_i, ny_i = altaz_to_normalized_xy(alt_i, az_i, viewer_data.view_center)
                arc_points.append((nx_i, ny_i))
            for raw_frag in split_by_gaps(arc_points):
                clipped_frags = _clip_polyline_to_radius(raw_frag, clip_radius)
                for frag in clipped_frags:
                    if len(frag) < 2:
                        continue
                    poly = QPolygonF([QPointF(*normalized_to_screen_xy(nx, ny, geometry)) for nx, ny in frag])
                    painter.setPen(outline_pen)
                    painter.drawPolyline(poly)
                    painter.setPen(line_pen)
                    painter.drawPolyline(poly)
                    label_points.extend(poly)
        return label_points

    def _draw_one_asterism(asterism: Any, outline_pen: QPen, line_pen: QPen) -> List[QPointF]:
        return _draw_segments(asterism.segments(), outline_pen, line_pen)

    if draw_base:
        base_outline_pen = _make_pen(base_outline_color, 4.0 * width_scale)
        base_line_pen = _make_pen(base_line_color, 2.5 * width_scale)
        base_segments: set[Tuple[str, str]] = set()
        for asterism in ASTERISMS:
            for source_a, source_b in asterism.segments():
                if source_a == source_b:
                    continue
                base_segments.add(tuple(sorted((source_a, source_b))))
        _draw_segments(sorted(base_segments), base_outline_pen, base_line_pen)

    highlighted_asterism = None
    if draw_highlight and highlighted_object is not None:
        hovered_obj, _ = highlighted_object
        if isinstance(hovered_obj, dict):
            hovered_source_id = str(hovered_obj.get("source_id", "")).strip()
            if hovered_source_id:
                second_slot = int(datetime.now().timestamp()) // 3
                highlighted_asterism = pick_rotating_asterism(hovered_source_id, second_slot)

    label_points: List[QPointF] = []
    if highlighted_asterism is not None:
        highlight_outline_pen = _make_pen(highlight_outline_color, 3.2 * width_scale)
        highlight_line_pen = _make_pen(highlight_line_color, 2.0 * width_scale)
        label_points = _draw_one_asterism(highlighted_asterism, highlight_outline_pen, highlight_line_pen)

    if label_points:
        cx = sum(pt.x() for pt in label_points) / len(label_points)
        cy = sum(pt.y() for pt in label_points) / len(label_points)
        label_pos = QPointF(cx + 8.0, cy - 8.0)
        text_color, outline_text_color = get_text_style(preset)
        outline_width = get_text_outline_width(preset)
        if label_candidates is not None:
            label_candidates.append(
                {
                    "text": highlighted_asterism.name,
                    "pos": label_pos,
                    "text_color": text_color,
                    "outline_color": outline_text_color,
                    "outline_width": outline_width,
                    "priority": 30,
                }
            )
        else:
            draw_outlined_text(
                painter,
                highlighted_asterism.name,
                label_pos,
                text_font,
                text_color,
                outline_text_color,
                outline_width=outline_width,
            )
            if label_reservations is not None:
                label_reservations.append(_text_bounds_at_baseline(highlighted_asterism.name, text_font, label_pos))

    painter.restore()
