import math
from typing import Any, Callable, Dict, List, Optional, Tuple

from PySide6.QtCore import QPoint, QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen

from ..paths import TEXT_STYLES_BY_PRESET
from ..types import CelestialData, CelestialObject, PlanetBody, ScreenGeometry, ViewerData
from .background import format_overlay_info_lines
from .deep_sky_objects import _DSO_HOVER_SIZE_GAIN, _dso_ellipse_polygon


def draw_overlay_info(
    painter: QPainter,
    geometry: ScreenGeometry,
    celestial_data: CelestialData,
    viewer_data: ViewerData,
    vmag_limit: float,
    enlarge_moon: bool,
    highlighted_dso: Optional[Tuple[CelestialObject, QPointF]],
    highlighted_object: Optional[Tuple[CelestialObject, QPointF]],
    text_font: QFont,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    label_reservations: Optional[List[QRectF]] = None,
    *,
    mouse_pos: Optional[QPoint] = None,
    preset: str = "night",
    draw_static_info: bool = True,
    draw_hover_info: bool = True,
    get_text_style_func: Callable[[str], Tuple[QColor, QColor]],
    draw_outlined_text_func: Callable[..., None],
    text_bounds_at_baseline_func: Callable[[str, QFont, QPointF], QRectF],
) -> None:
    """
    Draw overlay text information on the screen.

    This includes the current time, location, view direction, magnitude limit,
    and information about any highlighted celestial object.
    """
    del enlarge_moon

    text_color, outline_color = get_text_style_func(preset)
    outline_width = 3.0

    font_metrics = QFontMetrics(text_font)
    line_spacing = font_metrics.lineSpacing()
    line_height = int(line_spacing * 1.2)
    line_x = line_spacing
    sample_bounds = font_metrics.tightBoundingRect("Ag")
    first_line_baseline_y = line_spacing - float(sample_bounds.top())
    line_y = first_line_baseline_y - line_height

    def print_line(message: str) -> None:
        nonlocal line_x, line_y
        line_y += line_height
        draw_outlined_text_func(
            painter,
            message,
            QPointF(line_x, line_y),
            text_font,
            text_color,
            outline_color,
            outline_width=outline_width,
        )

    static_lines = format_overlay_info_lines(celestial_data, viewer_data, vmag_limit)
    if draw_static_info and mouse_pos is not None:
        overlay_top_y = float(line_spacing)
        overlay_bottom_y = float(line_spacing + len(static_lines) * line_height)
        mouse_y = float(mouse_pos.y())
        if overlay_top_y <= mouse_y <= overlay_bottom_y:
            draw_static_info = False

    if draw_static_info:
        for line in static_lines:
            print_line(line)

    if draw_hover_info and highlighted_dso:
        dso_obj, _ = highlighted_dso
        major_arcmin = float(dso_obj.get("major_arcmin", 0.0))
        minor_arcmin = float(dso_obj.get("minor_arcmin", 0.0))
        pa_deg = float(dso_obj.get("pa_deg", 0.0))
        alt = float(dso_obj.get("alt", 0.0))
        az = float(dso_obj.get("az", 0.0))
        if preset in ("white", "day"):
            dso_label_color = QColor(*TEXT_STYLES_BY_PRESET["white"].text, 228)
        else:
            dso_label_color = QColor(110, 195, 255, 230)
        hover_poly = _dso_ellipse_polygon(
            alt_deg=alt,
            az_deg=az,
            major_arcmin=major_arcmin,
            minor_arcmin=minor_arcmin,
            pa_deg=pa_deg if math.isfinite(pa_deg) else 0.0,
            viewer_data=viewer_data,
            geometry=geometry,
            gain=_DSO_HOVER_SIZE_GAIN,
            samples=60,
        )
        dso_name = str(dso_obj.get("name", "")).strip()
        if dso_name:
            bounds = hover_poly.boundingRect()
            label_pos = QPointF(bounds.right() + 10.0, bounds.top() - 6.0)
            if label_candidates is not None:
                label_candidates.append(
                    {
                        "text": dso_name,
                        "pos": label_pos,
                        "text_color": dso_label_color,
                        "outline_color": outline_color,
                        "outline_width": outline_width,
                        "priority": 20,
                    }
                )
            else:
                draw_outlined_text_func(
                    painter,
                    dso_name,
                    label_pos,
                    text_font,
                    dso_label_color,
                    outline_color,
                    outline_width=outline_width,
                )
                if label_reservations is not None:
                    label_reservations.append(text_bounds_at_baseline_func(dso_name, text_font, label_pos))

    if draw_hover_info and highlighted_object:
        obj, pos = highlighted_object
        painter.setPen(QPen(text_color, 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(pos, 10, 10)

        if not isinstance(obj, PlanetBody):
            name = getattr(obj, "name", "") if hasattr(obj, "name") else obj.get("name", "")
            label_pos = QPointF(pos.x() + 15, pos.y() - 15)
            if label_candidates is not None:
                label_candidates.append(
                    {
                        "text": str(name),
                        "pos": label_pos,
                        "text_color": text_color,
                        "outline_color": outline_color,
                        "outline_width": outline_width,
                        "priority": 10,
                    }
                )
            else:
                draw_outlined_text_func(
                    painter,
                    name,
                    label_pos,
                    text_font,
                    text_color,
                    outline_color,
                    outline_width=outline_width,
                )
                if label_reservations is not None:
                    label_reservations.append(text_bounds_at_baseline_func(str(name), text_font, label_pos))
