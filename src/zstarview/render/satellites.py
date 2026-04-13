from typing import Any, Dict, List, Optional

from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor, QPainter

from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import FIELD_OF_VIEW_DEG
from ..satellite_constants import SATELLITE_OVERLAY_MARKER_COLOR_RGB, SATELLITE_OVERLAY_MARKER_MAX_ALPHA
from ..satellites.types import SatelliteOverlayPoint
from ..types import ScreenGeometry
from .geometry import normalized_to_screen_xy
from .guides import draw_gauge_cross
from .text import blend_color_toward_white, resolve_text_style


def draw_satellite_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    satellite_points: list[SatelliteOverlayPoint] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 1.0,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    preset: str = "night",
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> None:
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if not satellite_points or layer_opacity <= 0.0:
        return

    painter.save()
    marker_color = QColor(
        *SATELLITE_OVERLAY_MARKER_COLOR_RGB,
        max(0, min(255, int(round(SATELLITE_OVERLAY_MARKER_MAX_ALPHA * layer_opacity)))),
    )
    label_style = resolve_text_style(preset, painter.font(), opacity=layer_opacity)
    label_color = blend_color_toward_white(QColor(*SATELLITE_OVERLAY_MARKER_COLOR_RGB), 0.1)
    label_color.setAlpha(label_style.text_color.alpha())
    label_style = type(label_style)(
        font=label_style.font,
        text_color=label_color,
        outline_color=label_style.outline_color,
        outline_width=label_style.outline_width,
    )
    for point in satellite_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        if not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pos = QPointF(float(px), float(py))
        draw_gauge_cross(
            painter,
            marker_color,
            pos,
            scale=float(point.marker_scale),
            pen_width=1.0,
        )
        if bool(point.show_label):
            label_text = str(point.satellite_name).strip()
            if label_text and label_candidates is not None:
                label_candidates.append(
                    {
                        "text": label_text,
                        "pos": QPointF(pos.x() + 10.0, pos.y() - 8.0),
                        "style": label_style,
                        "priority": 42,
                        "hide_on_overlap": True,
                    }
                )
    painter.restore()
