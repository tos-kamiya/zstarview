from typing import Any, Dict, List, Optional

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF

from ..aircraft.types import AircraftOverlayPoint
from ..aircraft_constants import AIRCRAFT_FADE_START_SECONDS, AIRCRAFT_OVERLAY_LINE_COLOR_RGB
from ..astro import altaz_to_normalized_xy, is_in_fov
from ..paths import FIELD_OF_VIEW_DEG
from ..types import ScreenGeometry
from .geometry import normalized_to_screen_xy
from .text import resolve_text_style

_AIRCRAFT_CALLSIGN_MAX_DISTANCE_KM = 10.0
_AIRCRAFT_MAX_DRAW_DISTANCE_KM = 50.0


def draw_aircraft_overlay(
    painter: QPainter,
    geometry: ScreenGeometry,
    aircraft_points: list[AircraftOverlayPoint] | None,
    view_center: tuple[float, float],
    *,
    opacity: float = 1.0,
    line_width_scale: float = 1.0,
    label_candidates: Optional[List[Dict[str, Any]]] = None,
    preset: str = "night",
    content_fov_deg: float = FIELD_OF_VIEW_DEG,
) -> None:
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if not aircraft_points or layer_opacity <= 0.0:
        return

    painter.save()
    width_scale = max(1.0, float(line_width_scale))
    line_color = QColor(*AIRCRAFT_OVERLAY_LINE_COLOR_RGB, 255)
    label_style = resolve_text_style(preset, painter.font(), opacity=layer_opacity)
    label_color = QColor(*AIRCRAFT_OVERLAY_LINE_COLOR_RGB)
    label_color.setAlpha(label_style.text_color.alpha())
    label_style = type(label_style)(
        font=label_style.font,
        text_color=label_color,
        outline_color=label_style.outline_color,
        outline_width=label_style.outline_width,
    )
    line_pen = QPen(line_color, 1.0, Qt.PenStyle.SolidLine)
    line_pen.setCosmetic(True)
    line_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    line_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(line_pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    for point in aircraft_points:
        alt = float(point.alt_deg)
        az = float(point.az_deg)
        distance_km = float(point.distance_km)
        if distance_km > _AIRCRAFT_MAX_DRAW_DISTANCE_KM:
            continue
        if alt <= 0.0 or not is_in_fov(alt, az, view_center, fov_deg=content_fov_deg):
            continue
        min_line_alpha = int(round(30.0 * layer_opacity))
        line_alpha = max(
            min_line_alpha,
            min(255, int(round(255.0 * float(point.alpha_scale) * layer_opacity))),
        )
        line_color.setAlpha(line_alpha)
        line_pen.setColor(line_color)
        line_pen.setWidthF(_aircraft_line_width_px(distance_km, width_scale=width_scale))
        painter.setPen(line_pen)
        trail_points = tuple(
            (float(sample_alt_deg), float(sample_az_deg))
            for sample_alt_deg, sample_az_deg in point.trail_alt_az_points
        )
        if any(
            is_in_fov(sample_alt_deg, sample_az_deg, view_center, fov_deg=content_fov_deg)
            for sample_alt_deg, sample_az_deg in trail_points
        ):
            polyline_points: list[QPointF] = []
            for sample_alt_deg, sample_az_deg in trail_points:
                sample_nx, sample_ny = altaz_to_normalized_xy(sample_alt_deg, sample_az_deg, view_center)
                sample_px, sample_py = normalized_to_screen_xy(sample_nx, sample_ny, geometry)
                polyline_points.append(QPointF(float(sample_px), float(sample_py)))
            if len(polyline_points) >= 2:
                painter.drawPolyline(QPolygonF(polyline_points))
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        px, py = normalized_to_screen_xy(nx, ny, geometry)
        pos = QPointF(float(px), float(py))
        callsign = (point.callsign or "").strip()
        if (
            callsign
            and float(point.distance_km) <= _AIRCRAFT_CALLSIGN_MAX_DISTANCE_KM
            and float(point.age_seconds) <= AIRCRAFT_FADE_START_SECONDS
            and label_candidates is not None
        ):
            label_candidates.append(
                {
                    "text": callsign,
                    "pos": QPointF(pos.x() + 8.0, pos.y() - 6.0),
                    "style": label_style,
                    "priority": 45,
                    "hide_on_overlap": True,
                }
            )
    painter.restore()


def _aircraft_line_width_px(distance_km: float, *, width_scale: float = 1.0) -> float:
    d = max(0.0, float(distance_km))
    scale = max(1.0, float(width_scale))
    aircraft_scale = 2.4 * scale
    if d <= 1.0:
        return 3.0 * aircraft_scale
    if d <= 3.0:
        return 2.2 * aircraft_scale
    if d <= 5.0:
        return 1.6 * aircraft_scale
    if d <= 10.0:
        return 1.0 * aircraft_scale
    if d <= 20.0:
        return 0.8 * aircraft_scale
    return 0.6 * aircraft_scale
