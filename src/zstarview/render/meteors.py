from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone

import astropy.time
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..meteors.types import MeteorTrail
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy

METEOR_CORE_COLOR = (255, 255, 255)
METEOR_GLOW_COLOR = (255, 220, 120)
METEOR_LABEL_COLOR = (255, 238, 188)
METEOR_FULL_OPACITY_AGE = timedelta(hours=24)
METEOR_FADE_SPAN = timedelta(hours=72)
METEOR_MIN_OPACITY = 0.3
METEOR_AGE_LABEL_PIXEL_SIZE = 9


def meteor_age_opacity(beginning_utc: datetime, display_time_utc: datetime) -> float:
    beginning = _utc(beginning_utc)
    display = _utc(display_time_utc)
    age = display - beginning
    if age < timedelta(0):
        return 0.0
    fade_age = age - METEOR_FULL_OPACITY_AGE
    if fade_age <= timedelta(0):
        return 1.0
    fade_span_seconds = METEOR_FADE_SPAN.total_seconds()
    faded = 1.0 - (1.0 - METEOR_MIN_OPACITY) * (fade_age.total_seconds() / fade_span_seconds)
    return max(METEOR_MIN_OPACITY, faded)


def meteor_age_label(beginning_utc: datetime, display_time_utc: datetime) -> str:
    """Return the compact signed hour age label for a meteor trail."""
    age_hours = int((_utc(display_time_utc) - _utc(beginning_utc)).total_seconds() // 3600)
    return f"{-age_hours:+d}h"


def draw_meteor_trails(painter: QPainter, geometry: ScreenGeometry, *,
                       viewer_data: ViewerData, trails: tuple[MeteorTrail, ...] | None,
                       time_obj: astropy.time.Time | None,
                       opacity: float = 1.0) -> None:
    if not trails or time_obj is None or opacity <= 0.0:
        return
    display_time_utc = time_obj.to_datetime(timezone=timezone.utc)
    painter.save()
    try:
        for trail in trails:
            alpha = meteor_age_opacity(trail.beginning_utc, display_time_utc) * opacity
            if alpha <= 0.0:
                continue
            points = []
            for alt, az in (
                (trail.begin_alt_deg, trail.begin_az_deg),
                (trail.end_alt_deg, trail.end_az_deg),
            ):
                nx, ny = altaz_to_normalized_xy(float(alt), float(az), viewer_data.view_center,
                                                edge_fov_deg=viewer_data.edge_fov_deg)
                x, y = normalized_to_screen_xy(nx, ny, geometry)
                points.append(QPointF(float(x), float(y)))
            _draw_meteor_trail_shape(
                painter,
                points[0],
                points[1],
                color=QColor(*METEOR_GLOW_COLOR, int(round(255 * alpha * 0.4))),
                start_half_width=0.4,
                peak_half_width=1.25,
                end_half_width=0.4,
            )
            _draw_meteor_trail_shape(
                painter,
                points[0],
                points[1],
                color=QColor(*METEOR_CORE_COLOR, int(round(255 * alpha * 0.8))),
                start_half_width=0.18,
                peak_half_width=0.72,
                end_half_width=0.18,
            )
            label_color = QColor(*METEOR_LABEL_COLOR, int(round(255 * alpha)))
            label_font = QFont("Sans Serif")
            label_font.setPixelSize(METEOR_AGE_LABEL_PIXEL_SIZE)
            painter.setFont(label_font)
            label_position = QPointF(points[0].x() + 3.0, points[0].y() - 3.0)
            label_text = meteor_age_label(trail.beginning_utc, display_time_utc)
            painter.setPen(label_color)
            painter.drawText(label_position, label_text)
    finally:
        painter.restore()


def _utc(value: datetime) -> datetime:
    return value.replace(tzinfo=timezone.utc) if value.tzinfo is None else value.astimezone(timezone.utc)


def _draw_meteor_trail_shape(
    painter: QPainter,
    beginning: QPointF,
    ending: QPointF,
    *,
    color: QColor,
    start_half_width: float,
    peak_half_width: float,
    end_half_width: float,
) -> None:
    dx = float(ending.x()) - float(beginning.x())
    dy = float(ending.y()) - float(beginning.y())
    length = math.hypot(dx, dy)
    if length <= 0.0:
        return

    normal_x = -dy / length
    normal_y = dx / length
    peak_fraction = 0.8
    peak_x = float(beginning.x()) + dx * peak_fraction
    peak_y = float(beginning.y()) + dy * peak_fraction

    polygon = QPolygonF(
        [
            QPointF(
                float(beginning.x()) + normal_x * start_half_width,
                float(beginning.y()) + normal_y * start_half_width,
            ),
            QPointF(peak_x + normal_x * peak_half_width, peak_y + normal_y * peak_half_width),
            QPointF(
                float(ending.x()) + normal_x * end_half_width,
                float(ending.y()) + normal_y * end_half_width,
            ),
            QPointF(peak_x - normal_x * peak_half_width, peak_y - normal_y * peak_half_width),
        ]
    )
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(color)
    painter.drawPolygon(polygon)
