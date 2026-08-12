from __future__ import annotations

from datetime import datetime, timedelta, timezone

import astropy.time
from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor, QFont, QPainter, QPen

from ..astro import altaz_to_normalized_xy
from ..meteors.types import MeteorTrail
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy

METEOR_TRAIL_COLOR = (221, 245, 168)
METEOR_FULL_OPACITY_AGE = timedelta(hours=18)
METEOR_MAX_AGE = timedelta(hours=24)
METEOR_AGE_LABEL_PIXEL_SIZE = 9


def meteor_age_opacity(beginning_utc: datetime, display_time_utc: datetime) -> float:
    beginning = _utc(beginning_utc)
    display = _utc(display_time_utc)
    age = display - beginning
    if age < timedelta(0) or age > METEOR_MAX_AGE:
        return 0.0
    if age <= METEOR_FULL_OPACITY_AGE:
        return 1.0
    return max(0.0, (METEOR_MAX_AGE - age) / (METEOR_MAX_AGE - METEOR_FULL_OPACITY_AGE))


def meteor_age_label(beginning_utc: datetime, display_time_utc: datetime) -> str:
    """Return the compact signed hour age label for a meteor trail."""
    age_hours = int((_utc(display_time_utc) - _utc(beginning_utc)).total_seconds() // 3600)
    return f"{-age_hours:+d}h"


def draw_meteor_trails(painter: QPainter, geometry: ScreenGeometry, *,
                       viewer_data: ViewerData, trails: tuple[MeteorTrail, ...] | None,
                       time_obj: astropy.time.Time | None,
                       fade_reference_utc: datetime | None = None,
                       opacity: float = 1.0) -> None:
    if not trails or time_obj is None or opacity <= 0.0:
        return
    display_time_utc = time_obj.to_datetime(timezone=timezone.utc)
    fade_reference = fade_reference_utc or display_time_utc
    painter.save()
    try:
        for trail in trails:
            alpha = meteor_age_opacity(trail.beginning_utc, fade_reference) * opacity
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
            painter.setPen(QPen(QColor(*METEOR_TRAIL_COLOR, int(round(255 * alpha))), 1.2))
            painter.drawLine(points[0], points[1])
            label_color = QColor(*METEOR_TRAIL_COLOR, int(round(180 * alpha)))
            painter.setPen(label_color)
            label_font = QFont("Sans Serif")
            label_font.setPixelSize(METEOR_AGE_LABEL_PIXEL_SIZE)
            painter.setFont(label_font)
            painter.drawText(
                QPointF(points[0].x() + 3.0, points[0].y() - 3.0),
                meteor_age_label(trail.beginning_utc, display_time_utc),
            )
    finally:
        painter.restore()


def _utc(value: datetime) -> datetime:
    return value.replace(tzinfo=timezone.utc) if value.tzinfo is None else value.astimezone(timezone.utc)
