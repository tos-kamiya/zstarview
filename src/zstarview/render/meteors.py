from __future__ import annotations

from datetime import datetime, timedelta, timezone

import astropy.time
from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor, QPainter, QPen

from ..astro import altaz_to_normalized_xy
from ..meteors.types import MeteorTrail
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy

METEOR_TRAIL_COLOR = (221, 245, 168)
METEOR_FULL_OPACITY_AGE = timedelta(hours=18)
METEOR_MAX_AGE = timedelta(hours=24)


def meteor_age_opacity(beginning_utc: datetime, display_time_utc: datetime) -> float:
    beginning = _utc(beginning_utc)
    display = _utc(display_time_utc)
    age = display - beginning
    if age < timedelta(0) or age > METEOR_MAX_AGE:
        return 0.0
    if age <= METEOR_FULL_OPACITY_AGE:
        return 1.0
    return max(0.0, (METEOR_MAX_AGE - age) / (METEOR_MAX_AGE - METEOR_FULL_OPACITY_AGE))


def draw_meteor_trails(painter: QPainter, geometry: ScreenGeometry, *,
                       viewer_data: ViewerData, trails: tuple[MeteorTrail, ...] | None,
                       time_obj: astropy.time.Time | None,
                       fade_reference_utc: datetime | None = None,
                       opacity: float = 1.0) -> None:
    if not trails or time_obj is None or opacity <= 0.0:
        return
    fade_reference = fade_reference_utc or time_obj.to_datetime(timezone=timezone.utc)
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
    finally:
        painter.restore()


def _utc(value: datetime) -> datetime:
    return value.replace(tzinfo=timezone.utc) if value.tzinfo is None else value.astimezone(timezone.utc)
