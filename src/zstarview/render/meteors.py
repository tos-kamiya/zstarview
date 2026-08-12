from __future__ import annotations

from datetime import datetime, timedelta, timezone

import astropy.time
import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor, QPainter, QPen

from ..astro import altaz_to_normalized_xy, apply_icrs_to_altaz_matrix, build_icrs_to_altaz_matrix
from ..meteors.types import CelestialMeteorTrail
from ..types import ScreenGeometry, ViewerData
from .geometry import normalized_to_screen_xy

METEOR_TRAIL_COLOR = (240, 72, 72)
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
                       viewer_data: ViewerData, trails: tuple[CelestialMeteorTrail, ...] | None,
                       time_obj: astropy.time.Time | None, opacity: float = 1.0) -> None:
    if not trails or time_obj is None or opacity <= 0.0:
        return
    location = EarthLocation(lat=viewer_data.lat_deg * u.deg, lon=viewer_data.lon_deg * u.deg,
                             height=viewer_data.observer_height_m * u.m)
    matrix = build_icrs_to_altaz_matrix(time_obj, location)
    display_time = time_obj.to_datetime(timezone=timezone.utc)
    painter.save()
    try:
        for trail in trails:
            alpha = meteor_age_opacity(trail.beginning_utc, display_time) * opacity
            if alpha <= 0.0:
                continue
            ra = np.radians([trail.begin_ra_deg, trail.end_ra_deg])
            dec = np.radians([trail.begin_dec_deg, trail.end_dec_deg])
            vectors = np.column_stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)))
            altitudes, azimuths = apply_icrs_to_altaz_matrix(vectors, matrix)
            points = []
            for alt, az in zip(altitudes, azimuths, strict=True):
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
