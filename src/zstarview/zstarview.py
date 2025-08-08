# -*- coding: utf-8 -*-
import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import math
import os.path
from pathlib import Path
import sys
import threading
import time
from typing import Any, cast, List, Dict, Tuple, Union, Optional
from zoneinfo import ZoneInfo

from PyQt5.QtWidgets import QApplication, QSizeGrip, QSplashScreen, QMainWindow
from PyQt5.QtCore import Qt, QPoint, QPointF, QRectF, QSize, QTimer, pyqtSignal
from PyQt5.QtGui import QKeyEvent, QMouseEvent, QPaintEvent
from PyQt5.QtGui import (
    QBrush,
    QColor,
    QFont,
    QFontDatabase,
    QIcon,
    QImage,
    QPainter,
    QPen,
    QPixmap,
    QPolygonF,
    QRadialGradient,
    QResizeEvent,
)

from appdirs import user_cache_dir
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, GeocentricTrueEcliptic
from astropy.time import Time
import astropy.units as u
from skyfield.api import Loader, Topos
import skyfield.api

import numpy as np
from PIL import Image

# --- Helper Functions ---
cache_path = Path(user_cache_dir(appname="zstarview", appauthor="tos-kamiya"))
cache_path.mkdir(parents=True, exist_ok=True)
starfield_load = Loader(str(cache_path))


def pil2qpixmap(img: Image.Image) -> QPixmap:
    """Converts a PIL Image to a PyQt5 QPixmap."""
    arr = np.array(img.convert("RGBA"))
    h, w, ch = arr.shape
    bytes_per_line = ch * w
    qimg = QImage(arr.data.tobytes(), w, h, bytes_per_line, QImage.Format.Format_RGBA8888)
    return QPixmap.fromImage(qimg)


# --- Constants and Data Loading ---
_dir = os.path.dirname(os.path.abspath(__file__))

EMOJI_FONT_PATH = os.path.join(_dir, "data", "Noto_Sans_Symbols", "NotoSansSymbols-VariableFont_wght.ttf")
CITY_COORD_FILE = os.path.join(_dir, "data", "cities1000.txt")
STARS_CSV_FILE = os.path.join(_dir, "data", "stars.csv")
APP_ICON_FILE = os.path.join(_dir, "data", "app_icon.ico")

TEXT_COLOR = QColor(120, 120, 120)
TEXT_FONT_SIZE = 14

HORIZON_LINE_COLOR = QColor(40, 50, 40)
CELESTIAL_EQUATOR_COLOR = QColor(40, 40, 40)
ECLIPTIC_COLOR = QColor(80, 60, 0)

FIELD_OF_VIEW_DEG = 240
ANGLE_BELOW_HORIZON = 2

DIRECTIONS = {
    "N": 0,
    "E": 90,
    "S": 180,
    "W": 270,
    "NE": 45,
    "SE": 135,
    "SW": 225,
    "NW": 315,
}

PLANET_SYMBOLS = {
    "sun": "☀",
    "moon": "🌛",
    "mercury": "☿",
    "venus": "♀",
    "mars": "♂",
    "jupiter barycenter": "♃",
    "saturn barycenter": "♄",
}


@dataclass
class StarData:
    name: str
    alt: float
    az: float
    vmag: float
    bv: float


@dataclass
class PlanetBody:
    name: str
    alt: float
    az: float
    symbol: str
    phase_angle: Optional[float] = None  # moon only


@dataclass
class SkyData:
    """Container for all calculated sky data for a specific time and location."""

    location: Tuple[float, float]  # (lat, lon)
    time: Time
    planets: List[PlanetBody]
    stars: List[StarData]
    celestial_equator_points: List[Tuple[float, float]]
    ecliptic_points: List[Tuple[float, float]]
    horizon_points: List[Tuple[float, float]]
    timezone_name: str
    city_name: str
    view_center: Tuple[float, float] = (90.0, 180.0)


def load_city_coords(filename: str) -> Dict[str, Tuple[float, float, str]]:
    """Loads city coordinates and timezone from the data file."""
    city_table = {}
    with open(filename, encoding="utf-8") as f:
        for line in f:
            cols = line.strip().split("\t")
            # Check if the row has enough columns including timezone information
            if len(cols) < 18:
                continue
            name = cols[1]
            lat = float(cols[4])
            lon = float(cols[5])
            country = cols[8]
            timezone_name = cols[17]  # Get timezone information
            key = f"{country.lower()}/{name.lower()}"
            city_table[key] = (lat, lon, timezone_name)
    return city_table


def bv_to_qcolor(bv: float) -> QColor:
    """Converts a B-V color index to a QColor."""
    if bv < 0.0:
        return QColor(170, 191, 255)
    elif bv < 0.3:
        return QColor(202, 215, 255)
    elif bv < 0.6:
        return QColor(248, 247, 255)
    elif bv < 1.0:
        return QColor(255, 210, 161)
    else:
        return QColor(255, 204, 111)


def generate_moon_phase_image(
    size: int,
    sun_dir_3d: np.ndarray,
    view_dir_3d: np.ndarray,
    moon_color: Tuple[int, int, int] = (230, 230, 230),
    dark_color: Tuple[int, int, int] = (30, 30, 30),
    earthshine_factor: float = 0.15,
) -> Image.Image:
    """Generates a spherical image of the moon phase as seen by the observer, including earthshine."""
    cx, cy = size // 2, size // 2
    r = size // 2
    img_array = np.zeros((size, size, 3), dtype=np.uint8)
    for y in range(size):
        for x in range(size):
            dx = (x - cx) / r
            dy = (y - cy) / r
            if dx**2 + dy**2 > 1:
                continue
            dz = np.sqrt(1 - dx**2 - dy**2)
            surf = np.array([dx, dy, dz])
            surf /= np.linalg.norm(surf)
            view = np.dot(surf, view_dir_3d)
            if view < 0:
                continue  # Don't draw the back side

            light = np.dot(surf, sun_dir_3d)
            if light > 0:
                c = np.array(moon_color) * light
            else:
                # Earthshine (make the dark side slightly brighter)
                c = np.array(dark_color) * earthshine_factor
            c = np.clip(c, 0, 255)
            img_array[y, x] = c
    return Image.fromarray(img_array)


def altaz_to_normalized_xy(alt: float, az: float, view_center: Tuple[float, float]) -> Tuple[float, float]:
    """Converts altitude and azimuth coordinates to custom screen coordinates relative to a view center."""
    center_alt, center_az = view_center

    # Convert degrees to radians
    alt1 = math.radians(center_alt)
    az1 = math.radians(center_az)
    alt2 = math.radians(alt)
    az2 = math.radians(az)

    # Angular distance using spherical trigonometry
    cos_theta = math.sin(alt1) * math.sin(alt2) + math.cos(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    theta = math.acos(cos_theta)  # [radian]

    # Map angular distance to a normalized radius (0 at center, 1 at 90 degrees difference)
    r = theta / (math.pi / 2)
    # If theta exceeds pi/2 (object is on the unseen opposite side), r > 1

    # Calculate the direction of the star relative to the view center
    # For simplicity, using the star's azimuth (can be adjusted based on projection type)

    # Projection direction
    dx = math.cos(alt2) * math.sin(az2 - az1)
    dy = math.cos(alt1) * math.sin(alt2) - math.sin(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    # Normalize
    length = math.hypot(dx, dy)
    if length != 0:
        dx /= length
        dy /= length
    nx = r * dx
    ny = -r * dy  # Adjust y-axis sign according to display system

    return (nx, ny)


def load_star_catalog(filename: str) -> List[Dict[str, Any]]:
    """Loads the star catalog from a CSV file."""
    star_catalog = []
    with open(filename, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                ra = float(row["RA"]) * 15
                name = row["Name"]
                dec = float(row["Dec"])
                vmag = float(row["Vmag"])
                bv = float(row["B-V"])
                coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
                star_catalog.append({"name": name, "coord": coord, "vmag": vmag, "bv": bv})
            except Exception:
                continue
    return star_catalog


def calculate_visible_stars(star_catalog: List[Dict[str, Any]], lat: float, lon: float, time_obj: Time, view_center: Tuple[float, float]) -> Tuple[List[StarData], EarthLocation]:
    """Calculates the positions of stars based on the given time and location."""
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    visible_stars = []
    for i, star in enumerate(star_catalog):
        coord = cast(SkyCoord, star["coord"])
        altaz = coord.transform_to(AltAz(obstime=time_obj, location=location))
        if altaz.alt.deg > -ANGLE_BELOW_HORIZON and is_in_fov(altaz.alt.deg, altaz.az.deg, view_center):
            visible_stars.append(StarData(name=star["name"], alt=altaz.alt.deg, az=altaz.az.deg, vmag=star["vmag"], bv=star["bv"]))
        if (i + 1) % 500 == 0:
            time.sleep(0)
    return (visible_stars, location)


def calculate_visible_planets(lat: float, lon: float, astropy_time: Time, view_center: Tuple[float, float]) -> List[PlanetBody]:
    """Calculates and returns a list of visible planets."""
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(astropy_time)
    planets = starfield_load("de421.bsp")
    observer = planets["earth"] + Topos(latitude_degrees=lat, longitude_degrees=lon)

    visible_bodies = []
    for name, symbol in PLANET_SYMBOLS.items():
        planet = planets[name]
        astrometric = observer.at(t).observe(planet).apparent()
        alt, az, _ = astrometric.altaz()
        if alt.degrees > -ANGLE_BELOW_HORIZON and is_in_fov(alt.degrees, az.degrees, view_center):
            phase_angle = None
            if name == "moon":
                phase_angle = calculate_moon_phase_angle(observer, t, planets)
            visible_bodies.append(
                PlanetBody(
                    name=name,
                    alt=alt.degrees,
                    az=az.degrees,
                    symbol=symbol,
                    phase_angle=phase_angle,
                )
            )
    return visible_bodies


def calculate_moon_phase_angle(observer: Topos, t: skyfield.timelib.Time, planets: Any) -> float:
    """Calculates the phase angle of the Moon."""
    moon = planets["moon"]
    sun = planets["sun"]
    e_to_moon = observer.at(t).observe(moon).apparent()
    e_to_sun = observer.at(t).observe(sun).apparent()
    return e_to_moon.separation_from(e_to_sun).degrees


def is_in_fov(alt: float, az: float, view_center: Tuple[float, float]) -> bool:
    """Checks if a celestial object is within the field of view."""
    center_alt, center_az = view_center
    alt1, az1 = math.radians(center_alt), math.radians(center_az)
    alt2, az2 = math.radians(alt), math.radians(az)
    cos_theta = math.sin(alt1) * math.sin(alt2) + math.cos(alt1) * math.cos(alt2) * math.cos(az2 - az1)
    theta = math.acos(min(1.0, max(-1.0, cos_theta)))
    return math.degrees(theta) <= FIELD_OF_VIEW_DEG / 2


def calculate_horizon_points(location: EarthLocation, time: Time, view_center: Tuple[float, float]) -> List[Tuple[float, float]]:
    """Calculates points along the horizon for drawing."""
    points = []
    alt = 0.0  # Horizon
    for az in range(0, 360 + 5, 5):
        if not is_in_fov(alt, az, view_center):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        points.append((nx, ny))
    return points


def calculate_celestial_equator_points(location: EarthLocation, time: Time, view_center: Tuple[float, float]) -> List[Tuple[float, float]]:
    """Calculates points along the celestial equator for drawing."""
    a = 5
    points = []
    for ra_deg in range(0, 360 + 5, 5):
        coord = SkyCoord(ra=ra_deg * u.deg, dec=0 * u.deg, frame="icrs")
        altaz = coord.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > -a and is_in_fov(altaz.alt.deg, altaz.az.deg, view_center):
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg, view_center)
            points.append((nx, ny))
    return points


def calculate_ecliptic_points(location: EarthLocation, time: Time, view_center: Tuple[float, float]) -> List[Tuple[float, float]]:
    """Calculates points along the ecliptic for drawing."""
    points = []
    for lon_deg in range(0, 360 + 5, 5):
        ecl = SkyCoord(
            lon=lon_deg * u.deg,
            lat=0 * u.deg,
            frame=GeocentricTrueEcliptic(obstime=time),
        )
        icrs = ecl.transform_to("icrs")
        altaz = icrs.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > -ANGLE_BELOW_HORIZON and is_in_fov(altaz.alt.deg, altaz.az.deg, view_center):
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg, view_center)
            points.append((nx, ny))
    return points


def calculate_sun_angle_on_moon(moon_altaz: Tuple[float, float], sun_altaz: Tuple[float, float]) -> float:
    """Calculates the angle of the sun relative to the moon on the screen."""
    m_alt, m_az = moon_altaz
    s_alt, s_az = sun_altaz

    # Relative direction of the sun on the celestial sphere with respect to the moon
    # Azimuth difference (azimuth: North=0, East=90, South=180, West=270)
    # In screen coordinates, y+ is North (up), x- is East (left)

    # Which direction is the sun relative to the moon?
    d_az = math.radians(s_az - m_az)
    d_alt = math.radians(s_alt - m_alt)

    # Relative sun direction vector on screen (x-left=East, y-up=North)
    dx = -math.sin(d_az) * math.cos(math.radians(s_alt))
    dy = math.sin(d_alt)

    # Angle on screen from (dx, dy) with North as up (North=0, East=90, South=180, West=270, clockwise positive)
    angle = math.atan2(dx, dy)
    return angle


def normalized_to_screen_xy(nx: float, ny: float, center: QPoint, radius: int) -> QPointF:
    """Converts normalized coordinates to screen coordinates."""
    return QPointF(center.x() + nx * radius, center.y() + ny * radius)


def find_highlighted_object(sky_data: Optional[SkyData], mouse_pos: QPoint, center: QPoint, radius: int) -> Optional[Tuple[Union[StarData, PlanetBody], QPointF]]:
    """Finds the celestial object closest to the mouse cursor, within FOV."""
    min_dist = 30**2  # squared pixels
    highlighted_object = None

    if not sky_data:
        return None

    all_objects = sky_data.stars + sky_data.planets
    for obj in all_objects:
        # Check if within field of view
        if not is_in_fov(obj.alt, obj.az, sky_data.view_center):
            continue
        nx, ny = altaz_to_normalized_xy(obj.alt, obj.az, sky_data.view_center)
        pos = normalized_to_screen_xy(nx, ny, center, radius)
        dist_sq = (mouse_pos.x() - pos.x()) ** 2 + (mouse_pos.y() - pos.y()) ** 2
        if dist_sq < min_dist:
            min_dist = dist_sq
            highlighted_object = (obj, pos)
    return highlighted_object


def draw_radial_background(painter: QPainter, rect: QRectF, center: QPoint, radius: int):
    assert radius >= 10
    fov_middle = 90 + (FIELD_OF_VIEW_DEG / 2 - 90) / 2
    r90 = float(radius * (90 / 90))  # = radius
    r_fov = float(radius * (fov_middle / 90))
    r_max = float(r_fov * 1.5)
    step_px = 0.5

    def pos(r):
        return max(0.0, min(1.0, r / r_max))

    def col(r, s):
        return QColor(0, 0, 0, max(0, 255 - (s + int(150 * (r - r90) / r_max))))

    g = QRadialGradient(center, r_max)

    g.setColorAt(pos(0), col(r90, 0))
    g.setColorAt(pos(r90), col(r90, 0))

    g.setColorAt(pos(r90 + step_px), col(r90, 10))
    g.setColorAt(pos(r_fov), col(r_fov, 10))

    g.setColorAt(pos(r_fov + step_px), col(r_fov, 20))
    g.setColorAt(1.0, col(r_max, 20))

    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(g)
    painter.drawRect(rect)
    painter.restore()


def split_by_gaps(points: List[Tuple[float, float]]) -> List[List[Tuple[float, float]]]:
    """Splits a list of points into fragments based on gaps between consecutive points."""

    def dist(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
        return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

    fragments = [[]]
    for p in points:
        if not fragments[-1] or dist(p, fragments[-1][-1]) < 0.3:
            fragments[-1].append(p)
        else:
            fragments.append([p])
    return fragments


def draw_sky_reference_lines(painter: QPainter, center: QPoint, radius: int, sky_data: SkyData):
    """Draws celestial lines (equator, ecliptic, horizon) on the screen."""
    point_list_pen_styles = [
        (sky_data.celestial_equator_points, (CELESTIAL_EQUATOR_COLOR, 2, Qt.PenStyle.DashLine)),
        (sky_data.ecliptic_points, (ECLIPTIC_COLOR, 2, Qt.PenStyle.DotLine)),
        (sky_data.horizon_points, (HORIZON_LINE_COLOR, 2, Qt.PenStyle.SolidLine)),
    ]
    for points, pen_style in point_list_pen_styles:
        for frag in split_by_gaps(points):
            if len(frag) >= 2:
                points = [normalized_to_screen_xy(nx, ny, center, radius) for nx, ny in frag]
                poly = QPolygonF(points)
                painter.setPen(QPen(*pen_style))
                painter.drawPolyline(poly)


def draw_stars(painter: QPainter, center: QPoint, radius: int, sky_data: SkyData, star_base_radius: float):
    """Draws the stars on the screen."""

    def mag_to_size(vmag: float) -> float:
        return max(0.1, star_base_radius * 10 ** (-0.2 * vmag)) * radius / 500

    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    for star in sky_data.stars:
        if not is_in_fov(star.alt, star.az, sky_data.view_center):
            continue
        pos = normalized_to_screen_xy(*altaz_to_normalized_xy(star.alt, star.az, sky_data.view_center), center, radius)
        color = bv_to_qcolor(star.bv)
        siz = mag_to_size(star.vmag)

        # For stars smaller than 4.0, draw as a minimum size (2x2) square
        if siz < 4.0:
            alpha_value = min(1.0, max(0.1, siz / 4.0))
            color.setAlphaF(alpha_value)

            painter.fillRect(QRectF(pos.x() - 1, pos.y() - 1, 2, 2), color)
        else:  # For larger stars (size >= 4.0), draw as a gradient circle
            r = math.sqrt(siz)
            gradient = QRadialGradient(pos, r)
            gradient.setColorAt(0, color)

            color_transparent = QColor(color)
            color_transparent.setAlpha(0)
            gradient.setColorAt(1, color_transparent)

            painter.setBrush(QBrush(gradient))
            painter.setPen(Qt.PenStyle.NoPen)
            painter.drawEllipse(pos, r, r)

    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)  # Reset mode


def draw_gauge_cross(painter: QPainter, color: QColor, center: QPointF):
    """Draws a cross gauge symbol at the given center point."""
    cross_outer_len, cross_inner_len = 15, 4
    x, y = center.x(), center.y()
    painter.setPen(QPen(color, 1))
    painter.drawLine(QPointF(x - cross_outer_len, y), QPointF(x - cross_inner_len, y))
    painter.drawLine(QPointF(x + cross_inner_len, y), QPointF(x + cross_outer_len, y))
    painter.drawLine(QPointF(x, y - cross_outer_len), QPointF(x, y - cross_inner_len))
    painter.drawLine(QPointF(x, y + cross_inner_len), QPointF(x, y + cross_outer_len))


def draw_moon(
    painter: QPainter,
    center: QPointF,
    radius: float,
    phase_angle_deg: float,
    sun_altaz: Optional[Tuple[float, float]] = None,
    moon_altaz: Optional[Tuple[float, float]] = None,
    opacity: float = 1.0,
):
    """Draws the moon with its phase and earthshine."""
    img_size = int(radius * 2)
    if img_size < 5:
        img_size = 5

    view_dir = np.array([0, 0, 1])
    phase_angle_rad = math.radians(phase_angle_deg)
    sun_dir = np.array([np.sin(phase_angle_rad), 0, -np.cos(phase_angle_rad)])
    sun_dir /= np.linalg.norm(sun_dir)
    moon_img_pil = generate_moon_phase_image(img_size, sun_dir, view_dir)

    rotate_deg = 0
    if sun_altaz is not None and moon_altaz is not None:
        angle = calculate_sun_angle_on_moon(moon_altaz, sun_altaz)
        rotate_deg = -math.degrees(angle) - 90

    moon_img_pil = moon_img_pil.rotate(rotate_deg, resample=Image.Resampling.BICUBIC, expand=False)

    pixmap = pil2qpixmap(moon_img_pil)
    target_rect = QRectF(center.x() - img_size / 2, center.y() - img_size / 2, img_size, img_size)
    painter.save()
    painter.setOpacity(opacity)
    painter.drawPixmap(target_rect, pixmap, QRectF(0, 0, img_size, img_size))
    painter.restore()


def draw_planets(painter: QPainter, center: QPoint, radius: int, sky_data: SkyData, enlarge_moon: bool, emoji_font: QFont):
    """Draws the planets, Sun, and Moon on the screen."""
    sun_altaz: Optional[Tuple[float, float]] = None
    moon_altaz: Optional[Tuple[float, float]] = None
    for body in sky_data.planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        if body.name == "moon":
            moon_altaz = (body.alt, body.az)

    for body in sky_data.planets:
        pos = normalized_to_screen_xy(*altaz_to_normalized_xy(body.alt, body.az, sky_data.view_center), center, radius)
        if body.name == "sun":
            draw_gauge_cross(painter, TEXT_COLOR, pos)
        elif body.name == "moon":
            moon_radius = 0.5 * (1 if not enlarge_moon else 3) / 2 * (radius / 90.0)
            draw_moon(
                painter,
                pos,
                moon_radius,
                body.phase_angle if body.phase_angle is not None else 0.0,
                sun_altaz=sun_altaz,
                moon_altaz=moon_altaz,
                opacity=1.0 if not enlarge_moon else 0.7,
            )
            draw_gauge_cross(painter, TEXT_COLOR, pos)
        else:
            painter.setFont(emoji_font)
            painter.setPen(TEXT_COLOR)
            painter.drawText(pos, body.symbol)


def draw_direction_labels(painter: QPainter, center: QPoint, radius: int, view_center: Tuple[float, float], text_font: QFont):
    """Draws direction labels (N, E, S, W) on the screen."""
    painter.setPen(TEXT_COLOR)
    painter.setFont(text_font)
    alt = 0
    for label, az in DIRECTIONS.items():
        if not is_in_fov(alt, az, view_center):
            continue
        nx, ny = altaz_to_normalized_xy(alt, az, view_center)
        pos = normalized_to_screen_xy(nx, ny, center, radius)
        painter.drawText(pos, label)


def draw_overlay_info(
    painter: QPainter,
    sky_data: SkyData,
    highlighted_object: Optional[Tuple[Union[StarData, PlanetBody], QPointF]],
    text_font: QFont,
):
    """Draws time information, city name, and the label for the highlighted object."""

    utc_time = sky_data.time
    tz_name = sky_data.timezone_name
    time_text = ""
    try:
        local_tz = ZoneInfo(tz_name)
        local_dt = utc_time.to_datetime(timezone=local_tz)
        time_text = local_dt.strftime("%Y-%m-%d %H:%M:%S %Z")
    except Exception:
        time_text = utc_time.to_datetime().strftime("%Y-%m-%d %H:%M:%S UTC")

    painter.setPen(QColor(180, 180, 180))
    painter.setFont(text_font)
    painter.drawText(QPoint(10, 20), time_text)

    city_name_text = sky_data.city_name.title()
    painter.drawText(QPoint(10, 40), city_name_text)

    if highlighted_object:
        obj, pos = highlighted_object
        painter.setPen(QPen(TEXT_COLOR, 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(pos, 10, 10)

        name = obj.name or ""
        painter.setPen(TEXT_COLOR)
        painter.drawText(QPointF(pos.x() + 15, pos.y() - 15), str(name))


def get_screen_geometry(width: int, height: int, alt: float) -> Tuple[QPoint, int]:
    """Calculates the center and radius for drawing based on window size and view altitude."""
    margin_x = 10
    margin_y = 10
    radius = (width - margin_x * 2) // 2
    ud = 90
    dd = alt
    center = QPoint(radius + margin_x, int((height - margin_y * 2) * ud / (ud + dd)) + margin_y)
    return center, radius


class SkyWindow(QMainWindow):
    data_updated = pyqtSignal(object)
    initial_data_loaded = pyqtSignal()

    def __init__(
        self,
        city_name: str,
        city_data: Tuple[float, float, str],
        star_catalog: List[Dict[str, Any]],
        delta_t: timedelta,
        enlarge_moon: bool,
        star_base_radius: float,
        view_center: Tuple[float, float],
    ):
        """Initializes the SkyWindow."""
        super().__init__()
        self.setWindowIcon(QIcon(APP_ICON_FILE))

        self.city_name = city_name
        self.lat, self.lon, self.tz_name = city_data
        self.star_catalog = star_catalog
        self.delta_t = delta_t
        self.enlarge_moon = enlarge_moon
        self.star_base_radius = star_base_radius
        self.view_center = view_center

        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setWindowFlags(Qt.FramelessWindowHint)
        self.setWindowTitle(f"Zenith Star View - {self.city_name.title()}")
        self.setGeometry(100, 100, 800, 800)

        # Size grip
        self.size_grip = QSizeGrip(self)
        self.size_grip.setFixedSize(24, 24)
        self.size_grip.raise_()

        # Fonts for drawing
        font_id = QFontDatabase.addApplicationFont(EMOJI_FONT_PATH)
        font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
        self.emoji_font = QFont(font_family, TEXT_FONT_SIZE + 4)
        self.text_font = QFont("Arial", TEXT_FONT_SIZE)

        self.setMouseTracking(True)
        self.setMinimumSize(400, 400)
        self.sky_data: Optional[SkyData] = None
        self.mouse_pos: Optional[QPoint] = None

        self.data_updated.connect(self.on_data_updated)
        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.start_background_update)

        self.start_background_update(is_initial_load=True)

    def resizeEvent(self, event: QResizeEvent):
        """Handles the resize event for the window."""
        grip_size = self.size_grip.size()
        self.size_grip.move(self.width() - grip_size.width(), self.height() - grip_size.height())
        super().resizeEvent(event)

    def mousePressEvent(self, event: QMouseEvent):
        """Handles mouse press events for window dragging."""
        if event.button() == Qt.LeftButton:
            self._drag_active = True
            self._drag_pos = event.globalPos() - self.frameGeometry().topLeft()
            event.accept()

    def leaveEvent(self, event):
        self.mouse_pos = None
        self.update()
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent):
        """Handles mouse move events for window dragging and object highlighting."""
        if getattr(self, "_drag_active", False) and event.buttons() & Qt.LeftButton:
            self.move(event.globalPos() - self._drag_pos)
            event.accept()
        else:
            self.mouse_pos = event.pos()
            self.update()
            event.accept()

    def mouseReleaseEvent(self, event: QMouseEvent):
        """Handles mouse release events."""
        self._drag_active = False
        event.accept()

    def set_star_base_radius(self, star_base_radius: float):
        """Sets the base radius for stars."""
        self.star_base_radius = star_base_radius

    def set_enlarge_moon(self, enlarge_moon: bool):
        """Sets whether the moon should be enlarged."""
        self.enlarge_moon = enlarge_moon

    def set_sky_data(self, data: SkyData):
        """Sets the sky data and triggers a repaint."""
        self.sky_data = data
        self.update()

    def paintEvent(self, event: QPaintEvent):
        """Handles the paint event for drawing the sky."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        if not self.sky_data:
            painter.setPen(Qt.GlobalColor.white)
            painter.setFont(QFont("Arial", 16))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Loading sky data...")
            return

        alt = self.sky_data.view_center[0]
        center, radius = get_screen_geometry(self.width(), self.height(), alt)

        painter.setCompositionMode(QPainter.CompositionMode_Clear)
        painter.fillRect(self.rect(), Qt.transparent)

        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)

        draw_radial_background(painter, self.rect(), center, radius)

        draw_sky_reference_lines(painter, center, radius, self.sky_data)
        draw_direction_labels(painter, center, radius, self.sky_data.view_center, self.text_font)

        draw_stars(painter, center, radius, self.sky_data, self.star_base_radius)
        draw_planets(painter, center, radius, self.sky_data, self.enlarge_moon, self.emoji_font)

        highlighted_object = None
        if self.mouse_pos is not None:
            highlighted_object = find_highlighted_object(self.sky_data, self.mouse_pos, center, radius)
        draw_overlay_info(painter, self.sky_data, highlighted_object, self.text_font)

    def on_data_updated(self, sky_data: SkyData):
        """Slot to receive updated sky data and trigger a repaint."""
        self.set_sky_data(sky_data)
        if not self.update_timer.isActive():
            self.update_timer.start(5 * 60 * 1000)
            self.initial_data_loaded.emit()

    def update_sky_data_in_background(self):
        """Updates sky data in a background thread and emits data_updated signal."""
        try:
            now = datetime.now(timezone.utc) + self.delta_t
            time_obj = Time(now)
            stars, loc = calculate_visible_stars(self.star_catalog, self.lat, self.lon, time_obj, self.view_center)
            planets = calculate_visible_planets(self.lat, self.lon, time_obj, self.view_center)
            celestial_equator_points = calculate_celestial_equator_points(loc, time_obj, self.view_center)
            ecliptic_points = calculate_ecliptic_points(loc, time_obj, self.view_center)
            horizon_points = calculate_horizon_points(loc, time_obj, self.view_center)
            sky_data = SkyData(
                (self.lat, self.lon),
                time_obj,
                planets,
                stars,
                celestial_equator_points,
                ecliptic_points,
                horizon_points,
                timezone_name=self.tz_name,
                city_name=self.city_name,
                view_center=self.view_center,
            )
            self.data_updated.emit(sky_data)
        except Exception as e:
            print(f"Error in background update thread: {e}", file=sys.stderr)
            import traceback

            traceback.print_exc()

    def start_background_update(self, is_initial_load: bool = False):
        """Starts a background thread to update sky data."""
        if is_initial_load:
            print("Performing initial data load...")
        else:
            print("Updating sky data...")
        thread = threading.Thread(target=self.update_sky_data_in_background)
        thread.daemon = True
        thread.start()

    def keyPressEvent(self, event: QKeyEvent):
        """Handles key press events."""
        if event and event.key() == Qt.Key.Key_F11:
            if self.isFullScreen():
                self.showNormal()
            else:
                self.showFullScreen()
        elif event and event.key() == Qt.Key.Key_Escape:
            if self.isFullScreen():
                self.showNormal()
        elif event and event.key() == Qt.Key.Key_Q:
            QApplication.quit()
        else:
            super().keyPressEvent(event)


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    parser.add_argument("city", type=str, nargs="?", default="Tokyo", help="City name (default: Tokyo)")
    parser.add_argument("-H", "--hours", type=float, default=0, help="Number of hours to add to current time (default: 0)")
    parser.add_argument("-D", "--days", type=float, default=0, help="Number of days to add to current time (default: 0)")
    parser.add_argument(
        "-m",
        "--enlarge-moon",
        action="store_true",
        help="Show the moon in 3x size.",
    )
    parser.add_argument("-s", "--star-base-radius", type=float, default=15.0, help="Base size of stars (default: 15.0)")
    parser.add_argument("-Z", "--view-center-az", type=float, default=180.0, help="Viewing azimuth angle [deg] (0=N, 90=E, 180=S, 270=W; default=180)")
    parser.add_argument("-A", "--view-center-alt", type=float, default=90.0, help="Viewing altitude angle [deg] (90=zenith, 0=horizon; default=90)")
    return parser.parse_args()


def main():
    """Main entry point for the star sky visualizer."""
    app = QApplication(sys.argv)

    splash = QSplashScreen(QPixmap(400, 200), Qt.WindowType.WindowStaysOnTopHint)
    splash.show()

    def show_splash_message(message: str, color: QColor):
        splash.showMessage(message, Qt.AlignmentFlag.AlignCenter, color)
        app.processEvents()

    try:
        city_table = load_city_coords(CITY_COORD_FILE)
    except FileNotFoundError:
        show_splash_message("Error: cities1000.txt not found.", Qt.GlobalColor.red)
        time.sleep(3)
        return

    app.setWindowIcon(QIcon(APP_ICON_FILE))

    args = parse_args()
    city = args.city
    city = city.lower()
    if city not in city_table:
        candidate_cities = [c for c in city_table.keys() if c.endswith("/" + city)]
        if not candidate_cities:
            print(f"Unknown city: {city}")
            return
        elif len(candidate_cities) > 1:
            print(f"Specify explicit country name: {candidate_cities}")
            return
        else:
            city = candidate_cities[0]

    delta_hours = args.hours
    delta_days = args.days

    delta_t = timedelta(days=delta_days, hours=delta_hours)
    view_center = (args.view_center_alt, args.view_center_az)

    # Show a splash screen while loading initial data
    pixmap = QPixmap(400, 200)
    pixmap.fill(Qt.GlobalColor.black)
    splash.setPixmap(pixmap)
    show_splash_message("Loading city and star data...", Qt.GlobalColor.white)

    try:
        star_catalog = load_star_catalog(STARS_CSV_FILE)
    except FileNotFoundError:
        show_splash_message("Error: stars.csv not found.", Qt.GlobalColor.red)
        time.sleep(3)
        return

    show_splash_message(f"Calculating sky for {city.title()}...", Qt.GlobalColor.white)

    lat, lon, tz_name = city_table[city]
    main_win = SkyWindow(
        city,
        (lat, lon, tz_name),
        star_catalog,
        delta_t,
        enlarge_moon=args.enlarge_moon,
        star_base_radius=args.star_base_radius,
        view_center=view_center,
    )

    # When the initial data is loaded, show the main window and close the splash screen
    main_win.initial_data_loaded.connect(main_win.show)
    main_win.initial_data_loaded.connect(splash.close)

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
