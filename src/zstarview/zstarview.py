# -*- coding: utf-8 -*-
import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import math
import os.path
import sys
import threading
import time
from typing import Any, cast
from zoneinfo import ZoneInfo

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QSplashScreen
from PyQt5.QtGui import (
    QPainter,
    QColor,
    QFont,
    QPen,
    QBrush,
    QPolygonF,
    QPixmap,
    QFontDatabase,
    QRadialGradient,
    QImage,
)
from PyQt5.QtCore import Qt, QPoint, QPointF, QRectF, QTimer, pyqtSignal, QObject

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, GeocentricTrueEcliptic
from astropy.time import Time
import astropy.units as u
from skyfield.api import Topos
import skyfield.api

import numpy as np
from PIL import Image


def pil2qpixmap(img: Image.Image) -> QPixmap:
    arr = np.array(img.convert("RGBA"))
    h, w, ch = arr.shape
    bytes_per_line = ch * w
    qimg = QImage(arr.data.tobytes(), w, h, bytes_per_line, QImage.Format.Format_RGBA8888)
    return QPixmap.fromImage(qimg)


# --- Constants and Data Loading (mostly unchanged) ---
_dir = os.path.dirname(os.path.abspath(__file__))

EMOJI_FONT_PATH = os.path.join(_dir, "data", "Noto_Sans_Symbols", "NotoSansSymbols-VariableFont_wght.ttf")
CITY_COORD_FILE = os.path.join(_dir, "data", "cities1000.txt")
STARS_CSV_FILE = os.path.join(_dir, "data", "stars.csv")

TEXT_COLOR = QColor(120, 120, 120)
TEXT_FONT_SIZE = 14

HORIZON_LINE_COLOR = QColor(40, 40, 40)
CELESTIAL_EQUATOR_COLOR = QColor(40, 40, 40)
ECLIPTIC_COLOR = QColor(80, 60, 0)

VIEWER_LOCATION_HEIGHT = 1.0  # km


PLANET_SYMBOLS = {
    "sun": "☀",
    "moon": "🌛",
    "mercury": "☿",
    "venus": "♀",
    "mars": "♂",
    "jupiter barycenter": "♃",
    "saturn barycenter": "♄",
}


def load_city_coords(filename: str) -> dict[str, tuple[float, float, str]]:
    """
    Loads city coordinates and timezone from the data file.

    Returns:
        dict[str, tuple[float, float, str]]: A dict mapping city key to (latitude, longitude, timezone_name).
    """
    city_table = {}
    with open(filename, encoding="utf-8") as f:
        for line in f:
            cols = line.strip().split("	")
            # タイムゾーン情報を含む列があるかチェック
            if len(cols) < 18:
                continue
            name = cols[1]
            lat = float(cols[4])
            lon = float(cols[5])
            country = cols[8]
            timezone_name = cols[17]  # タイムゾーン情報を取得
            key = f"{country.lower()}/{name.lower()}"
            city_table[key] = (lat, lon, timezone_name)
    return city_table


def bv_to_color(bv: float) -> QColor:
    """Bv to color."""
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


def render_moon_phase_img(
    size: int, sun_dir_3d: np.ndarray, view_dir_3d: np.ndarray, moon_color=(230, 230, 230), dark_color=(30, 30, 30)
) -> Image.Image:
    """
    観測者から見た月相の球体画像を生成
    """
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
            surf_norm = np.array([dx, dy, dz])
            surf_norm /= np.linalg.norm(surf_norm)
            # 観測方向へ
            view = np.dot(surf_norm, view_dir_3d)
            if view < 0:
                continue  # 裏側は描かない
            # 太陽光
            light = np.dot(surf_norm, sun_dir_3d)
            if light > 0:
                c = np.array(moon_color) * light
                c = np.clip(c, 0, 255)
            else:
                c = np.array(dark_color)
            img_array[y, x] = c
    return Image.fromarray(img_array)


def angle_below_horizon(h_km: float, R: float = 6371) -> float:
    ratio = R / (R + h_km)
    theta_rad = math.acos(ratio)
    return math.degrees(theta_rad)


def altaz_to_normalized_xy(alt: float, az: float) -> tuple[float, float]:
    """Altaz to normalized xy."""
    alt = float(alt)
    az = math.radians(float(az))
    r_norm = (90 - alt) / 90.0
    nx = -r_norm * math.sin(az)
    ny = -r_norm * math.cos(az)
    return (nx, ny)


def load_star_catalog(filename: str) -> list[dict[str, object]]:
    """Loads data: load star catalog."""
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


def update_star_positions(
    star_catalog: list[dict[str, object]], lat: float, lon: float, time_obj: Time
) -> tuple[list[dict[str, object]], EarthLocation]:
    """Updates update star positions."""
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    a = angle_below_horizon(VIEWER_LOCATION_HEIGHT)

    visible_stars = []
    for i, star in enumerate(star_catalog):
        coord = cast(SkyCoord, star["coord"])
        altaz = coord.transform_to(AltAz(obstime=time_obj, location=location))
        if altaz.alt.deg > -a:
            visible_stars.append(
                StarData(name=star["name"], alt=altaz.alt.deg, az=altaz.az.deg, vmag=star["vmag"], bv=star["bv"])
            )
        if (i + 1) % 500 == 0:
            time.sleep(0)
    return (visible_stars, location)


def get_visible_planets(lat: float, lon: float, astropy_time: Time) -> list[dict[str, object]]:
    """Gets visible planets using a given astropy Time object."""
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(astropy_time)
    planets = skyfield.api.load("de421.bsp")
    observer = planets["earth"] + Topos(latitude_degrees=lat, longitude_degrees=lon)
    a = angle_below_horizon(VIEWER_LOCATION_HEIGHT)

    visible_bodies = []
    for name, symbol in PLANET_SYMBOLS.items():
        planet = planets[name]
        astrometric = observer.at(t).observe(planet).apparent()
        alt, az, _ = astrometric.altaz()
        if alt.degrees > -a:
            phase_angle = None
            if name == "moon":
                phase_angle = moon_phase_angle(observer, t, planets)
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


def moon_phase_angle(observer: Topos, t: skyfield.timelib.Time, planets: Any) -> float:
    """Moon phase angle."""
    moon = planets["moon"]
    sun = planets["sun"]
    e_to_moon = observer.at(t).observe(moon).apparent()
    e_to_sun = observer.at(t).observe(sun).apparent()
    return e_to_moon.separation_from(e_to_sun).degrees


def rotate_points_at_max_gap(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Rotates points to avoid a large azimuth jump in the middle of the line."""
    if len(points) < 2:
        return points
    max_gap = 0
    max_index = -1
    for i, (p1, p2) in enumerate(zip(points, points[1:] + points[:1])):
        delta = (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2
        if delta > max_gap:
            max_gap, max_index = delta, i

    # Make the 'max_index' item the last one of the list
    if max_index == -1 or max_index == len(points) - 1:
        return points[:]
    else:
        return points[max_index + 1 :] + points[: max_index + 1]


def calculate_celestial_equator_points_norm(location: EarthLocation, time: Time) -> list[tuple[float, float]]:
    """Calculates calculate celestial equator points norm."""
    points = []
    for ra_deg in range(0, 360, 5):
        coord = SkyCoord(ra=ra_deg * u.deg, dec=0 * u.deg, frame="icrs")
        altaz = coord.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > -2:
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg)
            points.append((nx, ny))
    return rotate_points_at_max_gap(points)


def calculate_ecliptic_points_norm(location: EarthLocation, time: Time) -> list[tuple[float, float]]:
    """Calculates calculate ecliptic points norm."""
    points = []
    for lon_deg in range(0, 360, 5):
        ecl = SkyCoord(
            lon=lon_deg * u.deg,
            lat=0 * u.deg,
            frame=GeocentricTrueEcliptic(obstime=time),
        )
        icrs = ecl.transform_to("icrs")
        altaz = icrs.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > -2:
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg)
            points.append((nx, ny))
    return rotate_points_at_max_gap(points)


def calc_sun_angle_on_moon(moon_altaz: tuple[float, float], sun_altaz: tuple[float, float]) -> float:
    """
    画面（北上、東左）上での太陽方向角度を返す（ラジアン、y上方向=0、時計回り正）
    """
    m_alt, m_az = moon_altaz
    s_alt, s_az = sun_altaz

    # 月の位置を中心とした天球座標上で太陽の方向ベクトル
    # 方位角の差分（azは北=0,東=90,南=180,西=270）
    # 画面座標系ではy+が北（上）、x-が東（左）

    # 太陽が月に対してどちらの方位にあるか
    d_az = math.radians(s_az - m_az)
    d_alt = math.radians(s_alt - m_alt)

    # 画面上の相対的な太陽方向ベクトル（x左=東、y上=北）
    dx = -math.sin(d_az) * math.cos(math.radians(s_alt))
    dy = math.sin(d_alt)

    # (dx, dy)から北を上にした画面での角度（北=0, 東=90, 南=180, 西=270、時計回り正）
    angle = math.atan2(dx, dy)
    return angle


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
    phase_angle: float | None = None  # moonのみ


@dataclass
class SkyData:
    """Container for all calculated sky data for a specific time and location."""

    location: tuple[float, float]  # (lat, lon)
    time: Time
    planets: list[dict[str, PlanetBody]]
    stars: list[dict[str, StarData]]
    celestial_equator_points: list[tuple[float, float]]
    ecliptic_points: list[tuple[float, float]]
    timezone_name: str
    city_name: str


# --- Drawing Functions ---


def to_screen_xy(nx: float, ny: float, center: QPoint, radius: int) -> QPointF:
    """Converts normalized coordinates to screen coordinates."""
    return QPointF(center.x() + nx * radius, center.y() + ny * radius)


def find_highlighted_object(
    sky_data: SkyData | None, mouse_pos: QPoint, center: QPoint, radius: int
) -> tuple[dict[str, object], QPointF] | None:
    """Finds the celestial object closest to the mouse cursor."""
    min_dist = 30**2  # Use squared distance to avoid sqrt
    highlighted_object = None

    if not sky_data:
        return None

    all_objects = sky_data.stars + sky_data.planets
    for obj in all_objects:
        pos = to_screen_xy(*altaz_to_normalized_xy(obj.alt, obj.az), center, radius)
        dist_sq = (mouse_pos.x() - pos.x()) ** 2 + (mouse_pos.y() - pos.y()) ** 2
        if dist_sq < min_dist:
            min_dist = dist_sq
            highlighted_object = (obj, pos)
    return highlighted_object


def draw_horizon(painter: QPainter, center: QPoint, radius: int, text_font: QFont):
    """Draws the horizon circle and direction labels."""
    painter.setPen(QPen(HORIZON_LINE_COLOR, 1))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawEllipse(center, radius, radius)

    painter.setPen(TEXT_COLOR)
    painter.setFont(text_font)
    directions = {
        "N": 0,
        "E": 90,
        "S": 180,
        "W": 270,
        "NE": 45,
        "SE": 135,
        "SW": 225,
        "NW": 315,
    }
    for label, angle in directions.items():
        nx, ny = altaz_to_normalized_xy(-3, angle)  # Slightly outside the circle
        pos = to_screen_xy(nx, ny, center, radius)
        painter.drawText(pos, label)


def draw_celestial_lines(painter: QPainter, center: QPoint, radius: int, sky_data: SkyData):
    """Draws the celestial equator and ecliptic lines."""
    if len(sky_data.celestial_equator_points) >= 2:
        points = [to_screen_xy(nx, ny, center, radius) for nx, ny in sky_data.celestial_equator_points]
        poly = QPolygonF(points)
        painter.setPen(QPen(CELESTIAL_EQUATOR_COLOR, 2, Qt.PenStyle.DashLine))
        painter.drawPolyline(poly)

    if len(sky_data.ecliptic_points) >= 2:
        points = [to_screen_xy(nx, ny, center, radius) for nx, ny in sky_data.ecliptic_points]
        poly = QPolygonF(points)
        painter.setPen(QPen(ECLIPTIC_COLOR, 1, Qt.PenStyle.SolidLine))
        painter.drawPolyline(poly)


def draw_stars(painter: QPainter, center: QPoint, radius: int, sky_data: SkyData, star_base_radius: float):
    """Draws the stars."""
    def mag_to_size(vmag: float) -> float:
        return max(0.1, star_base_radius * 10 ** (-0.2 * vmag)) * radius / 500

    painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
    for star in sky_data.stars:
        pos = to_screen_xy(*altaz_to_normalized_xy(star.alt, star.az), center, radius)
        color = bv_to_color(star.bv)
        siz = mag_to_size(star.vmag)

        # サイズが4.0未満の星は、最小サイズ(2x2)で描画する
        if siz < 4.0:
            alpha_value = min(1.0, max(0.1, siz/4.0))
            color.setAlphaF(alpha_value)

            painter.fillRect(QRectF(pos.x() - 1, pos.y() - 1, 2, 2), color)
        else:  # サイズが4.0以上の大きな星は、グラデーション付きの円で描画
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


def draw_cross_gauge(painter: QPainter, color: QColor, center: QPointF):
    """Draws a cross gauge symbol."""
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
    sun_altaz: tuple[float, float] | None = None,
    moon_altaz: tuple[float, float] | None = None,
    opacity: float = 1.0,
):
    img_size = int(radius * 2)
    if img_size < 5:
        img_size = 5

    view_dir = np.array([0, 0, 1])
    phase_angle_rad = math.radians(phase_angle_deg)
    sun_dir = np.array([np.sin(phase_angle_rad), 0, -np.cos(phase_angle_rad)])
    sun_dir /= np.linalg.norm(sun_dir)
    moon_img_pil = render_moon_phase_img(img_size, sun_dir, view_dir)

    rotate_deg = 0
    if sun_altaz is not None and moon_altaz is not None:
        angle = calc_sun_angle_on_moon(moon_altaz, sun_altaz)
        rotate_deg = -math.degrees(angle) - 90

    moon_img_pil = moon_img_pil.rotate(rotate_deg, resample=Image.Resampling.BICUBIC, expand=False)

    pixmap = pil2qpixmap(moon_img_pil)
    target_rect = QRectF(center.x() - img_size / 2, center.y() - img_size / 2, img_size, img_size)
    painter.save()
    painter.setOpacity(opacity)
    painter.drawPixmap(target_rect, pixmap, QRectF(0, 0, img_size, img_size))
    painter.restore()


def draw_planets(
    painter: QPainter, center: QPoint, radius: int, sky_data: SkyData, enlarge_moon: bool, emoji_font: QFont
):
    """Draws the planets, Sun, and Moon."""
    sun_altaz = None
    moon_altaz = None
    for body in sky_data.planets:
        if body.name == "sun":
            sun_altaz = (body.alt, body.az)
        if body.name == "moon":
            moon_altaz = (body.alt, body.az)

    for body in sky_data.planets:
        pos = to_screen_xy(*altaz_to_normalized_xy(body.alt, body.az), center, radius)
        if body.name == "sun":
            draw_cross_gauge(painter, TEXT_COLOR, pos)
        elif body.name == "moon":
            moon_radius = (0.5 if not enlarge_moon else 1.5) / 2 * (radius / 90.0)
            draw_moon(
                painter,
                pos,
                moon_radius,
                body.phase_angle if body.phase_angle is not None else 0.0,
                sun_altaz=sun_altaz,
                moon_altaz=moon_altaz,
                opacity=1.0 if not enlarge_moon else 0.7,
            )
            draw_cross_gauge(painter, TEXT_COLOR, pos)
        else:
            painter.setFont(emoji_font)
            painter.setPen(TEXT_COLOR)
            painter.drawText(pos, body.symbol)


def draw_overlay_text(
    painter: QPainter,
    sky_data: SkyData,
    highlighted_object: tuple[dict[str, object], QPointF] | None,
    text_font: QFont,
):
    """Draws time info and the label for the highlighted object."""

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

    city_name_text = sky_data.city_name.replace("/", " - ").title()
    painter.drawText(QPoint(10, 40), city_name_text)

    if highlighted_object:
        obj, pos = highlighted_object
        painter.setPen(QPen(TEXT_COLOR, 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(pos, 10, 10)

        name = obj.name or ""
        painter.setPen(TEXT_COLOR)
        painter.drawText(QPointF(pos.x() + 15, pos.y() - 15), str(name))


class SkyWidget(QWidget):
    """The main widget for drawing the sky chart."""

    def __init__(self, parent: QObject | None = None):
        super().__init__(parent)
        self.enlarge_moon = False
        self.star_base_radius = 7.0
        self.sky_data: SkyData | None = None
        self.mouse_pos = QPoint()

        # Load emoji font
        font_id = QFontDatabase.addApplicationFont(EMOJI_FONT_PATH)
        font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
        self.emoji_font = QFont(font_family, TEXT_FONT_SIZE + 4)
        self.text_font = QFont("Arial", TEXT_FONT_SIZE)

        self.setMouseTracking(True)
        self.setMinimumSize(400, 400)

    def set_star_base_radius(self, star_base_radius: float):
        self.star_base_radius = star_base_radius

    def set_enlarge_moon(self, enlarge_moon: bool):
        self.enlarge_moon = enlarge_moon

    def set_sky_data(self, data: SkyData):
        """Receives new sky data and triggers a repaint."""
        self.sky_data = data
        self.update()  # Schedule a repaint

    def get_screen_geometry(self) -> tuple[QPoint, int]:
        """Calculates the center and radius of the sky circle."""
        width = self.width()
        height = self.height()
        center = QPoint(width // 2, height // 2)
        s1, s2 = (width, height) if width > height else (height, width)
        radius = max(50, int(s1 * 0.7 + s2 * 0.3) // 2 - 10)
        return center, radius

    def paintEvent(self, event: QObject | None):
        """The main drawing method, called whenever the widget needs to be repainted."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), Qt.GlobalColor.black)

        if not self.sky_data:
            painter.setPen(Qt.GlobalColor.white)
            painter.setFont(QFont("Arial", 16))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Loading sky data...")
            return

        center, radius = self.get_screen_geometry()

        # Draw horizon and direction labels
        draw_horizon(painter, center, radius, self.text_font)

        # Draw celestial lines
        draw_celestial_lines(painter, center, radius, self.sky_data)

        # Determine highlighted object
        highlighted_object = find_highlighted_object(self.sky_data, self.mouse_pos, center, radius)

        # Draw stars and planets
        draw_stars(painter, center, radius, self.sky_data, self.star_base_radius)
        draw_planets(painter, center, radius, self.sky_data, self.enlarge_moon, self.emoji_font)

        # Draw info text and highlighted object label
        draw_overlay_text(painter, self.sky_data, highlighted_object, self.text_font)

    def mouseMoveEvent(self, event: QObject | None):
        """Tracks the mouse position and triggers a repaint to update the highlight."""
        if event:
            self.mouse_pos = event.pos()
        self.update()


class MainWindow(QMainWindow):
    """The main application window."""

    data_updated = pyqtSignal(object)
    initial_data_loaded = pyqtSignal()

    def __init__(
        self,
        city_name: str,
        city_data: tuple[float, float, str],
        star_catalog: list[dict[str, object]],
        delta_t: timedelta,
        enlarge_moon: bool,
        star_base_radius: float,
    ):
        super().__init__()
        self.city_name = city_name
        self.lat, self.lon, self.tz_name = city_data
        self.star_catalog = star_catalog
        self.delta_t = delta_t
        self.enlarge_moon = enlarge_moon
        self.star_base_radius = star_base_radius

        self.setWindowTitle(f"Zenith Star View - {self.city_name.replace('/', ' - ').title()}")
        self.setGeometry(100, 100, 800, 800)

        self.sky_widget = SkyWidget(self)
        self.sky_widget.set_enlarge_moon(self.enlarge_moon)
        self.sky_widget.set_star_base_radius(self.star_base_radius)
        self.setCentralWidget(self.sky_widget)

        self.data_updated.connect(self.on_data_updated)

        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.start_background_update)

        self.start_background_update(is_initial_load=True)

    def on_data_updated(self, sky_data: SkyData):
        """Slot to handle new data from the background thread."""
        self.sky_widget.set_sky_data(sky_data)
        if not self.update_timer.isActive():
            self.update_timer.start(5 * 60 * 1000)
            self.initial_data_loaded.emit()

    def update_sky_data_in_background(self):
        """Performs the heavy calculations for sky data."""
        try:
            now = datetime.now(timezone.utc) + self.delta_t
            time_obj = Time(now)
            stars, loc = update_star_positions(self.star_catalog, self.lat, self.lon, time_obj)
            planets = get_visible_planets(self.lat, self.lon, time_obj)
            celestial_equator_points = calculate_celestial_equator_points_norm(loc, time_obj)
            ecliptic_points = calculate_ecliptic_points_norm(loc, time_obj)

            sky_data = SkyData(
                (self.lat, self.lon),
                time_obj,
                planets,
                stars,
                celestial_equator_points,
                ecliptic_points,
                timezone_name=self.tz_name,
                city_name=self.city_name,
            )
            self.data_updated.emit(sky_data)
        except Exception as e:
            print(f"Error in background update thread: {e}", file=sys.stderr)
            import traceback

            traceback.print_exc()

    def start_background_update(self, is_initial_load: bool = False):
        """Starts a new thread to calculate sky data."""
        if is_initial_load:
            print("Performing initial data load...")
        else:
            print("Updating sky data...")
        thread = threading.Thread(target=self.update_sky_data_in_background)
        thread.daemon = True
        thread.start()

    def keyPressEvent(self, event: QObject | None):
        """Handles key presses for fullscreen toggle."""
        if event and event.key() == Qt.Key.Key_F11:
            if self.isFullScreen():
                self.showNormal()
            else:
                self.showFullScreen()
        elif event and event.key() == Qt.Key.Key_Escape:
            if self.isFullScreen():
                self.showNormal()
        else:
            super().keyPressEvent(event)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    parser.add_argument("city", type=str, default="Tokyo", help="City name (default: Tokyo)")
    parser.add_argument(
        "-H", "--hours", type=float, default=0, help="Number of hours to add to current time (default: 0)"
    )
    parser.add_argument(
        "-D", "--days", type=float, default=0, help="Number of days to add to current time (default: 0)"
    )
    parser.add_argument(
        "-m",
        "--enlarge-moon",
        action="store_true",
        help="Show the moon in 3x size.",
    )
    parser.add_argument(
        "-s",
        "--star-base-radius",
        type=float,
        default=10.0,
        help="Base size of stars (default: 10.0)"
    )
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
    main_win = MainWindow(
        city, 
        (lat, lon, tz_name), 
        star_catalog, 
        delta_t, 
        enlarge_moon=args.enlarge_moon,
        star_base_radius=args.star_base_radius)

    # When the initial data is loaded, show the main window and close the splash screen
    main_win.initial_data_loaded.connect(main_win.show)
    main_win.initial_data_loaded.connect(splash.close)

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
