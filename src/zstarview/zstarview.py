# -*- coding: utf-8 -*-
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import math
from multiprocessing import Process, Pipe
import os.path
import pickle
import sys
import threading
import time
from zoneinfo import ZoneInfo

# --- PyQt5 Imports ---
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
)
from PyQt5.QtCore import Qt, QPoint, QPointF, QRectF, QSize, QTimer, pyqtSignal

# --- Astropy and Skyfield Imports (unchanged) ---
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, GeocentricTrueEcliptic
from astropy.time import Time
import astropy.units as u
from skyfield.api import Topos
import skyfield.api

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
            cols = line.strip().split("\t")
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


def mag_to_radius(vmag: float) -> float:
    """Mag to radius."""
    base_radius = 7.0
    return max(0.1, base_radius * 10 ** (-0.2 * vmag))


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
    visible_stars = []
    for i, star in enumerate(star_catalog):
        altaz = star["coord"].transform_to(AltAz(obstime=time_obj, location=location))
        if altaz.alt.deg > 0:
            visible_stars.append(
                {
                    "alt": altaz.alt.deg,
                    "az": altaz.az.deg,
                    "vmag": star["vmag"],
                    "bv": star["bv"],
                    "name": star["name"],
                }
            )
        if (i + 1) % 500 == 0:
            time.sleep(0)
    return (visible_stars, location)


def estimate_magnitude(body_name, observer, t: Time, planets: dict) -> float | None:
    """Estimate magnitude."""
    sun = planets["sun"]
    planet = planets[body_name]
    obs_to_planet = observer.at(t).observe(planet).apparent()
    r = planet.at(t).observe(sun).apparent().distance().au
    delta = obs_to_planet.distance().au
    if body_name == "moon":  # Skyfield provides moon magnitude
        return planet.magnitude
    phase_angle = obs_to_planet.separation_from(planet.at(t).observe(sun).apparent())
    i = phase_angle.degrees
    mag = 5 * math.log10(r * delta)
    if body_name == "mercury":
        mag += -0.613 + 6.328e-2 * i - 1.738e-3 * i**2 + 2.956e-5 * i**3 - 2.814e-7 * i**4 + 1.058e-9 * i**5
    elif body_name == "venus":
        mag += -4.47 + 1.03e-2 * i + 3.69e-4 * i**2 - 2.81e-6 * i**3 + 8.94e-9 * i**4
    elif body_name == "mars":
        mag += -1.52 + 0.016 * i
    elif body_name == "jupiter barycenter":
        mag += -9.40 + 0.005 * i
    elif body_name == "saturn barycenter":
        # Simplified, ring effects ignored
        mag += -8.88
    else:
        return None
    return mag


def get_visible_planets(lat: float, lon: float, astropy_time: Time) -> list[dict[str, object]]:
    """Gets visible planets using a given astropy Time object."""
    ts = skyfield.api.load.timescale()
    t = ts.from_astropy(astropy_time)
    planets = skyfield.api.load("de421.bsp")
    observer = planets["earth"] + Topos(latitude_degrees=lat, longitude_degrees=lon)

    visible_bodies = []
    for name, symbol in PLANET_SYMBOLS.items():
        planet = planets[name]
        astrometric = observer.at(t).observe(planet).apparent()
        alt, az, _ = astrometric.altaz()
        if alt.degrees > 0:
            mag = estimate_magnitude(name, observer, t, planets)
            body = {
                "alt": alt.degrees,
                "az": az.degrees,
                "symbol": symbol,
                "mag": mag,
                "name": name,
            }
            if name == "moon":
                body["phase_angle"] = moon_phase_angle(observer, t, planets)
            visible_bodies.append(body)
    return visible_bodies


def moon_phase_angle(observer, t: Time, planets: list[dict[str, object]]) -> float:
    """Moon phase angle."""
    moon = planets["moon"]
    sun = planets["sun"]
    e_to_moon = observer.at(t).observe(moon).apparent()
    e_to_sun = observer.at(t).observe(sun).apparent()
    return e_to_moon.separation_from(e_to_sun).degrees


def rotate_points_at_max_azimuth_gap(
    points: list[tuple[float, float]], azimuths: list[float]
) -> list[tuple[float, float]]:
    """Rotates points to avoid a large azimuth jump in the middle of the line."""
    if len(points) < 2:
        return points
    max_gap, max_index = -1, 0
    for i in range(1, len(azimuths)):
        delta = abs(azimuths[i] - azimuths[i - 1])
        delta = min(delta, 360 - delta)
        if delta > max_gap:
            max_gap, max_index = delta, i
    return points[max_index:] + points[:max_index]


def calculate_celestial_equator_points_norm(location: EarthLocation, time: Time) -> list[tuple[float, float]]:
    """Calculates calculate celestial equator points norm."""
    points, azimuths = [], []
    for ra_deg in range(0, 360, 5):
        coord = SkyCoord(ra=ra_deg * u.deg, dec=0 * u.deg, frame="icrs")
        altaz = coord.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > 0:
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg)
            points.append((nx, ny))
            azimuths.append(altaz.az.deg)
    return rotate_points_at_max_azimuth_gap(points, azimuths)


def calculate_ecliptic_points_norm(location: EarthLocation, time: Time) -> list[tuple[float, float]]:
    """Calculates calculate ecliptic points norm."""
    points, azimuths = [], []
    for lon_deg in range(0, 360, 5):
        ecl = SkyCoord(
            lon=lon_deg * u.deg,
            lat=0 * u.deg,
            frame=GeocentricTrueEcliptic(obstime=time),
        )
        icrs = ecl.transform_to("icrs")
        altaz = icrs.transform_to(AltAz(obstime=time, location=location))
        if altaz.alt.deg > 0:
            nx, ny = altaz_to_normalized_xy(altaz.alt.deg, altaz.az.deg)
            points.append((nx, ny))
            azimuths.append(altaz.az.deg)
    return rotate_points_at_max_azimuth_gap(points, azimuths)


@dataclass
class SkyData:
    """Container for all calculated sky data for a specific time and location."""

    location: tuple[float, float]  # (lat, lon)
    time: Time
    planets: list[dict[str, object]]
    stars: list[dict[str, object]]
    celestial_equator_points: list[tuple[float, float]]
    ecliptic_points: list[tuple[float, float]]
    timezone_name: str
    city_name: str


# --- PyQt5 UI Classes ---


class SkyWidget(QWidget):
    """The main widget for drawing the sky chart."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.sky_data = None
        self.highlighted_object = None
        self.mouse_pos = QPoint()

        # Load emoji font
        font_id = QFontDatabase.addApplicationFont(EMOJI_FONT_PATH)
        font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
        self.emoji_font = QFont(font_family, TEXT_FONT_SIZE + 4)
        self.text_font = QFont("Arial", TEXT_FONT_SIZE)

        self.setMouseTracking(True)
        self.setMinimumSize(400, 400)

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

    def to_screen_xy(self, nx: float, ny: float, center: QPoint, radius: int) -> QPointF:
        """Converts normalized coordinates to screen coordinates."""
        return QPointF(center.x() + nx * radius, center.y() + ny * radius)

    def paintEvent(self, event):
        """The main drawing method, called whenever the widget needs to be repainted."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.fillRect(self.rect(), Qt.black)

        if not self.sky_data:
            painter.setPen(Qt.white)
            painter.setFont(QFont("Arial", 16))
            painter.drawText(self.rect(), Qt.AlignCenter, "Loading sky data...")
            return

        center, radius = self.get_screen_geometry()

        # Draw horizon and direction labels
        self.draw_horizon(painter, center, radius)

        # Draw celestial lines
        self.draw_celestial_lines(painter, center, radius)

        # Determine highlighted object
        self.update_highlight(center, radius)

        # Draw stars and planets
        self.draw_stars(painter, center, radius)
        self.draw_planets(painter, center, radius)

        # Draw info text and highlighted object label
        self.draw_overlay_text(painter, center, radius)

    def draw_horizon(self, painter: QPainter, center: QPoint, radius: int):
        """Draws the horizon circle and direction labels."""
        painter.setPen(QPen(HORIZON_LINE_COLOR, 1))
        painter.setBrush(Qt.NoBrush)
        painter.drawEllipse(center, radius, radius)

        painter.setPen(TEXT_COLOR)
        painter.setFont(self.text_font)
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
            pos = self.to_screen_xy(nx, ny, center, radius)
            painter.drawText(pos, label)

    def draw_celestial_lines(self, painter: QPainter, center: QPoint, radius: int):
        """Draws the celestial equator and ecliptic lines."""
        if len(self.sky_data.celestial_equator_points) >= 2:
            points = [self.to_screen_xy(nx, ny, center, radius) for nx, ny in self.sky_data.celestial_equator_points]
            poly = QPolygonF(points)
            painter.setPen(QPen(CELESTIAL_EQUATOR_COLOR, 1, Qt.DashLine))
            painter.drawPolyline(poly)

        if len(self.sky_data.ecliptic_points) >= 2:
            points = [self.to_screen_xy(nx, ny, center, radius) for nx, ny in self.sky_data.ecliptic_points]
            poly = QPolygonF(points)
            painter.setPen(QPen(ECLIPTIC_COLOR, 1, Qt.SolidLine))
            painter.drawPolyline(poly)

    def draw_stars(self, painter: QPainter, center: QPoint, radius: int):
        """Draws the stars."""
        painter.setCompositionMode(QPainter.CompositionMode_Plus)
        for star in self.sky_data.stars:
            pos = self.to_screen_xy(*altaz_to_normalized_xy(star["alt"], star["az"]), center, radius)
            color = bv_to_color(star["bv"])
            rad = mag_to_radius(star["vmag"]) * radius / 2000

            # 半径が1.0ピクセル未満の星は、最小サイズ(2x2)で描画する
            if rad < 1.0:
                alpha_value = min(1.0, max(0.1, rad * 1.5))
                color.setAlphaF(alpha_value)

                painter.fillRect(QRectF(pos.x() - 1, pos.y() - 1, 2, 2), color)
            else:  # 半径が1.0ピクセル以上の大きな星は、グラデーション付きの円で描画
                gradient = QRadialGradient(pos, rad)
                gradient.setColorAt(0, color)

                color_transparent = QColor(color)
                color_transparent.setAlpha(0)
                gradient.setColorAt(1, color_transparent)

                painter.setBrush(QBrush(gradient))
                painter.setPen(Qt.NoPen)
                painter.drawEllipse(pos, rad, rad)

        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)  # Reset mode

    def draw_planets(self, painter: QPainter, center: QPoint, radius: int):
        """Draws the planets, Sun, and Moon."""
        for body in self.sky_data.planets:
            pos = self.to_screen_xy(*altaz_to_normalized_xy(body["alt"], body["az"]), center, radius)
            name = body.get("name")

            if name == "sun":
                self.draw_cross_gauge(painter, TEXT_COLOR, pos)
            elif name == "moon":
                moon_radius = 0.5 / 2 * (radius / 90.0)
                self.draw_moon(painter, pos, moon_radius, body["phase_angle"])
                self.draw_cross_gauge(painter, TEXT_COLOR, pos)
            else:
                painter.setFont(self.emoji_font)
                painter.setPen(TEXT_COLOR)
                painter.drawText(pos, body["symbol"])

    def draw_moon(self, painter: QPainter, center: QPointF, radius: float, phase_angle_deg: float):
        """Draws the Moon with its correct phase."""
        painter.setPen(Qt.NoPen)
        painter.setBrush(QColor(230, 230, 230))
        painter.drawEllipse(center, radius, radius)  # Draw the illuminated part

        # Draw the dark part (shadow)
        phase_rad = math.radians(phase_angle_deg - 90)  # Adjust for calculation
        x_offset = radius * math.cos(phase_rad)

        painter.setBrush(Qt.black)
        rect = QRectF(center.x() - x_offset, center.y() - radius, 2 * x_offset, 2 * radius)
        painter.drawChord(
            QRectF(center.x() - radius, center.y() - radius, 2 * radius, 2 * radius),
            90 * 16,
            -180 * 16,
        )

    def draw_cross_gauge(self, painter: QPainter, color: QColor, center: QPointF):
        """Draws a cross gauge symbol."""
        cross_outer_len, cross_inner_len = 15, 4
        x, y = center.x(), center.y()
        painter.setPen(QPen(color, 1))
        painter.drawLine(QPointF(x - cross_outer_len, y), QPointF(x - cross_inner_len, y))
        painter.drawLine(QPointF(x + cross_inner_len, y), QPointF(x + cross_outer_len, y))
        painter.drawLine(QPointF(x, y - cross_outer_len), QPointF(x, y - cross_inner_len))
        painter.drawLine(QPointF(x, y + cross_inner_len), QPointF(x, y + cross_outer_len))

    def update_highlight(self, center: QPoint, radius: int):
        """Finds the celestial object closest to the mouse cursor."""
        min_dist = 30**2  # Use squared distance to avoid sqrt
        self.highlighted_object = None

        all_objects = self.sky_data.stars + self.sky_data.planets
        for obj in all_objects:
            pos = self.to_screen_xy(*altaz_to_normalized_xy(obj["alt"], obj["az"]), center, radius)
            dist_sq = (self.mouse_pos.x() - pos.x()) ** 2 + (self.mouse_pos.y() - pos.y()) ** 2
            if dist_sq < min_dist:
                min_dist = dist_sq
                self.highlighted_object = (obj, pos)

    def draw_overlay_text(self, painter: QPainter, center: QPoint, radius: int):
        """Draws time info and the label for the highlighted object."""

        utc_time = self.sky_data.time
        tz_name = self.sky_data.timezone_name
        time_text = ""
        try:
            local_tz = ZoneInfo(tz_name)
            local_dt = utc_time.to_datetime(timezone=local_tz)
            time_text = local_dt.strftime("%Y-%m-%d %H:%M:%S %Z")
        except Exception:
            time_text = utc_time.to_datetime().strftime("%Y-%m-%d %H:%M:%S UTC")

        painter.setPen(QColor(180, 180, 180))
        painter.setFont(self.text_font)
        painter.drawText(QPoint(10, 20), time_text)

        city_name_text = self.sky_data.city_name.replace("/", " - ").title()
        painter.drawText(QPoint(10, 40), city_name_text)

        if self.highlighted_object:
            obj, pos = self.highlighted_object
            painter.setPen(QPen(TEXT_COLOR, 2))
            painter.setBrush(Qt.NoBrush)
            painter.drawEllipse(pos, 10, 10)

            name = obj.get("name") or ""
            painter.setPen(TEXT_COLOR)
            painter.drawText(QPointF(pos.x() + 15, pos.y() - 15), name)

    def mouseMoveEvent(self, event):
        """Tracks the mouse position and triggers a repaint to update the highlight."""
        self.mouse_pos = event.pos()
        self.update()


class MainWindow(QMainWindow):
    """The main application window."""

    data_updated = pyqtSignal(object)
    initial_data_loaded = pyqtSignal()

    def __init__(self, city_name: str, city_data: tuple, star_catalog: list):
        super().__init__()
        self.city_name = city_name
        # 緯度、経度、タイムゾーン名をアンパック
        self.lat, self.lon, self.tz_name = city_data
        self.star_catalog = star_catalog

        self.setWindowTitle(f"Zenith Star View - {self.city_name.replace('/', ' - ').title()}")
        self.setGeometry(100, 100, 800, 800)

        self.sky_widget = SkyWidget(self)
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
            now = datetime.now(timezone.utc)
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

    def start_background_update(self, is_initial_load=False):
        """Starts a new thread to calculate sky data."""
        if is_initial_load:
            print("Performing initial data load...")
        else:
            print("Updating sky data...")
        thread = threading.Thread(target=self.update_sky_data_in_background)
        thread.daemon = True
        thread.start()

    def keyPressEvent(self, event):
        """Handles key presses for fullscreen toggle."""
        if event.key() == Qt.Key_F11:
            if self.isFullScreen():
                self.showNormal()
            else:
                self.showFullScreen()
        elif event.key() == Qt.Key_Escape:
            if self.isFullScreen():
                self.showNormal()
        else:
            super().keyPressEvent(event)


def main():
    """Main entry point for the star sky visualizer."""
    app = QApplication(sys.argv)

    # Show a splash screen while loading initial data
    pixmap = QPixmap(400, 200)
    pixmap.fill(Qt.black)
    splash = QSplashScreen(pixmap, Qt.WindowStaysOnTopHint)
    splash.show()
    splash.showMessage("Loading city and star data...", Qt.AlignCenter, Qt.white)
    app.processEvents()

    try:
        city_table = load_city_coords(CITY_COORD_FILE)
    except FileNotFoundError:
        splash.showMessage("Error: cities1000.txt not found.", Qt.AlignCenter, Qt.red)
        time.sleep(3)
        return

    city = "Tokyo"
    if len(sys.argv) >= 2:
        city = sys.argv[1]
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

    try:
        star_catalog = load_star_catalog(STARS_CSV_FILE)
    except FileNotFoundError:
        splash.showMessage("Error: stars.csv not found.", Qt.AlignCenter, Qt.red)
        time.sleep(3)
        return

    splash.showMessage(f"Calculating sky for {city.title()}...", Qt.AlignCenter, Qt.white)
    app.processEvents()

    # lat, lon = city_table[city] の行を修正
    lat, lon, tz_name = city_table[city]
    main_win = MainWindow(city, (lat, lon, tz_name), star_catalog)

    # When the initial data is loaded, show the main window and close the splash screen
    main_win.initial_data_loaded.connect(lambda: (main_win.show(), splash.finish(main_win)))

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
