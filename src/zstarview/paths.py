import os
import os.path
from PyQt5.QtGui import QColor


# Application identifiers
APP_ID = "zstarview"
APP_AUTHOR = "tos-kamiya"


# Base directory of this package
_dir = os.path.dirname(os.path.abspath(__file__))

# Data file paths
EMOJI_FONT_PATH = os.path.join(_dir, "data", "Noto_Sans_Symbols", "NotoSansSymbols-VariableFont_wght.ttf")
CITY_COORD_FILE = os.path.join(_dir, "data", "cities1000.txt")
STARS_CSV_FILE = os.path.join(_dir, "data", "stars.csv")
APP_ICON_FILE = os.path.join(_dir, "data", "app_icon.ico")


# UI constants
TEXT_COLOR = QColor(120, 120, 120)
TEXT_FONT_SIZE = 14

HORIZON_LINE_COLOR = QColor(40, 50, 40)
CELESTIAL_EQUATOR_COLOR = QColor(40, 40, 40)
ECLIPTIC_COLOR = QColor(80, 60, 0)


# Rendering / FOV
FIELD_OF_VIEW_DEG = 240
ANGLE_BELOW_HORIZON = 2


# Direction labels (16-point compass rose)
DIRECTIONS = {
    "N": 0.0,
    "NNE": 22.5,
    "NE": 45.0,
    "ENE": 67.5,
    "E": 90.0,
    "ESE": 112.5,
    "SE": 135.0,
    "SSE": 157.5,
    "S": 180.0,
    "SSW": 202.5,
    "SW": 225.0,
    "WSW": 247.5,
    "W": 270.0,
    "WNW": 292.5,
    "NW": 315.0,
    "NNW": 337.5,
}


# Symbols for bodies
PLANET_SYMBOLS = {
    "sun": "☀",
    "moon": "🌛",
    "mercury": "☿",
    "venus": "♀",
    "mars": "♂",
    "jupiter barycenter": "♃",
    "saturn barycenter": "♄",
}
