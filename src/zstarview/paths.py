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


# Direction labels
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
