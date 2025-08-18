import os
import os.path
from PySide6.QtGui import QColor


# Application identifiers
APP_ID = "zstarview"
APP_AUTHOR = "tos-kamiya"


# Base directory of this package
_dir = os.path.dirname(os.path.abspath(__file__))

# Data file paths
TEXT_FONT_PATH = os.path.join(_dir, "data", "Noto_Sans", "NotoSans-VariableFont_wdth,wght.ttf")
EMOJI_FONT_PATH = os.path.join(_dir, "data", "Noto_Sans_Symbols", "NotoSansSymbols-VariableFont_wght.ttf")
CITY_COORD_FILE = os.path.join(_dir, "data", "cities1000.txt")
STARS_CSV_FILE = os.path.join(_dir, "data", "stars.csv")
APP_ICON_FILE = os.path.join(_dir, "data", "icon-256.png")

# Window UI
GUI_MENU_TEXT_COLOR = "#787878"
GUI_BUTTON_SIZE = 30
WINDOW_WIDTH = 400
WINDOW_HEIGHT = 400

# UI constants
TEXT_COLOR = QColor(120, 120, 120)
TEXT_FONT_SIZE = 10
EMOJI_FONT_SIZE = 18

HORIZON_LINE_COLOR = QColor(45, 70, 45)
CELESTIAL_EQUATOR_COLOR = QColor(60, 60, 60)
ECLIPTIC_COLOR = QColor(90, 65, 0)


# Rendering / FOV
FIELD_OF_VIEW_DEG = 235
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


# Planet labels/ids
PLANET_SYMBOLS = {
    "sun": "☀",
    "moon": "🌛",
    "mercury": "☿",
    "venus": "♀",
    "mars": "♂",
    "jupiter": "♃",
    "saturn": "♄",
    "uranus": "♅",
    "neptune": "♆",
    "pluto": "♇",
    "ceres": "⚳",
    "pallas": "⚴",
    "vesta": "⚵",
}

PLANET_IDS = {
    "sun": 10,
    "moon": 301,
    "mercury": 199,
    "venus": 299,
    "mars": 4,  # BARYCENTER
    "jupiter": 5,  # BARYCENTER
    "saturn": 6,  # BARYCENTER
    "uranus": 7,  # BARYCENTER
    "neptune": 8,  # BARYCENTER
    "pluto": 9,  # BARYCENTER
    "ceres": 1,
    "pallas": 2,
    "vesta": 4,
}

