from dataclasses import dataclass
import os
import os.path

from appdirs import user_cache_dir, user_log_dir

# Base directory of this package
_dir = os.path.dirname(os.path.abspath(__file__))

# Application identifiers
APP_ID = "zstarview"
APP_AUTHOR = "tos-kamiya"
# Human‑readable display name
APP_DISPLAY_NAME = "Zenith Star View"

# Data file paths
TEXT_FONT_PATH = os.path.join(_dir, "data", "Noto_Sans", "NotoSans-VariableFont_wdth,wght.ttf")
CITY_COORD_FILE = os.path.join(_dir, "data", "cities1000.txt")
CITY_ADMIN1_CODES_FILE = os.path.join(_dir, "data", "admin1CodesASCII.txt")
STARS_CSV_FILE = os.path.join(_dir, "data", "stars.csv")
APP_ICON_FILE = os.path.join(_dir, "data", "icon-256.png")
CACHE_PATH = user_cache_dir(appname=APP_ID, appauthor=APP_AUTHOR)
LOG_PATH = user_log_dir(appname=APP_ID, appauthor=APP_AUTHOR)

# Window UI
GUI_MENU_TEXT_COLOR = (128, 128, 128)
GUI_BUTTON_SIZE = 30
WINDOW_WIDTH = 600
WINDOW_HEIGHT = 600

# UI constants
TEXT_COLOR = (180, 180, 180)
TEXT_FONT_SIZE = 11
STATUS_LINE_COLOR = (190, 190, 160)
STATUS_LINE_FONT_SIZE = 8

HORIZON_LINE_COLOR = (72, 127, 71)
CELESTIAL_EQUATOR_COLOR = (139, 139, 136)
ECLIPTIC_COLOR = (236, 173, 2)

SKY_UPDATE_INTERVAL = 3 * 60  # seconds
CLOUD_UPDATE_INTERVAL = 10 * 60  # seconds

CLOUD_SHELL_KM = 6371.0 + 5.0  # 5km above Earth's surface

# Rendering / FOV
FIELD_OF_VIEW_DEG = 115
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
}


@dataclass(frozen=True, slots=True)
class HatchConfig:
    tile_w_px: int
    tile_h_px: int
    line_px: int
    strength: int


CLOUD_HATCH_DEFAULT = HatchConfig(20, 19, 8, 255)

# Skyfield ephemeris kernel filename
EPHEMERIS_FILENAME = "de440s.bsp"
