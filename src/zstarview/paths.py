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
STARS_CSV_FILE = os.path.join(_dir, "data", "stars", "stars_base.csv")
DSO_CSV_FILE = os.path.join(_dir, "data", "dso.csv")
TOWER_VIEWPOINTS_FILE = os.path.join(_dir, "data", "viewpoints", "tower_viewpoints.json")
MOUNTAIN_VIEWPOINTS_FILE = os.path.join(_dir, "data", "viewpoints", "mountain_viewpoints.json")
PLATEAU_DERIVED_ROOT_DIR = os.path.join(_dir, "data", "plateau_derived")
APP_ICON_FILE = os.path.join(_dir, "data", "icon-256.png")
CACHE_PATH = user_cache_dir(appname=APP_ID, appauthor=APP_AUTHOR)
LOG_PATH = user_log_dir(appname=APP_ID, appauthor=APP_AUTHOR)
OVERTURE_DERIVED_ROOT_DIR = os.path.join(CACHE_PATH, "overture_buildings")

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
LIGHT_LABEL_COLOR = (246, 249, 255)
TEXT_OUTLINE_COLOR_NIGHT_RGBA = (0, 0, 0, 76)
TEXT_STYLES_BY_PRESET = {
    "night": {
        "text": TEXT_COLOR,
        "outline": TEXT_OUTLINE_COLOR_NIGHT_RGBA,
    },
    "white": {
        "text": (44, 112, 196),
        "outline": (36, 80, 140, 122),
    },
    "day": {
        "text": (34, 106, 192),
        "outline": (30, 74, 136, 116),
    },
    "black": {
        "text": LIGHT_LABEL_COLOR,
        "outline": (2, 2, 3, 236),
    },
}
STATUS_LINE_STYLES_BY_PRESET = {
    "night": {
        "text": STATUS_LINE_COLOR,
        "outline": TEXT_OUTLINE_COLOR_NIGHT_RGBA,
    },
    "white": {
        "text": (42, 102, 184),
        "outline": (34, 72, 128, 126),
    },
    "day": {
        "text": (34, 96, 178),
        "outline": (28, 68, 122, 120),
    },
    "black": {
        "text": (255, 220, 220),
        "outline": (2, 2, 3, 236),
    },
}

HORIZON_LINE_COLOR = (72, 127, 71)
TERRAIN_HORIZON_LINE_COLOR = (93, 76, 33)
URBAN_DEBUG_LAYER_LINE_COLOR = LIGHT_LABEL_COLOR
CELESTIAL_EQUATOR_COLOR = (139, 139, 136)
ECLIPTIC_COLOR = (236, 173, 2)

CLOUD_UPDATE_INTERVAL = 10 * 60  # seconds

CLOUD_SHELL_KM = 6371.0 + 5.0  # 5km above Earth's surface

# Rendering / FOV
FIELD_OF_VIEW_DEG = 100
STAR_FIELD_OF_VIEW_DEG = 100
BACKGROUND_FIELD_OF_VIEW_DEG1 = 90
BACKGROUND_FIELD_OF_VIEW_DEG2 = 105
ASTERISM_CLIP_FIELD_OF_VIEW_DEG = BACKGROUND_FIELD_OF_VIEW_DEG2

# Direction labels (8-point compass rose)
DIRECTIONS = {
    "N": 0.0,
    "NE": 45.0,
    "E": 90.0,
    "SE": 135.0,
    "S": 180.0,
    "SW": 225.0,
    "W": 270.0,
    "NW": 315.0,
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
CLOUD_MISSING_HATCH_DEFAULT = HatchConfig(20, 19, 8, 160)
CLOUD_MISSING_TINT_RGBA = (255, 220, 80, 45)

# Skyfield ephemeris kernel filename
EPHEMERIS_FILENAME = "de440s.bsp"
