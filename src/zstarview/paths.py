from dataclasses import dataclass
import os
import os.path

from appdirs import user_cache_dir, user_log_dir

# Base directory of this package
_dir = os.path.dirname(os.path.abspath(__file__))

# Application identifiers
APP_ID = "zstarview"
APP_AUTHOR = "tos-kamiya"
# Human-readable display name
APP_DISPLAY_NAME = "Zenith Star View"

# Data file paths
TEXT_FONT_PATH = os.path.join(_dir, "data", "Noto_Sans", "NotoSans-VariableFont_wdth,wght.ttf")
CITY_COORD_FILE = os.path.join(_dir, "data", "cities1000.txt")
CITY_ADMIN1_CODES_FILE = os.path.join(_dir, "data", "admin1CodesASCII.txt")
STARS_CSV_FILE = os.path.join(_dir, "data", "stars", "stars_base.csv")
DSO_CSV_FILE = os.path.join(_dir, "data", "dso.csv")
TOWER_VIEWPOINTS_FILE = os.path.join(_dir, "data", "viewpoints", "tower_viewpoints.json")
MOUNTAIN_VIEWPOINTS_FILE = os.path.join(_dir, "data", "viewpoints", "mountain_viewpoints.json")
APP_ICON_FILE = os.path.join(_dir, "data", "icon-256.png")
SKYSCRAPER_TILES_FILE = os.path.join(_dir, "data", "skyscraper_tiles_z14.json")
CACHE_PATH = user_cache_dir(appname=APP_ID, appauthor=APP_AUTHOR)
LOG_PATH = user_log_dir(appname=APP_ID, appauthor=APP_AUTHOR)
OVERTURE_DERIVED_ROOT_DIR = os.path.join(CACHE_PATH, "overture_buildings")
OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR = os.path.join(CACHE_PATH, "overture_skyscrapers")
AIRCRAFT_CACHE_ROOT_DIR = os.path.join(CACHE_PATH, "aircraft", "opensky")
SATELLITE_CACHE_ROOT_DIR = os.path.join(CACHE_PATH, "satellites", "celestrak")

# Window UI
GUI_MENU_TEXT_COLOR = (128, 128, 128)
GUI_BUTTON_SIZE = 30
WINDOW_WIDTH = 600
WINDOW_HEIGHT = 600

# UI constants
TEXT_FONT_SIZE = 11
STATUS_LINE_FONT_SIZE = 8


@dataclass(frozen=True, slots=True)
class TextStyle:
    text: tuple[int, ...]
    outline: tuple[int, ...]


NIGHT_TEXT_STYLE = TextStyle(
    text=(180, 180, 180),
    outline=(0, 0, 0, 76),
)
TRANSPARENT_TEXT_STYLE = TextStyle(
    text=(246, 249, 255),
    outline=(2, 2, 3, 212),
)
WHITE_TEXT_STYLE = TextStyle(
    text=(220, 132, 66),
    outline=(156, 92, 42, 132),
)
DAY_TEXT_STYLE = TextStyle(
    text=(220, 132, 66),
    outline=(148, 86, 38, 126),
)
BLACK_TEXT_STYLE = TextStyle(
    text=(246, 249, 255),
    outline=(2, 2, 3, 236),
)

NIGHT_STATUS_LINE_STYLE = TextStyle(
    text=(190, 190, 160),
    outline=(0, 0, 0, 76),
)
TRANSPARENT_STATUS_LINE_STYLE = TextStyle(
    text=(226, 228, 234),
    outline=(2, 2, 3, 208),
)
WHITE_STATUS_LINE_STYLE = TextStyle(
    text=(208, 124, 62),
    outline=(150, 88, 40, 136),
)
DAY_STATUS_LINE_STYLE = TextStyle(
    text=(208, 124, 62),
    outline=(142, 82, 36, 130),
)
BLACK_STATUS_LINE_STYLE = TextStyle(
    text=(255, 220, 220),
    outline=(2, 2, 3, 236),
)


TEXT_STYLES_BY_PRESET = {
    "night": NIGHT_TEXT_STYLE,
    "transparent": TRANSPARENT_TEXT_STYLE,
    "white": WHITE_TEXT_STYLE,
    "day": DAY_TEXT_STYLE,
    "black": BLACK_TEXT_STYLE,
}
STATUS_LINE_STYLES_BY_PRESET = {
    "night": NIGHT_STATUS_LINE_STYLE,
    "transparent": TRANSPARENT_STATUS_LINE_STYLE,
    "white": WHITE_STATUS_LINE_STYLE,
    "day": DAY_STATUS_LINE_STYLE,
    "black": BLACK_STATUS_LINE_STYLE,
}

HORIZON_LINE_COLOR = (72, 127, 71)
TERRAIN_HORIZON_LINE_COLOR = (93, 76, 33)
URBAN_OUTLINE_LAYER_LINE_COLOR = (246, 249, 255)
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
    "sun": "Sun",
    "moon": "Moon",
    "mercury": "Mercury",
    "venus": "Venus",
    "mars": "Mars",
    "jupiter": "Jupiter",
    "saturn": "Saturn",
    "uranus": "Uranus",
    "neptune": "Neptune",
    "pluto": "Pluto",
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
EPHEMERIS_FILENAME = "de442s.bsp"
EPHEMERIS_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de442s.bsp"
