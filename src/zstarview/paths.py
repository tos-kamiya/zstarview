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
EARTH_GUIDE_LAND_FILE = os.path.join(_dir, "data", "earth_guide_land_110m.json")
CACHE_PATH = user_cache_dir(appname=APP_ID, appauthor=APP_AUTHOR)
LOG_PATH = user_log_dir(appname=APP_ID, appauthor=APP_AUTHOR)
COPERNICUS_DEM_CACHE_DIR = os.path.join(CACHE_PATH, "copernicus-dem")
OVERTURE_DERIVED_ROOT_DIR = os.path.join(CACHE_PATH, "overture_buildings")
OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR = os.path.join(CACHE_PATH, "overture_skyscrapers")
AIRCRAFT_CACHE_ROOT_DIR = os.path.join(CACHE_PATH, "aircraft", "opensky")
SATELLITE_CACHE_ROOT_DIR = os.path.join(CACHE_PATH, "satellites", "celestrak")

# Window UI
GUI_MENU_TEXT_COLOR = (128, 128, 128)
GUI_BUTTON_SIZE = 30
WINDOW_WIDTH = 600
WINDOW_HEIGHT = 600

# Minimum observer altitude (degrees)
OBSERVER_MIN_ALT_DEG = -5.0

# UI constants
TEXT_FONT_SIZE = 11
STATUS_LINE_FONT_SIZE = 8

# Palette-driven overlay accents.
PALETTE_EARTH_GUIDE_RGB = (112, 99, 89)
PALETTE_AIRCRAFT_AND_SATELLITE_RGB = (206, 122, 240)
PALETTE_HORIZON_AND_LABEL_RGB = (206, 240, 122)
PALETTE_ASTERISM_RGB = (122, 226, 240)


@dataclass(frozen=True, slots=True)
class TextStyle:
    text: tuple[int, ...]
    outline: tuple[int, ...]
    outline_width: float = 3.0


@dataclass(frozen=True, slots=True)
class WindowBackgroundStyle:
    base_rgb: tuple[int, int, int]
    delta_rgb: tuple[int, int, int]
    outer_alpha: int
    edge_alpha: int
    inner_rgba: tuple[int, int, int, int]
    border_rgba: tuple[int, int, int, int]

    def average_alpha(self) -> int:
        boundary_alpha = int(round(self.outer_alpha * 0.7 + self.edge_alpha * 0.3))
        return int(round((self.inner_rgba[3] + self.inner_rgba[3] + boundary_alpha + self.edge_alpha) / 4.0))


@dataclass(frozen=True, slots=True)
class SplashStyle:
    gradient_rgb: tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
    frame_rgb: tuple[int, int, int]
    info_text_rgb: tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class ThemeStyle:
    text: TextStyle
    status_text: TextStyle
    window_background: WindowBackgroundStyle
    splash: SplashStyle


THEME_STYLES_BY_PRESET = {
    "night": ThemeStyle(
        text=TextStyle(
            text=(180, 180, 180),
            outline=(0, 0, 0, 76),
        ),
        status_text=TextStyle(
            text=(190, 190, 160),
            outline=(0, 0, 0, 76),
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(10, 12, 16),
            delta_rgb=(7, 9, 11),
            outer_alpha=200,
            edge_alpha=80,
            inner_rgba=(4, 4, 4, 255),
            border_rgba=(30, 34, 40, 45),
        ),
        splash=SplashStyle(
            gradient_rgb=((12, 14, 20), (8, 10, 14), (4, 6, 9)),
            frame_rgb=(70, 76, 92),
            info_text_rgb=(228, 236, 250),
        ),
    ),
    "white": ThemeStyle(
        text=TextStyle(
            text=(228, 158, 92),
            outline=(130, 74, 30, 120),
            outline_width=3.6,
        ),
        status_text=TextStyle(
            text=(218, 150, 86),
            outline=(126, 72, 30, 124),
            outline_width=3.6,
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(242, 245, 250),
            delta_rgb=(46, 48, 50),
            outer_alpha=255,
            edge_alpha=200,
            inner_rgba=(4, 4, 4, 255),
            border_rgba=(254, 254, 255, 112),
        ),
        splash=SplashStyle(
            gradient_rgb=((252, 252, 252), (234, 234, 234), (206, 206, 206)),
            frame_rgb=(158, 178, 206),
            info_text_rgb=(228, 158, 92),
        ),
    ),
    "day": ThemeStyle(
        text=TextStyle(
            text=(232, 142, 104),
            outline=(128, 72, 40, 114),
            outline_width=3.6,
        ),
        status_text=TextStyle(
            text=(222, 136, 98),
            outline=(120, 68, 38, 118),
            outline_width=3.6,
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(230, 242, 255),
            delta_rgb=(28, 34, 34),
            outer_alpha=200,
            edge_alpha=80,
            inner_rgba=(4, 4, 4, 255),
            border_rgba=(250, 252, 255, 35),
        ),
        splash=SplashStyle(
            gradient_rgb=((240, 248, 255), (226, 240, 252), (206, 228, 246)),
            frame_rgb=(158, 182, 206),
            info_text_rgb=(232, 142, 104),
        ),
    ),
    "black": ThemeStyle(
        text=TextStyle(
            text=(246, 249, 255),
            outline=(2, 2, 3, 236),
        ),
        status_text=TextStyle(
            text=(255, 220, 220),
            outline=(2, 2, 3, 236),
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(12, 12, 12),
            delta_rgb=(9, 9, 9),
            outer_alpha=255,
            edge_alpha=200,
            inner_rgba=(4, 4, 4, 255),
            border_rgba=(34, 34, 36, 128),
        ),
        splash=SplashStyle(
            gradient_rgb=((6, 6, 6), (3, 3, 3), (0, 0, 0)),
            frame_rgb=(56, 56, 64),
            info_text_rgb=(244, 248, 255),
        ),
    ),
}

TEXT_STYLES_BY_PRESET = {
    preset: theme.text for preset, theme in THEME_STYLES_BY_PRESET.items()
}


HORIZON_LINE_COLOR = PALETTE_HORIZON_AND_LABEL_RGB
TERRAIN_HORIZON_LINE_COLOR = PALETTE_EARTH_GUIDE_RGB
EARTH_GUIDE_LINE_COLOR = TERRAIN_HORIZON_LINE_COLOR
URBAN_OUTLINE_LAYER_LINE_COLOR = (214, 232, 255)
CELESTIAL_EQUATOR_COLOR = (139, 139, 136)
ECLIPTIC_COLOR = (236, 173, 2)

CLOUD_UPDATE_INTERVAL = 10 * 60  # seconds

CLOUD_SHELLS_KM = (6371.0 + 3.0, 6371.0 + 5.0, 6371.0 + 7.0)  # representative cloud shells above Earth's surface

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
