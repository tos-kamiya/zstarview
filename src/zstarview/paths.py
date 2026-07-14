import colorsys
import os
import os.path
from dataclasses import dataclass, replace

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
GEOSATELLITE_DATA_DIR = os.path.join(_dir, "data", "geosatellite")
GEOSATELLITE_EQDC_LONLAT_FILE = os.path.join(GEOSATELLITE_DATA_DIR, "eqdc_lonlat.npz")
GEOSATELLITE_GRAY_COMMON_MASK_FILE = os.path.join(GEOSATELLITE_DATA_DIR, "Europe-IR-gray-common-mask.png")
CACHE_PATH = user_cache_dir(appname=APP_ID, appauthor=APP_AUTHOR)
LOG_PATH = user_log_dir(appname=APP_ID, appauthor=APP_AUTHOR)
COPERNICUS_DEM_CACHE_DIR = os.path.join(CACHE_PATH, "copernicus-dem")
OVERTURE_DERIVED_ROOT_DIR = os.path.join(CACHE_PATH, "overture_buildings")
OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR = os.path.join(CACHE_PATH, "overture_skyscrapers")
AIRCRAFT_CACHE_ROOT_DIR = os.path.join(CACHE_PATH, "aircraft", "opensky")
SATELLITE_CACHE_ROOT_DIR = os.path.join(CACHE_PATH, "satellites", "celestrak")
NIGHT_LIGHTS_CACHE_DIR = os.path.join(CACHE_PATH, "night_lights")
GEOSATELLITE_CACHE_ROOT_DIR = os.path.join(CACHE_PATH, "geosatellite")
TROPICAL_CYCLONE_CACHE_DIR = os.path.join(CACHE_PATH, "tropical_cyclones")

# Window UI
GUI_BUTTON_SIZE = 30
WINDOW_WIDTH = 600
WINDOW_HEIGHT = 600

# Observer altitude clamp bounds (degrees).
OBSERVER_MIN_ALT_DEG = -45.0
OBSERVER_MAX_ALT_DEG = 90.0

# UI constants
OVERLAY_FONT_SIZE_DEFAULT = 11
OVERLAY_FONT_SIZE_MIN = 8
OVERLAY_FONT_SIZE_MAX = 18
STATUS_LINE_FONT_SIZE = 8
TROPICAL_CYCLONE_DEFAULT_OPACITY = 0.7
CLOUD_DEFAULT_OPACITY = 0.05
NIGHT_LIGHT_DEFAULT_OPACITY = 0.05
RIDGE_GLOW_DEFAULT_OPACITY = 0.04

# Palette-driven overlay accents.
PALETTE_NEVER_RISES_RGB = (192, 192, 192)
PALETTE_NEVER_RISES_GUIDE_RGB = PALETTE_NEVER_RISES_RGB
PALETTE_EARTH_GUIDE_RGB = (112, 99, 89)
PALETTE_TERRAIN_HORIZON_RGB = (216, 206, 192)
PALETTE_AIRCRAFT_AND_SATELLITE_RGB = (227, 108, 240)
PALETTE_HORIZON_AND_LABEL_RGB = (206, 240, 122)
PALETTE_ATLAS_DIRECTION_GUIDE_RGB = (142, 166, 84)
PALETTE_ASTERISM_RGB = (122, 226, 240)

TRANSPARENT_THEME_OPACITY_VALUES = tuple(range(10, 100, 10))
TRANSPARENT_THEME_ALIAS = "transparent"
ATLAS_THEME_PRESET = "atlas-white"
TRANSPARENT_THEME_DEFAULT_PRESET = "transparent-40"
TRANSPARENT_THEME_PRESETS = tuple(
    f"transparent-{opacity}" for opacity in TRANSPARENT_THEME_OPACITY_VALUES
)
THEME_PRESET_NAMES = (
    "night",
    "day",
    "white",
    "black",
    TRANSPARENT_THEME_ALIAS,
    *TRANSPARENT_THEME_PRESETS,
)


@dataclass(frozen=True, slots=True)
class TextStyle:
    foreground_rgb: tuple[int, ...]
    outline_rgba: tuple[int, ...]
    outline_width: float = 3.0

    @property
    def text(self) -> tuple[int, ...]:
        return self.foreground_rgb

    @property
    def outline(self) -> tuple[int, ...]:
        return self.outline_rgba


@dataclass(frozen=True, slots=True)
class WindowBackgroundStyle:
    base_rgb: tuple[int, int, int]
    delta_rgb: tuple[int, int, int]
    outer_alpha: int
    edge_alpha: int
    inner_rgba: tuple[int, int, int, int]
    border_rgba: tuple[int, int, int, int]
    flat_background: bool = False

    def average_alpha(self) -> int:
        boundary_alpha = int(round(self.outer_alpha * 0.7 + self.edge_alpha * 0.3))
        return int(round((self.inner_rgba[3] + self.inner_rgba[3] + boundary_alpha + self.edge_alpha) / 4.0))


@dataclass(frozen=True, slots=True)
class SplashStyle:
    gradient_rgb: tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
    frame_rgb: tuple[int, int, int]
    info_text_rgb: tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class WindowChromeStyle:
    menu_fill_rgba: tuple[int, int, int, int]
    menu_icon_rgba: tuple[int, int, int, int]
    menu_button_text_rgb: tuple[int, int, int]
    menu_hover_text_rgb: tuple[int, int, int] = (255, 255, 255)
    menu_pressed_text_rgb: tuple[int, int, int] = (255, 255, 255)
    menu_hover_bg_rgba: tuple[int, int, int, int] = (255, 255, 255, 26)
    menu_pressed_bg_rgba: tuple[int, int, int, int] = (255, 255, 255, 40)


@dataclass(frozen=True, slots=True)
class SkyDiscStyle:
    opacity: float


@dataclass(frozen=True, slots=True)
class GuideStyle:
    reference_rgb: tuple[int, int, int]
    horizon_rgb: tuple[int, int, int]
    equator_rgb: tuple[int, int, int]
    ecliptic_rgb: tuple[int, int, int]
    simple_reference_lines: bool = False
    label_rgb: tuple[int, int, int] | None = None
    grid_rgb: tuple[int, int, int] | None = None
    marker_width: float = 1.6
    ecliptic_dash_pattern: tuple[float, ...] | None = None
    reference_width: float = 0.7
    grid_width: float = 0.51
    grid_minor_width: float = 0.51
    reference_alpha: int = 230
    grid_alpha: int = 190


DEFAULT_GUIDE_STYLE = GuideStyle(
    reference_rgb=(206, 240, 122),
    horizon_rgb=(206, 240, 122),
    equator_rgb=(192, 192, 192),
    ecliptic_rgb=(236, 173, 2),
)
ATLAS_GUIDE_STYLE = GuideStyle(
    reference_rgb=(24, 24, 24),
    horizon_rgb=(24, 24, 24),
    equator_rgb=(24, 24, 24),
    ecliptic_rgb=(24, 24, 24),
    simple_reference_lines=True,
    label_rgb=PALETTE_ATLAS_DIRECTION_GUIDE_RGB,
    grid_rgb=PALETTE_ATLAS_DIRECTION_GUIDE_RGB,
    marker_width=0.72,
    ecliptic_dash_pattern=(1.8, 3.2),
    reference_width=0.62,
    grid_width=0.46,
    grid_minor_width=0.42,
    reference_alpha=255,
    grid_alpha=225,
)


@dataclass(frozen=True, slots=True)
class OverlayLayerStyle:
    rgb: tuple[int, int, int]
    alpha_scale: float = 1.0
    width_scale: float = 1.0
    marker_width: float = 1.0
    outline_rgba: tuple[int, int, int, int] | None = None
    label_rgb: tuple[int, int, int] | None = None
    label_outline_rgba: tuple[int, int, int, int] | None = None
    line_alpha: int | None = None
    fill: bool = True


@dataclass(frozen=True, slots=True)
class CloudLayerStyle:
    rgb: tuple[int, int, int]
    alpha_scale: float = 1.0
    width_scale: float = 1.0
    missing_rgba: tuple[int, int, int, int] = (255, 220, 80, 45)


@dataclass(frozen=True, slots=True)
class OverlayStyles:
    aircraft: OverlayLayerStyle
    satellite: OverlayLayerStyle
    terrain_horizon: OverlayLayerStyle
    urban_outline: OverlayLayerStyle
    water: OverlayLayerStyle
    earth_guide: OverlayLayerStyle
    asterism: OverlayLayerStyle
    dso: OverlayLayerStyle
    cloud: CloudLayerStyle


DEFAULT_OVERLAY_STYLES = OverlayStyles(
    aircraft=OverlayLayerStyle(
        rgb=PALETTE_AIRCRAFT_AND_SATELLITE_RGB,
        label_rgb=PALETTE_AIRCRAFT_AND_SATELLITE_RGB,
    ),
    satellite=OverlayLayerStyle(
        rgb=PALETTE_AIRCRAFT_AND_SATELLITE_RGB,
        label_rgb=PALETTE_AIRCRAFT_AND_SATELLITE_RGB,
    ),
    terrain_horizon=OverlayLayerStyle(rgb=PALETTE_TERRAIN_HORIZON_RGB),
    urban_outline=OverlayLayerStyle(rgb=(214, 232, 255)),
    water=OverlayLayerStyle(rgb=(122, 218, 240)),
    earth_guide=OverlayLayerStyle(rgb=PALETTE_EARTH_GUIDE_RGB),
    asterism=OverlayLayerStyle(
        rgb=PALETTE_ASTERISM_RGB,
        label_rgb=(111, 207, 219),
    ),
    dso=OverlayLayerStyle(rgb=(122, 173, 240)),
    cloud=CloudLayerStyle(rgb=(255, 255, 255)),
)

ATLAS_OVERLAY_STYLES = OverlayStyles(
    aircraft=OverlayLayerStyle(
        rgb=(139, 24, 86),
        outline_rgba=(42, 24, 34, 80),
        label_rgb=(108, 20, 67),
        width_scale=0.85,
    ),
    satellite=OverlayLayerStyle(
        rgb=(139, 24, 86),
        outline_rgba=(42, 24, 34, 80),
        label_rgb=(108, 20, 67),
        width_scale=0.85,
        marker_width=2.0,
    ),
    terrain_horizon=OverlayLayerStyle(
        rgb=(82, 69, 56),
        outline_rgba=(42, 42, 42, 70),
        width_scale=0.82,
    ),
    urban_outline=OverlayLayerStyle(
        rgb=(60, 60, 60),
        outline_rgba=(32, 32, 32, 70),
        width_scale=0.78,
    ),
    water=OverlayLayerStyle(
        rgb=(28, 101, 115),
        outline_rgba=(14, 57, 68, 70),
        width_scale=0.82,
    ),
    earth_guide=OverlayLayerStyle(
        rgb=(76, 65, 55),
        outline_rgba=(38, 38, 38, 65),
        width_scale=0.82,
    ),
    asterism=OverlayLayerStyle(
        rgb=(35, 105, 120),
        label_rgb=(30, 85, 100),
        width_scale=0.66,
        line_alpha=105,
    ),
    dso=OverlayLayerStyle(
        rgb=(48, 97, 145),
        label_rgb=(45, 80, 120),
        width_scale=0.78,
        alpha_scale=1.2,
    ),
    cloud=CloudLayerStyle(
        rgb=(55, 78, 96),
        alpha_scale=0.88,
        width_scale=0.95,
        missing_rgba=(255, 220, 80, 45),
    ),
)


@dataclass(frozen=True, slots=True)
class ThemeStyle:
    text: TextStyle
    status_text: TextStyle
    window_background: WindowBackgroundStyle
    window_chrome: WindowChromeStyle
    sky_disc: SkyDiscStyle
    splash: SplashStyle
    star_visibility_boost: float = 1.0
    label_outline_suppressed: bool = False
    overlays: OverlayStyles = DEFAULT_OVERLAY_STYLES
    guide_style: GuideStyle = DEFAULT_GUIDE_STYLE


def _theme_background_luminance(base_rgb: tuple[int, int, int]) -> float:
    return (
        0.2126 * float(base_rgb[0])
        + 0.7152 * float(base_rgb[1])
        + 0.0722 * float(base_rgb[2])
    )


def _theme_chrome_icon_rgba(base_rgb: tuple[int, int, int]) -> tuple[int, int, int, int]:
    if _theme_background_luminance(base_rgb) >= 128.0:
        return (70, 70, 70, 220)
    return (210, 210, 210, 220)


def _theme_window_chrome(
    fill_rgba: tuple[int, int, int, int],
    base_rgb: tuple[int, int, int],
) -> WindowChromeStyle:
    icon_rgba = _theme_chrome_icon_rgba(base_rgb)
    return WindowChromeStyle(
        menu_fill_rgba=fill_rgba,
        menu_icon_rgba=icon_rgba,
        menu_button_text_rgb=icon_rgba[:3],
    )


def _theme_sky_disc(opacity: float) -> SkyDiscStyle:
    return SkyDiscStyle(opacity=float(opacity))


def _theme_with_sky_disc_opacity(theme: ThemeStyle, opacity: float) -> ThemeStyle:
    return replace(theme, sky_disc=_theme_sky_disc(opacity))


def _rgb_from_hsv(h_deg: float, s_pct: float, v_pct: float) -> tuple[int, int, int]:
    """Convert HSV values, given in degrees and percent, to an RGB tuple."""
    red, green, blue = colorsys.hsv_to_rgb(
        float(h_deg) / 360.0,
        float(s_pct) / 100.0,
        float(v_pct) / 100.0,
    )
    return (
        int(round(red * 255.0)),
        int(round(green * 255.0)),
        int(round(blue * 255.0)),
    )


THEME_STYLES_BY_PRESET = {
    "night": ThemeStyle(
        text=TextStyle(
            foreground_rgb=(180, 180, 180),
            outline_rgba=(0, 0, 0, 76),
        ),
        status_text=TextStyle(
            foreground_rgb=(190, 190, 160),
            outline_rgba=(0, 0, 0, 76),
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(10, 12, 16),
            delta_rgb=(7, 9, 11),
            outer_alpha=200,
            edge_alpha=120,
            inner_rgba=(4, 4, 4, 255),
            border_rgba=(30, 34, 40, 45),
        ),
        window_chrome=_theme_window_chrome((12, 14, 18, 112), (10, 12, 16)),
        sky_disc=_theme_sky_disc(1.0),
        splash=SplashStyle(
            gradient_rgb=((12, 14, 20), (8, 10, 14), (4, 6, 9)),
            frame_rgb=(70, 76, 92),
            info_text_rgb=(228, 236, 250),
        ),
        star_visibility_boost=1.0,
    ),
    "white": ThemeStyle(
        text=TextStyle(
            foreground_rgb=(214, 136, 103, 255),
            outline_rgba=(0, 0, 0, 76),
            outline_width=3.0,
        ),
        status_text=TextStyle(
            foreground_rgb=(190, 190, 160, 255),
            outline_rgba=(0, 0, 0, 76),
            outline_width=3.0,
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(240, 238, 237),
            delta_rgb=(38, 40, 42),
            outer_alpha=255,
            edge_alpha=214,
            inner_rgba=(4, 4, 4, 255),
            border_rgba=(254, 254, 255, 112),
        ),
        window_chrome=_theme_window_chrome((222, 220, 219, 100), (240, 238, 237)),
        sky_disc=_theme_sky_disc(1.0),
        splash=SplashStyle(
            gradient_rgb=((252, 252, 252), (234, 234, 234), (206, 206, 206)),
            frame_rgb=(158, 178, 206),
            info_text_rgb=(32, 32, 32),
        ),
        star_visibility_boost=1.12,
        label_outline_suppressed=True,
    ),
    ATLAS_THEME_PRESET: ThemeStyle(
        text=TextStyle(
            foreground_rgb=(24, 24, 24, 255),
            outline_rgba=(255, 255, 255, 220),
        ),
        status_text=TextStyle(
            foreground_rgb=(32, 32, 32, 255),
            outline_rgba=(255, 255, 255, 220),
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(255, 255, 255),
            delta_rgb=(0, 0, 0),
            outer_alpha=255,
            edge_alpha=255,
            inner_rgba=(255, 255, 255, 255),
            border_rgba=(160, 160, 160, 96),
            flat_background=True,
        ),
        window_chrome=_theme_window_chrome((244, 244, 244, 190), (255, 255, 255)),
        sky_disc=_theme_sky_disc(0.0),
        splash=SplashStyle(
            gradient_rgb=((255, 255, 255), (248, 248, 248), (232, 232, 232)),
            frame_rgb=(100, 100, 100),
            info_text_rgb=(24, 24, 24),
        ),
        star_visibility_boost=1.0,
        label_outline_suppressed=False,
        overlays=ATLAS_OVERLAY_STYLES,
        guide_style=ATLAS_GUIDE_STYLE,
    ),
    "day": ThemeStyle(
        text=TextStyle(
            foreground_rgb=_rgb_from_hsv(17.9, 51.9, 84.0),
            outline_rgba=(0, 0, 0, 76),
            outline_width=3.0,
        ),
        status_text=TextStyle(
            foreground_rgb=(190, 190, 160),
            outline_rgba=(0, 0, 0, 76),
            outline_width=3.0,
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(226, 223, 222),
            delta_rgb=(28, 34, 34),
            outer_alpha=200,
            edge_alpha=94,
            inner_rgba=(4, 4, 4, 255),
            border_rgba=(250, 252, 255, 35),
        ),
        window_chrome=_theme_window_chrome((228, 225, 224, 112), (226, 223, 222)),
        sky_disc=_theme_sky_disc(1.0),
        splash=SplashStyle(
            gradient_rgb=((240, 248, 255), (226, 240, 252), (206, 228, 246)),
            frame_rgb=(158, 182, 206),
            info_text_rgb=(32, 32, 32),
        ),
        star_visibility_boost=1.05,
        label_outline_suppressed=True,
    ),
    "black": ThemeStyle(
        text=TextStyle(
            foreground_rgb=(246, 249, 255),
            outline_rgba=(0, 0, 0, 76),
        ),
        status_text=TextStyle(
            foreground_rgb=(255, 220, 220),
            outline_rgba=(0, 0, 0, 76),
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(12, 12, 12),
            delta_rgb=(9, 9, 9),
            outer_alpha=255,
            edge_alpha=214,
            inner_rgba=(4, 4, 4, 255),
            border_rgba=(34, 34, 36, 128),
        ),
        window_chrome=_theme_window_chrome((34, 34, 34, 100), (12, 12, 12)),
        sky_disc=_theme_sky_disc(1.0),
        splash=SplashStyle(
            gradient_rgb=((6, 6, 6), (3, 3, 3), (0, 0, 0)),
            frame_rgb=(56, 56, 64),
            info_text_rgb=(244, 248, 255),
        ),
        star_visibility_boost=1.0,
    ),
    TRANSPARENT_THEME_ALIAS: ThemeStyle(
        text=TextStyle(
            foreground_rgb=(242, 245, 250),
            outline_rgba=(4, 4, 6, 214),
        ),
        status_text=TextStyle(
            foreground_rgb=(255, 224, 224),
            outline_rgba=(4, 4, 6, 214),
        ),
        window_background=WindowBackgroundStyle(
            base_rgb=(8, 8, 9),
            delta_rgb=(0, 0, 0),
            outer_alpha=80,
            edge_alpha=80,
            inner_rgba=(4, 4, 4, 80),
            border_rgba=(24, 24, 28, 60),
            flat_background=True,
        ),
        window_chrome=_theme_window_chrome((14, 14, 15, 96), (8, 8, 9)),
        sky_disc=_theme_sky_disc(0.4),
        splash=SplashStyle(
            gradient_rgb=((8, 8, 8), (4, 4, 4), (0, 0, 0)),
            frame_rgb=(56, 56, 64),
            info_text_rgb=(244, 248, 255),
        ),
        star_visibility_boost=1.0,
    ),
}

_transparent_base_theme = THEME_STYLES_BY_PRESET[TRANSPARENT_THEME_ALIAS]
for opacity_value in TRANSPARENT_THEME_OPACITY_VALUES:
    preset_name = f"transparent-{opacity_value}"
    THEME_STYLES_BY_PRESET[preset_name] = _theme_with_sky_disc_opacity(
        _transparent_base_theme,
        float(opacity_value) / 100.0,
    )

HORIZON_LINE_COLOR = PALETTE_HORIZON_AND_LABEL_RGB
TERRAIN_HORIZON_LINE_COLOR = PALETTE_TERRAIN_HORIZON_RGB
EARTH_GUIDE_LINE_COLOR = PALETTE_EARTH_GUIDE_RGB
URBAN_OUTLINE_LAYER_LINE_COLOR = (214, 232, 255)
CELESTIAL_EQUATOR_COLOR = PALETTE_NEVER_RISES_RGB
ECLIPTIC_COLOR = (236, 173, 2)
PALETTE_ASTERISM_LABEL_RGB = _rgb_from_hsv(187.1, 49.2, 86.0)

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
CLOUD_MISSING_TINT_RGBA = (255, 220, 80, 45)

# Skyfield ephemeris kernel filename
EPHEMERIS_FILENAME = "de442s.bsp"
EPHEMERIS_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de442s.bsp"
