import argparse
import sys
from typing import Sequence, Tuple, Union

from ..__about__ import __version__
from ..data.skyscraper_tiles import SKYSCRAPER_OUTER_RADIUS_KM
from ..paths import (
    CLOUD_MISSING_TINT_RGBA,
    DIRECTIONS,
    OVERLAY_FONT_SIZE_DEFAULT,
    OVERLAY_FONT_SIZE_MAX,
    OVERLAY_FONT_SIZE_MIN,
    THEME_PRESET_NAMES,
    TRANSPARENT_THEME_ALIAS,
    TRANSPARENT_THEME_DEFAULT_PRESET,
    TRANSPARENT_THEME_OPACITY_VALUES,
    TROPICAL_CYCLONE_DEFAULT_OPACITY,
    WINDOW_HEIGHT,
    WINDOW_WIDTH,
)

WindowGeometryArg = Union[str, Tuple[int, int, int, int]]
ImageSizeArg = Tuple[int, int]

_VMAG_MULTIPLIER_MIN = 10.0**0.2
_VMAG_MULTIPLIER_MAX = 10.0**0.4
SKY_OPACITY_DEFAULT = 0.1
_EDGE_FOV_MIN = 0.0
_EDGE_FOV_MAX = 135.0
_CONTENT_FOV_MIN = 90.0
_CONTENT_FOV_MAX = 135.0
_COMMITTED_VMAG_LIMIT_MAX = 10.5
_URBAN_OUTLINE_MAX_CANDIDATES_DEFAULT = 5000


def _parse_azimuth(value: str) -> float:
    """Parse azimuth given as degrees or compass points."""
    try:
        return float(value)
    except (TypeError, ValueError):
        pass

    compass = value.strip().upper()
    if compass in DIRECTIONS:
        return float(DIRECTIONS[compass])
    raise argparse.ArgumentTypeError(
        f"Invalid azimuth: {value!r}. Use degrees (e.g., 180) or compass (e.g., N, NE, E)."
    )


def _parse_theme(value: str) -> str:
    """Parse theme preset."""
    key = (value or "").strip().lower()
    allowed = {
        "night": "night",
        "day": "day",
        "white": "white",
        "black": "black",
        TRANSPARENT_THEME_ALIAS: TRANSPARENT_THEME_DEFAULT_PRESET,
        "translucent": TRANSPARENT_THEME_DEFAULT_PRESET,
    }
    if key in allowed:
        return allowed[key]
    if key.startswith(f"{TRANSPARENT_THEME_ALIAS}-"):
        suffix = key.removeprefix(f"{TRANSPARENT_THEME_ALIAS}-")
        if suffix.isdigit() and len(suffix) == 2:
            opacity = int(suffix)
            if opacity in TRANSPARENT_THEME_OPACITY_VALUES:
                return key
    raise argparse.ArgumentTypeError(
        f"Invalid theme: {value!r}. Use one of: {', '.join(THEME_PRESET_NAMES)}."
    )


def _parse_cloud_stripe(value: str) -> tuple[str, int, float]:
    """Parse cloud stripe style as 'mode[,count[,width_factor]]'."""
    text = (value or "").strip()
    parts = [p.strip() for p in text.split(",")]
    if len(parts) < 1 or len(parts) > 3:
        raise argparse.ArgumentTypeError(
            f"Invalid cloud stripe style: {value!r}. Use 'width[,count[,width]]' or 'alpha[,count[,width]]'."
        )
    mode = parts[0].lower()
    if mode not in {"width", "alpha"}:
        raise argparse.ArgumentTypeError(
            f"Invalid cloud stripe mode: {value!r}. Use 'width' or 'alpha'."
        )
    default_count = 50
    default_width = 0.85 if mode == "width" else 0.25
    try:
        count = default_count if len(parts) < 2 or parts[1] == "" else int(parts[1])
        width = default_width if len(parts) < 3 or parts[2] == "" else float(parts[2])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid cloud stripe style: {value!r}. Use numeric count/width after the mode."
        ) from exc
    if count < 0 or width < 0:
        raise argparse.ArgumentTypeError(
            f"Invalid cloud stripe style: {value!r}. count and width must be >= 0."
        )
    return mode, count, width


def _parse_vmag_brightness_multiplier(value: str) -> float:
    """Parse the brightness multiplier used per magnitude difference."""
    try:
        multiplier = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid magnitude brightness multiplier: {value!r}"
        ) from exc
    if not (_VMAG_MULTIPLIER_MIN <= multiplier <= _VMAG_MULTIPLIER_MAX):
        raise argparse.ArgumentTypeError(
            f"Value must be between {_VMAG_MULTIPLIER_MIN:.2f} and {_VMAG_MULTIPLIER_MAX:.2f}."
        )
    return multiplier


def _parse_positive_int(value: str) -> int:
    """Parse a strictly positive integer."""
    try:
        out = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid positive integer: {value!r}"
        ) from exc
    if out <= 0:
        raise argparse.ArgumentTypeError("Value must be > 0.")
    return out


def _parse_non_negative_int(value: str) -> int:
    """Parse an integer >= 0."""
    try:
        out = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid non-negative integer: {value!r}"
        ) from exc
    if out < 0:
        raise argparse.ArgumentTypeError("Value must be >= 0.")
    return out


def _parse_overlay_font_size(value: str) -> float:
    """Parse the overlay font size used for canvas text."""
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid overlay font size: {value!r}"
        ) from exc
    if not (OVERLAY_FONT_SIZE_MIN <= out <= OVERLAY_FONT_SIZE_MAX):
        raise argparse.ArgumentTypeError(
            f"Value must be between {OVERLAY_FONT_SIZE_MIN} and {OVERLAY_FONT_SIZE_MAX}."
        )
    return out


def _parse_non_negative_float(value: str) -> float:
    """Parse a float >= 0."""
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid non-negative float: {value!r}"
        ) from exc
    if out < 0.0:
        raise argparse.ArgumentTypeError("Value must be >= 0.")
    return out


def _parse_content_fov_deg(value: str) -> float:
    """Parse the shared overscan content FOV angle."""
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(f"Invalid content FOV: {value!r}") from exc
    if not (_CONTENT_FOV_MIN <= out <= _CONTENT_FOV_MAX):
        raise argparse.ArgumentTypeError(
            f"Value must be between {_CONTENT_FOV_MIN:.0f} and {_CONTENT_FOV_MAX:.0f} degrees."
        )
    return out


def _parse_edge_fov_deg(value: str) -> float:
    """Parse the projected edge FOV angle."""
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(f"Invalid edge FOV: {value!r}") from exc
    if not (_EDGE_FOV_MIN < out <= _EDGE_FOV_MAX):
        raise argparse.ArgumentTypeError(
            f"Value must be greater than {_EDGE_FOV_MIN:.0f} and at most {_EDGE_FOV_MAX:.0f} degrees."
        )
    return out


def _parse_bool(value: str) -> bool:
    """Parse boolean values for CLI options."""
    text = (value or "").strip().lower()
    if text in {"true", "t", "1", "yes", "y", "on"}:
        return True
    if text in {"false", "f", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(
        f"Invalid boolean value: {value!r}. Use true/false."
    )


def _parse_sky_disc_style(value: str) -> str:
    """Parse the sky-disc fill style."""
    text = (value or "").strip().lower()
    if text == "smooth":
        return text
    raise argparse.ArgumentTypeError(
        f"Invalid sky disc style: {value!r}. Use 'smooth'."
    )


def _parse_sky_disc_altaz_rings_mode(value: str) -> str:
    """Parse the sky-disc alt/az ring visibility mode."""
    text = (value or "").strip().lower()
    if text in {"off", "dimalt", "altaz"}:
        return text
    raise argparse.ArgumentTypeError(
        f"Invalid sky disc alt/az rings mode: {value!r}. Use 'off', 'dimalt', or 'altaz'."
    )


def _parse_window_geometry(value: str) -> WindowGeometryArg:
    """Parse window geometry as 'restore' or 'x,y,width,height'."""
    text = (value or "").strip()
    if text.lower() == "restore":
        return "restore"
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            f"Invalid window geometry: {value!r}. Use 'restore' or 'x,y,width,height'."
        )
    try:
        x, y, width, height = (int(p) for p in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid window geometry: {value!r}. Use integers for x,y,width,height."
        ) from exc
    return (x, y, width, height)


def _parse_window_frame(value: str) -> str:
    """Parse the window chrome mode."""
    text = (value or "").strip().lower()
    if text in {"frameless", "window"}:
        return text
    raise argparse.ArgumentTypeError(
        f"Invalid window frame: {value!r}. Use 'frameless' or 'window'."
    )


def _parse_image_size(value: str) -> ImageSizeArg:
    """Parse image size as 'width,height'."""
    text = (value or "").strip()
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"Invalid image size: {value!r}. Use 'width,height'."
        )
    try:
        width, height = (int(p) for p in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid image size: {value!r}. Use integers for width,height."
        ) from exc
    if width <= 0 or height <= 0:
        raise argparse.ArgumentTypeError("Image width and height must be > 0.")
    return (width, height)


def add_location_arguments(parser: argparse.ArgumentParser) -> None:
    """Add observing-location arguments."""
    parser.add_argument(
        "city",
        type=str,
        nargs="?",
        default="",
        help=(
            "Location name: city, lat;lon, lat,lon, @lat,lon, supported Google Maps URL, "
            "tower name, mountain name, or explicit t/NAME and m/NAME prefixes "
            "(default: same as the last run)"
        ),
    )
    parser.add_argument(
        "-p",
        "--place",
        type=str,
        default=None,
        metavar="QUERY",
        help=(
            "Search a place, station, or facility name via Nominatim and use the "
            "top candidate as the observing location."
        ),
    )
    parser.add_argument(
        "--place-countrycode",
        type=str,
        default=None,
        metavar="CODE",
        help=(
            "Restrict --place search to an ISO 3166-1 alpha-2 country code "
            "(for example: jp)."
        ),
    )
    parser.add_argument(
        "--place-lang",
        type=str,
        default="en",
        metavar="LANG",
        help="Accept-Language for --place Nominatim search (default: en).",
    )


def add_dataset_query_arguments(parser: argparse.ArgumentParser) -> None:
    """Add bundled-viewpoint dataset query arguments."""
    dataset_query_group = parser.add_mutually_exclusive_group()
    dataset_query_group.add_argument(
        "--list-viewpoints",
        type=str,
        default=None,
        metavar="{t,m}",
        help="List bundled viewpoint primary names for towers (t) or mountains (m) and exit.",
    )
    dataset_query_group.add_argument(
        "--list-viewpoint-names",
        type=str,
        default=None,
        metavar="{t,m}",
        help="List bundled viewpoint names including localized names for towers (t) or mountains (m) and exit.",
    )
    dataset_query_group.add_argument(
        "--show-viewpoint-json",
        type=str,
        default=None,
        metavar="NAME",
        help="Resolve a bundled viewpoint name and print its JSON metadata, then exit.",
    )


def add_time_arguments(parser: argparse.ArgumentParser) -> None:
    """Add time-selection arguments."""
    time_group = parser.add_argument_group("Time settings")
    time_group.add_argument(
        "-H",
        "--hours",
        type=float,
        default=0,
        help="Number of hours to add to current time (default: 0)",
    )
    time_group.add_argument(
        "-D",
        "--days",
        type=float,
        default=0,
        help="Number of days to add to current time (default: 0)",
    )
    time_group.add_argument(
        "--datetime",
        type=str,
        default=None,
        help=(
            "Set an absolute date and time in 'YYYY-MM-DD HH[:MM[:SS]] [TZ]' format "
            "(e.g., '2025-09-12 9', '2025-09-12 09:00', '2025-09-12 9:0:0 JST'). "
            "If TZ is omitted, the resolved location timezone is used. "
            "Overrides --hours and --days."
        ),
    )
    time_group.add_argument(
        "--timezone",
        type=str,
        default=None,
        metavar="TZ",
        help=(
            "Override the resolved location timezone for --datetime and on-screen time "
            "(for example: Asia/Tokyo, JST, UTC+9)."
        ),
        )


def add_search_arguments(
    parser: argparse.ArgumentParser, *, include_list: bool = True
) -> None:
    """Add object-search arguments."""
    parser.add_argument(
        "--search",
        type=str,
        default=None,
        metavar="QUERY",
        help=(
            "Search stars, asterisms, places, and JPL bodies at startup. "
            "Bare queries search label and id; label= and id= limit the search."
        ),
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help=(
            "List search candidates one per line instead of resolving a unique target."
            if include_list
            else argparse.SUPPRESS
        ),
    )
    parser.add_argument(
        "--search-keep-marker",
        action="store_true",
        help=(
            "Keep the selected search target as a persistent marker and label."
        ),
    )


def add_geo_satellite_argument(parser: argparse.ArgumentParser) -> None:
    """Add the Geo-satellite cloud data toggle."""
    parser.add_argument(
        "--geo-satellite",
        type=_parse_bool,
        default=False,
        metavar="true|false",
        help=(
            "Use Geo-satellite infrared cloud data instead of GOES/Himawari "
            "when the observer is inside the Europe workflow band (true/false)."
        ),
    )


def add_observing_arguments(parser: argparse._ActionsContainer) -> None:
    """Add observing location, time, and view-center arguments."""
    add_location_arguments(parser)
    parser.add_argument(
        "-H",
        "--hours",
        type=float,
        default=0,
        help="Number of hours to add to current time (default: 0)",
    )
    parser.add_argument(
        "-D",
        "--days",
        type=float,
        default=0,
        help="Number of days to add to current time (default: 0)",
    )
    parser.add_argument(
        "--datetime",
        type=str,
        default=None,
        help=(
            "Set an absolute date and time in 'YYYY-MM-DD HH[:MM[:SS]] [TZ]' format "
            "(e.g., '2025-09-12 9', '2025-09-12 09:00', '2025-09-12 9:0:0 JST'). "
            "If TZ is omitted, the resolved location timezone is used. "
            "Overrides --hours and --days."
        ),
    )
    parser.add_argument(
        "--timezone",
        type=str,
        default=None,
        metavar="TZ",
        help=(
            "Override the resolved location timezone for --datetime and on-screen time "
            "(for example: Asia/Tokyo, JST, UTC+9)."
        ),
    )
    parser.add_argument(
        "-V",
        "--vmag-limit",
        type=float,
        default=7.0,
        help=(
            "Limit stars to Vmag <= this value (default: 7.0). "
            f"Bundled catalogs clamp values above {_COMMITTED_VMAG_LIMIT_MAX:.1f}."
        ),
    )
    parser.add_argument(
        "--vmag-brightness-multiplier",
        type=_parse_vmag_brightness_multiplier,
        default=2.5,
        help="Brightness multiplier per magnitude step (allowed range: 1.58-2.512, default: 2.5; 2.512 is the classical Pogson value).",
    )
    parser.add_argument(
        "-m",
        "--enlarge-moon",
        action="store_true",
        help="Show the moon in 5x size.",
    )
    parser.add_argument(
        "--bright-bodies",
        choices=("outline", "fill"),
        default="outline",
        help="Bright bodies rendering mode: outline or fill (default: outline).",
    )
    parser.add_argument(
        "-s",
        "--star-base-radius",
        type=float,
        default=4.0,
        help="Base size of 2nd-magnitude stars (default: 4.0)",
    )
    parser.add_argument(
        "-w",
        "--expected-render-width",
        type=_parse_positive_int,
        default=WINDOW_WIDTH,
        help=(
            "Expected window width for full-resolution star rendering (default: 600). "
            "When the celestial disc width exceeds this, the star layer uses square-root scaling."
        ),
    )
    parser.add_argument(
        "-Z",
        "--view-center-az",
        type=_parse_azimuth,
        default=180.0,
        help=(
            "Viewing azimuth angle [deg or compass] "
            "(0=N, 90=E, 180=S, 270=W; also accepts N, NE, E, SE, S, SW, W, NW; default=180)"
        ),
    )
    parser.add_argument(
        "-A",
        "--view-center-alt",
        type=float,
        default=90.0,
        help="Viewing altitude angle [deg] (90=zenith, 0=horizon; default=90)",
    )
    parser.add_argument(
        "--edge-fov-deg",
        type=_parse_edge_fov_deg,
        default=95.0,
        help=(
            "Projected edge field of view in degrees (default: 95). "
            "This controls how angular distance maps to the window edge."
        ),
    )
    parser.add_argument(
        "--content-fov-deg",
        type=_parse_content_fov_deg,
        default=110.0,
        help=(
            "Shared overscan content FOV in degrees (allowed: 90-135, default: 110). "
            "The window edge remains fixed at 90 degrees from the view center."
        ),
    )
    parser.add_argument(
        "--height-add-m",
        "--observer-height-m",
        dest="observer_height_m",
        type=_parse_non_negative_float,
        default=None,
        help=(
            "Additional height above the active observation base in meters. "
            "Default: 1.7 for city/latlon/mountain and tower height + 1.7 for tower-name input. "
            "--observer-height-m remains as a compatibility alias."
        ),
    )
    parser.add_argument(
        "--use-building-top",
        action="store_true",
        help=(
            "If the resolved location lies inside a building footprint, use that building's "
            "highest top height as the observation base before adding observer eye height."
        ),
    )
    return


def add_sky_and_star_arguments(
    parser: argparse._ActionsContainer, *, include_sky_update_interval: bool = True
) -> None:
    """Add the sky and star rendering arguments."""
    parser.add_argument(
        "--sky-opacity",
        type=float,
        default=SKY_OPACITY_DEFAULT,
        help=(
            f"Opacity of the simulated sky-color disc (0.0 - 1.0, default: {SKY_OPACITY_DEFAULT}). "
            "Set to 0.0 to disable sky-color rendering."
        ),
    )
    parser.add_argument(
        "--sky-disc-style",
        type=_parse_sky_disc_style,
        default="smooth",
        metavar="{smooth}",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--sky-disc-altaz-rings",
        type=_parse_sky_disc_altaz_rings_mode,
        default="dimalt",
        metavar="{off,dimalt,altaz}",
        help=(
            "Draw the always-visible sky-disc alt/az grid overlay. "
            "Use 'dimalt' for subtle altitude rings or 'altaz' for the full grid."
        ),
    )
    parser.add_argument(
        "--sky-disc-altaz-rings-hover",
        type=_parse_sky_disc_altaz_rings_mode,
        default="altaz",
        metavar="{off,dimalt,altaz}",
        help=(
            "Draw the hover-time sky-disc alt/az grid overlay. "
            "Use 'dimalt' for subtle altitude rings or 'altaz' for the full grid."
        ),
    )
    parser.add_argument(
        "--show-dso-initial",
        type=_parse_bool,
        default=None,
        metavar="true|false",
        help="Whether to show DSO overlays at startup (true/false).",
    )
    parser.add_argument(
        "--show-asterisms-initial",
        type=_parse_bool,
        default=None,
        metavar="true|false",
        help="Whether to show asterism overlays at startup (true/false).",
    )
    if include_sky_update_interval:
        parser.add_argument(
            "-i",
            "--sky-update-interval",
            type=int,
            default=60,
            help="Interval for updating stars/sky-color disc in sec. (default: 60).",
        )


def add_overlay_arguments(parser: argparse._ActionsContainer) -> None:
    """Add overlay-related rendering arguments."""
    parser.add_argument(
        "-c",
        "--cloud-opacity",
        type=float,
        default=0.075,
        help=(
            "Opacity of the clouds (0.0 - 1.0, default: 0.075). "
            "Set to 0.0 to disable cloud rendering."
        ),
    )
    parser.add_argument(
        "--cloud-stripe",
        type=_parse_cloud_stripe,
        default=("width", 50, 0.85),
        metavar="MODE[,COUNT[,WIDTH]]",
        help=(
            "Cloud stripe style as 'mode[,count[,width]]' "
            "(count is scaled with the star render surface size; defaults: "
            "width -> width,50,0.85; alpha -> alpha,50,0.25). "
            "If either value is 0, cloud rendering is disabled."
        ),
    )
    parser.add_argument(
        "--cloud-missing-tint-opacity",
        type=float,
        default=float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0,
        help=(
            "Opacity for missing-cloud-data tint (0.0 - 1.0, default: "
            f"{float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0:.3f})."
        ),
    )
    parser.add_argument(
        "--tropical-cyclone-opacity",
        type=float,
        default=TROPICAL_CYCLONE_DEFAULT_OPACITY,
        help=(
            "Opacity of the tropical cyclone overlay (0.0 - 1.0, default: "
            f"{TROPICAL_CYCLONE_DEFAULT_OPACITY}). "
            "Set to 0.0 to disable cyclone API fetch and rendering for that run."
        ),
    )
    parser.add_argument(
        "-a",
        "--aircraft-opacity",
        type=float,
        default=0.4,
        help=(
            "Opacity of the aircraft overlay (0.0 - 1.0, default: 0.4). "
            "Set to 0.0 to disable aircraft queries and rendering."
        ),
    )
    parser.add_argument(
        "--satellite-opacity",
        type=float,
        default=0.7,
        help=(
            "Opacity of the artificial satellite overlay (0.0 - 1.0, default: 0.7). "
            "Set to 0.0 to disable satellite element fetch and rendering."
        ),
    )
    parser.add_argument(
        "--show-guidelines-initial",
        type=_parse_bool,
        default=None,
        metavar="true|false",
        help="Whether to show guideline overlays at startup (true/false).",
    )
    parser.add_argument(
        "-d",
        "--terrain-horizon-opacity",
        type=float,
        default=0.003,
        help=(
            "Opacity of the terrain horizon polyline (0.0 - 1.0, default: 0.003). "
            "Set to 0.0 to disable DEM download, terrain-horizon calculation, and drawing."
        ),
    )
    parser.add_argument(
        "-e",
        "--earth-guide-opacity",
        type=float,
        default=0.028,
        help=(
            "Opacity of the Earth guide line layer (0.0 - 1.0, default: 0.028). "
            "Set to 0.0 to disable Earth guide drawing and lock the GUI toggle off for that session."
        ),
    )
    parser.add_argument(
        "--ground-tint-opacity",
        type=float,
        default=0.04,
        help=(
            "Overlay opacity of the ground tint color below the geometric/terrain horizon "
            "(0.0 - 1.0, default: 0.04)."
        ),
    )
    parser.add_argument(
        "--water-surface-opacity",
        type=float,
        default=0.4,
        help=(
            "Opacity of the water surface layer (0.0 - 1.0, default: 0.4). "
            "Set to 0.0 to disable water surface rendering at startup."
        ),
    )
    parser.add_argument(
        "--night-light-opacity",
        type=float,
        default=0.04,
        help=(
            "Opacity of the street-light part of the night light overlay (0.0 - 1.0, default: 0.04). "
            "Set to 0.0 to disable street-light rendering and lock the GUI toggle off for that session."
        ),
    )
    parser.add_argument(
        "--sky-glow-opacity",
        type=float,
        default=0.2,
        help=(
            "Opacity of the sky glow drawn over the terrain horizon (0.0 - 1.0, default: 0.2). "
            "Set to 0.0 to disable the sky glow independently of the street-light layer."
        ),
    )
    parser.add_argument(
        "-u",
        "--urban-outline-opacity",
        type=float,
        default=0.2,
        help=(
            "Opacity of the urban outline overlay (0.0 - 1.0, default: 0.2). "
            "Set to 0.0 to disable urban outline rendering at startup."
        ),
    )
    parser.add_argument(
        "--urban-outline-feature-type",
        choices=("both", "building"),
        default="both",
        help=(
            "Urban outline Overture mode "
            "(default: both). Use building to skip building_part overlays."
        ),
    )
    parser.add_argument(
        "-r",
        "--urban-outline-radius-km",
        type=_parse_non_negative_float,
        default=2.5,
        help=(
            "Download and cache radius for the urban outline overlay in kilometers "
            "(default: 2.5). This value becomes part of the cache key."
        ),
    )
    parser.add_argument(
        "--urban-outline-skyscraper-radius-km",
        type=_parse_non_negative_float,
        default=SKYSCRAPER_OUTER_RADIUS_KM,
        help=(
            "Outer radius for the far-range skyscraper urban outline layer in kilometers "
            f"(default: {SKYSCRAPER_OUTER_RADIUS_KM:.1f}). "
            "Use 0 to disable skyscraper-tile lookup for that run. "
            "Otherwise it must be greater than or equal to --urban-outline-radius-km."
        ),
    )
    parser.add_argument(
        "-b",
        "--urban-outline-min-building-height-m",
        dest="urban_outline_min_height_m",
        type=_parse_non_negative_float,
        default=0.0,
        help=(
            "Deprecated. Minimum building height in meters for buildings included in the urban outline overlay "
            "(default: 0.0). Use --urban-outline-max-candidates for performance tuning."
        ),
    )
    parser.add_argument(
        "--urban-outline-max-candidates",
        type=_parse_non_negative_int,
        default=_URBAN_OUTLINE_MAX_CANDIDATES_DEFAULT,
        help=(
            "Maximum number of urban-outline ring candidates to keep before expensive sampling "
            f"(default: {_URBAN_OUTLINE_MAX_CANDIDATES_DEFAULT}). "
            "Set to 0 to disable the layer."
        ),
    )
    parser.add_argument(
        "--urban-outline-skyscraper-only",
        action="store_true",
        help=(
            "Validation mode: draw only the far-range skyscraper urban outline layer "
            "and skip the normal near-range outline layer."
        ),
    )


def add_general_arguments(
    parser: argparse._ActionsContainer,
    *,
    include_window_geometry: bool = True,
    include_window_frame: bool = True,
) -> None:
    """Add general-purpose CLI arguments."""
    if include_window_geometry:
        parser.add_argument(
            "--window-geometry",
            type=_parse_window_geometry,
            default=None,
            metavar="restore|X,Y,W,H",
            help=(
                "Window position and size. "
                "Use 'restore' to load the last saved geometry, "
                "or 'x,y,width,height' to set explicit values."
            ),
        )
    if include_window_frame:
        parser.add_argument(
            "--window-frame",
            type=_parse_window_frame,
            default="frameless",
            metavar="{frameless,window}",
            help=(
                "Window decoration mode. "
                "Use 'frameless' for the current borderless window or "
                "'window' for a standard titled OS window."
            ),
        )
    parser.add_argument(
        "-t",
        "--theme",
        type=_parse_theme,
        default="night",
        metavar="{night,day,white,black,transparent,transparent-10..90}",
        help="Theme preset for background and star contrast (default: night).",
    )
    parser.add_argument(
        "--observation-info",
        type=str,
        default="bottom",
        choices=("auto", "top", "bottom", "off"),
        metavar="auto|top|bottom|off",
        help=(
            "Observation info overlay mode at startup: auto (hover-avoid), "
            "top (fixed top), bottom (fixed bottom, default), or off (hidden)."
        ),
    )
    parser.add_argument(
        "--visibility-boost",
        type=float,
        default=1.0,
        help=(
            "Tiered opacity boost for hard-to-see layers (0.0 - 1.0 base, default: 1.0). "
            "Values above 1.0 increase supplemental layers more than small figure layers."
        ),
    )
    parser.add_argument(
        "--overlay-font-size",
        type=_parse_overlay_font_size,
        default=OVERLAY_FONT_SIZE_DEFAULT,
        metavar="POINTS",
        help=(
            "Base font size for window-drawn labels and HUD text only "
            f"({OVERLAY_FONT_SIZE_MIN}-{OVERLAY_FONT_SIZE_MAX}, decimals allowed, default: {OVERLAY_FONT_SIZE_DEFAULT}). "
            "This does not affect the status line, menus, dialogs, or standard Qt widgets."
        ),
    )
    parser.add_argument(
        "--clear-long-lived-cache",
        action="store_true",
        help=(
            "Delete cached DEM and urban-outline data before startup. "
            "This clears copernicus-dem, overture_buildings, and overture_skyscrapers."
        ),
    )
    parser.add_argument(
        "--close-on-startup-error",
        action="store_true",
        help=(
            "Close the app automatically with a non-zero exit code when startup "
            "fails instead of keeping the error log visible in the window."
        ),
    )


def add_export_image_arguments(parser: argparse.ArgumentParser) -> None:
    """Add grouped arguments for the headless export-image CLI."""
    observing_group = parser.add_argument_group("Observing Location and Time")
    search_group = parser.add_argument_group("Search Objects at startup")
    sky_group = parser.add_argument_group("Sky and Stars")
    overlay_group = parser.add_argument_group("Overlays")
    export_group = parser.add_argument_group("Export")
    general_group = parser.add_argument_group("General")

    add_observing_arguments(observing_group)
    add_search_arguments(search_group, include_list=True)
    add_sky_and_star_arguments(sky_group, include_sky_update_interval=False)
    add_overlay_arguments(overlay_group)
    add_geo_satellite_argument(export_group)
    add_general_arguments(
        general_group,
        include_window_geometry=False,
        include_window_frame=False,
    )

    export_group.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        metavar="PATH",
        help="Output image file path.",
    )
    export_group.add_argument(
        "--image-size",
        type=_parse_image_size,
        default=(WINDOW_WIDTH, WINDOW_HEIGHT),
        metavar="WIDTH,HEIGHT",
        help=(
            f"Output image size in pixels (default: {WINDOW_WIDTH},{WINDOW_HEIGHT})."
        ),
    )
    export_group.add_argument(
        "--layer-timeout-seconds",
        type=_parse_non_negative_float,
        default=90.0,
        metavar="SECONDS",
        help="Maximum time to wait for enabled external layers (default: 90).",
    )
    export_group.add_argument(
        "--allow-partial-data",
        action="store_true",
        help="Allow saving an image even when enabled external layers fail or time out.",
    )
    export_group.add_argument(
        "--include-direction-grid",
        action="store_true",
        help="Include the direction grid in exported images, with 30-degree major lines and 10-degree intersection crosses.",
    )
    export_group.add_argument(
        "--print-cache-dir",
        action="store_true",
        help="Print the cache root directory and exit.",
    )
    export_group.add_argument(
        "--sixel",
        action="store_true",
        help="Display the generated image in the terminal via img2sixel.",
    )


def add_render_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_window_geometry: bool = True,
    include_window_frame: bool = True,
    include_sky_update_interval: bool = True,
    include_startup_overlay_arguments: bool = True,
) -> None:
    """Add shared view and rendering arguments."""
    parser.add_argument(
        "-V",
        "--vmag-limit",
        type=float,
        default=7.0,
        help=(
            "Limit stars to Vmag <= this value (default: 7.0). "
            f"Bundled catalogs clamp values above {_COMMITTED_VMAG_LIMIT_MAX:.1f}."
        ),
    )
    parser.add_argument(
        "--vmag-brightness-multiplier",
        type=_parse_vmag_brightness_multiplier,
        default=2.5,
        help="Brightness multiplier per magnitude step (allowed range: 1.58-2.512, default: 2.5; 2.512 is the classical Pogson value).",
    )
    parser.add_argument(
        "-m",
        "--enlarge-moon",
        action="store_true",
        help="Show the moon in 5x size.",
    )
    parser.add_argument(
        "--bright-bodies",
        choices=("outline", "fill"),
        default="outline",
        help="Bright bodies rendering mode: outline or fill (default: outline).",
    )
    parser.add_argument(
        "-s",
        "--star-base-radius",
        type=float,
        default=4.0,
        help="Base size of 2nd-magnitude stars (default: 4.0)",
    )
    parser.add_argument(
        "-w",
        "--expected-render-width",
        type=_parse_positive_int,
        default=WINDOW_WIDTH,
        help=(
            "Expected window width for full-resolution star rendering (default: 600). "
            "When the celestial disc width exceeds this, the star layer uses square-root scaling."
        ),
    )
    parser.add_argument(
        "-Z",
        "--view-center-az",
        type=_parse_azimuth,
        default=180.0,
        help=(
            "Viewing azimuth angle [deg or compass] "
            "(0=N, 90=E, 180=S, 270=W; also accepts N, NE, E, SE, S, SW, W, NW; default=180)"
        ),
    )
    parser.add_argument(
        "-A",
        "--view-center-alt",
        type=float,
        default=90.0,
        help="Viewing altitude angle [deg] (90=zenith, 0=horizon; default=90)",
    )
    parser.add_argument(
        "--content-fov-deg",
        type=_parse_content_fov_deg,
        default=110.0,
        help=(
            "Shared overscan content FOV in degrees (allowed: 90-135, default: 110). "
            "The window edge remains fixed at 90 degrees from the view center."
        ),
    )
    parser.add_argument(
        "--height-add-m",
        "--observer-height-m",
        dest="observer_height_m",
        type=_parse_non_negative_float,
        default=None,
        help=(
            "Additional height above the active observation base in meters. "
            "Default: 1.7 for city/latlon/mountain and tower height + 1.7 for tower-name input. "
            "--observer-height-m remains as a compatibility alias."
        ),
    )
    parser.add_argument(
        "--use-building-top",
        action="store_true",
        help=(
            "If the resolved location lies inside a building footprint, use that building's "
            "highest top height as the observation base before adding observer eye height."
        ),
    )
    parser.add_argument(
        "--sky-opacity",
        type=float,
        default=SKY_OPACITY_DEFAULT,
        help=(
            f"Opacity of the simulated sky-color disc (0.0 - 1.0, default: {SKY_OPACITY_DEFAULT}). "
            "Set to 0.0 to disable sky-color rendering."
        ),
    )
    parser.add_argument(
        "-c",
        "--cloud-opacity",
        type=float,
        default=0.075,
        help=(
            "Opacity of the clouds (0.0 - 1.0, default: 0.075). "
            "Set to 0.0 to disable cloud rendering."
        ),
    )
    parser.add_argument(
        "--tropical-cyclone-opacity",
        type=float,
        default=TROPICAL_CYCLONE_DEFAULT_OPACITY,
        help=(
            "Opacity of the tropical cyclone overlay (0.0 - 1.0, default: "
            f"{TROPICAL_CYCLONE_DEFAULT_OPACITY}). "
            "Set to 0.0 to disable cyclone API fetch and rendering for that run."
        ),
    )
    parser.add_argument(
        "-a",
        "--aircraft-opacity",
        type=float,
        default=0.4,
        help=(
            "Opacity of the aircraft overlay (0.0 - 1.0, default: 0.4). "
            "Set to 0.0 to disable aircraft queries and rendering."
        ),
    )
    parser.add_argument(
        "--satellite-opacity",
        type=float,
        default=0.7,
        help=(
            "Opacity of the artificial satellite overlay (0.0 - 1.0, default: 0.7). "
            "Set to 0.0 to disable satellite element fetch and rendering."
        ),
    )
    parser.add_argument(
        "--terrain-horizon-opacity",
        type=float,
        default=0.003,
        help=(
            "Opacity of the terrain horizon polyline (0.0 - 1.0, default: 0.003). "
            "Set to 0.0 to disable DEM download, terrain-horizon calculation, and drawing."
        ),
    )
    parser.add_argument(
        "--earth-guide-opacity",
        type=float,
        default=0.028,
        help=(
            "Opacity of the Earth guide line layer (0.0 - 1.0, default: 0.028). "
            "Set to 0.0 to disable Earth guide drawing and lock the GUI toggle off for that session."
        ),
    )
    parser.add_argument(
        "--urban-outline-opacity",
        type=float,
        default=0.2,
        help=(
            "Opacity of the urban outline overlay (0.0 - 1.0, default: 0.2). "
            "Set to 0.0 to disable urban outline rendering at startup."
        ),
    )
    parser.add_argument(
        "--night-light-opacity",
        type=float,
        default=0.04,
        help=(
            "Opacity of the street-light part of the night light overlay (0.0 - 1.0, default: 0.04). "
            "Set to 0.0 to disable street-light rendering and lock the GUI toggle off for that session."
        ),
    )
    parser.add_argument(
        "--sky-glow-opacity",
        type=float,
        default=0.2,
        help=(
            "Opacity of the sky glow drawn over the terrain horizon (0.0 - 1.0, default: 0.2). "
            "Set to 0.0 to disable the sky glow independently of the street-light layer."
        ),
    )
    parser.add_argument(
        "-r",
        "--urban-outline-radius-km",
        type=_parse_non_negative_float,
        default=2.5,
        help=(
            "Download and cache radius for the urban outline overlay in kilometers "
            "(default: 2.5). This value becomes part of the cache key."
        ),
    )
    parser.add_argument(
        "--urban-outline-skyscraper-radius-km",
        type=_parse_non_negative_float,
        default=SKYSCRAPER_OUTER_RADIUS_KM,
        help=(
            "Outer radius for the far-range skyscraper urban outline layer in kilometers "
            f"(default: {SKYSCRAPER_OUTER_RADIUS_KM:.1f}). "
            "Use 0 to disable skyscraper-tile lookup for that run. "
            "Otherwise it must be greater than or equal to --urban-outline-radius-km."
        ),
    )
    parser.add_argument(
        "-b",
        "--urban-outline-min-building-height-m",
        dest="urban_outline_min_height_m",
        type=_parse_non_negative_float,
        default=0.0,
        help=(
            "Deprecated. Minimum building height in meters for buildings included in the urban outline overlay "
            "(default: 0.0). Use --urban-outline-max-candidates for performance tuning."
        ),
    )
    parser.add_argument(
        "--urban-outline-max-candidates",
        type=_parse_non_negative_int,
        default=_URBAN_OUTLINE_MAX_CANDIDATES_DEFAULT,
        help=(
            "Maximum number of urban-outline ring candidates to keep before expensive sampling "
            f"(default: {_URBAN_OUTLINE_MAX_CANDIDATES_DEFAULT}). "
            "Set to 0 to disable the layer."
        ),
    )
    parser.add_argument(
        "--urban-outline-feature-type",
        choices=("both", "building"),
        default="both",
        help=(
            "Urban outline Overture mode "
            "(default: both). Use building to skip building_part overlays."
        ),
    )
    parser.add_argument(
        "--urban-outline-skyscraper-only",
        action="store_true",
        help=(
            "Validation mode: draw only the far-range skyscraper urban outline layer "
            "and skip the normal near-range outline layer."
        ),
    )
    add_general_arguments(
        parser,
        include_window_geometry=include_window_geometry,
        include_window_frame=include_window_frame,
    )
    if include_sky_update_interval:
        parser.add_argument(
            "-i",
            "--sky-update-interval",
            type=int,
            default=60,
            help="Interval for updating stars/sky-color disc in sec. (default: 60).",
        )
    if include_startup_overlay_arguments:
        parser.add_argument(
            "--show-dso-initial",
            type=_parse_bool,
            default=None,
            metavar="true|false",
            help="Whether to show DSO overlays at startup (true/false).",
        )
        parser.add_argument(
            "--show-asterisms-initial",
            type=_parse_bool,
            default=None,
            metavar="true|false",
            help="Whether to show asterism overlays at startup (true/false).",
        )
        parser.add_argument(
            "--show-guidelines-initial",
            type=_parse_bool,
            default=None,
            metavar="true|false",
            help="Whether to show guideline overlays at startup (true/false).",
        )
        parser.add_argument(
            "--observation-info",
            type=str,
            default="bottom",
            choices=("auto", "top", "bottom", "off"),
            metavar="auto|top|bottom|off",
            help=(
                "Observation info overlay mode at startup: auto (hover-avoid), "
                "top (fixed top), bottom (fixed bottom, default), or off (hidden)."
            ),
        )


def add_main_arguments(parser: argparse.ArgumentParser) -> None:
    """Add main-CLI arguments grouped like the README."""
    observing_group = parser.add_argument_group("Observing Location and Time")
    search_group = parser.add_argument_group("Search Objects at startup")
    dataset_group = parser.add_argument_group("Viewpoint dataset queries")
    sky_group = parser.add_argument_group("Sky and Stars")
    overlay_group = parser.add_argument_group("Overlays")
    general_group = parser.add_argument_group("General")

    add_observing_arguments(observing_group)
    add_search_arguments(search_group, include_list=False)
    add_dataset_query_arguments(dataset_group)
    add_sky_and_star_arguments(sky_group)
    add_overlay_arguments(overlay_group)
    add_geo_satellite_argument(general_group)
    add_general_arguments(general_group)


def build_main_argument_parser() -> argparse.ArgumentParser:
    """Build the main zstarview argument parser."""
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    add_main_arguments(parser)
    return parser


def build_export_image_argument_parser() -> argparse.ArgumentParser:
    """Build the headless export-image CLI parser."""
    parser = argparse.ArgumentParser(
        description="Export a zstarview image and exit",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    add_export_image_arguments(parser)
    return parser


def _normalize_location_arguments(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    if args.place is not None:
        args.place = args.place.strip()
        if not args.place:
            parser.error("--place must not be empty")
    if args.place_countrycode is not None:
        args.place_countrycode = args.place_countrycode.strip().lower()
        if not args.place_countrycode:
            parser.error("--place-countrycode must not be empty")
    if args.place_lang is not None:
        args.place_lang = args.place_lang.strip()
        if not args.place_lang:
            parser.error("--place-lang must not be empty")
    if getattr(args, "timezone", None) is not None:
        args.timezone = args.timezone.strip()
        if not args.timezone:
            parser.error("--timezone must not be empty")


def _normalize_dataset_query_arguments(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    for option_name in ("list_viewpoints", "list_viewpoint_names"):
        value = getattr(args, option_name)
        if value is not None:
            normalized = value.strip().lower()
            if normalized not in {"t", "m"}:
                parser.error(f"--{option_name.replace('_', '-')} must be 't' or 'm'")
            setattr(args, option_name, normalized)


def _validate_dataset_query_compatibility(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    def has_non_default(name: str) -> bool:
        return hasattr(args, name) and getattr(args, name) != parser.get_default(name)

    dataset_cli_requested = bool(
        args.list_viewpoints or args.list_viewpoint_names or args.show_viewpoint_json
    )
    if dataset_cli_requested:
        if args.city:
            parser.error(
                "dataset query options cannot be used with the location argument"
            )
        if args.place:
            parser.error("dataset query options cannot be used with --place")

        incompatible_non_default = (
            has_non_default("place_countrycode")
            or has_non_default("place_lang")
            or has_non_default("hours")
            or has_non_default("days")
            or has_non_default("datetime")
            or has_non_default("timezone")
            or has_non_default("vmag_limit")
            or has_non_default("vmag_brightness_multiplier")
            or has_non_default("enlarge_moon")
            or has_non_default("bright_bodies")
            or has_non_default("star_base_radius")
            or has_non_default("expected_render_width")
            or has_non_default("window_geometry")
            or has_non_default("window_frame")
            or has_non_default("view_center_az")
            or has_non_default("view_center_alt")
            or has_non_default("content_fov_deg")
            or has_non_default("observer_height_m")
            or has_non_default("sky_opacity")
            or has_non_default("sky_disc_style")
            or has_non_default("sky_disc_altaz_rings")
            or has_non_default("sky_disc_altaz_rings_hover")
            or has_non_default("cloud_opacity")
            or has_non_default("aircraft_opacity")
            or has_non_default("satellite_opacity")
            or has_non_default("terrain_horizon_opacity")
            or has_non_default("earth_guide_opacity")
            or has_non_default("night_light_opacity")
            or has_non_default("sky_glow_opacity")
            or has_non_default("urban_outline_opacity")
            or has_non_default("water_surface_opacity")
            or has_non_default("urban_outline_radius_km")
            or has_non_default("urban_outline_min_height_m")
            or has_non_default("urban_outline_max_candidates")
            or has_non_default("urban_outline_feature_type")
            or has_non_default("urban_outline_skyscraper_only")
            or has_non_default("clear_long_lived_cache")
            or has_non_default("ground_tint_opacity")
            or has_non_default("visibility_boost")
            or has_non_default("cloud_stripe")
            or has_non_default("cloud_missing_tint_opacity")
            or has_non_default("geo_satellite")
            or has_non_default("sky_update_interval")
            or has_non_default("show_dso_initial")
            or has_non_default("show_asterisms_initial")
            or has_non_default("show_guidelines_initial")
            or has_non_default("include_direction_grid")
            or has_non_default("theme")
        )
        if incompatible_non_default:
            parser.error(
                "dataset query options cannot be used with time or rendering options"
            )


def _validate_location_argument_combinations(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    if args.city and args.place:
        parser.error("--place cannot be used with the location argument")
    if args.place is None and args.place_countrycode is not None:
        parser.error("--place-countrycode requires --place")
    if args.place is None and args.place_lang != "en":
        parser.error("--place-lang requires --place")


def _validate_urban_outline_argument_combinations(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    skyscraper_radius_km = getattr(args, "urban_outline_skyscraper_radius_km", None)
    base_radius_km = getattr(args, "urban_outline_radius_km", None)
    if skyscraper_radius_km is None or base_radius_km is None:
        return
    if float(skyscraper_radius_km) == 0.0:
        return
    if float(skyscraper_radius_km) < float(base_radius_km):
        parser.error(
            "--urban-outline-skyscraper-radius-km must be 0 or greater than or equal to "
            "--urban-outline-radius-km"
        )


def _validate_fov_relationship(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    edge_fov_deg = getattr(args, "edge_fov_deg", None)
    content_fov_deg = getattr(args, "content_fov_deg", None)
    if edge_fov_deg is None or content_fov_deg is None:
        return
    if float(edge_fov_deg) > float(content_fov_deg):
        parser.error("--edge-fov-deg must be less than or equal to --content-fov-deg")


def _validate_main_search_arguments(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    argv: Sequence[str] | None = None,
) -> None:
    if argv is not None:
        search_option_count = sum(
            1 for token in argv if token == "--search" or token.startswith("--search=")
        )
        if search_option_count > 1:
            parser.error("--search may be specified only once")
    if getattr(args, "list", False):
        parser.error("--list is only supported by zstarview-export-image")


def _normalize_vmag_limit(args: argparse.Namespace) -> None:
    if hasattr(args, "vmag_limit"):
        args.vmag_limit = min(float(args.vmag_limit), _COMMITTED_VMAG_LIMIT_MAX)


def _warn_deprecated_urban_outline_min_height_option(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> None:
    if not hasattr(args, "urban_outline_min_height_m"):
        return
    if float(args.urban_outline_min_height_m) == float(parser.get_default("urban_outline_min_height_m")):
        return
    print(
        "warning: --urban-outline-min-building-height-m is deprecated; use --urban-outline-max-candidates for performance tuning",
        file=sys.stderr,
    )


def _argv_has_option(argv: Sequence[str], *option_names: str) -> bool:
    for token in argv:
        for option_name in option_names:
            if token == option_name or token.startswith(f"{option_name}="):
                return True
            if option_name.startswith("-") and not option_name.startswith("--"):
                if token.startswith(option_name) and len(token) > len(option_name):
                    return True
    return False


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the main zstarview app."""
    parser = build_main_argument_parser()
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv)
    _normalize_location_arguments(parser, args)
    _normalize_dataset_query_arguments(parser, args)
    _normalize_vmag_limit(args)
    _validate_dataset_query_compatibility(parser, args)
    _validate_location_argument_combinations(parser, args)
    _validate_urban_outline_argument_combinations(parser, args)
    _warn_deprecated_urban_outline_min_height_option(parser, args)
    _validate_fov_relationship(parser, args)
    _validate_main_search_arguments(parser, args, raw_argv)
    args.view_center_alt_specified = _argv_has_option(
        raw_argv, "-A", "--view-center-alt"
    )
    args.view_center_az_specified = _argv_has_option(
        raw_argv, "-Z", "--view-center-az"
    )

    return args


def parse_export_image_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the export-image CLI."""
    parser = build_export_image_argument_parser()
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv)
    _normalize_location_arguments(parser, args)
    _normalize_vmag_limit(args)
    _validate_location_argument_combinations(parser, args)
    _validate_urban_outline_argument_combinations(parser, args)
    _warn_deprecated_urban_outline_min_height_option(parser, args)
    _validate_fov_relationship(parser, args)
    _validate_main_search_arguments(parser, args, raw_argv)
    if args.print_cache_dir:
        if args.output or args.sixel:
            parser.error("--print-cache-dir cannot be used with --output or --sixel")
        return args
    args.view_center_alt_specified = _argv_has_option(
        raw_argv, "-A", "--view-center-alt"
    )
    args.view_center_az_specified = _argv_has_option(
        raw_argv, "-Z", "--view-center-az"
    )
    if args.list and not getattr(args, "search", None):
        parser.error("--list requires --search")
    if not args.output and not args.sixel and not args.list:
        parser.error("either --output or --sixel is required")
    if args.output == "-" and args.sixel:
        parser.error("--output - cannot be used together with --sixel")
    return args
