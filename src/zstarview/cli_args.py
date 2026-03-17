import argparse
from typing import Sequence
from typing import Tuple, Union

from .paths import CLOUD_MISSING_TINT_RGBA, DIRECTIONS, WINDOW_WIDTH

WindowGeometryArg = Union[str, Tuple[int, int, int, int]]

_VMAG_MULTIPLIER_MIN = 10.0 ** 0.2
_VMAG_MULTIPLIER_MAX = 10.0 ** 0.4


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
    }
    if key in allowed:
        return allowed[key]
    raise argparse.ArgumentTypeError(
        f"Invalid theme: {value!r}. Use one of: night, day, white, black."
    )


def _parse_cloud_stripe(value: str) -> Tuple[int, float]:
    """Parse cloud stripe style as 'count,width_factor' (e.g. '50,0.2')."""
    text = (value or "").strip()
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"Invalid cloud stripe style: {value!r}. Use 'count,width' (e.g. 50,0.2)."
        )
    try:
        count = int(parts[0])
        width = float(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid cloud stripe style: {value!r}. Use numeric 'count,width'."
        ) from exc
    if count < 0 or width < 0:
        raise argparse.ArgumentTypeError(
            f"Invalid cloud stripe style: {value!r}. count and width must be >= 0."
        )
    return count, width


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
        raise argparse.ArgumentTypeError(f"Invalid positive integer: {value!r}") from exc
    if out <= 0:
        raise argparse.ArgumentTypeError("Value must be > 0.")
    return out


def _parse_non_negative_float(value: str) -> float:
    """Parse a float >= 0."""
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(f"Invalid non-negative float: {value!r}") from exc
    if out < 0.0:
        raise argparse.ArgumentTypeError("Value must be >= 0.")
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


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    parser.add_argument(
        "city",
        type=str,
        nargs="?",
        default="",
        help=(
            "Location name: city, lat;lon, tower name, mountain name, "
            "or explicit t/NAME and m/NAME prefixes (default: same as the last run)"
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
            "If TZ is omitted, UTC is assumed. Overrides --hours and --days."
        ),
    )

    parser.add_argument(
        "-V",
        "--vmag-limit",
        type=float,
        default=6.0,
        help="Limit stars to Vmag <= this value (default: 6.0). Use a larger number to show more stars.",
    )
    parser.add_argument(
        "--vmag-brightness-multiplier",
        type=_parse_vmag_brightness_multiplier,
        default=2.5,
        help="Brightness multiplier per magnitude step (allowed range: 1.58–2.512, default: 2.5; 2.512 is the classical Pogson value).",
    )
    parser.add_argument(
        "-m",
        "--enlarge-moon",
        action="store_true",
        help="Show the moon in 5x size.",
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
        "--observer-height-m",
        type=_parse_non_negative_float,
        default=None,
        help=(
            "Observer height above local ground in meters. "
            "Default: 1.7 for city/latlon/mountain and tower height + 1.7 for tower-name input."
        ),
    )

    parser.add_argument(
        "--sky-opacity",
        type=float,
        default=0.2,
        help=(
            "Opacity of the simulated sky-color disc (0.0 - 1.0, default: 0.2). "
            "Set to 0.0 to disable sky-color rendering."
        ),
    )
    parser.add_argument(
        "-c",
        "--cloud-opacity",
        type=float,
        default=0.2,
        help=(
            "Opacity of the clouds (0.0 - 1.0, default: 0.2). "
            "Set to 0.0 to disable cloud rendering."
        ),
    )
    parser.add_argument(
        "--terrain-horizon-opacity",
        type=float,
        default=0.05,
        help=(
            "Opacity of the terrain horizon polyline (0.0 - 1.0, default: 0.05). "
            "Set to 0.0 to disable DEM download, terrain-horizon calculation, and drawing."
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
        "-b",
        "--urban-outline-min-building-height-m",
        dest="urban_outline_min_height_m",
        type=_parse_non_negative_float,
        default=0.0,
        help=(
            "Minimum building height in meters for buildings included in the urban outline overlay "
            "(default: 0.0). This value becomes part of the cache key."
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
            "and skip the normal 0-2.5km outline layer."
        ),
    )
    parser.add_argument(
        "--ground-tint-opacity",
        type=float,
        default=0.1,
        help=(
            "Overlay opacity of the ground tint color below the geometric/terrain horizon "
            "(0.0 - 1.0, default: 0.1)."
        ),
    )
    parser.add_argument(
        "--cloud-stripe",
        type=_parse_cloud_stripe,
        default=(50, 0.2),
        metavar="COUNT,WIDTH",
        help=(
            "Cloud stripe style as 'count,width' (default: 50,0.2). "
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
        "-i",
        "--sky-update-interval",
        type=int,
        default=60,
        help="Interval for updating stars/sky-color disc in sec. (default: 60).",
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
    theme_default = "night"
    parser.add_argument(
        "-t",
        "--theme",
        type=_parse_theme,
        default=theme_default,
        metavar="{night,day,white,black}",
        help="Theme preset for background and star contrast (default: night).",
    )
    args = parser.parse_args(argv)
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
    for option_name in ("list_viewpoints", "list_viewpoint_names"):
        value = getattr(args, option_name)
        if value is not None:
            normalized = value.strip().lower()
            if normalized not in {"t", "m"}:
                parser.error(f"--{option_name.replace('_', '-')} must be 't' or 'm'")
            setattr(args, option_name, normalized)
    dataset_cli_requested = bool(
        args.list_viewpoints
        or args.list_viewpoint_names
        or args.show_viewpoint_json
    )
    if dataset_cli_requested:
        if args.city:
            parser.error("dataset query options cannot be used with the location argument")
        if args.place:
            parser.error("dataset query options cannot be used with --place")

        incompatible_non_default = (
            args.place_countrycode is not None
            or args.place_lang != "en"
            or args.hours != 0
            or args.days != 0
            or args.datetime is not None
            or args.vmag_limit != 6.0
            or args.vmag_brightness_multiplier != 2.5
            or args.enlarge_moon
            or args.star_base_radius != 4.0
            or args.expected_render_width != WINDOW_WIDTH
            or args.window_geometry is not None
            or args.view_center_az != 180.0
            or args.view_center_alt != 90.0
            or args.observer_height_m is not None
            or args.sky_opacity != 0.2
            or args.cloud_opacity != 0.2
            or args.terrain_horizon_opacity != 0.05
            or args.urban_outline_feature_type != "both"
            or args.ground_tint_opacity != 0.1
            or args.cloud_stripe != (50, 0.2)
            or args.cloud_missing_tint_opacity != float(CLOUD_MISSING_TINT_RGBA[3]) / 255.0
            or args.sky_update_interval != 60
            or args.show_dso_initial is not None
            or args.show_asterisms_initial is not None
            or args.theme != theme_default
        )
        if incompatible_non_default:
            parser.error("dataset query options cannot be used with time or rendering options")
    if args.city and args.place:
        parser.error("--place cannot be used with the location argument")
    if args.place is None and args.place_countrycode is not None:
        parser.error("--place-countrycode requires --place")
    if args.place is None and args.place_lang != "en":
        parser.error("--place-lang requires --place")

    return args
