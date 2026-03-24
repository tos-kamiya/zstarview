import argparse
from typing import Sequence
from typing import Tuple, Union

from ..paths import CLOUD_MISSING_TINT_RGBA, DIRECTIONS, WINDOW_HEIGHT, WINDOW_WIDTH

WindowGeometryArg = Union[str, Tuple[int, int, int, int]]
ImageSizeArg = Tuple[int, int]

_VMAG_MULTIPLIER_MIN = 10.0 ** 0.2
_VMAG_MULTIPLIER_MAX = 10.0 ** 0.4
_CONTENT_FOV_MIN = 90.0
_CONTENT_FOV_MAX = 127.0
_COMMITTED_VMAG_LIMIT_MAX = 10.5


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
        "transparent": "transparent",
    }
    if key in allowed:
        return allowed[key]
    raise argparse.ArgumentTypeError(
        f"Invalid theme: {value!r}. Use one of: night, day, white, black, transparent."
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
            "Location name: city, lat;lon, @lat,lon, supported Google Maps URL, "
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


def add_render_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_window_geometry: bool = True,
    include_sky_update_interval: bool = True,
    include_startup_overlay_arguments: bool = True,
) -> None:
    """Add shared view and rendering arguments."""
    parser.add_argument(
        "-V",
        "--vmag-limit",
        type=float,
        default=6.0,
        help=(
            "Limit stars to Vmag <= this value (default: 6.0). "
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
        default=100.0,
        help=(
            "Shared overscan content FOV in degrees (allowed: 90-127, default: 100). "
            "The window edge remains fixed at 90 degrees from the view center."
        ),
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
        default=0.15,
        help=(
            "Opacity of the simulated sky-color disc (0.0 - 1.0, default: 0.15). "
            "Set to 0.0 to disable sky-color rendering."
        ),
    )
    parser.add_argument(
        "-c",
        "--cloud-opacity",
        type=float,
        default=0.15,
        help=(
            "Opacity of the clouds (0.0 - 1.0, default: 0.15). "
            "Set to 0.0 to disable cloud rendering."
        ),
    )
    parser.add_argument(
        "-a",
        "--aircraft-opacity",
        type=float,
        default=0.4,
        help=(
            "Opacity of the aircraft overlay (0.0 - 1.0, default: 0.5). "
            "Set to 0.0 to disable aircraft queries and rendering."
        ),
    )
    parser.add_argument(
        "--satellite-opacity",
        type=float,
        default=0.5,
        help=(
            "Opacity of the artificial satellite overlay (0.0 - 1.0, default: 0.5). "
            "Set to 0.0 to disable satellite element fetch and rendering."
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
    theme_default = "night"
    parser.add_argument(
        "-t",
        "--theme",
        type=_parse_theme,
        default=theme_default,
        metavar="{night,day,white,black,transparent}",
        help="Theme preset for background and star contrast (default: night).",
    )


def build_main_argument_parser() -> argparse.ArgumentParser:
    """Build the main zstarview argument parser."""
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    add_location_arguments(parser)
    add_dataset_query_arguments(parser)
    add_time_arguments(parser)
    add_render_arguments(parser)
    return parser


def build_export_image_argument_parser() -> argparse.ArgumentParser:
    """Build the headless export-image CLI parser."""
    parser = argparse.ArgumentParser(
        description="Export a zstarview image and exit",
    )
    add_location_arguments(parser)
    add_time_arguments(parser)
    add_render_arguments(
        parser,
        include_window_geometry=False,
        include_sky_update_interval=False,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        metavar="PATH",
        help="Output image file path.",
    )
    parser.add_argument(
        "--image-size",
        type=_parse_image_size,
        default=(WINDOW_WIDTH, WINDOW_HEIGHT),
        metavar="WIDTH,HEIGHT",
        help=(
            "Output image size in pixels (default: "
            f"{WINDOW_WIDTH},{WINDOW_HEIGHT})."
        ),
    )
    parser.add_argument(
        "--layer-timeout-seconds",
        type=_parse_non_negative_float,
        default=90.0,
        metavar="SECONDS",
        help="Maximum time to wait for enabled external layers (default: 90).",
    )
    parser.add_argument(
        "--allow-partial-data",
        action="store_true",
        help="Allow saving an image even when enabled external layers fail or time out.",
    )
    parser.add_argument(
        "--sixel",
        action="store_true",
        help="Display the generated image in the terminal via img2sixel.",
    )
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
            has_non_default("place_countrycode")
            or has_non_default("place_lang")
            or has_non_default("hours")
            or has_non_default("days")
            or has_non_default("datetime")
            or has_non_default("timezone")
            or has_non_default("vmag_limit")
            or has_non_default("vmag_brightness_multiplier")
            or has_non_default("enlarge_moon")
            or has_non_default("star_base_radius")
            or has_non_default("expected_render_width")
            or has_non_default("window_geometry")
            or has_non_default("view_center_az")
            or has_non_default("view_center_alt")
            or has_non_default("content_fov_deg")
            or has_non_default("observer_height_m")
            or has_non_default("sky_opacity")
            or has_non_default("cloud_opacity")
            or has_non_default("aircraft_opacity")
            or has_non_default("satellite_opacity")
            or has_non_default("terrain_horizon_opacity")
            or has_non_default("urban_outline_opacity")
            or has_non_default("urban_outline_radius_km")
            or has_non_default("urban_outline_min_height_m")
            or has_non_default("urban_outline_feature_type")
            or has_non_default("urban_outline_skyscraper_only")
            or has_non_default("ground_tint_opacity")
            or has_non_default("cloud_stripe")
            or has_non_default("cloud_missing_tint_opacity")
            or has_non_default("sky_update_interval")
            or has_non_default("show_dso_initial")
            or has_non_default("show_asterisms_initial")
            or has_non_default("theme")
        )
        if incompatible_non_default:
            parser.error("dataset query options cannot be used with time or rendering options")


def _validate_location_argument_combinations(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    if args.city and args.place:
        parser.error("--place cannot be used with the location argument")
    if args.place is None and args.place_countrycode is not None:
        parser.error("--place-countrycode requires --place")
    if args.place is None and args.place_lang != "en":
        parser.error("--place-lang requires --place")


def _normalize_vmag_limit(args: argparse.Namespace) -> None:
    if hasattr(args, "vmag_limit"):
        args.vmag_limit = min(float(args.vmag_limit), _COMMITTED_VMAG_LIMIT_MAX)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the main zstarview app."""
    parser = build_main_argument_parser()
    args = parser.parse_args(argv)
    _normalize_location_arguments(parser, args)
    _normalize_dataset_query_arguments(parser, args)
    _normalize_vmag_limit(args)
    _validate_dataset_query_compatibility(parser, args)
    _validate_location_argument_combinations(parser, args)

    return args


def parse_export_image_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the export-image CLI."""
    parser = build_export_image_argument_parser()
    args = parser.parse_args(argv)
    _normalize_location_arguments(parser, args)
    _normalize_vmag_limit(args)
    _validate_location_argument_combinations(parser, args)
    if not args.output and not args.sixel:
        parser.error("either --output or --sixel is required")
    if args.output == "-" and args.sixel:
        parser.error("--output - cannot be used together with --sixel")
    return args
