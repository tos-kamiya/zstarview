import argparse
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


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    parser.add_argument(
        "city",
        type=str,
        nargs="?",
        default="",
        help="Location name: city, lat;lon, or tower name (default: same as the last run)",
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
            "Default: 1.7 for city/latlon, tower height for tower-name input."
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
        default=0.25,
        help=(
            "Opacity of the terrain horizon polyline (0.0 - 1.0, default: 0.25). "
            "Set to 0.0 to disable DEM download, terrain-horizon calculation, and drawing."
        ),
    )
    parser.add_argument(
        "--ground-tint-opacity",
        type=float,
        default=0.2,
        help=(
            "Overlay opacity of the ground tint color below the geometric/terrain horizon "
            "(0.0 - 1.0, default: 0.2)."
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
    parser.add_argument(
        "-t",
        "--theme",
        type=_parse_theme,
        default="night",
        metavar="{night,day,white,black}",
        help="Theme preset for background and star contrast (default: night).",
    )
    return parser.parse_args()
