import logging

logger = logging.getLogger(__name__)
logging.getLogger("satpy.readers.core.utils").setLevel(
    logging.WARNING
)  # suppress "[INFO] satpy.readers.core.utils: Using temp file for BZ2 decompression: /tmp/..."

import argparse
from datetime import datetime, timedelta, timezone
from pathlib import Path
import re
import sys
import threading
import time
from typing import Callable, List, Optional, Tuple

from PySide6.QtWidgets import QApplication, QSplashScreen
from PySide6.QtCore import Qt
from PySide6.QtGui import (
    QColor,
    QLinearGradient,
    QPainter,
    QIcon,
    QPixmap,
)
from PySide6.QtGui import QGuiApplication

import polars as pl

from .paths import (
    APP_ID,
    APP_AUTHOR,
    APP_DISPLAY_NAME,
    CACHE_PATH,
    LOG_PATH,
    CITY_COORD_FILE,
    CITY_ADMIN1_CODES_FILE,
    STARS_CSV_FILE,
    APP_ICON_FILE,
    DIRECTIONS,
)
from .__about__ import __version__
from .config import load_last_city, save_last_city
from .catalog import load_star_catalog
from .ui.window import SkyWindow
from .utils.resolve_city import CityRec, load_admin1_names, resolve_city, resolve_city_by_name, resolve_city_by_geonameid
from .utils.timezone_parser import parse_tz_string

# --- Helper Functions ---
_cache_path = Path(CACHE_PATH)
_cache_path.mkdir(parents=True, exist_ok=True)

SPLASH_INFO_COLOR = QColor(18, 29, 48)
SPLASH_WARN_COLOR = QColor(130, 82, 20)
SPLASH_ERROR_COLOR = QColor(146, 34, 34)


def _get_splash_palette(visual_preset: str) -> tuple[list[QColor], QColor, QColor]:
    """Return gradient colors, frame color, and default message color for splash."""
    if visual_preset == "night":
        return (
            [QColor(16, 20, 34), QColor(10, 14, 25), QColor(6, 9, 17)],
            QColor(66, 84, 118),
            QColor(214, 228, 255),
        )
    if visual_preset == "white":
        return (
            [QColor(244, 250, 255), QColor(224, 236, 250), QColor(196, 214, 238)],
            QColor(158, 178, 206),
            QColor(19, 31, 50),
        )
    # day
    return (
        [QColor(236, 244, 255), QColor(214, 227, 246), QColor(186, 204, 232)],
        QColor(148, 167, 194),
        QColor(18, 29, 48),
    )


def _build_splash_pixmap(visual_preset: str, width: int = 400, height: int = 200) -> QPixmap:
    """Create a splash background matching the selected visual preset."""
    grad_colors, frame_color, _ = _get_splash_palette(visual_preset)
    pixmap = QPixmap(width, height)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    grad = QLinearGradient(0, 0, width, height)
    grad.setColorAt(0.0, grad_colors[0])
    grad.setColorAt(0.55, grad_colors[1])
    grad.setColorAt(1.0, grad_colors[2])
    painter.fillRect(0, 0, width, height, grad)

    # Subtle frame to separate splash from desktop background.
    painter.setPen(frame_color)
    painter.drawRect(0, 0, width - 1, height - 1)
    painter.end()
    return pixmap


def _parse_azimuth(value: str) -> float:
    """Parse azimuth given as degrees or compass points.

    Examples:
      - "180" -> 180.0
      - "E" -> 90.0
      - "NE" -> 45.0
      - "WNW" -> 292.5
    """
    try:
        return float(value)
    except (TypeError, ValueError):
        pass

    compass = value.strip().upper()
    if compass in DIRECTIONS:
        return float(DIRECTIONS[compass])
    raise argparse.ArgumentTypeError(f"Invalid azimuth: {value!r}. Use degrees (e.g., 180) or compass (e.g., N, NE, E).")


def _parse_theme(value: str) -> str:
    """Parse theme preset."""
    key = (value or "").strip().lower()
    allowed = {
        "night": "night",
        "day": "day",
        "white": "white",
    }
    if key in allowed:
        return allowed[key]
    raise argparse.ArgumentTypeError(
        f"Invalid theme: {value!r}. Use one of: night, day, white."
    )


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    parser.add_argument("city", type=str, nargs="?", default="", help="City name (default: same as the last run)")
    time_group = parser.add_argument_group("Time settings")
    time_group.add_argument("-H", "--hours", type=float, default=0, help="Number of hours to add to current time (default: 0)")
    time_group.add_argument("-D", "--days", type=float, default=0, help="Number of days to add to current time (default: 0)")
    time_group.add_argument(
        "--datetime",
        type=str,
        default=None,
        help=("Set an absolute date and time in 'YYYY-MM-DD HH[:MM[:SS]] [TZ]' format "
            "(e.g., '2025-09-12 9', '2025-09-12 09:00', '2025-09-12 9:0:0 JST'). "
            "If TZ is omitted, UTC is assumed. Overrides --hours and --days."),
    )

    parser.add_argument(
        "-V",
        "--vmag-limit",
        type=float,
        default=6.0,
        help="Limit stars to Vmag <= this value (default: 6.0). Use a larger number to show more stars.",
    )
    parser.add_argument(
        "-m",
        "--enlarge-moon",
        action="store_true",
        help="Show the moon in 5x size.",
    )
    parser.add_argument("-s", "--star-base-radius", type=float, default=8.0, help="Base size of stars (default: 8.0)")
    parser.add_argument(
        "-Z",
        "--view-center-az",
        type=_parse_azimuth,
        default=180.0,
        help=("Viewing azimuth angle [deg or compass] " "(0=N, 90=E, 180=S, 270=W; also accepts N, NE, E, SE, S, SW, W, NW; default=180)"),
    )
    parser.add_argument(
        "-A",
        "--view-center-alt",
        type=float,
        default=90.0,
        help="Viewing altitude angle [deg] (90=zenith, 0=horizon; default=90)",
    )

    parser.add_argument(
        "--sky-opacity",
        type=float,
        default=0.2,
        help=("Opacity of the simulated sky-color disc (0.0 - 1.0, default: 0.2). " "Set to 0.0 to disable sky-color rendering."),
    )
    parser.add_argument(
        "-c",
        "--cloud-opacity",
        type=float,
        default=0.2,
        help=("Opacity of the clouds (0.0 - 1.0, default: 0.2). " "Set to 0.0 to disable cloud rendering."),
    )
    parser.add_argument(
        "-i",
        "--sky-update-interval",
        type=int,
        default=3 * 60,
        help=("Interval for updating stars/sky-color disc in sec. (default: 180)."),
    )
    parser.add_argument(
        "-t",
        "--theme",
        type=_parse_theme,
        default="night",
        metavar="{night,day,white}",
        help="Theme preset for background and star contrast (default: night).",
    )
    return parser.parse_args()


class SplashLogHandler(logging.Handler):
    """A temporary log handler to display logs on the splash screen."""

    def __init__(self, show_fn: Callable[[str, QColor], None], info_color: QColor):
        """
        Initializes the SplashLogHandler.

        Args:
            show_fn: A function that takes a message string and a QColor
                     and displays it on the splash screen.
        """
        super().__init__()
        self.show_fn = show_fn
        self._main_thread_id = threading.get_ident()
        self._info_color = info_color
        # Use a concise and visible format for the splash screen.
        self.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        """
        Formats and emits a log record to the splash screen.

        Args:
            record: The log record to process.
        """
        try:
            # QSplashScreen painting must run on the main thread.
            # Background worker logs are still handled by other log handlers.
            if threading.get_ident() != self._main_thread_id:
                return
            msg = self.format(record)
            # Color-code messages based on log level.
            color = (
                self._info_color if record.levelno < logging.WARNING else SPLASH_WARN_COLOR if record.levelno < logging.ERROR else SPLASH_ERROR_COLOR
            )
            self.show_fn(msg, color)
        except Exception:
            self.handleError(record)


def setup_root_logger() -> logging.Logger:
    """
    Configures and returns the root logger for the application.

    Sets up logging to both stderr and a file (`app.log`).

    Returns:
        The configured root logger instance.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stderr,
    )
    root_logger = logging.getLogger()

    log_dir = Path(LOG_PATH)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "app.log"

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    root_logger.addHandler(file_handler)

    logger.info("Logging to file: %s", log_path)

    return root_logger


def setup_app(app_name: str) -> QApplication:
    """
    Initializes and configures the QApplication instance.

    Args:
        app_name: The name of the application.

    Returns:
        The configured QApplication instance.
    """
    app = QApplication(sys.argv)
    QGuiApplication.setDesktopFileName(APP_ID)
    app.setApplicationName(app_name)
    app.setApplicationDisplayName(app_name)
    app.setOrganizationName(APP_AUTHOR)
    app.setApplicationVersion(__version__)
    app.setWindowIcon(QIcon(APP_ICON_FILE))
    return app


def setup_splash_and_attach_logger(
    app: QApplication,
    app_name: str,
    root_logger: logging.Logger,
    visual_preset: str,
) -> Tuple[QSplashScreen, SplashLogHandler]:
    """
    Creates a splash screen and attaches a log handler to it.

    This allows startup logs to be displayed on the splash screen.

    Args:
        app: The QApplication instance.
        app_name: The name of the application.
        root_logger: The root logger to attach the splash handler to.

    Returns:
        A tuple containing the QSplashScreen and SplashLogHandler instances.
    """
    splash = QSplashScreen(QPixmap(400, 200), Qt.WindowType.WindowStaysOnTopHint)
    _, _, info_color = _get_splash_palette(visual_preset)
    pixmap = _build_splash_pixmap(visual_preset, 400, 200)
    splash.setPixmap(pixmap)
    splash.show()

    def show_splash_message(message: str, color: QColor):
        splash.showMessage(f"{app_name} ver. {__version__}\n{message}", Qt.AlignmentFlag.AlignCenter, color)
        app.processEvents()

    splash_handler = SplashLogHandler(show_splash_message, info_color)

    root_logger.addHandler(splash_handler)

    return splash, splash_handler


class StartupAbortError(Exception):
    """Abort the startup sequence (handled by main to show splash for 3s)."""

    ...


def _startup_resolve_city(args_city: Optional[str]) -> CityRec:
    """
    Resolves the target city from arguments or last used city.
    Also handles direct latitude/longitude input like "35.68;139.76" or "N35.68;E139.76".

    Args:
        args_city: The city name or lat;lon string from command-line arguments.

    Returns:
        The resolved city record.

    Raises:
        StartupAbortError: If data cannot be loaded or the location cannot be resolved.
    """
    last_city = load_last_city()
    if not args_city:
        args_city = last_city or "Tokyo"

    # Handle direct latitude/longitude input
    if ";" in args_city:
        try:
            lat_str, lon_str = [s.strip() for s in args_city.split(";")]

            def _parse_coord(s: str, dirs: str) -> float:
                s_upper = s.strip().upper()

                # Extract direction letters (N, S, E, W) that appear in the string
                found = {ch for ch in s_upper if ch in "NSEW"}
                allowed = set(dirs)

                # Ensure any direction letters used are valid for the current axis
                if found and not found.issubset(allowed):
                    raise ValueError(f"Invalid direction in '{s}' (expected one of {sorted(allowed)}).")

                # Determine sign: S or W implies negative
                sign = -1.0 if (("S" in found) or ("W" in found)) else 1.0

                # Remove non-numeric characters (keep digits, dot, minus)
                val_str = re.sub(r"[^0-9.-]", "", s)
                if not val_str:
                    raise ValueError(f"No numeric value found in '{s}'")
                val = float(val_str)

                # Explicit negative number takes precedence over direction sign
                return val if val < 0 else val * sign

            lat = _parse_coord(lat_str, "NS")
            lon = _parse_coord(lon_str, "EW")

            if not (-90 <= lat <= 90):
                raise ValueError(f"Latitude out of range (-90 to 90): {lat}")
            if not (-180 <= lon <= 180):
                raise ValueError(f"Longitude out of range (-180 to 180): {lon}")

            logger.info(f"Parsed location: Lat={lat:.6f}, Lon={lon:.6f}, Timezone=UTC")

            name = f"Lat: {lat:.2f}, Lon: {lon:.2f}"
            return CityRec(
                geonameid=0,
                name=name,
                asciiname=name,
                lat=lat,
                lon=lon,
                cc="",
                admin1_code="",
                admin1_name=None,
                pop=0,
                tz="UTC",
            )
        except (ValueError, IndexError) as e:
            logger.error(f"Invalid latitude/longitude format: '{args_city}'. {e}")
            raise StartupAbortError()

    # --- City name resolution logic ---
    logger.info("Loading city data...")
    try:
        admin1_map = load_admin1_names(CITY_ADMIN1_CODES_FILE)
    except FileNotFoundError:
        logger.error("Fail to load admin1CodesASCII.txt.")
        raise StartupAbortError()

    recs: List[CityRec] = []
    try:
        if re.match(r"^\d+$", args_city):
            # If input is just a geonameid, resolve it directly
            rec = resolve_city_by_geonameid(int(args_city), CITY_COORD_FILE)
            if rec:
                recs.append(rec)
            else:
                logger.error(f"No city found for geonameid {args_city}")
                raise StartupAbortError()
        else:
            if not "/" in args_city:
                recs = resolve_city_by_name(args_city, CITY_COORD_FILE, admin1_map)
            else:
                recs = resolve_city(args_city, CITY_COORD_FILE, admin1_map)
            if recs:
                logger.info(f"Found {len(recs)} match(es) for '{args_city}':")
                for rec in recs:
                    logger.info(f"- {rec.cc}/{rec.name}, lat: {rec.lat:.6f}, lon: {rec.lon:.6f}, tz: {rec.tz}  (geonameid={rec.geonameid})")
                if len(recs) > 1:
                    logger.warning(f"Multiple matches found for '{args_city}'")
            else:
                logger.error(f"No match for '{args_city}'")
                raise StartupAbortError()
    except FileNotFoundError:
        logger.error("Fail to load cities1000.txt.")
        raise StartupAbortError()

    city = recs[0]  # Take the city with the highest population

    # Save city/coordinates for next launch
    if getattr(city, "geonameid", 0) == 0:
        # Coordinate-based input: save as "lat;lon" string
        city_str = f"{city.lat:.6f};{city.lon:.6f}"
    else:
        # Named city: save as "CC/Name"
        city_str = f"{city.cc}/{city.name}"
    save_last_city(city_str)
    logger.info(f"City: {city_str}")

    return city


def _parse_flexible_time(time_str: str) -> Tuple[int, int, int]:
    """
    Parse time string that may omit minutes/seconds.
    Accepts:
      - "H" or "HH"
      - "H:M" or "HH:MM"
      - "H:M:S" or "HH:MM:SS"
    Returns (hour, minute, second)
    """
    m = re.fullmatch(r"\s*(\d{1,2})(?::(\d{1,2}))?(?::(\d{1,2}))?\s*", time_str)
    if not m:
        raise ValueError(f"Invalid time: {time_str!r}. Use HH, HH:MM, or HH:MM:SS.")
    h = int(m.group(1))
    mi = int(m.group(2)) if m.group(2) is not None else 0
    s = int(m.group(3)) if m.group(3) is not None else 0

    if not (0 <= h <= 23):
        raise ValueError(f"Hour out of range (0-23): {h}")
    if not (0 <= mi <= 59):
        raise ValueError(f"Minute out of range (0-59): {mi}")
    if not (0 <= s <= 59):
        raise ValueError(f"Second out of range (0-59): {s}")
    return h, mi, s



# TODO: This function name is confusingly similar to _startup_parse_time_arguments; consider renaming.
def _startup_parse_time_arguments(args_datetime: Optional[str], args_days: int, args_hours: int) -> timedelta:
    """
    Parses time-related arguments and returns a timedelta from now.

    Supports --datetime in 'YYYY-MM-DD HH[:MM[:SS]] [TZ]' (TZ optional).
    """
    if not args_datetime:
        delta_t = timedelta(days=args_days, hours=args_hours)
        return delta_t

    if args_hours != 0 or args_days != 0:
        logger.error("Invalid option: --datetime cannot be used with --hours or --days.")
        raise StartupAbortError()

    try:
        parts = args_datetime.split()
        if len(parts) < 2:
            raise ValueError("Missing time part. Use 'YYYY-MM-DD HH[:MM[:SS]] [TZ]'.")

        date_str = parts[0]
        time_str = parts[1]
        tz_str = parts[2] if len(parts) >= 3 else None

        # parse date (strict)
        try:
            date_only = datetime.strptime(date_str, "%Y-%m-%d")
        except ValueError:
            raise ValueError(f"Invalid date: {date_str!r}. Use YYYY-MM-DD.")

        # parse flexible time
        h, mi, s = _parse_flexible_time(time_str)

        dt_naive = date_only.replace(hour=h, minute=mi, second=s)

        if tz_str:
            try:
                tz = parse_tz_string(tz_str)
                dt_local = dt_naive.replace(tzinfo=tz)
                target_time_utc = dt_local.astimezone(timezone.utc)
            except Exception as e:
                logger.error(f"Invalid timezone '{tz_str}'. {e}")
                raise StartupAbortError()
        else:
            # If no timezone specified, assume UTC
            target_time_utc = dt_naive.replace(tzinfo=timezone.utc)

        now_utc = datetime.now(timezone.utc)
        delta_t = target_time_utc - now_utc

    except ValueError as e:
        logger.error(f"{e} Input was: {args_datetime!r}")
        raise StartupAbortError()

    return delta_t


def _startup_load_stars(args_vmag_limit: Optional[float]) -> pl.DataFrame:
    """
    Loads the star catalog from the source file.

    Args:
        args_vmag_limit: The magnitude limit for visible stars.

    Returns:
        A Polars DataFrame containing the star catalog.

    Raises:
        StartupAbortError: If the star catalog file cannot be found.
    """
    logger.info("Loading city and star data...")

    try:
        star_catalog = load_star_catalog(STARS_CSV_FILE, vmag_threshold=args_vmag_limit)
    except FileNotFoundError:
        logger.error(f"Fail to load file: {STARS_CSV_FILE}")
        raise StartupAbortError()

    limit_str = args_vmag_limit if args_vmag_limit is not None else "no limit"
    logger.info(f"Loaded {len(star_catalog)} stars (Vmag ≤ {limit_str})")
    return star_catalog


def main() -> None:
    """Main entry point for the star sky visualizer."""

    app_name = APP_DISPLAY_NAME
    app = setup_app(app_name)
    args = parse_args()

    root_logger = setup_root_logger()
    logger.info(f"{APP_DISPLAY_NAME} starting...")

    splash, splash_handler = setup_splash_and_attach_logger(app, app_name, root_logger, args.theme)

    try:
        city = _startup_resolve_city(args.city)
        delta_t = _startup_parse_time_arguments(args.datetime, args.days, args.hours)
        star_catalog = _startup_load_stars(args.vmag_limit)
    except StartupAbortError:
        time.sleep(3)
        return

    view_center = (args.view_center_alt, args.view_center_az)
    view_center = (min(90.0, max(0.0, view_center[0])), view_center[1] % 360)
    star_base_radius = max(2.0, args.star_base_radius)
    sky_opacity = min(1.0, max(0.0, args.sky_opacity))
    cloud_opacity = min(1.0, max(0.0, args.cloud_opacity))
    sky_update_interval = max(1, args.sky_update_interval)
    visual_preset = args.theme
    star_visibility_boost = 1.12 if visual_preset == "white" else 1.05 if visual_preset == "day" else 1.0

    city_str = f"{city.cc}/{city.name}"
    main_win = SkyWindow(
        city_str,
        (city.lat, city.lon, city.tz),
        view_center,
        star_catalog,
        delta_t,
        sky_disc_alpha=sky_opacity,
        cloud_disc_alpha=cloud_opacity,
        enlarge_moon=args.enlarge_moon,
        star_base_radius=star_base_radius,
        vmag_limit=args.vmag_limit,
        sky_update_interval=sky_update_interval,
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
    )

    def _on_initial_loaded():
        """
        Handles the signal that initial data has been loaded.

        Shows the main window, closes the splash screen, and removes the
        temporary splash screen log handler.
        """
        main_win.show()
        splash.close()
        root_logger.removeHandler(splash_handler)

    main_win.initial_data_loaded.connect(_on_initial_loaded)

    sys.exit(app.exec())
