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
import time
from typing import Callable, Optional, Tuple

from PySide6.QtWidgets import QApplication, QSplashScreen
from PySide6.QtCore import Qt
from PySide6.QtGui import (
    QColor,
    QIcon,
    QPixmap,
)
from PySide6.QtGui import QGuiApplication

import pandas as pl

from .paths import (
    APP_ID,
    APP_AUTHOR,
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
        help="Set an absolute date and time in 'YYYY-MM-DD HH:MM:SS [TZ]' format. If TZ is omitted, UTC is assumed. Overrides --hours and --days.",
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
    return parser.parse_args()


class SplashLogHandler(logging.Handler):
    """A temporary log handler to display logs on the splash screen."""

    def __init__(self, show_fn: Callable[[str, QColor], None]):
        """
        Initializes the SplashLogHandler.

        Args:
            show_fn: A function that takes a message string and a QColor
                     and displays it on the splash screen.
        """
        super().__init__()
        self.show_fn = show_fn
        # Use a concise and visible format for the splash screen.
        self.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        """
        Formats and emits a log record to the splash screen.

        Args:
            record: The log record to process.
        """
        try:
            msg = self.format(record)
            # Color-code messages based on log level.
            color = (
                Qt.GlobalColor.white if record.levelno < logging.WARNING else QColor(255, 200, 120) if record.levelno < logging.ERROR else QColor(255, 100, 100)
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
    pixmap = QPixmap(400, 200)
    pixmap.fill(Qt.GlobalColor.black)
    splash.setPixmap(pixmap)
    splash.show()

    def show_splash_message(message: str, color: QColor):
        splash.showMessage(f"{app_name} ver. {__version__}\n{message}", Qt.AlignmentFlag.AlignCenter, color)
        app.processEvents()

    splash_handler = SplashLogHandler(show_splash_message)

    root_logger.addHandler(splash_handler)

    return splash, splash_handler


class StartupAbortError(Exception):
    """Abort the startup sequence (handled by main to show splash for 3s)."""

    ...


def _startup_resolve_city(args_city: Optional[str]) -> CityRec:
    """
    Resolves the target city from arguments or last used city.

    Args:
        args_city: The city name from command-line arguments.

    Returns:
        The resolved city record.

    Raises:
        StartupAbortError: If city data cannot be loaded or the city cannot be found.
    """
    last_city = load_last_city()
    if not args_city:
        args_city = last_city or "Tokyo"

    # load city data
    logger.info("Loading city data...")
    try:
        admin1_map = load_admin1_names(CITY_ADMIN1_CODES_FILE)
    except FileNotFoundError:
        logger.error("Fail to load admin1CodesASCII.txt.")
        raise StartupAbortError()

    recs = []
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

    city_str = f"{city.cc}/{city.name}"
    logger.info(f"City: {city_str}")

    return city


def _startup_parse_time_arguments(args_datetime: Optional[str], args_days: int, args_hours: int) -> timedelta:
    """
    Parses time-related arguments and returns a timedelta from now.

    Args:
        args_datetime: The absolute datetime string.
        args_days: The number of days to add.
        args_hours: The number of hours to add.

    Returns:
        A timedelta object representing the offset from the current UTC time.

    Raises:
        StartupAbortError: If the arguments are invalid.
    """
    if not args_datetime:
        delta_t = timedelta(days=args_days, hours=args_hours)
        return delta_t

    if args_hours != 0 or args_days != 0:
        logger.error("Invalid option: --datetime cannot be used with --hours or --days.")
        raise StartupAbortError()

    try:
        # Split the datetime string to check for a timezone suffix
        parts = args_datetime.split(" ")
        dt_str_naive = " ".join(parts[:2])  # YYYY-MM-DD HH:MM:SS
        tz_str = None
        if len(parts) > 2:  # If there's a third part, it might be the timezone
            tz_str = parts[2]

        dt_naive = datetime.strptime(dt_str_naive, "%Y-%m-%d %H:%M:%S")

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

    except ValueError:
        logger.error(f"Invalid datetime format: {args_datetime}. Use 'YYYY-MM-DD HH:MM:SS [TZ]'.")
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

    app_name = "Zenith Star View"
    app = setup_app(app_name)
    args = parse_args()

    root_logger = setup_root_logger()
    logger.info("Zenith Star View starting...")

    splash, splash_handler = setup_splash_and_attach_logger(app, app_name, root_logger)

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

    save_last_city(city_str)

    sys.exit(app.exec())
