import argparse
from datetime import datetime, timedelta, timezone
from pathlib import Path
import sys
import time

from PySide6.QtWidgets import QApplication, QSplashScreen
from PySide6.QtCore import Qt
from PySide6.QtGui import (
    QColor,
    QIcon,
    QPixmap,
)
from PySide6.QtGui import QGuiApplication

from appdirs import user_cache_dir

from .paths import (
    APP_ID,
    APP_AUTHOR,
    CITY_COORD_FILE,
    STARS_CSV_FILE,
    APP_ICON_FILE,
    DIRECTIONS,
)
from .__about__ import __version__
from .config import load_last_city, save_last_city
from .catalog import load_city_coords, load_star_catalog
from .ui.window import SkyWindow
from .utils.timezone_parser import parse_tz_string

# --- Helper Functions ---
cache_path = Path(user_cache_dir(appname=APP_ID, appauthor=APP_AUTHOR))
cache_path.mkdir(parents=True, exist_ok=True)


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
    raise argparse.ArgumentTypeError(
        f"Invalid azimuth: {value!r}. Use degrees (e.g., 180) or compass (e.g., N, NE, E)."
    )


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    parser.add_argument("city", type=str, nargs="?", default="", help="City name (default: same as the last run)")
    time_group = parser.add_argument_group("Time settings")
    time_group.add_argument(
        "-H", "--hours", type=float, default=0, help="Number of hours to add to current time (default: 0)"
    )
    time_group.add_argument(
        "-D", "--days", type=float, default=0, help="Number of days to add to current time (default: 0)"
    )
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
        "--sky-opacity",
        type=float,
        default=0.25,
        help=(
            "Opacity of the simulated sky-color disc (0.0 - 1.0, default: 0.25). "
            "Set to 0.0 to disable sky-color rendering."
        ),
    )
    return parser.parse_args()


def main():
    """Main entry point for the star sky visualizer."""
    app = QApplication(sys.argv)
    app_name = "Zenith Star View"
    # Ensure GNOME/Wayland associates the window with our desktop file
    QGuiApplication.setDesktopFileName(APP_ID)
    app.setApplicationName(app_name)
    app.setApplicationDisplayName(app_name)
    app.setOrganizationName(APP_AUTHOR)
    app.setApplicationVersion(__version__)
    app.setWindowIcon(QIcon(APP_ICON_FILE))

    splash = QSplashScreen(QPixmap(400, 200), Qt.WindowType.WindowStaysOnTopHint)
    splash.show()

    def show_splash_message(message: str, color: QColor):
        splash.showMessage(f"{app_name} ver. {__version__}\n{message}", Qt.AlignmentFlag.AlignCenter, color)
        app.processEvents()

    try:
        city_table = load_city_coords(CITY_COORD_FILE)
    except FileNotFoundError:
        show_splash_message("Error: cities1000.txt not found.", Qt.GlobalColor.red)
        time.sleep(3)
        return

    last_city = load_last_city()

    args = parse_args()
    city = args.city or last_city or "Tokyo"
    city = city.lower()
    if city not in city_table:
        candidate_cities = [c for c in city_table.keys() if c.endswith("/" + city)]
        if not candidate_cities:
            candidate_cities = [c for c in city_table.keys() if city in c]
        if not candidate_cities:
            print(f"Unknown city: {city}")
            return
        elif len(candidate_cities) > 1:
            print(f"Specify explicit country name: {candidate_cities}")
            return
        else:
            city = candidate_cities[0]

    print(f"City: {city}")

    sky_opacity = min(1.0, max(0.0, args.sky_opacity))
    star_base_radius = max(2.0, args.star_base_radius)

    # --- Determine the time for calculation ---
    if args.datetime:
        if args.hours != 0 or args.days != 0:
            print("Error: --datetime cannot be used with --hours or --days.", file=sys.stderr)
            return

        try:
            # Split the datetime string to check for a timezone suffix
            parts = args.datetime.split(" ")
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
                    print(f"Error: Invalid timezone '{tz_str}'. {e}", file=sys.stderr)
                    return
            else:
                # If no timezone specified, assume UTC
                target_time_utc = dt_naive.replace(tzinfo=timezone.utc)

            now_utc = datetime.now(timezone.utc)
            delta_t = target_time_utc - now_utc

        except ValueError:
            print(
                f"Error: Invalid datetime format: {args.datetime}. Use 'YYYY-MM-DD HH:MM:SS [TZ]'.",
                file=sys.stderr,
            )
            return
    else:
        delta_t = timedelta(days=args.days, hours=args.hours)

    view_center = (args.view_center_alt, args.view_center_az)
    view_center = (min(90.0, max(0.0, view_center[0])), view_center[1] % 360)

    # Show a splash screen while loading initial data
    pixmap = QPixmap(400, 200)
    pixmap.fill(Qt.GlobalColor.black)
    splash.setPixmap(pixmap)
    show_splash_message("Loading city and star data...", Qt.GlobalColor.white)

    try:
        star_catalog = load_star_catalog(STARS_CSV_FILE, vmag_threshold=args.vmag_limit)
    except FileNotFoundError:
        show_splash_message(f"Error: star data file not found: {STARS_CSV_FILE}", Qt.GlobalColor.red)
        time.sleep(3)
        return

    limit_str = args.vmag_limit if args.vmag_limit is not None else "no limit"
    print(f"Loaded {len(star_catalog)} stars (Vmag ≤ {limit_str})")

    show_splash_message(f"Calculating sky for {city.title()}...", Qt.GlobalColor.white)

    lat, lon, tz_name = city_table[city]
    main_win = SkyWindow(
        city,
        (lat, lon, tz_name),
        view_center,
        star_catalog,
        delta_t,
        sky_disc_alpha=sky_opacity,
        enlarge_moon=args.enlarge_moon,
        star_base_radius=star_base_radius,
        vmag_limit=args.vmag_limit,
    )

    # When the initial data is loaded, show the main window and close the splash screen
    main_win.initial_data_loaded.connect(main_win.show)
    main_win.initial_data_loaded.connect(splash.close)

    save_last_city(city)

    sys.exit(app.exec())
