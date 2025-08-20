import argparse
from datetime import datetime, timedelta, timezone
from pathlib import Path
import re
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
    CITY_ADMIN1_CODES_FILE,
    STARS_CSV_FILE,
    APP_ICON_FILE,
    DIRECTIONS,
)
from .__about__ import __version__
from .config import load_last_city, save_last_city
from .catalog import load_star_catalog
from .ui.window import SkyWindow
from .utils.resolve_city import load_admin1_names, resolve_city, resolve_city_by_name, resolve_city_by_geonameid
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

    bright_red = QColor(255, 100, 100)

    print("Loading city data...", file=sys.stderr)
    try:
        admin1_map = load_admin1_names(CITY_ADMIN1_CODES_FILE)
    except FileNotFoundError:
        print("Error: admin1CodesASCII.txt not found.", file=sys.stderr)
        return

    last_city = load_last_city()

    args = parse_args()
    city = args.city or last_city or "Tokyo"

    recs = []
    try:
        if re.match(r"^\d+$", args.city):
            # If input is just a geonameid, resolve it directly
            rec = resolve_city_by_geonameid(int(args.city), CITY_COORD_FILE)
            if rec:
                recs.append(rec)
            else:
                print(f"Error: No city found for geonameid {args.city}", file=sys.stderr)
                return
        else:
            if not "/" in args.city:
                recs = resolve_city_by_name(args.city, CITY_COORD_FILE, admin1_map)
            else:
                recs = resolve_city(args.city, CITY_COORD_FILE, admin1_map)
            if recs:
                print(f"Found {len(recs)} match(es) for '{args.city}':")
                for rec in recs:
                    print(f"- {rec.cc}/{rec.name}, lat: {rec.lat:.6f}, lon: {rec.lon:.6f}, tz: {rec.tz}  (geonameid={rec.geonameid})")
                if len(recs) > 1:
                    print(f"Warning: Multiple matches found for '{args.city}'", file=sys.stderr)
            else:
                print(f"Error: No match for '{args.city}'", file=sys.stderr)
                return
    except FileNotFoundError:
        print("Error: cities1000.txt not found.", file=sys.stderr)
        return

    city = recs[0]  # Take the city with the highest population
    city_str = f"{city.cc}/{city.name}"
    print(f"City: {city_str}")

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
        show_splash_message(f"Error: star data file not found: {STARS_CSV_FILE}", bright_red)
        time.sleep(3)
        return

    limit_str = args.vmag_limit if args.vmag_limit is not None else "no limit"
    print(f"Loaded {len(star_catalog)} stars (Vmag ≤ {limit_str})")

    show_splash_message(f"Calculating sky for {city_str}...", Qt.GlobalColor.white)

    main_win = SkyWindow(
        city_str,
        (city.lat, city.lon, city.tz),
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

    save_last_city(city_str)

    sys.exit(app.exec())
