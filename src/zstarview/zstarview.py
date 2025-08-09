# -*- coding: utf-8 -*-
import argparse
from datetime import timedelta
from pathlib import Path
import sys
import time

from PyQt5.QtWidgets import QApplication, QSplashScreen
from PyQt5.QtCore import Qt
from PyQt5.QtGui import (
    QColor,
    QIcon,
    QPixmap,
)

from appdirs import user_cache_dir

from .paths import (
    APP_ID,
    APP_AUTHOR,
    CITY_COORD_FILE,
    STARS_CSV_FILE,
    APP_ICON_FILE,
)
from .config import load_last_city, save_last_city
from .catalog import load_city_coords, load_star_catalog
from .ui.window import SkyWindow

# --- Helper Functions ---
cache_path = Path(user_cache_dir(appname=APP_ID, appauthor=APP_AUTHOR))
cache_path.mkdir(parents=True, exist_ok=True)


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Star sky visualizer")
    parser.add_argument("city", type=str, nargs="?", default="", help="City name (default: same as the last run)")
    parser.add_argument("-H", "--hours", type=float, default=0, help="Number of hours to add to current time (default: 0)")
    parser.add_argument("-D", "--days", type=float, default=0, help="Number of days to add to current time (default: 0)")
    parser.add_argument(
        "-V", "--vmag-limit",
        type=float,
        default=7.0,
        help="Limit stars to Vmag <= this value (default: 7.0). Use a larger number to show more stars."
    )
    parser.add_argument(
        "-m",
        "--enlarge-moon",
        action="store_true",
        help="Show the moon in 3x size.",
    )
    parser.add_argument("-s", "--star-base-radius", type=float, default=15.0, help="Base size of stars (default: 15.0)")
    parser.add_argument("-Z", "--view-center-az", type=float, default=180.0, help="Viewing azimuth angle [deg] (0=N, 90=E, 180=S, 270=W; default=180)")
    parser.add_argument("-A", "--view-center-alt", type=float, default=90.0, help="Viewing altitude angle [deg] (90=zenith, 0=horizon; default=90)")
    return parser.parse_args()


def main():
    """Main entry point for the star sky visualizer."""
    app = QApplication(sys.argv)

    splash = QSplashScreen(QPixmap(400, 200), Qt.WindowType.WindowStaysOnTopHint)
    splash.show()

    def show_splash_message(message: str, color: QColor):
        splash.showMessage(message, Qt.AlignmentFlag.AlignCenter, color)
        app.processEvents()

    try:
        city_table = load_city_coords(CITY_COORD_FILE)
    except FileNotFoundError:
        show_splash_message("Error: cities1000.txt not found.", Qt.GlobalColor.red)
        time.sleep(3)
        return

    app.setWindowIcon(QIcon(APP_ICON_FILE))

    last_city = load_last_city()

    args = parse_args()
    city = args.city or last_city or "Tokyo"
    city = city.lower()
    if city not in city_table:
        candidate_cities = [c for c in city_table.keys() if c.endswith("/" + city)]
        if not candidate_cities:
            print(f"Unknown city: {city}")
            return
        elif len(candidate_cities) > 1:
            print(f"Specify explicit country name: {candidate_cities}")
            return
        else:
            city = candidate_cities[0]

    print(f"City: {city}")

    delta_hours = args.hours
    delta_days = args.days

    delta_t = timedelta(days=delta_days, hours=delta_hours)
    view_center = (args.view_center_alt, args.view_center_az)

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
        star_catalog,
        delta_t,
        enlarge_moon=args.enlarge_moon,
        star_base_radius=args.star_base_radius,
        view_center=view_center,
    )

    # When the initial data is loaded, show the main window and close the splash screen
    main_win.initial_data_loaded.connect(main_win.show)
    main_win.initial_data_loaded.connect(splash.close)

    save_last_city(city)

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
