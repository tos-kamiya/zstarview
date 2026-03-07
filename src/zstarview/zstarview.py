import sys
import time
import logging
import math

from .cli_args import (
    WindowGeometryArg,
    _parse_theme,
    _parse_window_geometry,
    parse_args,
)
from .paths import APP_DISPLAY_NAME
from .startup import (
    StartupAbortError,
    _format_splash_location,
    _startup_load_dso,
    _startup_load_stars,
    _startup_parse_time_arguments,
    _startup_resolve_city,
    setup_root_logger,
)

logger = logging.getLogger(__name__)


def main() -> None:
    """Main entry point for the star sky visualizer."""
    from .splash import setup_app, setup_splash_and_attach_logger
    from .ui.window import SkyWindow

    app_name = APP_DISPLAY_NAME
    app = setup_app(app_name)
    app.setQuitOnLastWindowClosed(False)
    args = parse_args()

    root_logger = setup_root_logger()
    logger.info(f"{APP_DISPLAY_NAME} starting...")

    splash, splash_handler, set_splash_context = setup_splash_and_attach_logger(app, app_name, root_logger, args.theme)

    try:
        city = _startup_resolve_city(args.city)
        set_splash_context(_format_splash_location(city))
        delta_t = _startup_parse_time_arguments(args.datetime, args.days, args.hours)
        star_catalog = _startup_load_stars(args.vmag_limit)
        dso_catalog = _startup_load_dso()
    except StartupAbortError:
        time.sleep(3)
        return

    view_center = (args.view_center_alt, args.view_center_az)
    view_center = (min(90.0, max(0.0, view_center[0])), view_center[1] % 360)
    star_base_radius = max(2.0, args.star_base_radius)
    sky_opacity = min(1.0, max(0.0, args.sky_opacity))
    cloud_opacity = min(1.0, max(0.0, args.cloud_opacity))
    cloud_missing_tint_opacity = min(1.0, max(0.0, args.cloud_missing_tint_opacity))
    cloud_stripe_count, cloud_stripe_width = args.cloud_stripe
    if cloud_stripe_count == 0 or cloud_stripe_width == 0.0:
        cloud_opacity = 0.0
    sky_update_interval = max(1, args.sky_update_interval)
    visual_preset = args.theme
    star_visibility_boost = 1.12 if visual_preset == "white" else 1.05 if visual_preset == "day" else 1.0
    vmag_brightness_scale = -math.log10(args.vmag_brightness_multiplier)

    city_str = f"{city.cc}/{city.name}"
    main_win = SkyWindow(
        city_str,
        (city.lat, city.lon, city.tz),
        view_center,
        star_catalog,
        dso_catalog=dso_catalog,
        delta_t=delta_t,
        sky_disc_alpha=sky_opacity,
        cloud_disc_alpha=cloud_opacity,
        enlarge_moon=args.enlarge_moon,
        star_base_radius=star_base_radius,
        vmag_limit=args.vmag_limit,
        sky_update_interval=sky_update_interval,
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
        vmag_brightness_scale=vmag_brightness_scale,
        cloud_stripe_style=(cloud_stripe_count, cloud_stripe_width),
        cloud_missing_tint_opacity=cloud_missing_tint_opacity,
        star_render_expected_width=args.expected_render_width,
        window_geometry_arg=args.window_geometry,
        show_dso_initial=args.show_dso_initial,
        show_asterisms_initial=args.show_asterisms_initial,
    )

    def _on_initial_loaded():
        """
        Handles the signal that initial data has been loaded.

        Shows the main window, closes the splash screen, and removes the
        temporary splash screen log handler.
        """
        main_win.show()
        app.setQuitOnLastWindowClosed(True)
        splash.close()
        root_logger.removeHandler(splash_handler)

    main_win.initial_data_loaded.connect(_on_initial_loaded)

    sys.exit(app.exec())
