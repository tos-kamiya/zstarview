import sys
import time
import logging
import math
import json

from .cli_args import (
    _parse_theme,
    _parse_window_geometry,
    parse_args,
)
from .mountain_viewpoints import (
    list_mountain_all_names,
    list_mountain_primary_names,
    mountain_viewpoint_to_dict,
    resolve_mountain_viewpoint,
)
from .ui.window_inputs import (
    prepare_window_catalogs,
    prepare_window_runtime_options,
    prepare_window_user_options,
    prepare_window_viewer_data,
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
from .tower_viewpoints import (
    list_tower_all_names,
    list_tower_primary_names,
    resolve_tower_viewpoint,
    tower_viewpoint_to_dict,
)

logger = logging.getLogger(__name__)


def _handle_dataset_query_cli(args: object) -> int | None:
    if getattr(args, "list_towers", False):
        for name in list_tower_primary_names():
            print(name)
        return 0
    if getattr(args, "list_tower_names", False):
        for name in list_tower_all_names():
            print(name)
        return 0
    tower_name = getattr(args, "show_tower_json", None)
    if tower_name:
        tower = resolve_tower_viewpoint(tower_name)
        if tower is None:
            print(f"No tower found for {tower_name!r}", file=sys.stderr)
            return 1
        print(json.dumps(tower_viewpoint_to_dict(tower), ensure_ascii=False, indent=2, sort_keys=True))
        return 0
    if getattr(args, "list_mountains", False):
        for name in list_mountain_primary_names():
            print(name)
        return 0
    if getattr(args, "list_mountain_names", False):
        for name in list_mountain_all_names():
            print(name)
        return 0
    mountain_name = getattr(args, "show_mountain_json", None)
    if mountain_name:
        mountain = resolve_mountain_viewpoint(mountain_name)
        if mountain is None:
            print(f"No mountain found for {mountain_name!r}", file=sys.stderr)
            return 1
        print(json.dumps(mountain_viewpoint_to_dict(mountain), ensure_ascii=False, indent=2, sort_keys=True))
        return 0
    return None


def main() -> None:
    """Main entry point for the star sky visualizer."""
    args = parse_args()
    cli_exit_code = _handle_dataset_query_cli(args)
    if cli_exit_code is not None:
        raise SystemExit(cli_exit_code)

    from .splash import setup_app, setup_splash_and_attach_logger
    from .ui.window import SkyWindow

    app_name = APP_DISPLAY_NAME
    app = setup_app(app_name)
    app.setQuitOnLastWindowClosed(False)

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
    cloud_stripe_count, cloud_stripe_width = args.cloud_stripe
    visual_preset = args.theme
    star_visibility_boost = 1.12 if visual_preset == "white" else 1.05 if visual_preset == "day" else 1.0
    vmag_brightness_scale = -math.log10(args.vmag_brightness_multiplier)
    catalogs = prepare_window_catalogs(
        star_catalog,
        dso_catalog=dso_catalog,
        vmag_brightness_scale=vmag_brightness_scale,
    )

    viewer_data = prepare_window_viewer_data(
        city.display_name,
        (city.lat, city.lon, city.tz),
        view_center,
        observer_height_m=city.observer_height_m if args.observer_height_m is None else args.observer_height_m,
    )
    user_options = prepare_window_user_options(
        sky_disc_alpha=args.sky_opacity,
        cloud_disc_alpha=0.0 if cloud_stripe_count == 0 or cloud_stripe_width == 0.0 else args.cloud_opacity,
        terrain_horizon_opacity=args.terrain_horizon_opacity,
        ground_tint_opacity=args.ground_tint_opacity,
        enlarge_moon=args.enlarge_moon,
        star_base_radius=args.star_base_radius,
        vmag_limit=args.vmag_limit,
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
        show_dso_initial=args.show_dso_initial,
        show_asterisms_initial=args.show_asterisms_initial,
        terrain_horizon_gui_allowed=args.terrain_horizon_opacity > 0.0,
    )
    runtime_options = prepare_window_runtime_options(
        delta_t=delta_t,
        sky_update_interval=args.sky_update_interval,
        cloud_stripe_style=(cloud_stripe_count, cloud_stripe_width),
        cloud_missing_tint_opacity=args.cloud_missing_tint_opacity,
        star_render_expected_width=args.expected_render_width,
        window_geometry_arg=args.window_geometry,
    )
    main_win = SkyWindow(
        viewer_data,
        catalogs,
        user_options=user_options,
        runtime_options=runtime_options,
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
