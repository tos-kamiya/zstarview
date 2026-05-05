import sys
import time
import logging
import math
import json
from dataclasses import replace

from PySide6.QtCore import QTimer

from ..astro import load_ephemeris
from ..cache_maintenance import LongLivedCacheClearCooldownError, clear_long_lived_cache
from ..catalog import load_dso_catalog, load_star_catalog
from ..location_resolver import (
    format_splash_location,
    find_exact_viewpoint_matches,
    list_mountain_all_names,
    list_mountain_primary_names,
    list_tower_all_names,
    list_tower_primary_names,
    load_mountain_viewpoints,
    load_tower_viewpoints,
    LocationResolveError,
    mountain_viewpoint_to_dict,
    prefixed_viewpoint_name,
    resolve_mountain_viewpoint,
    resolve_launch_location,
    resolve_tower_viewpoint,
    split_prefixed_viewpoint,
    tower_viewpoint_to_dict,
)
from ..launch_time import (
    LaunchSetupError,
    parse_launch_time_arguments,
)
from ..logging_utils import setup_root_logger
from ..paths import (
    APP_DISPLAY_NAME,
    DSO_CSV_FILE,
    EPHEMERIS_FILENAME,
    THEME_STYLES_BY_PRESET,
    STARS_CSV_FILE,
)
from ..gui.window_inputs import (
    prepare_window_catalogs,
    prepare_window_runtime_options,
    prepare_window_user_options,
    prepare_window_viewer_data,
)
from ..overlay_time import target_time_utc_from_delta
from ..search.jpl import search_jpl_targets
from ..search.satellites import search_satellite_targets
from ..search.resolver import resolve_search_targets
from ..cli.args import parse_args

logger = logging.getLogger(__name__)


def _load_star_catalog_for_launch(vmag_limit: float | None):
    logger.info("Loading city and star data...")
    try:
        star_catalog = load_star_catalog(STARS_CSV_FILE, vmag_threshold=vmag_limit)
    except FileNotFoundError as exc:
        logger.error("Fail to load file: %s", STARS_CSV_FILE)
        raise LaunchSetupError() from exc

    limit_str = vmag_limit if vmag_limit is not None else "no limit"
    logger.info("Loaded %d stars (Vmag <= %s)", len(star_catalog), limit_str)
    return star_catalog


def _load_dso_catalog_for_launch():
    try:
        dso_catalog = load_dso_catalog(DSO_CSV_FILE)
    except FileNotFoundError:
        logger.info("DSO catalog not found: %s (shape overlay disabled)", DSO_CSV_FILE)
        return None
    except Exception:
        logger.warning("Failed to load DSO catalog: %s", DSO_CSV_FILE, exc_info=True)
        return None
    logger.info("Loaded %d DSO rows", len(dso_catalog))
    return dso_catalog


def _verify_ephemeris_for_launch() -> None:
    logger.info("Checking ephemeris cache...")
    try:
        load_ephemeris()
    except OSError as exc:
        logger.error(
            "Failed to load ephemeris %s. The app cannot continue until this file "
            "is available in the cache. %s",
            EPHEMERIS_FILENAME,
            exc,
        )
        raise LaunchSetupError() from exc
    except Exception as exc:
        logger.error("Unexpected ephemeris load failure for %s: %s", EPHEMERIS_FILENAME, exc)
        raise LaunchSetupError() from exc
    logger.info("Ephemeris ready: %s", EPHEMERIS_FILENAME)


def _handle_dataset_query_cli(args: object) -> int | None:
    list_kind = getattr(args, "list_viewpoints", None)
    if list_kind:
        if list_kind == "t":
            for name in list_tower_primary_names():
                print(prefixed_viewpoint_name("tower", name))
        else:
            for name in list_mountain_primary_names():
                print(prefixed_viewpoint_name("mountain", name))
        return 0
    list_names_kind = getattr(args, "list_viewpoint_names", None)
    if list_names_kind:
        if list_names_kind == "t":
            for name in list_tower_all_names():
                print(prefixed_viewpoint_name("tower", name))
        else:
            for name in list_mountain_all_names():
                print(prefixed_viewpoint_name("mountain", name))
        return 0
    viewpoint_name = getattr(args, "show_viewpoint_json", None)
    if viewpoint_name:
        explicit = split_prefixed_viewpoint(viewpoint_name)
        if explicit is not None:
            kind, name = explicit
            if kind == "tower":
                tower = resolve_tower_viewpoint(name)
                if tower is None:
                    print(f"No tower found for {viewpoint_name!r}", file=sys.stderr)
                    return 1
                print(
                    json.dumps(
                        tower_viewpoint_to_dict(tower),
                        ensure_ascii=False,
                        indent=2,
                        sort_keys=True,
                    )
                )
                return 0
            mountain = resolve_mountain_viewpoint(name)
            if mountain is None:
                print(f"No mountain found for {viewpoint_name!r}", file=sys.stderr)
                return 1
            print(
                json.dumps(
                    mountain_viewpoint_to_dict(mountain),
                    ensure_ascii=False,
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0

        exact_towers = find_exact_viewpoint_matches(viewpoint_name, load_tower_viewpoints())
        exact_mountains = find_exact_viewpoint_matches(viewpoint_name, load_mountain_viewpoints())
        if exact_towers and exact_mountains:
            candidates = [
                *(prefixed_viewpoint_name("tower", tower.name) for tower in exact_towers),
                *(prefixed_viewpoint_name("mountain", mountain.name) for mountain in exact_mountains),
            ]
            print(
                f"Ambiguous viewpoint name {viewpoint_name!r}. Matches:\n"
                + "\n".join(f"- {candidate}" for candidate in candidates),
                file=sys.stderr,
            )
            return 1

        tower = resolve_tower_viewpoint(viewpoint_name)
        if tower is not None:
            print(json.dumps(tower_viewpoint_to_dict(tower), ensure_ascii=False, indent=2, sort_keys=True))
            return 0
        mountain = resolve_mountain_viewpoint(viewpoint_name)
        if mountain is not None:
            print(json.dumps(mountain_viewpoint_to_dict(mountain), ensure_ascii=False, indent=2, sort_keys=True))
            return 0
        print(f"No viewpoint found for {viewpoint_name!r}", file=sys.stderr)
        return 1
    return None


def main() -> None:
    """Main entry point for the star sky visualizer."""
    args = parse_args()
    cli_exit_code = _handle_dataset_query_cli(args)
    if cli_exit_code is not None:
        raise SystemExit(cli_exit_code)

    from ..splash import setup_app, setup_splash_and_attach_logger
    from ..gui.window import FramelessSkyWindow, StandardSkyWindow

    app_name = APP_DISPLAY_NAME
    app = setup_app(app_name)
    app.setQuitOnLastWindowClosed(False)

    root_logger = setup_root_logger()
    logger.info(f"{APP_DISPLAY_NAME} starting...")

    theme = THEME_STYLES_BY_PRESET.get(args.theme, THEME_STYLES_BY_PRESET["night"])
    splash, splash_handler, set_splash_context = setup_splash_and_attach_logger(
        app,
        app_name,
        root_logger,
        theme,
    )
    if getattr(args, "clear_long_lived_cache", False):
        try:
            logger.info("Clearing long-lived cache on user request...")
            clear_long_lived_cache()
        except LongLivedCacheClearCooldownError as exc:
            logger.error("%s", exc)
            time.sleep(3)
            splash.close()
            root_logger.removeHandler(splash_handler)
            return

    try:
        city = resolve_launch_location(
            args.city,
            place_query=args.place,
            place_countrycode=args.place_countrycode,
            place_lang=args.place_lang,
            use_building_top=bool(getattr(args, "use_building_top", False)),
        )
        if args.timezone is not None:
            city = replace(city, tz=args.timezone)
        set_splash_context(format_splash_location(city))
        delta_t = parse_launch_time_arguments(
            args.datetime,
            args.days,
            args.hours,
            timezone_name=city.tz,
            timezone_override=args.timezone,
        )
        star_catalog = _load_star_catalog_for_launch(args.vmag_limit)
        dso_catalog = _load_dso_catalog_for_launch()
        _verify_ephemeris_for_launch()
    except (LocationResolveError, LaunchSetupError):
        time.sleep(3)
        splash.close()
        root_logger.removeHandler(splash_handler)
        return

    view_center = (args.view_center_alt, args.view_center_az)
    view_center = (min(90.0, max(-5.0, view_center[0])), view_center[1] % 360)
    cloud_stripe_mode, cloud_stripe_count, cloud_stripe_width = args.cloud_stripe
    visual_preset = args.theme
    star_visibility_boost = theme.star_visibility_boost
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
        edge_fov_deg=args.edge_fov_deg,
        content_fov_deg=args.content_fov_deg,
        observer_height_m=city.observer_height_m if args.observer_height_m is None else args.observer_height_m,
        ground_elevation_m=city.ground_elevation_m,
        location_height_label=city.location_height_label,
        location_height_m=city.location_height_m,
        show_observer_height=args.observer_height_m is not None,
    )
    user_options = prepare_window_user_options(
        sky_disc_alpha=args.sky_opacity,
        sky_disc_style=args.sky_disc_style,
        cloud_disc_alpha=0.0 if cloud_stripe_count == 0 or cloud_stripe_width == 0.0 else args.cloud_opacity,
        satellite_opacity=args.satellite_opacity,
        aircraft_opacity=args.aircraft_opacity,
        terrain_horizon_opacity=args.terrain_horizon_opacity,
        earth_guide_opacity=args.earth_guide_opacity,
        urban_outline_opacity=args.urban_outline_opacity,
        ground_tint_opacity=args.ground_tint_opacity,
        enlarge_moon=args.enlarge_moon,
        bright_bodies_mode=args.bright_bodies,
        star_base_radius=args.star_base_radius,
        vmag_limit=args.vmag_limit,
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
        visibility_boost=args.visibility_boost,
        show_dso_initial=args.show_dso_initial,
        show_asterisms_initial=args.show_asterisms_initial,
        show_guidelines_initial=args.show_guidelines_initial,
        observation_info_mode=args.observation_info,
        sky_disc_gui_allowed=args.sky_opacity > 0.0,
        cloud_gui_allowed=args.cloud_opacity > 0.0,
        satellite_gui_allowed=args.satellite_opacity > 0.0,
        aircraft_gui_allowed=args.aircraft_opacity > 0.0,
        terrain_horizon_gui_allowed=args.terrain_horizon_opacity > 0.0,
        earth_guide_gui_allowed=args.earth_guide_opacity > 0.0,
        urban_outline_gui_allowed=args.urban_outline_opacity > 0.0,
    )
    runtime_options = prepare_window_runtime_options(
        delta_t=delta_t,
        sky_update_interval=args.sky_update_interval,
        urban_outline_radius_km=args.urban_outline_radius_km,
        urban_outline_skyscraper_radius_km=args.urban_outline_skyscraper_radius_km,
        urban_outline_min_height_m=args.urban_outline_min_height_m,
        urban_outline_feature_type=args.urban_outline_feature_type,
        urban_outline_skyscraper_only=args.urban_outline_skyscraper_only,
        cloud_stripe_style=(cloud_stripe_count, cloud_stripe_width),
        cloud_stripe_mode=cloud_stripe_mode,
        cloud_missing_tint_opacity=args.cloud_missing_tint_opacity,
        visibility_boost=args.visibility_boost,
        star_render_expected_width=args.expected_render_width,
        content_fov_deg=args.content_fov_deg,
        window_geometry_arg=args.window_geometry,
        window_frame_mode=args.window_frame,
    )
    window_cls = FramelessSkyWindow if args.window_frame == "frameless" else StandardSkyWindow
    main_win = window_cls(
        viewer_data,
        catalogs,
        user_options=user_options,
        runtime_options=runtime_options,
    )
    main_win._search_view_center_alt_specified = bool(
        getattr(args, "view_center_alt_specified", False)
    )
    main_win._search_view_center_az_specified = bool(
        getattr(args, "view_center_az_specified", False)
    )

    startup_search = str(getattr(args, "search", "") or "").strip()
    startup_target_time_utc = target_time_utc_from_delta(delta_t) if startup_search else None

    def _run_startup_search() -> None:
        if not startup_search or startup_target_time_utc is None:
            return
        try:
            resolution = resolve_search_targets(
                startup_search,
                catalogs.named_stars_search_all,
                satellite_search_callback=lambda query: search_satellite_targets(
                    query,
                    target_time_utc=startup_target_time_utc,
                ),
                jpl_search_callback=lambda query: search_jpl_targets(
                    query,
                    target_time_utc=startup_target_time_utc,
                ),
            )
        except Exception as exc:
            logger.error("Startup search failed: %s", exc)
            return
        if len(resolution.candidates) == 1 and resolution.selected_target is not None:
            target = resolution.selected_target
            if bool(getattr(args, "search_keep_marker", False)):
                target = replace(target, persistent_keep_marker=True)
            main_win._jump_to_search_target(target)

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
        if startup_search:
            QTimer.singleShot(0, _run_startup_search)

    main_win.initial_data_loaded.connect(_on_initial_loaded)

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
