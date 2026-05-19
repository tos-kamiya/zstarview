import faulthandler
import json
import logging
import math
import sys
from dataclasses import replace
from datetime import timedelta

from PySide6.QtCore import QObject, QTimer, Signal

from ..astro import load_ephemeris
from ..cache_maintenance import LongLivedCacheClearCooldownError, clear_long_lived_cache
from ..catalog import load_dso_catalog, load_star_catalog
from ..cli.args import parse_args
from ..gui.window_inputs import (
    prepare_window_catalogs,
    prepare_window_runtime_options,
    prepare_window_user_options,
    prepare_window_viewer_data,
)
from ..launch_time import (
    LaunchSetupError,
    parse_launch_time_arguments,
)
from ..location_resolver import (
    LocationResolveError,
    find_exact_viewpoint_matches,
    list_mountain_all_names,
    list_mountain_primary_names,
    list_tower_all_names,
    list_tower_primary_names,
    load_mountain_viewpoints,
    load_tower_viewpoints,
    mountain_viewpoint_to_dict,
    prefixed_viewpoint_name,
    resolve_launch_location,
    resolve_mountain_viewpoint,
    resolve_tower_viewpoint,
    split_prefixed_viewpoint,
    tower_viewpoint_to_dict,
)
from ..logging_utils import setup_root_logger
from ..overlay_time import target_time_utc_from_delta
from ..paths import (
    APP_DISPLAY_NAME,
    DSO_CSV_FILE,
    EPHEMERIS_FILENAME,
    OBSERVER_MAX_ALT_DEG,
    OBSERVER_MIN_ALT_DEG,
    STARS_CSV_FILE,
    THEME_STYLES_BY_PRESET,
)
from ..search.jpl import search_jpl_targets
from ..search.resolver import resolve_search_targets
from ..search.satellites import search_satellite_targets
from ..startup_log import BufferedStartupLogHandler
from ..types import ViewerData
from .worker_pool import submit_gui_work

logger = logging.getLogger(__name__)


def _enable_faulthandler() -> None:
    try:
        faulthandler.enable(all_threads=True)
    except Exception:
        pass


class _StartupBootstrap(QObject):
    """Resolve launch-time state off the UI thread."""

    finished = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        *,
        args: object,
        catalogs: object,
        view_center: tuple[float, float],
    ) -> None:
        super().__init__()
        self._args = args
        self._catalogs = catalogs
        self._view_center = (float(view_center[0]), float(view_center[1]))
        self._started = False

    def start(self) -> None:
        if self._started:
            return
        self._started = True
        submit_gui_work(self._run)

    def _run(self) -> None:
        try:
            startup_search = str(getattr(self._args, "search", "") or "").strip()
            city = resolve_launch_location(
                getattr(self._args, "city", None),
                place_query=getattr(self._args, "place", None),
                place_countrycode=getattr(self._args, "place_countrycode", None),
                place_lang=getattr(self._args, "place_lang", "en"),
                use_building_top=bool(getattr(self._args, "use_building_top", False)),
            )
            timezone_override = getattr(self._args, "timezone", None)
            if timezone_override is not None:
                city = replace(city, tz=timezone_override)
            delta_t = parse_launch_time_arguments(
                getattr(self._args, "datetime", None),
                int(getattr(self._args, "days", 0)),
                int(getattr(self._args, "hours", 0)),
                timezone_name=city.tz,
                timezone_override=timezone_override,
            )
            startup_search_target = None
            if startup_search:
                startup_target_time_utc = target_time_utc_from_delta(delta_t)
                try:
                    resolution = resolve_search_targets(
                        startup_search,
                        self._catalogs.named_stars_search_all,
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
                else:
                    if len(resolution.candidates) == 1 and resolution.selected_target is not None:
                        startup_search_target = resolution.selected_target
                        if bool(getattr(self._args, "search_keep_marker", False)):
                            startup_search_target = replace(
                                startup_search_target,
                                persistent_keep_marker=True,
                            )
            viewer_data = prepare_window_viewer_data(
                city.display_name,
                (city.lat, city.lon, city.tz),
                self._view_center,
                edge_fov_deg=float(getattr(self._args, "edge_fov_deg", 95.0)),
                content_fov_deg=float(getattr(self._args, "content_fov_deg", 110.0)),
                ground_elevation_m=city.ground_elevation_m,
                location_height_label=city.location_height_label,
                location_height_m=city.location_height_m,
                height_add_m=float(
                    getattr(self._args, "observer_height_m", None)
                    if getattr(self._args, "observer_height_m", None) is not None
                    else getattr(city, "height_add_m", 1.7)
                ),
            )
            self.finished.emit(
                {
                    "viewer_data": viewer_data,
                    "delta_t": delta_t,
                    "startup_search_target": startup_search_target,
                }
            )
        except (LocationResolveError, LaunchSetupError) as exc:
            exc_text = str(exc).strip()
            if exc_text:
                logger.error("Startup failed: %s: %s", exc.__class__.__name__, exc_text)
            else:
                logger.error("Startup failed: %s", exc.__class__.__name__)
            self.failed.emit(str(exc))
        except Exception as exc:
            logger.error("Startup failed: %s", exc, exc_info=True)
            self.failed.emit(str(exc))


class _StartupRevealGate:
    """Track when startup overlay can be hidden."""

    def __init__(self, *, requires_terrain: bool) -> None:
        self._initial_loaded = False
        self._terrain_resolved = not bool(requires_terrain)

    def mark_initial_loaded(self) -> bool:
        self._initial_loaded = True
        return self.is_ready()

    def mark_terrain_resolved(self) -> bool:
        self._terrain_resolved = True
        return self.is_ready()

    def is_ready(self) -> bool:
        return self._initial_loaded and self._terrain_resolved


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
    _enable_faulthandler()
    args = parse_args()
    cli_exit_code = _handle_dataset_query_cli(args)
    if cli_exit_code is not None:
        raise SystemExit(cli_exit_code)

    from ..gui.window import SkyWindow, StandardSkyWindow
    from ..splash import setup_app

    app_name = APP_DISPLAY_NAME
    app = setup_app(app_name)

    root_logger = setup_root_logger()
    startup_log_handler = BufferedStartupLogHandler()
    root_logger.addHandler(startup_log_handler)
    logger.info(f"{APP_DISPLAY_NAME} starting...")

    try:
        _verify_ephemeris_for_launch()
    except LaunchSetupError:
        raise SystemExit(1)

    theme = THEME_STYLES_BY_PRESET.get(args.theme, THEME_STYLES_BY_PRESET["night"])
    star_catalog = _load_star_catalog_for_launch(args.vmag_limit)
    dso_catalog = _load_dso_catalog_for_launch()
    cloud_stripe_mode, cloud_stripe_count, cloud_stripe_width = args.cloud_stripe
    visual_preset = args.theme
    star_visibility_boost = theme.star_visibility_boost
    vmag_brightness_scale = -math.log10(args.vmag_brightness_multiplier)
    catalogs = prepare_window_catalogs(
        star_catalog,
        dso_catalog=dso_catalog,
        vmag_brightness_scale=vmag_brightness_scale,
    )
    view_center = (args.view_center_alt, args.view_center_az)
    view_center = (
        min(OBSERVER_MAX_ALT_DEG, max(OBSERVER_MIN_ALT_DEG, view_center[0])),
        view_center[1] % 360,
    )
    user_options = prepare_window_user_options(
        sky_disc_alpha=args.sky_opacity,
        sky_disc_style=args.sky_disc_style,
        sky_disc_altaz_rings=args.sky_disc_altaz_rings,
        sky_disc_altaz_rings_hover=args.sky_disc_altaz_rings_hover,
        night_light_opacity=args.night_light_opacity,
        cloud_disc_alpha=0.0 if cloud_stripe_count == 0 or cloud_stripe_width == 0.0 else args.cloud_opacity,
        satellite_opacity=args.satellite_opacity,
        aircraft_opacity=args.aircraft_opacity,
        terrain_horizon_opacity=args.terrain_horizon_opacity,
        earth_guide_opacity=args.earth_guide_opacity,
        urban_outline_opacity=args.urban_outline_opacity,
        water_overlay_opacity=args.water_surface_opacity,
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
        night_light_gui_allowed=args.night_light_opacity > 0.0,
        urban_outline_gui_allowed=args.urban_outline_opacity > 0.0,
    )
    runtime_options = prepare_window_runtime_options(
        delta_t=timedelta(0),
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
    window_cls = SkyWindow if args.window_frame == "frameless" else StandardSkyWindow
    main_win = window_cls(
        ViewerData(
            location=(0.0, 0.0),
            timezone_name="UTC",
            city_name="Loading...",
            view_center=view_center,
            edge_fov_deg=args.edge_fov_deg,
            content_fov_deg=args.content_fov_deg,
            observer_height_m=1.7 if args.observer_height_m is None else args.observer_height_m,
            height_add_m=1.7 if args.observer_height_m is None else args.observer_height_m,
            ground_elevation_m=0.0,
            location_height_label=None,
            location_height_m=0.0,
        ),
        catalogs,
        user_options=user_options,
        runtime_options=runtime_options,
        defer_initial_load=True,
    )
    main_win._search_view_center_alt_specified = bool(
        getattr(args, "view_center_alt_specified", False)
    )
    main_win._search_view_center_az_specified = bool(
        getattr(args, "view_center_az_specified", False)
    )
    startup_log_overlay = main_win._ensure_startup_log_overlay()
    startup_log_handler.set_consumer(startup_log_overlay.append_line)
    startup_log_overlay.show()
    startup_log_overlay.raise_()
    main_win._raise_overlay_widgets()
    main_win.show()
    app.setQuitOnLastWindowClosed(True)
    app.processEvents()

    close_on_startup_error = bool(getattr(args, "close_on_startup_error", False))
    startup_handler_attached = True
    pending_startup_search_target = None

    def _detach_startup_logging(*, hide_overlay: bool) -> None:
        nonlocal startup_handler_attached
        if not startup_handler_attached:
            return
        startup_log_handler.set_consumer(None)
        if hide_overlay:
            main_win._hide_startup_log_overlay()
        root_logger.removeHandler(startup_log_handler)
        startup_handler_attached = False

    startup_error = False
    if getattr(args, "clear_long_lived_cache", False):
        try:
            logger.info("Clearing long-lived cache on user request...")
            clear_long_lived_cache()
        except LongLivedCacheClearCooldownError as exc:
            logger.error("Startup failed: %s", exc)
            startup_error = True

    if startup_error:
        main_win._raise_overlay_widgets()
        if close_on_startup_error:
            QTimer.singleShot(0, lambda: app.exit(1))
    else:
        startup_bootstrap = _StartupBootstrap(
            args=args,
            catalogs=catalogs,
            view_center=view_center,
        )

        def _on_initial_loaded() -> None:
            nonlocal pending_startup_search_target
            if not startup_reveal_gate.mark_initial_loaded():
                return
            _detach_startup_logging(hide_overlay=True)
            if pending_startup_search_target is not None:
                target = pending_startup_search_target
                pending_startup_search_target = None
                QTimer.singleShot(
                    0,
                    lambda: main_win._jump_to_search_target(target),
                )

        def _on_startup_terrain_resolved(_payload: object) -> None:
            nonlocal pending_startup_search_target
            if not startup_reveal_gate.mark_terrain_resolved():
                return
            _detach_startup_logging(hide_overlay=True)
            if pending_startup_search_target is not None:
                target = pending_startup_search_target
                pending_startup_search_target = None
                QTimer.singleShot(
                    0,
                    lambda: main_win._jump_to_search_target(target),
                )

        def _on_startup_ready(payload: object) -> None:
            nonlocal pending_startup_search_target
            data = dict(payload)
            viewer_data = data["viewer_data"]
            delta_t = data["delta_t"]
            pending_startup_search_target = data.get("startup_search_target")
            main_win.apply_startup_delta_t(delta_t)
            main_win.apply_startup_viewer_data(viewer_data)
            main_win.start_initial_data_load()

        def _on_startup_failed(_message: str) -> None:
            if close_on_startup_error:
                QTimer.singleShot(0, lambda: app.exit(1))

        startup_reveal_gate = _StartupRevealGate(
            requires_terrain=main_win.terrain_horizon_opacity > 0.0
        )
        main_win.initial_data_loaded.connect(_on_initial_loaded)
        if main_win.terrain_horizon_opacity > 0.0:
            main_win._terrain_horizon_controller.terrain_ready.connect(
                _on_startup_terrain_resolved
            )
            main_win._terrain_horizon_controller.terrain_failed.connect(
                _on_startup_terrain_resolved
            )
        startup_bootstrap.finished.connect(_on_startup_ready)
        startup_bootstrap.failed.connect(_on_startup_failed)
        QTimer.singleShot(0, startup_bootstrap.start)

    exit_code = app.exec()
    _detach_startup_logging(hide_overlay=False)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
