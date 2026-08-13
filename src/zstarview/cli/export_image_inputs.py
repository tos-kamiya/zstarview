"""Build export-image catalogs, viewer data, and launch options."""

from __future__ import annotations

import math
from dataclasses import replace
from datetime import datetime, timezone

from ..astro import load_ephemeris
from ..catalog import load_dso_catalog, load_star_catalog
from ..gui.window_inputs import (
    PreparedWindowCatalogs,
    SkyWindowRuntimeOptions,
    SkyWindowUserOptions,
    prepare_window_user_options,
)
from ..launch_time import LaunchSetupError
from ..location_resolver import (
    LocationResolveError,
    ResolvedLocation,
)
from ..overlay_time import overlay_availability_for_delta
from ..paths import (
    DSO_CSV_FILE,
    EPHEMERIS_FILENAME,
    OBSERVER_MAX_ALT_DEG,
    OBSERVER_MIN_ALT_DEG,
    STARS_CSV_FILE,
    THEME_STYLES_BY_PRESET,
)
from ..render.molecular_cloud_overlay import is_molecular_cloud_cache_available
from ..search.jpl import resolve_jpl_target_state_vector, search_jpl_targets
from ..search.models import SearchJumpTarget
from ..search.resolver import compute_search_target_altaz, resolve_search_targets
from ..search.satellites import resolve_satellite_target_altaz, search_satellite_targets
from ..types import ViewerData
from .export_image_output import (
    _clamp_view_center_altitude,
    _format_search_candidate_line,
    _format_search_failure_message,
    _search_view_center_for_target,
)
from .export_image_support import (
    DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL,
    OPEN_METEO_CONSENT_REQUIRED_MESSAGE,
    host,
    logger,
)


def _require_open_meteo_consent_for_export(precipitation_opacity: float) -> None:
    if (
        float(precipitation_opacity) > 0.0
        and not host().open_meteo_noncommercial_terms_accepted()
    ):
        print(OPEN_METEO_CONSENT_REQUIRED_MESSAGE, file=host().sys.stderr)
        raise SystemExit(1)

def _load_star_catalog_for_export(vmag_limit: float | None):
    logger.info("Loading city and star data...")
    try:
        star_catalog = load_star_catalog(STARS_CSV_FILE, vmag_threshold=vmag_limit)
    except FileNotFoundError as exc:
        logger.error("Fail to load file: %s", STARS_CSV_FILE)
        raise LaunchSetupError() from exc

    limit_str = vmag_limit if vmag_limit is not None else "no limit"
    logger.info("Loaded %d stars (Vmag <= %s)", len(star_catalog), limit_str)
    return star_catalog

def _load_dso_catalog_for_export():
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

def _verify_ephemeris_for_export() -> None:
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
        logger.error(
            "Unexpected ephemeris load failure for %s: %s", EPHEMERIS_FILENAME, exc
        )
        raise LaunchSetupError() from exc
    logger.info("Ephemeris ready: %s", EPHEMERIS_FILENAME)

def _build_window_inputs_from_args(
    args: object,
) -> tuple[
    PreparedWindowCatalogs,
    ViewerData,
    SkyWindowUserOptions,
    SkyWindowRuntimeOptions,
    SearchJumpTarget | None,
    ResolvedLocation | None,
]:
    try:
        city = host().resolve_launch_location(
            args.city,
            place_query=args.place,
            place_countrycode=args.place_countrycode,
            place_lang=args.place_lang,
            use_building_top=bool(getattr(args, "use_building_top", False)),
        )
    except LocationResolveError as exc:
        raise LaunchSetupError() from exc
    timezone_arg = args.timezone
    if timezone_arg is not None:
        city = replace(city, tz=timezone_arg)
    delta_t = host().parse_launch_time_arguments(
        args.datetime,
        args.days,
        args.hours,
        timezone_name=city.tz,
        timezone_override=timezone_arg,
    )
    overlay_availability = overlay_availability_for_delta(delta_t)
    star_catalog = host()._load_star_catalog_for_export(args.vmag_limit)
    view_center = (
        args.view_center_alt,
        args.view_center_az,
    )
    view_center = (
        min(OBSERVER_MAX_ALT_DEG, max(OBSERVER_MIN_ALT_DEG, view_center[0])),
        view_center[1] % 360.0,
    )
    cloud_stripe_mode, cloud_stripe_count, cloud_stripe_width = args.cloud_stripe
    visual_preset = args.theme
    star_visibility_boost = THEME_STYLES_BY_PRESET.get(
        visual_preset, THEME_STYLES_BY_PRESET["night"]
    ).star_visibility_boost
    vmag_brightness_scale = -math.log10(args.vmag_brightness_multiplier)

    catalogs = host().prepare_window_catalogs(
        star_catalog,
        dso_catalog=None,
        vmag_brightness_scale=vmag_brightness_scale,
    )
    viewer_data = host().prepare_window_viewer_data(
        city.display_name,
        (city.lat, city.lon, city.tz),
        view_center,
        edge_fov_deg=args.edge_fov_deg,
        content_fov_deg=args.content_fov_deg,
        ground_elevation_m=city.ground_elevation_m,
        location_height_label=city.location_height_label,
        location_height_m=city.location_height_m,
        height_add_m=(
            getattr(city, "height_add_m", 1.7)
            if args.observer_height_m is None
            else args.observer_height_m
        ),
    )

    search_query = str(args.search or "").strip()
    search_overlay_target: SearchJumpTarget | None = None
    if search_query:
        fixed_alt = bool(args.view_center_alt_specified)
        fixed_az = bool(args.view_center_az_specified)
        target_time_utc = datetime.now(timezone.utc) + delta_t
        try:
            resolution = resolve_search_targets(
                search_query,
                catalogs.named_stars_search_all,
                satellite_search_callback=lambda query: search_satellite_targets(
                    query,
                    target_time_utc=target_time_utc,
                ),
                jpl_search_callback=lambda query: search_jpl_targets(
                    query,
                    target_time_utc=target_time_utc,
                ),
            )
        except Exception as exc:
            host().sys.stderr.write(f"{exc}\n")
            raise SystemExit(1) from exc
        candidates = resolution.candidates
        lines = [_format_search_candidate_line(target) for target in candidates]
        if args.list:
            if lines:
                host().sys.stdout.write("\n".join(lines) + "\n")
                host().sys.stdout.flush()
            return
        if len(candidates) != 1 or resolution.selected_target is None:
            host().sys.stderr.write(
                _format_search_failure_message(search_query, len(candidates)) + "\n"
            )
            if lines:
                host().sys.stderr.write("\n".join(lines) + "\n")
                host().sys.stderr.flush()
            raise SystemExit(1)
        target = resolution.selected_target
        if target.kind in {"jpl_small_body", "jpl_body"} and (
            target.horizons_epoch_utc is None
            or target.horizons_position_km is None
            or target.horizons_velocity_km_s is None
        ):
            state_vector = resolve_jpl_target_state_vector(
                target,
                target_time_utc=target_time_utc,
            )
            if state_vector is not None:
                horizons_epoch_utc, horizons_position_km, horizons_velocity_km_s = (
                    state_vector
                )
                target = replace(
                    target,
                    horizons_epoch_utc=horizons_epoch_utc,
                    horizons_position_km=horizons_position_km,
                    horizons_velocity_km_s=horizons_velocity_km_s,
                )
        altaz = compute_search_target_altaz(
            target,
            observer_lat=float(viewer_data.lat_deg),
            observer_lon=float(viewer_data.lon_deg),
            observer_height_m=float(viewer_data.observer_height_m),
            target_time_utc=target_time_utc,
            satellite_altaz_resolver=lambda satellite_target: (
                resolve_satellite_target_altaz(
                    satellite_target,
                    observer_lat=float(viewer_data.lat_deg),
                    observer_lon=float(viewer_data.lon_deg),
                    observer_height_m=float(viewer_data.observer_height_m),
                    target_time_utc=target_time_utc,
                )
            ),
        )
        if altaz is not None:
            target_alt = _clamp_view_center_altitude(float(altaz[0]))
            target_az = float(altaz[1]) % 360.0
            viewer_data = replace(
                viewer_data,
                view_center=_search_view_center_for_target(
                    base_view_center=view_center,
                    target_alt_deg=target_alt,
                    target_az_deg=target_az,
                    fixed_alt=fixed_alt,
                    fixed_az=fixed_az,
                ),
            )
            search_overlay_target = replace(
                target,
                alt_deg=target_alt,
                az_deg=target_az,
            )

    dso_catalog = host()._load_dso_catalog_for_export()
    host()._verify_ephemeris_for_export()

    catalogs = host().prepare_window_catalogs(
        star_catalog,
        dso_catalog=dso_catalog,
        vmag_brightness_scale=vmag_brightness_scale,
    )
    user_options = prepare_window_user_options(
        sky_disc_alpha=args.sky_opacity,
        sky_disc_style=args.sky_disc_style,
        sky_disc_altaz_rings=args.sky_disc_altaz_rings,
        sky_disc_altaz_rings_hover=args.sky_disc_altaz_rings_hover,
        night_light_opacity=args.night_light_opacity,
        road_light_opacity=float(getattr(args, "road_light_opacity", 0.0)),
        akari_ir_bands_opacity=float(
            getattr(args, "akari_ir_bands_opacity", 0.10)
        ),
        ridge_glow_opacity=args.ridge_glow_opacity,
        precipitation_opacity=(
            float(getattr(args, "precipitation_opacity", 0.0))
            if overlay_availability.precipitation
            else 0.0
        ),
        cloud_disc_alpha=(
            0.0
            if (not overlay_availability.cloud)
            or cloud_stripe_count == 0
            or cloud_stripe_width == 0.0
            else args.cloud_opacity
        ),
        geo_satellite=bool(args.geo_satellite),
        satellite_opacity=(
            args.satellite_opacity if overlay_availability.satellite else 0.0
        ),
        aircraft_opacity=(
            args.aircraft_opacity if overlay_availability.aircraft else 0.0
        ),
        tropical_cyclone_opacity=(
            args.tropical_cyclone_opacity
            if overlay_availability.tropical_cyclone
            else 0.0
        ),
        terrain_horizon_opacity=args.terrain_horizon_opacity,
        earth_guide_opacity=args.earth_guide_opacity,
        urban_outline_opacity=args.urban_outline_opacity,
        water_overlay_opacity=args.water_surface_opacity,
        ground_tint_opacity=args.ground_tint_opacity,
        overlay_font_size=args.overlay_font_size,
        enlarge_moon=bool(args.enlarge_moon),
        bright_bodies_mode=str(args.bright_bodies),
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
        cloud_gui_allowed=overlay_availability.cloud and args.cloud_opacity > 0.0,
        satellite_gui_allowed=overlay_availability.satellite
        and args.satellite_opacity > 0.0,
        aircraft_gui_allowed=overlay_availability.aircraft
        and args.aircraft_opacity > 0.0,
        tropical_cyclone_gui_allowed=(
            overlay_availability.tropical_cyclone
            and args.tropical_cyclone_opacity > 0.0
        ),
        terrain_horizon_gui_allowed=args.terrain_horizon_opacity > 0.0,
        earth_guide_gui_allowed=args.earth_guide_opacity > 0.0,
        night_light_gui_allowed=args.night_light_opacity > 0.0,
        road_light_gui_allowed=float(getattr(args, "road_light_opacity", 0.0)) > 0.0,
        akari_ir_bands_gui_allowed=(
            is_molecular_cloud_cache_available()
            and float(getattr(args, "akari_ir_bands_opacity", 0.10)) > 0.0
        ),
        urban_outline_gui_allowed=args.urban_outline_opacity > 0.0,
    )
    runtime_options = host().prepare_window_runtime_options(
        delta_t=delta_t,
        sky_update_interval=DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL,
        urban_outline_radius_km=args.urban_outline_radius_km,
        urban_outline_skyscraper_radius_km=args.urban_outline_skyscraper_radius_km,
        urban_outline_min_height_m=args.urban_outline_min_height_m,
        urban_outline_max_candidates=args.urban_outline_max_candidates,
        road_light_max_candidates=getattr(args, "road_light_max_candidates", 5000),
        urban_outline_feature_type=args.urban_outline_feature_type,
        urban_outline_skyscraper_only=bool(args.urban_outline_skyscraper_only),
        urban_outline_download_timeout_seconds=args.urban_outline_download_timeout_seconds,
        cloud_stripe_style=(cloud_stripe_count, cloud_stripe_width),
        cloud_stripe_mode=cloud_stripe_mode,
        cloud_missing_tint_opacity=args.cloud_missing_tint_opacity,
        visibility_boost=args.visibility_boost,
        star_render_expected_width=args.expected_render_width,
        content_fov_deg=args.content_fov_deg,
        window_geometry_arg=None,
        window_frame_mode="frameless",
    )
    place_location = city if args.place is not None else None
    return (
        catalogs,
        viewer_data,
        user_options,
        runtime_options,
        search_overlay_target,
        place_location,
    )
