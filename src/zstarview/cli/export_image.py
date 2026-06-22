from __future__ import annotations

import logging
import math
import os
import select
import threading
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, TypedDict

import numpy as np
from PySide6.QtCore import QBuffer, QByteArray, QIODevice, QPoint, QRect
from PySide6.QtGui import QFont, QFontDatabase, QImage, QPainter

try:
    import termios
except ImportError:  # pragma: no cover - non-Unix fallback
    termios = None  # type: ignore[assignment]

from ..aircraft import (
    build_observer_bbox,
    fetch_cached_opensky_states,
)
from ..astro import load_ephemeris
from ..cache_maintenance import LongLivedCacheClearCooldownError, clear_long_lived_cache
from ..catalog import load_dso_catalog, load_star_catalog
from ..clouddisc import CloudDisc, CloudDiscConfig, VisibilityError
from ..clouddisc.altaz_grid import CloudAltAzGrid
from ..data.import_overture_buildings import (
    derive_dataset_name,
    import_overture_buildings,
    import_overture_buildings_for_bbox,
    is_derived_dataset_stale,
    resolve_overture_release_for_cache_root,
)
from ..data.skyscraper_tiles import (
    SKYSCRAPER_TILES_FILE,
    select_skyscraper_seed_tiles_for_viewer,
    skyscraper_tile_derived_dir,
)
from ..gui.composite import SkyCompositorCache, build_cloud_amount_field_from_rgba
from ..gui.sky_worker import compute_sky_snapshot
from ..gui.window_inputs import (
    PreparedWindowCatalogs,
    SkyWindowRuntimeOptions,
    SkyWindowUserOptions,
    prepare_window_catalogs,
    prepare_window_runtime_options,
    prepare_window_user_options,
    prepare_window_viewer_data,
)
from ..launch_time import (
    LaunchSetupError,
    parse_launch_time_arguments,
)
from ..location_resolver import LocationResolveError, resolve_launch_location
from ..logging_utils import setup_root_logger
from ..night_lights import compute_night_light_glow_profile
from ..night_lights import is_night_light_enabled
from ..overlay_time import classify_target_time, overlay_availability_for_delta
from ..paths import (
    APP_DISPLAY_NAME,
    CACHE_PATH,
    CLOUD_MISSING_TINT_RGBA,
    CLOUD_SHELLS_KM,
    DSO_CSV_FILE,
    EPHEMERIS_FILENAME,
    OBSERVER_MAX_ALT_DEG,
    OBSERVER_MIN_ALT_DEG,
    OVERTURE_DERIVED_ROOT_DIR,
    OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR,
    STARS_CSV_FILE,
    STATUS_LINE_FONT_SIZE,
    TEXT_FONT_PATH,
    OVERLAY_FONT_SIZE_DEFAULT,
    THEME_STYLES_BY_PRESET,
    ThemeStyle,
)
from ..geosatellite.pipeline import is_within_europe_band, run_geo_satellite_pipeline
from ..geosatellite.projection import render_gray_image_to_cloud_rgba
from ..render import background as render_background
from ..render import geometry as render_geometry
from ..render import guides as render_guides
from ..render import text as render_text
from ..render.pipeline import (
    FrameContext,
    RenderHudState,
    RenderSceneData,
    RenderStyle,
    compute_star_render_upscale_factor,
    render_base_scene_into_painter,
)
from ..render.search_overlay import draw_search_target_overlay
from ..satellite_constants import SATELLITE_HORIZONS_CACHE_KEY, SATELLITE_ISS_CACHE_KEY
from ..satellites import resolve_satellite_elements_for_time
from ..search.jpl import resolve_jpl_target_state_vector, search_jpl_targets
from ..search.models import SearchJumpTarget
from ..search.resolver import compute_search_target_altaz, resolve_search_targets
from ..search.satellites import resolve_satellite_target_altaz, search_satellite_targets
from ..splash import setup_app
from ..terrain import (
    EARTH_MEAN_RADIUS_M,
    DEFAULT_TERRAIN_DISTANCE_SAMPLE_STEP_M,
    WGS84_GEOD,
    GeoTiffDem,
    ObserverLocation,
    build_distance_samples,
    build_download_bbox,
    compute_flat_ground_horizon_layers,
    compute_horizon_layers,
    fetch_copernicus_dem,
    reduce_profile_to_altaz,
    sample_ground_elevation,
)
from ..types import CelestialData, UrbanOutlinePolyline, ViewerData
from ..urban_outline_layer import resolve_urban_outline_layer_for_viewer
from ..water_mask_interface import (
    WaterSurfaceBandStats,
    sample_water_surface_interface_points_with_stats,
)
from ..water_overlay import (
    DEFAULT_WATER_RADIUS_KM,
    DEFAULT_WATER_OVERPASS_ENDPOINT,
    DEFAULT_WATER_USER_AGENT,
    WaterOverlayPoint,
    WaterPolygonFootprint,
    fetch_water_overlay_footprints,
    resolve_water_scan_radius_km,
    resolve_water_surface_azimuth_step_deg,
    sample_water_overlay_points,
    simplify_water_footprints_for_observer,
)
from ..gui.water_overlay_cache import (
    WATER_OVERLAY_CACHE_RETENTION_SECONDS,
    WaterOverlayCacheSnapshot,
    load_water_overlay_cache,
    save_water_overlay_cache,
    water_overlay_cache_is_recent,
    water_overlay_cache_scope_key,
)
from .args import parse_export_image_args

logger = logging.getLogger(__name__)

sample_water_overlay_points_for_observer = sample_water_overlay_points


class TerrainHorizonPayload(TypedDict):
    profile_altaz: list[tuple[float, float]]
    profile_distances_m: list[float]
    secondary_ridges_altaz_layers: list[list[tuple[float, float]]]
    secondary_ridges_distances_m_layers: list[list[float]]
    sample_distances_m: np.ndarray
    sample_terrain_elevation_m: np.ndarray


DEFAULT_CLOUD_ALT_MIN_DEG = 1.0
DEFAULT_CLOUD_FOV_OVERSCAN_DEG = 2.0
DEFAULT_CLOUD_BASE_SIZE = 256
DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL = 60
DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG = 0.5
DEFAULT_EXPORT_IMAGE_TERRAIN_SAMPLE_STEP_M = DEFAULT_TERRAIN_DISTANCE_SAMPLE_STEP_M


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


def _deadline_after(timeout_seconds: float) -> float | None:
    if float(timeout_seconds) <= 0.0:
        return None
    return time.monotonic() + float(timeout_seconds)


def _remaining_timeout_seconds(deadline: float | None) -> float | None:
    if deadline is None:
        return None
    return max(0.0, deadline - time.monotonic())


def _timed_out(deadline: float | None) -> bool:
    remaining = _remaining_timeout_seconds(deadline)
    return remaining is not None and remaining <= 0.0


def _water_overlay_band_stats_text(stats: WaterSurfaceBandStats) -> str:
    return (
        f"{stats.band_name} tiles={int(stats.loaded_tile_count)} "
        f"raw={int(stats.raw_point_count)} "
        f"collapsed={int(stats.collapsed_point_count)} "
        f"visible={int(stats.visible_point_count)}"
    )


def _water_overlay_band_counts(
    points: list[WaterOverlayPoint] | tuple[WaterOverlayPoint, ...],
) -> tuple[int, int, int]:
    counts = Counter(str(point.water_category).strip().lower() for point in points)
    return (
        int(counts.get("sea-125", 0)),
        int(counts.get("sea-250", 0)),
        int(counts.get("sea-500", 0)),
    )


def _build_window_inputs_from_args(
    args: object,
) -> tuple[
    PreparedWindowCatalogs,
    ViewerData,
    SkyWindowUserOptions,
    SkyWindowRuntimeOptions,
    SearchJumpTarget | None,
]:
    try:
        city = resolve_launch_location(
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
    delta_t = parse_launch_time_arguments(
        args.datetime,
        args.days,
        args.hours,
        timezone_name=city.tz,
        timezone_override=timezone_arg,
    )
    overlay_availability = overlay_availability_for_delta(delta_t)
    star_catalog = _load_star_catalog_for_export(args.vmag_limit)
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

    catalogs = prepare_window_catalogs(
        star_catalog,
        dso_catalog=None,
        vmag_brightness_scale=vmag_brightness_scale,
    )
    viewer_data = prepare_window_viewer_data(
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
            sys.stderr.write(f"{exc}\n")
            raise SystemExit(1) from exc
        candidates = resolution.candidates
        lines = [_format_search_candidate_line(target) for target in candidates]
        if args.list:
            if lines:
                sys.stdout.write("\n".join(lines) + "\n")
                sys.stdout.flush()
            return
        if len(candidates) != 1 or resolution.selected_target is None:
            sys.stderr.write(
                _format_search_failure_message(search_query, len(candidates)) + "\n"
            )
            if lines:
                sys.stderr.write("\n".join(lines) + "\n")
                sys.stderr.flush()
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

    dso_catalog = _load_dso_catalog_for_export()
    _verify_ephemeris_for_export()

    catalogs = prepare_window_catalogs(
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
        ridge_glow_opacity=args.ridge_glow_opacity,
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
        cloud_altaz_grid=bool(getattr(args, "cloud_altaz_grid", False)),
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
        urban_outline_gui_allowed=args.urban_outline_opacity > 0.0,
    )
    runtime_options = prepare_window_runtime_options(
        delta_t=delta_t,
        sky_update_interval=DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL,
        urban_outline_radius_km=args.urban_outline_radius_km,
        urban_outline_skyscraper_radius_km=args.urban_outline_skyscraper_radius_km,
        urban_outline_min_height_m=args.urban_outline_min_height_m,
        urban_outline_max_candidates=args.urban_outline_max_candidates,
        urban_outline_feature_type=args.urban_outline_feature_type,
        urban_outline_skyscraper_only=bool(args.urban_outline_skyscraper_only),
        cloud_stripe_style=(cloud_stripe_count, cloud_stripe_width),
        cloud_stripe_mode=cloud_stripe_mode,
        cloud_missing_tint_opacity=args.cloud_missing_tint_opacity,
        visibility_boost=args.visibility_boost,
        star_render_expected_width=args.expected_render_width,
        content_fov_deg=args.content_fov_deg,
        window_geometry_arg=None,
        window_frame_mode="frameless",
    )
    return catalogs, viewer_data, user_options, runtime_options, search_overlay_target


def _load_fonts(
    overlay_font_size: float = float(OVERLAY_FONT_SIZE_DEFAULT),
) -> tuple[QFont, QFont]:
    text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
    text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
    text_font = QFont(text_font_family)
    text_font.setPointSizeF(float(overlay_font_size))
    status_line_font = QFont(text_font_family)
    status_line_font.setPointSizeF(float(STATUS_LINE_FONT_SIZE))
    return (text_font, status_line_font)


def _build_compositor(
    runtime_options: SkyWindowRuntimeOptions, user_options: SkyWindowUserOptions
) -> SkyCompositorCache:
    target_stripes, width_factor = runtime_options.cloud_stripe_style
    missing_tint_alpha = int(round(255.0 * runtime_options.cloud_missing_tint_opacity))
    missing_tint_rgba = (
        int(CLOUD_MISSING_TINT_RGBA[0]),
        int(CLOUD_MISSING_TINT_RGBA[1]),
        int(CLOUD_MISSING_TINT_RGBA[2]),
        missing_tint_alpha,
    )
    return SkyCompositorCache(
        cloud_target_stripes=int(target_stripes),
        cloud_stripe_width_factor=float(width_factor),
        cloud_stripe_mode=runtime_options.cloud_stripe_mode,
        missing_tint_rgba=missing_tint_rgba,
        ground_tint_opacity=user_options.ground_tint_opacity,
    )


def _abort_export_without_partial_data() -> None:
    logger.error(
        "Export aborted because partial data is not allowed. Re-run with "
        "--allow-partial-data to continue and still output an image when "
        "terrain, urban, water, aircraft, or satellite data cannot be downloaded."
    )
    raise SystemExit(1)


def _start_background_task(
    *,
    name: str,
    target: Callable[[], object],
) -> tuple[threading.Thread, threading.Event, dict[str, object]]:
    result: dict[str, object] = {}
    done = threading.Event()

    def _runner() -> None:
        try:
            result["value"] = target()
        except Exception as exc:  # pragma: no cover - exercised through caller
            result["error"] = exc
        finally:
            done.set()

    thread = threading.Thread(target=_runner, name=name, daemon=True)
    thread.start()
    return thread, done, result


def _fetch_cloud_layer(
    *,
    viewer_data: ViewerData,
    user_options: SkyWindowUserOptions,
    deadline: float | None,
) -> tuple[np.ndarray | None, np.ndarray | None, object | None, float | None, CloudAltAzGrid | None]:
    if user_options.cloud_disc_alpha <= 0.0:
        return (None, None, None, None, None)
    if _timed_out(deadline):
        raise TimeoutError("cloud timed out")

    requested_geo_satellite = bool(getattr(user_options, "geo_satellite", False))
    within_geo_satellite_band = is_within_europe_band(
        float(viewer_data.lat_deg),
        float(viewer_data.lon_deg),
    )
    if requested_geo_satellite and within_geo_satellite_band:
        logger.info("Geo-sat + Downloading")
        result = run_geo_satellite_pipeline(
            observer_lat=float(viewer_data.lat_deg),
            observer_lon=float(viewer_data.lon_deg),
            alt=float(viewer_data.view_alt_deg),
            az=float(viewer_data.view_az_deg),
            fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
        )
        logger.info("Calculating initial cloud image...")
        download_result = result.download
        captured_at_utc = getattr(download_result, "captured_at_utc", None) or getattr(
            download_result,
            "fetched_at_utc",
        )
        logger.info(
            "Geo-sat + %s",
            captured_at_utc.astimezone(timezone.utc).isoformat(),
        )
        cloud_rgba = render_gray_image_to_cloud_rgba(result.disc_gray)
        cloud_amount_field = build_cloud_amount_field_from_rgba(cloud_rgba)
        missing_mask = None
        cloud_coverage_ratio = float(
            np.count_nonzero(cloud_rgba[..., 3]) / max(1, cloud_rgba[..., 3].size)
        )
        logger.info("Geo-sat + Projecting")
        return (cloud_rgba, missing_mask, cloud_amount_field, cloud_coverage_ratio, None)

    if within_geo_satellite_band and not requested_geo_satellite:
        logger.warning(
            "Cloud rendering is unavailable for this Europe-band location without --geo-satellite true; skipping the cloud layer."
        )
        return (None, None, None, None, None)

    clouddisc = CloudDisc(
        CloudDiscConfig(
            cache_dir=CACHE_PATH,
            sat_priority=("AUTO",),
            bt_warm_k=310.0,
            bt_cold_k=190.0,
            alt_min_deg=DEFAULT_CLOUD_ALT_MIN_DEG,
            search_back_minutes=120,
            use_altaz_grid=bool(getattr(user_options, "cloud_altaz_grid", False)),
        )
    )
    try:
        source = clouddisc.fetch_source(
            lat=float(viewer_data.lat_deg),
            lon=float(viewer_data.lon_deg),
        )
    except VisibilityError as exc:
        logger.warning("Cloud rendering is unavailable for this location: %s", exc)
        return (None, None, None, None, None)

    if clouddisc.cfg.use_altaz_grid:
        logger.info("Building alt/az cloud grid...")
        source.altaz_grid = clouddisc.build_altaz_grid_from_source(
            source=source,
            lat=float(viewer_data.lat_deg),
            lon=float(viewer_data.lon_deg),
            cloud_shells_km=CLOUD_SHELLS_KM,
        )
        logger.info("Alt/az cloud grid ready.")

    logger.info("Calculating initial cloud image...")

    altaz_grid = getattr(source, "altaz_grid", None)
    if isinstance(altaz_grid, CloudAltAzGrid) and clouddisc.cfg.use_altaz_grid:
        from ..clouddisc.altaz_render import render_altaz_grid_circles, render_altaz_missing_mask
        cloud_rgba = render_altaz_grid_circles(
            altaz_grid,
            width=DEFAULT_CLOUD_BASE_SIZE * 2 + 1,
            height=DEFAULT_CLOUD_BASE_SIZE * 2 + 1,
            center_alt_deg=float(viewer_data.view_alt_deg),
            center_az_deg=float(viewer_data.view_az_deg),
            edge_fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
            mask_fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
        )
        missing_mask = render_altaz_missing_mask(
            altaz_grid,
            width=DEFAULT_CLOUD_BASE_SIZE * 2 + 1,
            height=DEFAULT_CLOUD_BASE_SIZE * 2 + 1,
            center_alt_deg=float(viewer_data.view_alt_deg),
            center_az_deg=float(viewer_data.view_az_deg),
            edge_fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
            mask_fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
        )
        missing_mask_alpha = np.where(missing_mask > 0, 255, 0).astype(np.uint8)
        cloud_amount_field = None
        _coverage_ratio = altaz_grid.coverage_ratio
    else:
        cloud_rgba, _meta, missing_mask, _coverage_ratio = (
            clouddisc.render_from_source_with_coverage(
                source=source,
                lat=float(viewer_data.lat_deg),
                lon=float(viewer_data.lon_deg),
                alt=float(viewer_data.view_alt_deg),
                az=float(viewer_data.view_az_deg),
                radius_px=DEFAULT_CLOUD_BASE_SIZE,
                edge_fov_deg=float(viewer_data.edge_fov_deg)
                + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
                mask_fov_deg=float(viewer_data.edge_fov_deg)
                + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
                cloud_shells_km=CLOUD_SHELLS_KM,
            )
        )
        missing_mask_alpha = np.where(missing_mask > 0, 255, 0).astype(np.uint8)
        cloud_amount_field = build_cloud_amount_field_from_rgba(cloud_rgba)
        altaz_grid = None
    return (cloud_rgba, missing_mask_alpha, cloud_amount_field, float(_coverage_ratio), altaz_grid)


def _start_cloud_layer_fetch(
    *,
    viewer_data: ViewerData,
    user_options: SkyWindowUserOptions,
    deadline: float | None,
) -> tuple[threading.Thread, threading.Event, dict[str, object]]:
    return _start_background_task(
        name="zstarview-export-cloud",
        target=lambda: _fetch_cloud_layer(
            viewer_data=viewer_data,
            user_options=user_options,
            deadline=deadline,
        ),
    )


def _await_background_task_result(
    *,
    label: str,
    thread: threading.Thread,
    done: threading.Event,
    state: dict[str, object],
    deadline: float | None,
    layer_failures: list[str],
    allow_partial_data: bool,
) -> dict[str, object] | None:
    remaining_timeout = _remaining_timeout_seconds(deadline)
    done.wait(timeout=remaining_timeout)
    if not done.is_set():
        logger.warning("Export layer unavailable: %s (timeout)", label)
        layer_failures.append(label)
        if not allow_partial_data:
            _abort_export_without_partial_data()
        return None

    thread.join(timeout=0.0)
    layer_error = state.get("error")
    if isinstance(layer_error, Exception):
        logger.warning("Export layer unavailable: %s (%s)", label, layer_error)
        layer_failures.append(label)
        if not allow_partial_data:
            _abort_export_without_partial_data()
        return None
    return state


def _fetch_terrain_horizon_layer(
    *,
    viewer_data: ViewerData,
    deadline: float | None,
) -> TerrainHorizonPayload:
    if _timed_out(deadline):
        raise TimeoutError("terrain timed out")
    try:
        download = fetch_copernicus_dem(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            max_distance_km=128.0,
            margin_km=10.0,
            cache_dir=Path(CACHE_PATH) / "copernicus-dem",
        )
    except RuntimeError as exc:
        if (
            str(exc)
            != "No Copernicus DEM tiles were downloaded for the requested area."
        ):
            raise
        observer = ObserverLocation(
            latitude_deg=float(viewer_data.lat_deg),
            longitude_deg=float(viewer_data.lon_deg),
            observer_ground_m=0.0,
            observer_eye_m=float(viewer_data.observer_height_m),
        )
        layers = compute_flat_ground_horizon_layers(
            geod=WGS84_GEOD,
            observer=observer,
            azimuth_step_deg=DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG,
            distance_samples_m=build_distance_samples(
                128.0,
                DEFAULT_EXPORT_IMAGE_TERRAIN_SAMPLE_STEP_M,
            ),
            earth_radius_m=EARTH_MEAN_RADIUS_M,
            refraction_coefficient=0.13,
        )
        return {
            "profile_altaz": reduce_profile_to_altaz(layers.main_profile),
            "profile_distances_m": [
                float(point.distance_m) for point in layers.main_profile
            ],
            "secondary_ridges_altaz_layers": [
                reduce_profile_to_altaz(layer) for layer in layers.secondary_layers
            ],
            "secondary_ridges_distances_m_layers": [
                [float(point.distance_m) for point in layer]
                for layer in layers.secondary_layers
            ],
            "sample_distances_m": layers.sample_distances_m,
            "sample_terrain_elevation_m": layers.sample_terrain_elevation_m,
        }
    dem = GeoTiffDem(download.paths, default_elevation_m=0.0)
    try:
        bbox = build_download_bbox(
            lat_deg=float(viewer_data.lat_deg),
            lon_deg=float(viewer_data.lon_deg),
            radius_km=138.0,
        )
        dem_grid = dem.build_grid(bbox)
        ground_m = sample_ground_elevation(
            dem_grid,
            latitude_deg=float(viewer_data.lat_deg),
            longitude_deg=float(viewer_data.lon_deg),
            dem_resampling="bilinear",
        )
        observer = ObserverLocation(
            latitude_deg=float(viewer_data.lat_deg),
            longitude_deg=float(viewer_data.lon_deg),
            observer_ground_m=ground_m,
            observer_eye_m=float(viewer_data.observer_height_m),
        )
        layers = compute_horizon_layers(
            dem_grid=dem_grid,
            geod=WGS84_GEOD,
            observer=observer,
            azimuth_step_deg=DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG,
            distance_samples_m=build_distance_samples(
                128.0,
                DEFAULT_EXPORT_IMAGE_TERRAIN_SAMPLE_STEP_M,
            ),
            dem_resampling="bilinear",
            earth_radius_m=EARTH_MEAN_RADIUS_M,
            refraction_coefficient=0.13,
        )
    finally:
        dem.close()
    return {
        "profile_altaz": reduce_profile_to_altaz(layers.main_profile),
        "profile_distances_m": [
            float(point.distance_m) for point in layers.main_profile
        ],
        "secondary_ridges_altaz_layers": [
            reduce_profile_to_altaz(layer) for layer in layers.secondary_layers
        ],
        "secondary_ridges_distances_m_layers": [
            [float(point.distance_m) for point in layer]
            for layer in layers.secondary_layers
        ],
        "sample_distances_m": layers.sample_distances_m,
        "sample_terrain_elevation_m": layers.sample_terrain_elevation_m,
    }


def _build_water_target_ground_sampler(
    *,
    viewer_data: ViewerData,
    scan_radius_km: float,
    deadline: float | None,
) -> Callable[[float, float], float] | None:
    if _timed_out(deadline):
        raise TimeoutError("water timed out")
    try:
        download = fetch_copernicus_dem(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            max_distance_km=scan_radius_km,
            margin_km=10.0,
            cache_dir=Path(CACHE_PATH) / "copernicus-dem",
        )
    except Exception:
        return None

    dem = GeoTiffDem(download.paths, default_elevation_m=0.0)
    try:
        bbox = build_download_bbox(
            lat_deg=float(viewer_data.lat_deg),
            lon_deg=float(viewer_data.lon_deg),
            radius_km=scan_radius_km + 10.0,
        )
        dem_grid = dem.build_grid(bbox)
    except Exception:
        dem.close()
        return None

    dem.close()

    def sampler(latitude_deg: float, longitude_deg: float) -> float:
        return sample_ground_elevation(
            dem_grid,
            latitude_deg=float(latitude_deg),
            longitude_deg=float(longitude_deg),
            dem_resampling="bilinear",
        )

    return sampler


def _load_or_fetch_water_overlay_footprints(
    *,
    viewer_data: ViewerData,
    scan_radius_km: float,
    deadline: float | None,
) -> tuple[WaterPolygonFootprint, ...]:
    now_utc = datetime.now(timezone.utc)
    scope_key = water_overlay_cache_scope_key(
        observer_lat_deg=float(viewer_data.lat_deg),
        observer_lon_deg=float(viewer_data.lon_deg),
        radius_km=float(scan_radius_km),
    )
    snapshot = load_water_overlay_cache(
        scope_key,
        observer_lat_deg=float(viewer_data.lat_deg),
        observer_lon_deg=float(viewer_data.lon_deg),
    )
    if snapshot is not None and water_overlay_cache_is_recent(
        snapshot,
        now_utc=now_utc,
        max_age_seconds=WATER_OVERLAY_CACHE_RETENTION_SECONDS,
    ):
        logger.info(
            "Water overlay cache hit: scope=%s footprints=%d",
            scope_key,
            len(snapshot.footprints),
        )
        return snapshot.footprints

    remaining_s = _remaining_timeout_seconds(deadline)
    timeout_s = 60.0 if remaining_s is None else max(0.1, min(60.0, remaining_s))
    try:
        footprints = fetch_water_overlay_footprints(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            max_distance_km=float(scan_radius_km),
            endpoint=DEFAULT_WATER_OVERPASS_ENDPOINT,
            user_agent=DEFAULT_WATER_USER_AGENT,
            timeout_s=timeout_s,
        )
        footprints = simplify_water_footprints_for_observer(
            footprints,
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
        )
        fresh_snapshot = WaterOverlayCacheSnapshot(
            footprints=footprints,
            water_polygon_count=len(footprints),
            fetched_at_utc=now_utc,
        )
        save_water_overlay_cache(scope_key, fresh_snapshot)
        logger.info(
            "Water overlay cache miss: scope=%s footprints=%d",
            scope_key,
            len(footprints),
        )
        return footprints
    except Exception:
        if snapshot is not None and snapshot.footprints:
            logger.info(
                "Water overlay cache fallback: scope=%s footprints=%d",
                scope_key,
                len(snapshot.footprints),
            )
            return snapshot.footprints
        raise


def _fetch_water_overlay_dots_layer(
    *,
    viewer_data: ViewerData,
    surface_size_px: tuple[int, int],
    deadline: float | None,
    target_ground_sampler: Callable[[float, float], float] | None = None,
) -> list[WaterOverlayPoint] | None:
    if _timed_out(deadline):
        raise TimeoutError("water timed out")

    observer_ground_m = float(viewer_data.ground_elevation_m or 0.0)
    scan_radius_km = resolve_water_scan_radius_km(
        float(viewer_data.observer_height_m) + observer_ground_m,
        minimum_distance_km=DEFAULT_WATER_RADIUS_KM,
    )
    azimuth_step_deg = resolve_water_surface_azimuth_step_deg(*surface_size_px)
    sea_dots, band_stats = sample_water_surface_interface_points_with_stats(
        observer_lat_deg=float(viewer_data.lat_deg),
        observer_lon_deg=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m) + observer_ground_m,
        max_distance_km=scan_radius_km,
        azimuth_step_deg=azimuth_step_deg,
    )
    water_footprints = _load_or_fetch_water_overlay_footprints(
        viewer_data=viewer_data,
        scan_radius_km=scan_radius_km,
        deadline=deadline,
    )
    inland_dots = sample_water_overlay_points_for_observer(
        footprints=water_footprints,
        observer_lat_deg=float(viewer_data.lat_deg),
        observer_lon_deg=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m) + observer_ground_m,
        fallback_surface_height_m=float(observer_ground_m),
        target_ground_elevation_m_sampler=target_ground_sampler,
        max_distance_km=scan_radius_km,
        azimuth_step_deg=azimuth_step_deg,
        front_hemisphere_view_center=tuple(
            float(value) for value in viewer_data.view_center
        ),
        front_hemisphere_fov_deg=float(viewer_data.content_fov_deg),
    )
    water_dots = tuple(sea_dots) + tuple(inland_dots)
    nearest_distance_km = min(
        (float(dot.distance_km) for dot in water_dots), default=None
    )
    band_100_count, band_250_count, band_500_count = _water_overlay_band_counts(
        water_dots
    )
    for band_stat in band_stats:
        logger.info("Water band stats: %s", _water_overlay_band_stats_text(band_stat))
    if nearest_distance_km is None:
        logger.info(
            "Water mask dots: 0 visible, nearest sea dot n/a, bands: 125m=%d 250m=%d 500m=%d",
            band_100_count,
            band_250_count,
            band_500_count,
        )
    else:
        logger.info(
            "Water mask dots: %d visible, nearest sea dot %.3f km, bands: 125m=%d 250m=%d 500m=%d",
            len(water_dots),
            nearest_distance_km,
            band_100_count,
            band_250_count,
            band_500_count,
        )
    return list(water_dots) if water_dots else None


def _required_feature_types(feature_type: str) -> tuple[str, ...]:
    return ("building", "building_part") if feature_type == "both" else (feature_type,)


def _fetch_urban_outline_layer(
    *,
    viewer_data: ViewerData,
    runtime_options: SkyWindowRuntimeOptions,
    deadline: float | None,
) -> list[UrbanOutlinePolyline] | None:
    if _timed_out(deadline):
        raise TimeoutError("urban timed out")
    current_overture_release = resolve_overture_release_for_cache_root(
        cache_root_dir=Path(CACHE_PATH),
        now_utc=datetime.now(timezone.utc),
    )
    derived_root_dir = Path(OVERTURE_DERIVED_ROOT_DIR)
    required_feature_types = (
        ()
        if runtime_options.urban_outline_skyscraper_only
        else _required_feature_types(runtime_options.urban_outline_feature_type)
    )
    required_dirs: list[Path] = []
    for overture_feature_type in required_feature_types:
        dataset_name = (
            Path(derived_root_dir)
            / derive_dataset_name(
                float(viewer_data.lat_deg),
                float(viewer_data.lon_deg),
                float(runtime_options.urban_outline_radius_km),
                overture_feature_type,
                float(runtime_options.urban_outline_min_height_m),
            )
            / "bldg"
        )
        required_dirs.append(dataset_name)
        if dataset_name.exists() and not is_derived_dataset_stale(
            dataset_name,
            now_utc=datetime.now(timezone.utc),
            expected_overture_release=current_overture_release,
        ):
            continue
        import_overture_buildings(
            lat_deg=float(viewer_data.lat_deg),
            lon_deg=float(viewer_data.lon_deg),
            radius_km=float(runtime_options.urban_outline_radius_km),
            derived_root_dir=derived_root_dir,
            min_building_height_m=float(runtime_options.urban_outline_min_height_m),
            feature_type=overture_feature_type,
            fmt="geojsonseq",
            overturemaps_bin="overturemaps",
            dataset_name=dataset_name.parent.name,
            keep_download=None,
            no_stac=False,
            overture_release=current_overture_release,
            skip_release_lookup=True,
        )

    outlines = None
    if required_dirs:
        outlines = resolve_urban_outline_layer_for_viewer(
            viewer_data,
            derived_root_dir=derived_root_dir,
            derived_dirs=tuple(required_dirs),
            max_candidates=int(runtime_options.urban_outline_max_candidates),
            front_hemisphere_view_center=tuple(
                float(value) for value in viewer_data.view_center
            ),
            front_hemisphere_fov_deg=float(viewer_data.content_fov_deg),
        )

    skyscraper_tiles = ()
    if float(runtime_options.urban_outline_skyscraper_radius_km) > 0.0:
        skyscraper_tiles = select_skyscraper_seed_tiles_for_viewer(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            inner_radius_km=float(runtime_options.urban_outline_radius_km),
            outer_radius_km=float(runtime_options.urban_outline_skyscraper_radius_km),
            seed_file=Path(SKYSCRAPER_TILES_FILE),
        )
    skyscraper_derived_root = Path(OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR)
    skyscraper_dirs: list[Path] = []
    for tile in skyscraper_tiles:
        derived_dir = skyscraper_tile_derived_dir(
            tile, derived_root_dir=skyscraper_derived_root
        )
        skyscraper_dirs.append(derived_dir)
        if derived_dir.exists() and not is_derived_dataset_stale(
            derived_dir,
            now_utc=datetime.now(timezone.utc),
            expected_overture_release=current_overture_release,
        ):
            continue
        import_overture_buildings_for_bbox(
            bbox=(
                tile.envelope.min_lon_deg,
                tile.envelope.min_lat_deg,
                tile.envelope.max_lon_deg,
                tile.envelope.max_lat_deg,
            ),
            derived_root_dir=skyscraper_derived_root,
            min_building_height_m=max(
                150.0, float(runtime_options.urban_outline_min_height_m)
            ),
            feature_type="building",
            fmt="geojsonseq",
            overturemaps_bin="overturemaps",
            dataset_name=tile.cache_key,
            keep_download=None,
            no_stac=False,
            overture_release=current_overture_release,
            skip_release_lookup=True,
        )
    if skyscraper_dirs:
        skyscraper_outlines = resolve_urban_outline_layer_for_viewer(
            viewer_data,
            derived_root_dir=skyscraper_derived_root,
            derived_dirs=tuple(skyscraper_dirs),
            radius_km=float(runtime_options.urban_outline_skyscraper_radius_km),
            min_distance_km=float(runtime_options.urban_outline_radius_km),
            min_height_m=max(150.0, float(runtime_options.urban_outline_min_height_m)),
            max_candidates=int(runtime_options.urban_outline_max_candidates),
            front_hemisphere_view_center=tuple(
                float(value) for value in viewer_data.view_center
            ),
            front_hemisphere_fov_deg=float(viewer_data.content_fov_deg),
        )
        outlines = _merge_outline_layers(outlines, skyscraper_outlines)
    return outlines


def _merge_outline_layers(
    base_outlines: list[UrbanOutlinePolyline] | None,
    extra_outlines: list[UrbanOutlinePolyline] | None,
) -> list[UrbanOutlinePolyline] | None:
    merged: list[UrbanOutlinePolyline] = []
    if base_outlines:
        merged.extend(base_outlines)
    if extra_outlines:
        merged.extend(
            UrbanOutlinePolyline(
                points=list(outline.points),
                height_m=float(outline.height_m),
                distance_km=float(outline.distance_km),
                source="skyscraper",
            )
            for outline in extra_outlines
        )
    return merged or None


def _fetch_aircraft_snapshots(
    *,
    viewer_data: ViewerData,
    deadline: float | None,
) -> list[object] | None:
    if _timed_out(deadline):
        raise TimeoutError("aircraft timed out")
    remaining = _remaining_timeout_seconds(deadline)
    timeout_s = 20.0 if remaining is None else max(0.1, min(20.0, remaining))
    bbox = build_observer_bbox(float(viewer_data.lat_deg), float(viewer_data.lon_deg))
    fetched = fetch_cached_opensky_states(bbox, timeout_s=timeout_s)
    logger.info("Aircraft source: %s", fetched.source)
    return list(fetched.snapshots)


def _fetch_satellite_records_by_group(
    *,
    viewer_data: ViewerData,
    target_time_utc,
    deadline: float | None,
    enabled_groups: tuple[str, ...] = (
        SATELLITE_ISS_CACHE_KEY,
        SATELLITE_HORIZONS_CACHE_KEY,
    ),
) -> dict[str, list[dict[str, object]]] | None:
    if _timed_out(deadline):
        raise TimeoutError("satellites timed out")
    remaining = _remaining_timeout_seconds(deadline)
    timeout_s = 20.0 if remaining is None else max(0.1, min(20.0, remaining))
    records_by_group: dict[str, list[dict[str, object]]] = {}
    time_mode = classify_target_time(target_time_utc)
    for group_key in enabled_groups:
        fetched = resolve_satellite_elements_for_time(
            group_key,
            target_time_utc=target_time_utc,
            time_mode=time_mode,
            timeout_s=timeout_s,
            observer_lat=float(viewer_data.lat_deg),
            observer_lon=float(viewer_data.lon_deg),
            observer_height_m=float(viewer_data.observer_height_m),
        )
        logger.info("Satellite source [%s]: %s", group_key, fetched.source)
        records_by_group[group_key] = list(fetched.records)
    return records_by_group


def _render_image(
    *,
    image_size: tuple[int, int],
    scene: RenderSceneData,
    style: RenderStyle,
    compositor: SkyCompositorCache,
    draw_direction_grid: bool = False,
    search_overlay_target: SearchJumpTarget | None = None,
) -> QImage:
    width, height = image_size
    image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
    try:
        geometry = render_geometry.get_screen_geometry(
            width, height, scene.viewer.view_alt_deg
        )
        frame = FrameContext(
            viewer=scene.viewer,
            time_obj=scene.time_obj,
            geometry=geometry,
            viewport_rect=QRect(0, 0, width, height),
        )
        label_candidates: list[dict[str, object]] = []
        render_base_scene_into_painter(
            painter,
            frame=frame,
            scene=scene,
            style=style,
            hud=RenderHudState(
                mouse_pos=QPoint(),
                overlay_info_bottom_left=False,
                viewport_interaction_mode=False,
                viewport_interaction_stars=None,
                status_message=None,
            ),
            compositor=compositor,
            label_candidates=label_candidates,
            draw_labels=False,
        )
        if draw_direction_grid:
            render_guides.draw_direction_grid_overlay(
                painter,
                frame.geometry,
                frame.viewer,
                (width, height),
            )
        if search_overlay_target is not None:
            draw_search_target_overlay(
                painter,
                frame.geometry,
                search_overlay_target,
                viewer_data=frame.viewer,
                theme=style.theme,
                text_font=style.text_font,
                draw_marker=True,
                marker_scale=compute_star_render_upscale_factor(
                    frame.geometry.radius * 2,
                    style.star_render_expected_width,
                ),
                label_candidates=label_candidates,
            )
        if label_candidates:
            render_text._draw_label_candidates(
                painter, label_candidates, style.text_font
            )
    finally:
        painter.end()
    return image


def _build_render_style(
    *,
    text_font: QFont,
    status_line_font: QFont,
    catalogs: PreparedWindowCatalogs,
    user_options: SkyWindowUserOptions,
    runtime_options: SkyWindowRuntimeOptions,
    theme: ThemeStyle,
) -> RenderStyle:
    show_dso = catalogs.dso_catalog_np is not None
    if user_options.show_dso_initial is not None:
        show_dso = (
            bool(user_options.show_dso_initial) and catalogs.dso_catalog_np is not None
        )
    show_asterisms = (
        True
        if user_options.show_asterisms_initial is None
        else bool(user_options.show_asterisms_initial)
    )
    show_guidelines = (
        True
        if user_options.show_guidelines_initial is None
        else bool(user_options.show_guidelines_initial)
    )
    return RenderStyle(
        theme=theme,
        visual_preset=user_options.visual_preset,
        text_font=text_font,
        status_line_font=status_line_font,
        show_background_gradient=False,
        show_custom_window_frame=False,
        show_observation_info=False,
        show_dso=show_dso,
        show_asterisms=show_asterisms,
        show_guidelines=show_guidelines,
        enlarge_moon=bool(user_options.enlarge_moon),
        bright_bodies_mode=str(user_options.bright_bodies_mode),
        star_base_radius=float(user_options.star_base_radius),
        star_visibility_boost=float(user_options.star_visibility_boost),
        asterism_visibility_boost=float(user_options.asterism_visibility_boost),
        earth_guide_visibility_boost=float(user_options.earth_guide_visibility_boost),
        vmag_limit=float(user_options.vmag_limit),
        sky_disc_altaz_rings=str(user_options.sky_disc_altaz_rings),
        sky_disc_altaz_rings_hover=str(user_options.sky_disc_altaz_rings_hover),
        cloud_disc_alpha=float(user_options.cloud_disc_alpha),
        satellite_opacity=float(user_options.satellite_opacity),
        terrain_horizon_opacity=float(user_options.terrain_horizon_opacity),
        earth_guide_opacity=float(user_options.earth_guide_opacity),
        night_light_opacity=float(user_options.night_light_opacity),
        urban_outline_opacity=float(user_options.urban_outline_opacity),
        show_urban_outline_layer=float(user_options.urban_outline_opacity) > 0.0,
        aircraft_opacity=float(user_options.aircraft_opacity),
        star_render_expected_width=int(runtime_options.star_render_expected_width),
    )


def _require_img2sixel_binary() -> str:
    executable = shutil.which("img2sixel")
    if executable:
        return executable
    logger.error("--sixel was requested, but 'img2sixel' was not found in PATH.")
    raise SystemExit(1)


def _response_indicates_sixel_support(response: bytes) -> bool:
    if not response.startswith(b"\x1b[") or not response.endswith(b"c"):
        return False
    for token in response[2:-1].split(b";"):
        token = token.lstrip(b"?")
        if token == b"4":
            return True
    return False


def _require_sixel_terminal_support(timeout_seconds: float = 0.25) -> None:
    if os.environ.get("TERM", "").startswith("yaft"):
        return
    if termios is None:
        logger.error(
            "--sixel was requested, but terminal control is unavailable on this platform."
        )
        raise SystemExit(1)
    tty_fd: int | None = None
    old_attrs = None
    response = bytearray()
    try:
        tty_fd = os.open("/dev/tty", os.O_RDWR | os.O_NOCTTY)
        old_attrs = termios.tcgetattr(tty_fd)
        new_attrs = termios.tcgetattr(tty_fd)
        new_attrs[3] &= ~(termios.ECHO | termios.ICANON)
        new_attrs[6][termios.VMIN] = 0
        new_attrs[6][termios.VTIME] = 0
        termios.tcsetattr(tty_fd, termios.TCSANOW, new_attrs)
        os.write(tty_fd, b"\x1b[c")
        deadline = time.monotonic() + timeout_seconds
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            ready, _, _ = select.select([tty_fd], [], [], remaining)
            if not ready:
                break
            chunk = os.read(tty_fd, 1024)
            if not chunk:
                break
            response.extend(chunk)
            if b"c" in chunk:
                break
    except OSError:
        logger.error(
            "--sixel was requested, but the terminal does not report SIXEL graphics support."
        )
        raise SystemExit(1)
    finally:
        if tty_fd is not None and old_attrs is not None:
            try:
                termios.tcsetattr(tty_fd, termios.TCSANOW, old_attrs)
            except OSError:
                pass
        if tty_fd is not None:
            try:
                os.close(tty_fd)
            except OSError:
                pass

    if not _response_indicates_sixel_support(bytes(response)):
        logger.error(
            "--sixel was requested, but the terminal does not report SIXEL graphics support."
        )
        raise SystemExit(1)


def _encode_image_as_png_bytes(image: QImage) -> bytes:
    ba = QByteArray()
    buf = QBuffer(ba)
    if not buf.open(QIODevice.OpenModeFlag.WriteOnly):
        logger.error("Failed to open in-memory buffer for PNG encoding.")
        raise SystemExit(1)
    try:
        if not image.save(buf, "PNG"):
            logger.error("Failed to encode image as PNG for SIXEL output.")
            raise SystemExit(1)
    finally:
        buf.close()
    return bytes(ba.data())


def _write_sixel_to_stdout(image: QImage, *, img2sixel_bin: str) -> bool:
    png_bytes = _encode_image_as_png_bytes(image)
    try:
        proc = subprocess.run(
            [img2sixel_bin, "-"],
            input=png_bytes,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except OSError as exc:
        logger.warning("SIXEL output failed: %s", exc)
        return False
    if proc.returncode != 0:
        stderr_text = proc.stderr.decode("utf-8", errors="replace").strip()
        if stderr_text:
            logger.warning("SIXEL output failed: %s", stderr_text)
        else:
            logger.warning("SIXEL output failed with exit status %s.", proc.returncode)
        return False
    sys.stdout.buffer.write(proc.stdout)
    sys.stdout.buffer.flush()
    return True


def _write_png_to_stdout(image: QImage) -> bool:
    png_bytes = _encode_image_as_png_bytes(image)
    try:
        sys.stdout.buffer.write(png_bytes)
        sys.stdout.buffer.flush()
    except OSError as exc:
        logger.error("Failed to write PNG image to stdout: %s", exc)
        return False
    return True


def _write_export_overlay_summary_to_stderr(
    *,
    viewer_data: ViewerData,
    celestial_data: CelestialData,
    vmag_limit: float,
    cloud_coverage_ratio: float | None = None,
    search_overlay_target: SearchJumpTarget | None = None,
) -> None:
    lines = render_background.format_overlay_info_lines(
        celestial_data,
        viewer_data,
        vmag_limit,
        include_vmag_limit=True,
    )
    if cloud_coverage_ratio is None:
        lines.append("Cloud coverage n/a")
    else:
        lines.append(f"Cloud coverage {float(cloud_coverage_ratio) * 100.0:.1f}%")
    if search_overlay_target is not None:
        lines.append(_format_search_target_line(search_overlay_target))
    sys.stderr.write("\n".join(lines) + "\n")
    sys.stderr.flush()


def _format_search_target_line(target: SearchJumpTarget) -> str:
    parts = [
        f"Search target label={target.label}",
        f"id={target.object_key or target.command or target.label}",
        f"kind={target.kind}",
    ]
    if target.alt_deg is not None:
        parts.append(f"alt={float(target.alt_deg):.1f} deg")
    if target.az_deg is not None:
        parts.append(f"az={float(target.az_deg):.1f} deg")
    return " | ".join(parts)


def _format_search_candidate_line(target: SearchJumpTarget) -> str:
    parts = [
        f"label={target.label}",
        f"id={target.object_key or target.command or target.label}",
        f"kind={target.kind}",
    ]
    return " | ".join(parts)


def _format_search_failure_message(query: str, candidate_count: int) -> str:
    if candidate_count <= 0:
        return f"No search candidates found for {query!r}"
    return f"Search query {query!r} matched {candidate_count} candidates"


def _clamp_view_center_altitude(alt_deg: float) -> float:
    return max(
        float(OBSERVER_MIN_ALT_DEG), min(float(OBSERVER_MAX_ALT_DEG), float(alt_deg))
    )


def _search_view_center_for_target(
    *,
    base_view_center: tuple[float, float],
    target_alt_deg: float,
    target_az_deg: float,
    fixed_alt: bool,
    fixed_az: bool,
) -> tuple[float, float]:
    view_center_alt = (
        float(base_view_center[0])
        if fixed_alt
        else _clamp_view_center_altitude(target_alt_deg)
    )
    view_center_az = (
        float(base_view_center[1]) if fixed_az else float(target_az_deg) % 360.0
    )
    return view_center_alt, view_center_az


def main() -> None:
    args = parse_export_image_args()
    if args.print_cache_dir:
        print(CACHE_PATH)
        return
    setup_root_logger()
    logger.info("%s export-image starting...", APP_DISPLAY_NAME)
    if args.clear_long_lived_cache:
        try:
            logger.info("Clearing long-lived cache on user request...")
            clear_long_lived_cache()
        except LongLivedCacheClearCooldownError as exc:
            logger.error("%s", exc)
            raise SystemExit(1)
    wants_sixel = bool(args.sixel)
    img2sixel_bin = _require_img2sixel_binary() if wants_sixel else None
    if wants_sixel:
        _require_sixel_terminal_support()

    try:
        catalogs, viewer_data, user_options, runtime_options, search_overlay_target = (
            _build_window_inputs_from_args(args)
        )
    except LaunchSetupError:
        raise SystemExit(1)

    app = setup_app(f"{APP_DISPLAY_NAME} Export Image")
    app.setQuitOnLastWindowClosed(False)

    text_font, status_line_font = _load_fonts(user_options.overlay_font_size)
    compositor = _build_compositor(runtime_options, user_options)
    output_arg = args.output
    output_path = None if output_arg in {None, "-"} else Path(output_arg).expanduser()
    image_size = tuple(int(v) for v in args.image_size)
    layer_timeout_seconds = float(args.layer_timeout_seconds)
    allow_partial_data = bool(args.allow_partial_data)

    use_lod6_catalog = float(user_options.vmag_limit) <= 6.0
    star_catalog = catalogs.star_catalog_np
    star_subset_indices = (
        catalogs.star_catalog_lod6_indices if use_lod6_catalog else None
    )
    star_vmag_limit = None if use_lod6_catalog else float(user_options.vmag_limit)
    theme = THEME_STYLES_BY_PRESET.get(
        user_options.visual_preset, THEME_STYLES_BY_PRESET["night"]
    )
    ephemeris = load_ephemeris()
    sky_payload = compute_sky_snapshot(
        ephemeris=ephemeris,
        viewer_data=viewer_data,
        geometry=render_geometry.get_screen_geometry(
            max(2, int(image_size[0])),
            max(2, int(image_size[1])),
            viewer_data.view_alt_deg,
        ),
        star_catalog=star_catalog,
        dso_catalog=catalogs.dso_catalog_np,
        star_vmag_limit=star_vmag_limit,
        star_subset_indices=star_subset_indices,
        delta_t=runtime_options.delta_t,
        sky_disc_alpha=float(user_options.sky_disc_alpha),
        theme=theme,
        star_catalog_meta=catalogs.star_catalog_meta,
        image_size=(int(image_size[0]), int(image_size[1])),
        render_generation=0,
    )
    celestial_data = sky_payload["celestial"]
    sky_disc_image = sky_payload["sky_disc"]
    logger.info("Initial sky data ready.")

    layer_failures: list[str] = []
    cloud_image = None
    cloud_missing_mask = None
    cloud_amount_field = None
    cloud_altaz_grid = None
    cloud_coverage_ratio: float | None = None
    cloud_fetch_thread: threading.Thread | None = None
    cloud_fetch_done: threading.Event | None = None
    cloud_fetch_state: dict[str, object] = {}
    cloud_deadline: float | None = None
    use_geo_satellite = bool(
        user_options.geo_satellite
        and is_within_europe_band(
            float(viewer_data.lat_deg), float(viewer_data.lon_deg)
        )
    )
    if user_options.cloud_disc_alpha > 0.0:
        logger.info("Fetching initial cloud data...")
        cloud_deadline = _deadline_after(layer_timeout_seconds)
        (
            cloud_fetch_thread,
            cloud_fetch_done,
            cloud_fetch_state,
        ) = _start_cloud_layer_fetch(
            viewer_data=viewer_data,
            user_options=user_options,
            deadline=cloud_deadline,
        )

    aircraft_snapshots = None
    aircraft_fetch_thread: threading.Thread | None = None
    aircraft_fetch_done: threading.Event | None = None
    aircraft_fetch_state: dict[str, object] = {}
    if user_options.aircraft_opacity > 0.0:
        logger.info("Fetching initial aircraft state...")
        aircraft_deadline = _deadline_after(layer_timeout_seconds)
        (
            aircraft_fetch_thread,
            aircraft_fetch_done,
            aircraft_fetch_state,
        ) = _start_background_task(
            name="zstarview-export-aircraft",
            target=lambda: _fetch_aircraft_snapshots(
                viewer_data=viewer_data,
                deadline=aircraft_deadline,
            ),
        )

    satellite_records_by_group = None
    satellite_fetch_thread: threading.Thread | None = None
    satellite_fetch_done: threading.Event | None = None
    satellite_fetch_state: dict[str, object] = {}
    if user_options.satellite_opacity > 0.0:
        logger.info("Fetching initial satellite data...")
        satellite_deadline = _deadline_after(layer_timeout_seconds)
        (
            satellite_fetch_thread,
            satellite_fetch_done,
            satellite_fetch_state,
        ) = _start_background_task(
            name="zstarview-export-satellite",
            target=lambda: _fetch_satellite_records_by_group(
                viewer_data=viewer_data,
                target_time_utc=celestial_data.time.to_datetime(timezone=timezone.utc),
                deadline=satellite_deadline,
                enabled_groups=(
                    SATELLITE_ISS_CACHE_KEY,
                    SATELLITE_HORIZONS_CACHE_KEY,
                ),
            ),
        )

    terrain_horizon_profile = None
    terrain_horizon_profile_distances_m = None
    terrain_secondary_ridges_altaz_layers = None
    terrain_secondary_ridges_distances_m_layers = None
    terrain_horizon_payload: TerrainHorizonPayload | None = None
    terrain_fetch_thread: threading.Thread | None = None
    terrain_fetch_done: threading.Event | None = None
    terrain_fetch_state: dict[str, object] = {}
    if user_options.terrain_horizon_opacity > 0.0:
        logger.info("Fetching initial terrain horizon data...")
        terrain_deadline = _deadline_after(layer_timeout_seconds)
        (
            terrain_fetch_thread,
            terrain_fetch_done,
            terrain_fetch_state,
        ) = _start_background_task(
            name="zstarview-export-terrain",
            target=lambda: _fetch_terrain_horizon_layer(
                viewer_data=viewer_data,
                deadline=terrain_deadline,
            ),
        )

    if terrain_fetch_thread is not None and terrain_fetch_done is not None:
        terrain_state = _await_background_task_result(
            label="terrain",
            thread=terrain_fetch_thread,
            done=terrain_fetch_done,
            state=terrain_fetch_state,
            deadline=terrain_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        if terrain_state is not None:
            terrain_value = terrain_state.get("value")
            if isinstance(terrain_value, dict):
                terrain_horizon_payload = terrain_value
                terrain_horizon_profile = terrain_horizon_payload["profile_altaz"]
                terrain_horizon_profile_distances_m = terrain_horizon_payload[
                    "profile_distances_m"
                ]
                terrain_secondary_ridges_altaz_layers = terrain_horizon_payload[
                    "secondary_ridges_altaz_layers"
                ]
                terrain_secondary_ridges_distances_m_layers = terrain_horizon_payload[
                    "secondary_ridges_distances_m_layers"
                ]
                logger.info("Initial terrain horizon data ready.")

    urban_outlines = None
    urban_fetch_thread: threading.Thread | None = None
    urban_fetch_done: threading.Event | None = None
    urban_fetch_state: dict[str, object] = {}
    if user_options.urban_outline_opacity > 0.0:
        logger.info("Fetching initial urban outline data...")
        urban_deadline = _deadline_after(layer_timeout_seconds)
        (
            urban_fetch_thread,
            urban_fetch_done,
            urban_fetch_state,
        ) = _start_background_task(
            name="zstarview-export-urban",
            target=lambda: _fetch_urban_outline_layer(
                viewer_data=viewer_data,
                runtime_options=runtime_options,
                deadline=urban_deadline,
            ),
        )

    night_light_glow_profile = None
    night_light_fetch_thread: threading.Thread | None = None
    night_light_fetch_done: threading.Event | None = None
    night_light_fetch_state: dict[str, object] = {}
    night_light_deadline: float | None = None
    sun_alt_deg = None
    for body in celestial_data.planets:
        if body.name == "sun":
            sun_alt_deg = float(body.alt)
            break
    if sun_alt_deg is not None and is_night_light_enabled(float(sun_alt_deg)):
        logger.info("Calculating initial night light alpha grid...")
        night_light_deadline = _deadline_after(layer_timeout_seconds)
        night_light_fetch_thread, night_light_fetch_done, night_light_fetch_state = (
            _start_background_task(
                name="zstarview-export-night-light",
                target=lambda: compute_night_light_glow_profile(
                    observer_lat_deg=float(viewer_data.lat_deg),
                    observer_lon_deg=float(viewer_data.lon_deg),
                    sun_alt_deg=float(sun_alt_deg),
                    terrain_profile_altaz=terrain_horizon_profile,
                    terrain_profile_distances_m=terrain_horizon_profile_distances_m,
                    terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
                    terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
                    terrain_sample_distances_m=terrain_horizon_payload.get(
                        "sample_distances_m"
                    )
                    if terrain_horizon_payload is not None
                    else None,
                    terrain_sample_terrain_elevation_m=terrain_horizon_payload.get(
                        "sample_terrain_elevation_m"
                    )
                    if terrain_horizon_payload is not None
                    else None,
                    include_night_light_tiles=float(
                        getattr(user_options, "night_light_opacity", 0.0)
                    )
                    > 0.0,
                ),
            )
        )

    water_overlay_dots = None
    water_fetch_thread: threading.Thread | None = None
    water_fetch_done: threading.Event | None = None
    water_fetch_state: dict[str, object] = {}
    water_overlay_opacity = float(user_options.water_overlay_opacity)
    if water_overlay_opacity > 0.0:
        logger.info("Fetching initial water surface mask...")
        water_deadline = _deadline_after(layer_timeout_seconds)
        # Match the GUI water path: use the terrain-controller ground height
        # already present in ViewerData, but do not refit inland-water
        # heights against a fresh DEM inside export-image.
        (
            water_fetch_thread,
            water_fetch_done,
            water_fetch_state,
        ) = _start_background_task(
            name="zstarview-export-water",
            target=lambda: _fetch_water_overlay_dots_layer(
                viewer_data=viewer_data,
                surface_size_px=tuple(int(value) for value in args.image_size),
                deadline=water_deadline,
                target_ground_sampler=None,
            ),
        )

    if aircraft_fetch_thread is not None and aircraft_fetch_done is not None:
        aircraft_state = _await_background_task_result(
            label="aircraft",
            thread=aircraft_fetch_thread,
            done=aircraft_fetch_done,
            state=aircraft_fetch_state,
            deadline=aircraft_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        if aircraft_state is not None:
            aircraft_value = aircraft_state.get("value")
            if aircraft_value is not None:
                aircraft_snapshots = list(aircraft_value)
            logger.info("Initial aircraft state ready.")

    if satellite_fetch_thread is not None and satellite_fetch_done is not None:
        satellite_state = _await_background_task_result(
            label="satellites",
            thread=satellite_fetch_thread,
            done=satellite_fetch_done,
            state=satellite_fetch_state,
            deadline=satellite_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        if satellite_state is not None:
            satellite_records_by_group = satellite_state.get("value")
            logger.info("Initial satellite data ready.")

    if urban_fetch_thread is not None and urban_fetch_done is not None:
        urban_state = _await_background_task_result(
            label="urban",
            thread=urban_fetch_thread,
            done=urban_fetch_done,
            state=urban_fetch_state,
            deadline=urban_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        if urban_state is not None:
            urban_outlines = urban_state.get("value")
            logger.info("Initial urban outline data ready.")

    if night_light_fetch_thread is not None and night_light_fetch_done is not None:
        night_light_state = _await_background_task_result(
            label="night lights",
            thread=night_light_fetch_thread,
            done=night_light_fetch_done,
            state=night_light_fetch_state,
            deadline=night_light_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        if night_light_state is not None:
            night_light_glow_profile = night_light_state.get("value")
            logger.info("Night light alpha grid computed.")

    if water_fetch_thread is not None and water_fetch_done is not None:
        water_state = _await_background_task_result(
            label="water",
            thread=water_fetch_thread,
            done=water_fetch_done,
            state=water_fetch_state,
            deadline=water_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        if water_state is not None:
            water_overlay_dots = water_state.get("value")
            logger.info("Initial water surface mask ready.")

    if cloud_fetch_thread is not None and cloud_fetch_done is not None:
        cloud_state = _await_background_task_result(
            label="Geo-sat" if use_geo_satellite else "cloud",
            thread=cloud_fetch_thread,
            done=cloud_fetch_done,
            state=cloud_fetch_state,
            deadline=cloud_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        cloud_value = None if cloud_state is None else cloud_state.get("value")
        if isinstance(cloud_value, tuple) and len(cloud_value) == 4:
            (
                cloud_image,
                cloud_missing_mask,
                cloud_amount_field,
                cloud_coverage_ratio,
            ) = cloud_value
            logger.info("Initial cloud data ready.")
        elif isinstance(cloud_value, tuple) and len(cloud_value) == 5:
            (
                cloud_image,
                cloud_missing_mask,
                cloud_amount_field,
                cloud_coverage_ratio,
                cloud_altaz_grid,
            ) = cloud_value
            logger.info("Initial cloud data ready.")

    if layer_failures and not allow_partial_data:
        _abort_export_without_partial_data()

    style = _build_render_style(
        text_font=text_font,
        status_line_font=status_line_font,
        catalogs=catalogs,
        user_options=user_options,
        runtime_options=runtime_options,
        theme=theme,
    )
    scene = RenderSceneData(
        viewer=viewer_data,
        celestial_data=celestial_data,
        sky_disc_image=sky_disc_image,
        cloud_image=cloud_image,
        cloud_missing_mask=cloud_missing_mask,
        cloud_amount_field=cloud_amount_field,
        cloud_altaz_grid=cloud_altaz_grid if isinstance(cloud_altaz_grid, CloudAltAzGrid) else None,
        terrain_horizon_profile=terrain_horizon_profile,
        terrain_horizon_profile_distances_m=terrain_horizon_profile_distances_m,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
        urban_outlines=urban_outlines,
        water_overlay_dots=water_overlay_dots,
        satellite_records_by_group=satellite_records_by_group,
        aircraft_snapshots=aircraft_snapshots,
        night_light_glow_profile=night_light_glow_profile,
    )
    image = _render_image(
        image_size=image_size,
        scene=scene,
        style=style,
        compositor=compositor,
        draw_direction_grid=bool(args.include_direction_grid),
        search_overlay_target=search_overlay_target,
    )
    saved_output = False
    if output_arg == "-":
        if not _write_png_to_stdout(image):
            raise SystemExit(1)
        saved_output = True
    elif output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if not image.save(str(output_path), "PNG"):
            logger.error("Failed to save image: %s", output_path)
            raise SystemExit(1)
        saved_output = True
        logger.info("Saved image: %s", output_path)

    if wants_sixel:
        _write_export_overlay_summary_to_stderr(
            viewer_data=viewer_data,
            celestial_data=celestial_data,
            vmag_limit=float(style.vmag_limit),
            cloud_coverage_ratio=cloud_coverage_ratio,
            search_overlay_target=search_overlay_target,
        )
        assert img2sixel_bin is not None
        sixel_ok = _write_sixel_to_stdout(image, img2sixel_bin=img2sixel_bin)
        if not sixel_ok and not saved_output:
            raise SystemExit(1)
        return

    _write_export_overlay_summary_to_stderr(
        viewer_data=viewer_data,
        celestial_data=celestial_data,
        vmag_limit=float(style.vmag_limit),
        cloud_coverage_ratio=cloud_coverage_ratio,
        search_overlay_target=search_overlay_target,
    )


if __name__ == "__main__":
    main()
