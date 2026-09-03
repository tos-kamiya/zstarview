"""Headless one-shot sky-image export.

``main`` stays here so existing tests can patch names on this module. Layer
fetching, launch-input construction, rendering, and output helpers live in
sibling ``export_image_*.py`` modules and are re-exported from this facade.
"""

from __future__ import annotations

# ruff: noqa: F401
import json
import logging
import math
import os
import select
import shutil
import subprocess
import sys
import threading
import time
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import TypedDict

import numpy as np
from pyproj import Transformer
from PySide6.QtCore import QBuffer, QByteArray, QIODevice, QPoint, QRect
from PySide6.QtGui import QFont, QFontDatabase, QImage, QImageWriter, QPainter

try:
    import termios
except ImportError:  # pragma: no cover - non-Unix fallback
    termios = None  # type: ignore[assignment]

from ..__about__ import __version__
from ..aircraft import (
    build_observer_bbox,
    fetch_cached_opensky_states,
)
from ..astro import load_ephemeris
from ..cache_maintenance import LongLivedCacheClearCooldownError, clear_long_lived_cache
from ..catalog import load_dso_catalog, load_star_catalog
from ..clouddisc import CloudDisc, CloudDiscConfig, VisibilityError
from ..clouddisc.altaz_grid import CloudAltAzGrid
from ..coastline_tiles import PREVIEW_RADIUS_KM, load_coastline_overlay_polylines
from ..config import open_meteo_noncommercial_terms_accepted
from ..data.building_source import select_prepared_building_source
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
from ..geosatellite.pipeline import is_within_europe_band, run_geo_satellite_pipeline
from ..geosatellite.projection import render_gray_image_to_cloud_rgba
from ..gui.composite import SkyCompositorCache
from ..gui.sky_worker import compute_sky_snapshot
from ..gui.water_overlay_cache import (
    WATER_OVERLAY_CACHE_RETENTION_SECONDS,
    WaterOverlayCacheSnapshot,
    load_water_overlay_cache,
    save_water_overlay_cache,
    water_overlay_cache_is_recent,
    water_overlay_cache_scope_key,
)
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
from ..location_resolver import (
    LocationResolveError,
    ResolvedLocation,
    resolve_launch_location,
)
from ..logging_utils import setup_root_logger
from ..night_lights import compute_night_light_glow_profile, is_night_light_enabled
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
    OVERLAY_FONT_SIZE_DEFAULT,
    OVERTURE_DERIVED_ROOT_DIR,
    OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR,
    STARS_CSV_FILE,
    STATUS_LINE_FONT_SIZE,
    TEXT_FONT_PATH,
    THEME_STYLES_BY_PRESET,
    ThemeStyle,
)
from ..precipitation import (
    PrecipitationRenderItem,
    fetch_open_meteo_precipitation,
    generate_precipitation_request_samples,
    project_precipitation_columns,
)
from ..render import background as render_background
from ..render import geometry as render_geometry
from ..render import guides as render_guides
from ..render import text as render_text
from ..render.molecular_cloud_overlay import is_molecular_cloud_cache_available
from ..render.pipeline import (
    FrameContext,
    RenderHudState,
    RenderSceneData,
    RenderStyle,
    compute_star_render_upscale_factor,
    render_base_scene_into_painter,
)
from ..render.search_overlay import draw_search_target_overlay
from ..road_night_lights import (
    ROAD_NIGHT_LIGHT_MAX_CANDIDATES,
    ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
    RoadNightLightPolyline,
    build_road_night_light_ground_sampler,
    clip_road_night_lights_to_annulus,
    load_or_fetch_road_night_lights_with_source,
    project_road_night_lights,
    select_road_night_light_way_candidates,
    simplify_road_night_light_way_for_observer,
)
from ..satellite_constants import SATELLITE_HORIZONS_CACHE_KEY, SATELLITE_ISS_CACHE_KEY
from ..satellites import resolve_satellite_elements_for_time
from ..search.jpl import resolve_jpl_target_state_vector, search_jpl_targets
from ..search.models import SearchJumpTarget
from ..search.resolver import compute_search_target_altaz, resolve_search_targets
from ..search.satellites import resolve_satellite_target_altaz, search_satellite_targets
from ..splash import setup_app
from ..terrain import (
    DEFAULT_TERRAIN_DISTANCE_SAMPLE_STEP_M,
    EARTH_MEAN_RADIUS_M,
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
    DEFAULT_WATER_BOUNDARY_RADIUS_KM,
    DEFAULT_WATER_OVERPASS_ENDPOINT,
    DEFAULT_WATER_RADIUS_KM,
    DEFAULT_WATER_USER_AGENT,
    WaterOverlayPoint,
    WaterPolygonFootprint,
    build_water_overlay_polylines,
    fetch_water_overlay_footprints,
    resolve_water_scan_radius_km,
    resolve_water_surface_azimuth_step_deg,
    sample_water_overlay_points,
    simplify_water_footprints_for_observer,
)
from ..water_surface_mesh import make_local_transformer
from .args import parse_export_image_args
from .export_image_inputs import (
    _build_window_inputs_from_args,
    _load_dso_catalog_for_export,
    _load_star_catalog_for_export,
    _require_open_meteo_consent_for_export,
    _verify_ephemeris_for_export,
)
from .export_image_layers import (
    _abort_export_without_partial_data,
    _await_background_task_result,
    _build_water_target_ground_sampler,
    _fetch_aircraft_snapshots,
    _fetch_cloud_layer,
    _fetch_precipitation_layer,
    _fetch_road_night_lights_layer,
    _fetch_satellite_records_by_group,
    _fetch_terrain_horizon_layer,
    _fetch_urban_outline_layer,
    _fetch_water_overlay_dots_layer,
    _fetch_water_overlay_layer,
    _load_or_fetch_water_overlay_footprints,
    _merge_outline_layers,
    _required_feature_types,
    _start_background_task,
    _start_cloud_layer_fetch,
)
from .export_image_output import (
    _build_export_image_metadata_payload,
    _clamp_view_center_altitude,
    _encode_image_as_png_bytes,
    _format_search_candidate_line,
    _format_search_failure_message,
    _format_search_target_line,
    _require_img2sixel_binary,
    _require_sixel_terminal_support,
    _resolved_location_metadata,
    _resolved_location_source,
    _response_indicates_sixel_support,
    _search_target_metadata,
    _search_view_center_for_target,
    _write_export_overlay_summary_to_stderr,
    _write_png_to_path,
    _write_png_to_stdout,
    _write_sixel_to_stdout,
)
from .export_image_render import (
    _build_compositor,
    _build_render_style,
    _load_fonts,
    _render_image,
)
from .export_image_support import (
    DEFAULT_CLOUD_ALT_MIN_DEG,
    DEFAULT_CLOUD_BASE_SIZE,
    DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
    DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL,
    DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG,
    DEFAULT_EXPORT_IMAGE_TERRAIN_SAMPLE_STEP_M,
    EXPORT_IMAGE_METADATA_SCHEMA,
    EXPORT_IMAGE_METADATA_TEXT_KEY,
    LOGGER_NAME,
    OPEN_METEO_CONSENT_REQUIRED_MESSAGE,
    TerrainHorizonPayload,
    _deadline_after,
    _remaining_timeout_seconds,
    _timed_out,
    _UrbanOutlineFetchResult,
    _water_overlay_band_counts,
    _water_overlay_band_stats_text,
    sample_water_overlay_points_for_observer,
)

logger = logging.getLogger(LOGGER_NAME)

def main() -> None:
    args = parse_export_image_args()
    if args.print_cache_dir:
        print(CACHE_PATH)
        return
    _require_open_meteo_consent_for_export(
        float(getattr(args, "precipitation_opacity", 0.0))
    )
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
        (
            catalogs,
            viewer_data,
            user_options,
            runtime_options,
            search_overlay_target,
            place_location,
        ) = _build_window_inputs_from_args(args)
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
            edge_fov_deg=viewer_data.edge_fov_deg,
            content_fov_deg=viewer_data.content_fov_deg,
        ),
        star_catalog=star_catalog,
        dso_catalog=catalogs.dso_catalog_np,
        star_vmag_limit=star_vmag_limit,
        star_subset_indices=star_subset_indices,
        star_data_policy="scenic_view_scoped",
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
    _cloud_image = None
    cloud_missing_mask = None
    _cloud_amount_field = None
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

    precipitation_columns = None
    precipitation_fetch_thread: threading.Thread | None = None
    precipitation_fetch_done: threading.Event | None = None
    precipitation_fetch_state: dict[str, object] = {}
    precipitation_deadline: float | None = None
    if float(user_options.precipitation_opacity) > 0.0:
        logger.info("Fetching initial precipitation forecast...")
        precipitation_deadline = _deadline_after(layer_timeout_seconds)
        (
            precipitation_fetch_thread,
            precipitation_fetch_done,
            precipitation_fetch_state,
        ) = _start_background_task(
            name="zstarview-export-precipitation",
            target=lambda: _fetch_precipitation_layer(
                viewer_data=viewer_data,
                target_time_utc=celestial_data.time.to_datetime(timezone=timezone.utc),
                deadline=precipitation_deadline,
            ),
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

    road_night_light_polylines = None
    road_fetch_thread: threading.Thread | None = None
    road_fetch_done: threading.Event | None = None
    road_fetch_state: dict[str, object] = {}
    road_deadline: float | None = None
    if float(user_options.road_light_opacity) > 0.0:
        logger.info("Fetching initial road night lights data...")
        road_deadline = _deadline_after(layer_timeout_seconds)
        (
            road_fetch_thread,
            road_fetch_done,
            road_fetch_state,
        ) = _start_background_task(
            name="zstarview-export-road-lights",
            target=lambda: _fetch_road_night_lights_layer(
                viewer_data=viewer_data,
                deadline=road_deadline,
                max_candidates=int(runtime_options.road_light_max_candidates),
            ),
        )

    urban_outlines = None
    urban_outline_source = None
    urban_outline_count = None
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
        terrain_ground_value = (
            terrain_horizon_payload.get("ground_elevation_m")
            if terrain_horizon_payload is not None
            else None
        )
        terrain_ground_elevation_m = (
            float(terrain_ground_value)
            if isinstance(terrain_ground_value, (int, float))
            else float(viewer_data.ground_elevation_m)
        )
        observer_elevation_m = max(
            0.0,
            terrain_ground_elevation_m + float(viewer_data.observer_height_m),
        )
        night_light_deadline = _deadline_after(layer_timeout_seconds)
        night_light_fetch_thread, night_light_fetch_done, night_light_fetch_state = (
            _start_background_task(
                name="zstarview-export-night-light",
                target=lambda: compute_night_light_glow_profile(
                    observer_lat_deg=float(viewer_data.lat_deg),
                    observer_lon_deg=float(viewer_data.lon_deg),
                    observer_elevation_m=observer_elevation_m,
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
                    include_night_light_tiles=float(user_options.night_light_opacity)
                    > 0.0,
                ),
            )
        )

    water_overlay_dots = None
    water_overlay_polylines = None
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
            target=lambda: _fetch_water_overlay_layer(
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
            urban_result = urban_state.get("value")
            if isinstance(urban_result, _UrbanOutlineFetchResult):
                urban_outlines = urban_result.outlines
                urban_outline_source = urban_result.source
                urban_outline_count = (
                    len(urban_outlines) if urban_outlines else None
                )
            logger.info("Initial urban outline data ready.")

    if road_fetch_thread is not None and road_fetch_done is not None:
        road_state = _await_background_task_result(
            label="road lights",
            thread=road_fetch_thread,
            done=road_fetch_done,
            state=road_fetch_state,
            deadline=road_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        if road_state is not None:
            road_value = road_state.get("value")
            if isinstance(road_value, list):
                road_night_light_polylines = road_value
            logger.info("Initial road night lights data ready.")

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
            water_value = water_state.get("value")
            if isinstance(water_value, dict):
                water_overlay_dots = water_value.get("dots")
                water_overlay_polylines = water_value.get("polylines")
            else:
                water_overlay_dots = water_value
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
                _cloud_image,
                cloud_missing_mask,
                _cloud_amount_field,
                cloud_coverage_ratio,
            ) = cloud_value
            logger.info("Initial cloud data ready.")
        elif isinstance(cloud_value, tuple) and len(cloud_value) == 5:
            (
                _cloud_image,
                cloud_missing_mask,
                _cloud_amount_field,
                cloud_coverage_ratio,
                cloud_altaz_grid,
            ) = cloud_value
            logger.info("Initial cloud data ready.")

    if (
        precipitation_fetch_thread is not None
        and precipitation_fetch_done is not None
    ):
        precipitation_state = _await_background_task_result(
            label="precipitation",
            thread=precipitation_fetch_thread,
            done=precipitation_fetch_done,
            state=precipitation_fetch_state,
            deadline=precipitation_deadline,
            layer_failures=layer_failures,
            allow_partial_data=allow_partial_data,
        )
        if precipitation_state is not None:
            precipitation_value = precipitation_state.get("value")
            if isinstance(precipitation_value, list):
                precipitation_columns = precipitation_value
            logger.info("Initial precipitation forecast ready.")
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
    if (
        float(user_options.precipitation_opacity) > 0.0
        and precipitation_columns is None
    ):
        style = replace(style, precipitation_opacity=0.0)
    scene = RenderSceneData(
        viewer=viewer_data,
        celestial_data=celestial_data,
        sky_disc_image=sky_disc_image,
        cloud_missing_mask=cloud_missing_mask,
        cloud_altaz_grid=cloud_altaz_grid
        if isinstance(cloud_altaz_grid, CloudAltAzGrid)
        else None,
        terrain_horizon_profile=terrain_horizon_profile,
        terrain_horizon_profile_distances_m=terrain_horizon_profile_distances_m,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
        urban_outlines=urban_outlines,
        water_overlay_dots=water_overlay_dots,
        water_overlay_polylines=water_overlay_polylines,
        road_night_light_polylines=road_night_light_polylines,
        precipitation_columns=precipitation_columns,
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
    metadata_payload = _build_export_image_metadata_payload(
        app_version=__version__,
        viewer_data=viewer_data,
        celestial_data=celestial_data,
        style=style,
        place_query=getattr(args, "place", None),
        place_location=place_location,
        search_overlay_target=search_overlay_target,
        cloud_coverage_ratio=cloud_coverage_ratio,
        urban_outline_source=urban_outline_source,
        urban_outline_count=urban_outline_count,
    )
    saved_output = False
    if output_arg == "-":
        if not _write_png_to_stdout(image, metadata_payload=metadata_payload):
            raise SystemExit(1)
        saved_output = True
    elif output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if not _write_png_to_path(
            image,
            output_path,
            metadata_payload=metadata_payload,
        ):
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
        sixel_ok = _write_sixel_to_stdout(
            image,
            img2sixel_bin=img2sixel_bin,
            metadata_payload=metadata_payload,
        )
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
