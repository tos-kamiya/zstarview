"""Fetch export-image overlay layers in background tasks."""

from __future__ import annotations

import threading
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pyproj import Transformer

from ..aircraft import build_observer_bbox, fetch_cached_opensky_states
from ..clouddisc import CloudDisc, CloudDiscConfig, VisibilityError
from ..clouddisc.altaz_grid import CloudAltAzGrid
from ..coastline_tiles import PREVIEW_RADIUS_KM, load_coastline_overlay_polylines
from ..data.import_overture_buildings import (
    derive_dataset_name,
    import_overture_buildings,
    import_overture_buildings_for_bbox,
    is_derived_dataset_stale,
)
from ..data.skyscraper_tiles import (
    SKYSCRAPER_TILES_FILE,
    skyscraper_tile_derived_dir,
)
from ..gui.water_overlay_cache import (
    WATER_OVERLAY_CACHE_RETENTION_SECONDS,
    WaterOverlayCacheSnapshot,
    load_water_overlay_cache,
    save_water_overlay_cache,
    water_overlay_cache_is_recent,
    water_overlay_cache_scope_key,
)
from ..gui.window_inputs import SkyWindowRuntimeOptions, SkyWindowUserOptions
from ..overlay_time import classify_target_time
from ..paths import (
    CACHE_PATH,
    CLOUD_SHELLS_KM,
    OVERTURE_DERIVED_ROOT_DIR,
    OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR,
)
from ..precipitation import (
    PrecipitationRenderItem,
    fetch_open_meteo_precipitation,
    generate_precipitation_request_samples,
    project_precipitation_columns,
)
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
from ..terrain import (
    EARTH_MEAN_RADIUS_M,
    WGS84_GEOD,
    GeoTiffDem,
    ObserverLocation,
    build_distance_samples,
    build_download_bbox,
    compute_flat_ground_horizon_layers,
    compute_horizon_layers,
    reduce_profile_to_altaz,
    sample_ground_elevation,
)
from ..types import UrbanOutlinePolyline, ViewerData
from ..water_overlay import (
    DEFAULT_WATER_BOUNDARY_RADIUS_KM,
    DEFAULT_WATER_OVERPASS_ENDPOINT,
    DEFAULT_WATER_RADIUS_KM,
    DEFAULT_WATER_USER_AGENT,
    WaterOverlayPoint,
    WaterPolygonFootprint,
    build_water_overlay_polylines,
    resolve_water_scan_radius_km,
    resolve_water_surface_azimuth_step_deg,
    simplify_water_footprints_for_observer,
)
from ..water_surface_mesh import make_local_transformer
from .export_image_support import (
    DEFAULT_CLOUD_ALT_MIN_DEG,
    DEFAULT_CLOUD_BASE_SIZE,
    DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
    DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG,
    DEFAULT_EXPORT_IMAGE_TERRAIN_SAMPLE_STEP_M,
    TerrainHorizonPayload,
    _remaining_timeout_seconds,
    _UrbanOutlineFetchResult,
    _water_overlay_band_counts,
    _water_overlay_band_stats_text,
    host,
    logger,
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
) -> tuple[
    np.ndarray | None,
    np.ndarray | None,
    object | None,
    float | None,
    CloudAltAzGrid | None,
]:
    if user_options.cloud_disc_alpha <= 0.0:
        return (None, None, None, None, None)
    if host()._timed_out(deadline):
        raise TimeoutError("cloud timed out")

    requested_geo_satellite = bool(getattr(user_options, "geo_satellite", False))
    within_geo_satellite_band = host().is_within_europe_band(
        float(viewer_data.lat_deg),
        float(viewer_data.lon_deg),
    )
    if requested_geo_satellite and within_geo_satellite_band:
        logger.info("Geo-sat + Downloading")
        result = host().run_geo_satellite_pipeline(
            observer_lat=float(viewer_data.lat_deg),
            observer_lon=float(viewer_data.lon_deg),
            alt=float(viewer_data.view_alt_deg),
            az=float(viewer_data.view_az_deg),
            fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
        )
        logger.info("Calculating initial cloud image...")
        download_result = result.download
        captured_at_utc = getattr(download_result, "captured_at_utc", None) or download_result.fetched_at_utc
        logger.info(
            "Geo-sat + %s",
            captured_at_utc.astimezone(timezone.utc).isoformat(),
        )
        cloud_rgba = host().render_gray_image_to_cloud_rgba(result.disc_gray)
        _cloud_amount_field = None
        missing_mask = None
        cloud_coverage_ratio = float(
            np.count_nonzero(cloud_rgba[..., 3]) / max(1, cloud_rgba[..., 3].size)
        )
        logger.info("Geo-sat + Projecting")
        return (
            cloud_rgba,
            missing_mask,
            _cloud_amount_field,
            cloud_coverage_ratio,
            result.altaz_grid,
        )

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

    logger.info("Building alt/az cloud grid...")
    source.altaz_grid = clouddisc.build_altaz_grid_from_source(
        source=source,
        lat=float(viewer_data.lat_deg),
        lon=float(viewer_data.lon_deg),
        cloud_shells_km=CLOUD_SHELLS_KM,
    )
    logger.info("Alt/az cloud grid ready.")

    logger.info("Calculating initial cloud image...")

    from ..clouddisc.altaz_render import (
        render_altaz_grid_circles,
        render_altaz_missing_mask,
    )

    cloud_rgba = render_altaz_grid_circles(
        source.altaz_grid,
        width=DEFAULT_CLOUD_BASE_SIZE * 2 + 1,
        height=DEFAULT_CLOUD_BASE_SIZE * 2 + 1,
        center_alt_deg=float(viewer_data.view_alt_deg),
        center_az_deg=float(viewer_data.view_az_deg),
        edge_fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
        mask_fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
    )
    missing_mask = render_altaz_missing_mask(
        source.altaz_grid,
        width=DEFAULT_CLOUD_BASE_SIZE * 2 + 1,
        height=DEFAULT_CLOUD_BASE_SIZE * 2 + 1,
        center_alt_deg=float(viewer_data.view_alt_deg),
        center_az_deg=float(viewer_data.view_az_deg),
        edge_fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
        mask_fov_deg=float(viewer_data.edge_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
    )
    missing_mask_alpha = np.where(missing_mask > 0, 255, 0).astype(np.uint8)
    _coverage_ratio = source.altaz_grid.coverage_ratio
    return (
        cloud_rgba,
        missing_mask_alpha,
        None,
        float(_coverage_ratio),
        source.altaz_grid,
    )

def _start_cloud_layer_fetch(
    *,
    viewer_data: ViewerData,
    user_options: SkyWindowUserOptions,
    deadline: float | None,
) -> tuple[threading.Thread, threading.Event, dict[str, object]]:
    return _start_background_task(
        name="zstarview-export-cloud",
        target=lambda: host()._fetch_cloud_layer(
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
    if host()._timed_out(deadline):
        raise TimeoutError("terrain timed out")
    try:
        download = host().fetch_copernicus_dem(
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
    if host()._timed_out(deadline):
        raise TimeoutError("water timed out")
    try:
        download = host().fetch_copernicus_dem(
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

def _fetch_road_night_lights_layer(
    *,
    viewer_data: ViewerData,
    deadline: float | None,
    max_candidates: int = ROAD_NIGHT_LIGHT_MAX_CANDIDATES,
) -> list[RoadNightLightPolyline] | None:
    if max(0, int(max_candidates)) == 0:
        return []
    if host()._timed_out(deadline):
        raise TimeoutError("road lights timed out")
    snapshot, cache_hit = load_or_fetch_road_night_lights_with_source(
        observer_lat_deg=float(viewer_data.lat_deg),
        observer_lon_deg=float(viewer_data.lon_deg),
        radius_km=ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
    )
    forward_transformer = make_local_transformer(
        float(viewer_data.lat_deg), float(viewer_data.lon_deg)
    )
    inverse_transformer = Transformer.from_crs(
        forward_transformer.target_crs, "EPSG:4326", always_xy=True
    )
    selected = select_road_night_light_way_candidates(
        snapshot.ways,
        observer_lat_deg=float(viewer_data.lat_deg),
        observer_lon_deg=float(viewer_data.lon_deg),
        max_candidates=max_candidates,
        forward_transformer=forward_transformer,
    )
    simplified = tuple(
        simplify_road_night_light_way_for_observer(
            way,
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            forward_transformer=forward_transformer,
            inverse_transformer=inverse_transformer,
        )
        for way in selected
    )
    clipped = clip_road_night_lights_to_annulus(
        simplified,
        observer_lat_deg=float(viewer_data.lat_deg),
        observer_lon_deg=float(viewer_data.lon_deg),
        forward_transformer=forward_transformer,
        inverse_transformer=inverse_transformer,
    )
    polylines: list[RoadNightLightPolyline] = []
    if clipped:
        ground_sampler = build_road_night_light_ground_sampler(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            radius_km=ROAD_NIGHT_LIGHT_MAX_DISTANCE_KM,
        )
        polylines = list(
            project_road_night_lights(
                clipped,
                observer_lat_deg=float(viewer_data.lat_deg),
                observer_lon_deg=float(viewer_data.lon_deg),
                observer_height_m=float(viewer_data.observer_height_m),
                ground_elevation_m_sampler=ground_sampler,
                forward_transformer=forward_transformer,
                inverse_transformer=inverse_transformer,
            )
        )
    logger.info(
        "Road night lights ready: source=%s ways=%d polylines=%d",
        "cache" if cache_hit else "API",
        len(snapshot.ways),
        len(polylines),
    )
    return polylines or None

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
        footprints = host().fetch_water_overlay_footprints(
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
    if host()._timed_out(deadline):
        raise TimeoutError("water timed out")

    observer_ground_m = float(viewer_data.ground_elevation_m or 0.0)
    scan_radius_km = resolve_water_scan_radius_km(
        float(viewer_data.observer_height_m) + observer_ground_m,
        minimum_distance_km=DEFAULT_WATER_RADIUS_KM,
    )
    azimuth_step_deg = resolve_water_surface_azimuth_step_deg(*surface_size_px)
    sea_dots, band_stats = host().sample_water_surface_interface_points_with_stats(
        observer_lat_deg=float(viewer_data.lat_deg),
        observer_lon_deg=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m) + observer_ground_m,
        max_distance_km=scan_radius_km,
        azimuth_step_deg=azimuth_step_deg,
    )
    water_footprints = host()._load_or_fetch_water_overlay_footprints(
        viewer_data=viewer_data,
        scan_radius_km=scan_radius_km,
        deadline=deadline,
    )
    inland_dots = host().sample_water_overlay_points_for_observer(
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

def _fetch_water_overlay_layer(
    *,
    viewer_data: ViewerData,
    surface_size_px: tuple[int, int],
    deadline: float | None,
    target_ground_sampler: Callable[[float, float], float] | None = None,
) -> dict[str, object]:
    dots = host()._fetch_water_overlay_dots_layer(
        viewer_data=viewer_data,
        surface_size_px=surface_size_px,
        deadline=deadline,
        target_ground_sampler=target_ground_sampler,
    )
    if not dots:
        return {"dots": dots, "polylines": []}
    observer_ground_m = float(viewer_data.ground_elevation_m or 0.0)
    scan_radius_km = resolve_water_scan_radius_km(
        float(viewer_data.observer_height_m) + observer_ground_m,
        minimum_distance_km=DEFAULT_WATER_RADIUS_KM,
    )
    footprints = _load_or_fetch_water_overlay_footprints(
        viewer_data=viewer_data,
        scan_radius_km=scan_radius_km,
        deadline=deadline,
    )
    polylines = list(
        build_water_overlay_polylines(
            footprints,
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            observer_height_m=float(viewer_data.observer_height_m) + observer_ground_m,
            fallback_surface_height_m=observer_ground_m,
            target_ground_elevation_m_sampler=target_ground_sampler,
            max_distance_km=DEFAULT_WATER_BOUNDARY_RADIUS_KM,
        )
    )
    polylines = [polyline for polyline in polylines if polyline.water_category != "coastline"]
    polylines.extend(
        load_coastline_overlay_polylines(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            observer_height_m=float(viewer_data.observer_height_m) + observer_ground_m,
            max_distance_km=PREVIEW_RADIUS_KM,
            view_center=tuple(float(value) for value in viewer_data.view_center),
            fov_deg=float(viewer_data.content_fov_deg),
        )
    )
    return {"dots": dots, "polylines": list(polylines)}

def _required_feature_types(feature_type: str) -> tuple[str, ...]:
    return ("building", "building_part") if feature_type == "both" else (feature_type,)

def _fetch_urban_outline_layer(
    *,
    viewer_data: ViewerData,
    runtime_options: SkyWindowRuntimeOptions,
    deadline: float | None,
) -> _UrbanOutlineFetchResult:
    if host()._timed_out(deadline):
        raise TimeoutError("urban timed out")
    if not runtime_options.urban_outline_skyscraper_only:
        building_source = host().select_prepared_building_source(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            radius_km=float(runtime_options.urban_outline_radius_km),
        )
        if building_source.source == "plateau":
            outlines = host().resolve_urban_outline_layer_for_viewer(
                viewer_data,
                derived_dirs=building_source.derived_dirs,
                max_candidates=int(runtime_options.urban_outline_max_candidates),
                front_hemisphere_view_center=tuple(
                    float(value) for value in viewer_data.view_center
                ),
                front_hemisphere_fov_deg=float(viewer_data.content_fov_deg),
            )
            return _UrbanOutlineFetchResult(outlines=outlines, source="PLATEAU")
    current_overture_release = host().resolve_overture_release_for_cache_root(
        cache_root_dir=Path(CACHE_PATH),
        now_utc=datetime.now(timezone.utc),
    )
    derived_root_dir = Path(OVERTURE_DERIVED_ROOT_DIR)
    required_feature_types = (
        ()
        if runtime_options.urban_outline_skyscraper_only
        else host()._required_feature_types(runtime_options.urban_outline_feature_type)
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
        outlines = host().resolve_urban_outline_layer_for_viewer(
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
        skyscraper_tiles = host().select_skyscraper_seed_tiles_for_viewer(
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
        skyscraper_outlines = host().resolve_urban_outline_layer_for_viewer(
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
    return _UrbanOutlineFetchResult(outlines=outlines, source="Overture Maps")

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
    if host()._timed_out(deadline):
        raise TimeoutError("aircraft timed out")
    remaining = _remaining_timeout_seconds(deadline)
    timeout_s = 20.0 if remaining is None else max(0.1, min(20.0, remaining))
    bbox = build_observer_bbox(float(viewer_data.lat_deg), float(viewer_data.lon_deg))
    fetched = fetch_cached_opensky_states(
        bbox,
        timeout_s=timeout_s,
        enforce_global_rate_limit=False,
    )
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
    if host()._timed_out(deadline):
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

def _fetch_precipitation_layer(
    *, viewer_data: ViewerData, deadline: float | None
) -> list[PrecipitationRenderItem]:
    if host()._timed_out(deadline):
        raise TimeoutError("precipitation timed out")
    remaining = _remaining_timeout_seconds(deadline)
    timeout_seconds = 20.0 if remaining is None else max(0.1, min(20.0, remaining))
    samples = generate_precipitation_request_samples(
        float(viewer_data.lat_deg), float(viewer_data.lon_deg)
    )
    snapshot = fetch_open_meteo_precipitation(
        samples, timeout_seconds=timeout_seconds
    )
    return list(project_precipitation_columns(snapshot, viewer_data))
