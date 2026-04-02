from __future__ import annotations

from datetime import timezone
import logging
import math
import shutil
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from PySide6.QtCore import QBuffer, QByteArray, QIODevice, QPoint, QRect
from PySide6.QtGui import QFont, QFontDatabase, QImage, QPainter

from ..aircraft import build_observer_bbox, fetch_cached_opensky_states, project_aircraft_snapshots
from ..overlay_time import classify_target_time, overlay_availability_for_delta
from ..satellites import project_satellite_records, resolve_satellite_elements_for_time
from ..astro import _starfield_load
from ..cache_maintenance import LongLivedCacheClearCooldownError, clear_long_lived_cache
from ..catalog import load_dso_catalog, load_star_catalog
from ..clouddisc import CloudDisc, CloudDiscConfig
from ..data.import_overture_buildings import (
    derive_dataset_name,
    import_overture_buildings,
    import_overture_buildings_for_bbox,
)
from ..data.skyscraper_tiles import (
    SKYSCRAPER_OUTER_RADIUS_KM,
    SKYSCRAPER_TILES_FILE,
    select_skyscraper_seed_tiles_for_viewer,
    skyscraper_tile_derived_dir,
)
from ..launch_time import (
    LaunchSetupError,
    parse_launch_time_arguments,
)
from ..location_resolver import LocationResolveError, resolve_launch_location
from ..logging_utils import setup_root_logger
from ..paths import (
    APP_DISPLAY_NAME,
    CACHE_PATH,
    CLOUD_MISSING_TINT_RGBA,
    CLOUD_SHELLS_KM,
    DSO_CSV_FILE,
    EPHEMERIS_FILENAME,
    OVERTURE_DERIVED_ROOT_DIR,
    OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR,
    STARS_CSV_FILE,
    STATUS_LINE_FONT_SIZE,
    TEXT_FONT_PATH,
    TEXT_FONT_SIZE,
)
from ..render import background as render_background
from ..render import geometry as render_geometry
from ..render.pipeline import (
    RenderHudState,
    RenderSceneData,
    RenderStyle,
    render_base_scene_into_painter,
)
from ..satellite_constants import SATELLITE_ISS_CACHE_KEY
from ..splash import setup_app
from ..terrain import (
    EARTH_MEAN_RADIUS_M,
    GeoTiffDem,
    ObserverLocation,
    WGS84_GEOD,
    build_distance_samples,
    build_download_bbox,
    compute_horizon_profile,
    fetch_copernicus_dem,
    reduce_profile_to_altaz,
    sample_ground_elevation,
)
from ..types import CelestialData, UrbanOutlinePolyline, ViewerData
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
from ..urban_outline_layer import resolve_urban_outline_layer_for_viewer
from .args import parse_export_image_args

logger = logging.getLogger(__name__)

DEFAULT_CLOUD_ALT_MIN_DEG = 1.0
DEFAULT_CLOUD_FOV_OVERSCAN_DEG = 2.0
DEFAULT_CLOUD_BASE_SIZE = 256


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
        _starfield_load(EPHEMERIS_FILENAME)
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


def _build_window_inputs_from_args(
    args: object,
) -> tuple[PreparedWindowCatalogs, ViewerData, SkyWindowUserOptions, SkyWindowRuntimeOptions]:
    try:
        city = resolve_launch_location(
            getattr(args, "city", ""),
            place_query=getattr(args, "place", None),
            place_countrycode=getattr(args, "place_countrycode", None),
            place_lang=getattr(args, "place_lang", "en"),
            use_building_top=bool(getattr(args, "use_building_top", False)),
        )
    except LocationResolveError as exc:
        raise LaunchSetupError() from exc
    if getattr(args, "timezone", None) is not None:
        city = replace(city, tz=getattr(args, "timezone"))
    delta_t = parse_launch_time_arguments(
        getattr(args, "datetime", None),
        getattr(args, "days", 0),
        getattr(args, "hours", 0),
        timezone_name=city.tz,
        timezone_override=getattr(args, "timezone", None),
    )
    overlay_availability = overlay_availability_for_delta(delta_t)
    star_catalog = _load_star_catalog_for_export(getattr(args, "vmag_limit", 6.0))
    dso_catalog = _load_dso_catalog_for_export()
    _verify_ephemeris_for_export()

    view_center = (getattr(args, "view_center_alt", 90.0), getattr(args, "view_center_az", 180.0))
    view_center = (min(90.0, max(0.0, view_center[0])), view_center[1] % 360.0)
    cloud_stripe_mode, cloud_stripe_count, cloud_stripe_width = getattr(
        args,
        "cloud_stripe",
        ("width", 50, 0.85),
    )
    visual_preset = getattr(args, "theme", "night")
    star_visibility_boost = 1.12 if visual_preset == "white" else 1.05 if visual_preset == "day" else 1.0
    vmag_brightness_scale = -math.log10(getattr(args, "vmag_brightness_multiplier", 2.5))

    catalogs = prepare_window_catalogs(
        star_catalog,
        dso_catalog=dso_catalog,
        vmag_brightness_scale=vmag_brightness_scale,
    )
    viewer_data = prepare_window_viewer_data(
        city.display_name,
        (city.lat, city.lon, city.tz),
        view_center,
        content_fov_deg=getattr(args, "content_fov_deg", 100.0),
        observer_height_m=(
            city.observer_height_m
            if getattr(args, "observer_height_m", None) is None
            else getattr(args, "observer_height_m")
        ),
        location_height_label=city.location_height_label,
        location_height_m=city.location_height_m,
        show_observer_height=getattr(args, "observer_height_m", None) is not None,
    )
    user_options = prepare_window_user_options(
        sky_disc_alpha=getattr(args, "sky_opacity", 0.15),
        cloud_disc_alpha=(
            0.0
            if (not overlay_availability.cloud) or cloud_stripe_count == 0 or cloud_stripe_width == 0.0
            else getattr(args, "cloud_opacity", 0.075)
        ),
        satellite_opacity=(getattr(args, "satellite_opacity", 0.5) if overlay_availability.satellite else 0.0),
        aircraft_opacity=(getattr(args, "aircraft_opacity", 0.5) if overlay_availability.aircraft else 0.0),
        terrain_horizon_opacity=getattr(args, "terrain_horizon_opacity", 0.05),
        urban_outline_opacity=getattr(args, "urban_outline_opacity", 0.2),
        ground_tint_opacity=getattr(args, "ground_tint_opacity", 0.1),
        enlarge_moon=bool(getattr(args, "enlarge_moon", False)),
        star_base_radius=getattr(args, "star_base_radius", 4.0),
        vmag_limit=getattr(args, "vmag_limit", 6.0),
        visual_preset=visual_preset,
        star_visibility_boost=star_visibility_boost,
        show_dso_initial=getattr(args, "show_dso_initial", None),
        show_asterisms_initial=getattr(args, "show_asterisms_initial", None),
        show_guidelines_initial=getattr(args, "show_guidelines_initial", None),
        show_overlay_info_initial=getattr(args, "show_observation_info_initial", None),
        sky_disc_gui_allowed=getattr(args, "sky_opacity", 0.15) > 0.0,
        cloud_gui_allowed=overlay_availability.cloud and getattr(args, "cloud_opacity", 0.075) > 0.0,
        satellite_gui_allowed=overlay_availability.satellite and getattr(args, "satellite_opacity", 0.5) > 0.0,
        aircraft_gui_allowed=overlay_availability.aircraft and getattr(args, "aircraft_opacity", 0.5) > 0.0,
        terrain_horizon_gui_allowed=getattr(args, "terrain_horizon_opacity", 0.05) > 0.0,
        urban_outline_gui_allowed=getattr(args, "urban_outline_opacity", 0.2) > 0.0,
    )
    runtime_options = prepare_window_runtime_options(
        delta_t=delta_t,
        sky_update_interval=60,
        urban_outline_radius_km=getattr(args, "urban_outline_radius_km", 2.5),
        urban_outline_skyscraper_radius_km=getattr(args, "urban_outline_skyscraper_radius_km", SKYSCRAPER_OUTER_RADIUS_KM),
        urban_outline_min_height_m=getattr(args, "urban_outline_min_height_m", 0.0),
        urban_outline_feature_type=getattr(args, "urban_outline_feature_type", "both"),
        urban_outline_skyscraper_only=bool(getattr(args, "urban_outline_skyscraper_only", False)),
        cloud_stripe_style=(cloud_stripe_count, cloud_stripe_width),
        cloud_stripe_mode=cloud_stripe_mode,
        cloud_missing_tint_opacity=getattr(args, "cloud_missing_tint_opacity", 0.0),
        star_render_expected_width=getattr(args, "expected_render_width", 600),
        content_fov_deg=getattr(args, "content_fov_deg", 100.0),
        window_geometry_arg=None,
    )
    return catalogs, viewer_data, user_options, runtime_options


def _load_fonts() -> tuple[QFont, QFont]:
    text_font_id = QFontDatabase.addApplicationFont(TEXT_FONT_PATH)
    text_font_family = QFontDatabase.applicationFontFamilies(text_font_id)[0]
    return (
        QFont(text_font_family, TEXT_FONT_SIZE),
        QFont(text_font_family, STATUS_LINE_FONT_SIZE),
    )


def _build_compositor(runtime_options: SkyWindowRuntimeOptions, user_options: SkyWindowUserOptions) -> SkyCompositorCache:
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


def _fetch_cloud_layer(
    *,
    viewer_data: ViewerData,
    user_options: SkyWindowUserOptions,
    deadline: float | None,
) -> tuple[np.ndarray | None, np.ndarray | None, object | None]:
    if user_options.cloud_disc_alpha <= 0.0:
        return (None, None, None)
    if _timed_out(deadline):
        raise TimeoutError("cloud timed out")

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
    source = clouddisc.fetch_source(
        lat=float(viewer_data.lat_deg),
        lon=float(viewer_data.lon_deg),
    )
    cloud_rgba, _meta, missing_mask, _coverage_ratio = clouddisc.render_from_source_with_coverage(
        source=source,
        lat=float(viewer_data.lat_deg),
        lon=float(viewer_data.lon_deg),
        alt=float(viewer_data.view_alt_deg),
        az=float(viewer_data.view_az_deg),
        radius_px=DEFAULT_CLOUD_BASE_SIZE,
        edge_fov_deg=float(viewer_data.content_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
        mask_fov_deg=float(viewer_data.content_fov_deg) + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
        cloud_shells_km=CLOUD_SHELLS_KM,
    )
    if _timed_out(deadline):
        raise TimeoutError("cloud timed out")
    missing_mask_alpha = np.where(missing_mask > 0, 255, 0).astype(np.uint8)
    cloud_amount_field = build_cloud_amount_field_from_rgba(cloud_rgba)
    return (cloud_rgba, missing_mask_alpha, cloud_amount_field)


def _fetch_terrain_horizon_layer(
    *,
    viewer_data: ViewerData,
    deadline: float | None,
) -> list[tuple[float, float]] | None:
    if _timed_out(deadline):
        raise TimeoutError("terrain timed out")
    try:
        download = fetch_copernicus_dem(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            max_distance_km=120.0,
            margin_km=10.0,
            cache_dir=Path(CACHE_PATH) / "copernicus-dem",
        )
    except RuntimeError as exc:
        if str(exc) != "No Copernicus DEM tiles were downloaded for the requested area.":
            raise
        return []
    dem = GeoTiffDem(download.paths, default_elevation_m=0.0)
    try:
        bbox = build_download_bbox(
            lat_deg=float(viewer_data.lat_deg),
            lon_deg=float(viewer_data.lon_deg),
            radius_km=130.0,
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
        points = compute_horizon_profile(
            dem_grid=dem_grid,
            geod=WGS84_GEOD,
            observer=observer,
            azimuth_step_deg=1.0,
            distance_samples_m=build_distance_samples(120.0, 90.0),
            dem_resampling="bilinear",
            earth_radius_m=EARTH_MEAN_RADIUS_M,
            refraction_coefficient=0.13,
        )
    finally:
        dem.close()
    if _timed_out(deadline):
        raise TimeoutError("terrain timed out")
    return reduce_profile_to_altaz(points)


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
    derived_root_dir = Path(OVERTURE_DERIVED_ROOT_DIR)
    required_feature_types = () if runtime_options.urban_outline_skyscraper_only else _required_feature_types(
        runtime_options.urban_outline_feature_type
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
        if dataset_name.exists():
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
        )
        if _timed_out(deadline):
            raise TimeoutError("urban timed out")

    outlines = None
    if required_dirs:
        outlines = resolve_urban_outline_layer_for_viewer(
            viewer_data,
            derived_root_dir=derived_root_dir,
            derived_dirs=tuple(required_dirs),
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
        derived_dir = skyscraper_tile_derived_dir(tile, derived_root_dir=skyscraper_derived_root)
        skyscraper_dirs.append(derived_dir)
        if derived_dir.exists():
            continue
        import_overture_buildings_for_bbox(
            bbox=(
                tile.envelope.min_lon_deg,
                tile.envelope.min_lat_deg,
                tile.envelope.max_lon_deg,
                tile.envelope.max_lat_deg,
            ),
            derived_root_dir=skyscraper_derived_root,
            min_building_height_m=max(150.0, float(runtime_options.urban_outline_min_height_m)),
            feature_type="building",
            fmt="geojsonseq",
            overturemaps_bin="overturemaps",
            dataset_name=tile.cache_key,
            keep_download=None,
            no_stac=False,
        )
        if _timed_out(deadline):
            raise TimeoutError("urban timed out")
    if skyscraper_dirs:
        skyscraper_outlines = resolve_urban_outline_layer_for_viewer(
            viewer_data,
            derived_root_dir=skyscraper_derived_root,
            derived_dirs=tuple(skyscraper_dirs),
            radius_km=float(runtime_options.urban_outline_skyscraper_radius_km),
            min_distance_km=float(runtime_options.urban_outline_radius_km),
            min_height_m=max(150.0, float(runtime_options.urban_outline_min_height_m)),
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
                source="skyscraper",
            )
            for outline in extra_outlines
        )
    return merged or None


def _fetch_aircraft_overlay_points(
    *,
    viewer_data: ViewerData,
    celestial_time_obj: object,
    deadline: float | None,
) -> object | None:
    if _timed_out(deadline):
        raise TimeoutError("aircraft timed out")
    remaining = _remaining_timeout_seconds(deadline)
    timeout_s = 20.0 if remaining is None else max(0.1, min(20.0, remaining))
    bbox = build_observer_bbox(float(viewer_data.lat_deg), float(viewer_data.lon_deg))
    fetched = fetch_cached_opensky_states(bbox, timeout_s=timeout_s)
    snapshots = fetched.snapshots
    logger.info("Aircraft source: %s", fetched.source)
    if _timed_out(deadline):
        raise TimeoutError("aircraft timed out")
    return project_aircraft_snapshots(
        snapshots,
        observer_lat=float(viewer_data.lat_deg),
        observer_lon=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        time_obj=celestial_time_obj,
    )


def _fetch_satellite_overlay_points(
    *,
    viewer_data: ViewerData,
    celestial_time_obj: object,
    target_time_utc,
    deadline: float | None,
    enabled_groups: tuple[str, ...] = (SATELLITE_ISS_CACHE_KEY,),
) -> object | None:
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
        )
        logger.info("Satellite source [%s]: %s", group_key, fetched.source)
        records_by_group[group_key] = list(fetched.records)
    if _timed_out(deadline):
        raise TimeoutError("satellites timed out")
    return project_satellite_records(
        records_by_group,
        observer_lat=float(viewer_data.lat_deg),
        observer_lon=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        time_obj=celestial_time_obj,
    )


def _render_image(
    *,
    image_size: tuple[int, int],
    scene: RenderSceneData,
    style: RenderStyle,
    compositor: SkyCompositorCache,
) -> QImage:
    width, height = image_size
    image = QImage(width, height, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
    try:
        geometry = render_geometry.get_screen_geometry(width, height, scene.viewer.view_alt_deg)
        render_base_scene_into_painter(
            painter,
            geometry=geometry,
            viewport_rect=QRect(0, 0, width, height),
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
) -> RenderStyle:
    show_dso = catalogs.dso_catalog_np is not None
    if user_options.show_dso_initial is not None:
        show_dso = bool(user_options.show_dso_initial) and catalogs.dso_catalog_np is not None
    show_asterisms = True if user_options.show_asterisms_initial is None else bool(user_options.show_asterisms_initial)
    show_guidelines = True if user_options.show_guidelines_initial is None else bool(user_options.show_guidelines_initial)
    return RenderStyle(
        visual_preset=user_options.visual_preset,
        text_font=text_font,
        status_line_font=status_line_font,
        show_background_gradient=False,
        show_overlay_info=False,
        show_dso=show_dso,
        show_asterisms=show_asterisms,
        show_guidelines=show_guidelines,
        enlarge_moon=bool(user_options.enlarge_moon),
        star_base_radius=float(user_options.star_base_radius),
        star_visibility_boost=float(user_options.star_visibility_boost),
        vmag_limit=float(user_options.vmag_limit),
        cloud_disc_alpha=float(user_options.cloud_disc_alpha),
        satellite_opacity=float(user_options.satellite_opacity),
        terrain_horizon_opacity=float(user_options.terrain_horizon_opacity),
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
) -> None:
    lines = render_background.format_overlay_info_lines(
        celestial_data,
        viewer_data,
        vmag_limit,
        include_vmag_limit=True,
    )
    sys.stderr.write("\n".join(lines) + "\n")
    sys.stderr.flush()


def main() -> None:
    args = parse_export_image_args()
    if getattr(args, "print_cache_dir", False):
        print(CACHE_PATH)
        return
    setup_root_logger()
    logger.info("%s export-image starting...", APP_DISPLAY_NAME)
    if getattr(args, "clear_long_lived_cache", False):
        try:
            logger.info("Clearing long-lived cache on user request...")
            clear_long_lived_cache()
        except LongLivedCacheClearCooldownError as exc:
            logger.error("%s", exc)
            raise SystemExit(1)
    wants_sixel = bool(getattr(args, "sixel", False))
    img2sixel_bin = _require_img2sixel_binary() if wants_sixel else None

    try:
        catalogs, viewer_data, user_options, runtime_options = _build_window_inputs_from_args(args)
    except LaunchSetupError:
        raise SystemExit(1)

    app = setup_app(f"{APP_DISPLAY_NAME} Export Image")
    app.setQuitOnLastWindowClosed(False)

    text_font, status_line_font = _load_fonts()
    compositor = _build_compositor(runtime_options, user_options)
    output_arg = getattr(args, "output", None)
    output_path = None if output_arg in {None, "-"} else Path(output_arg).expanduser()
    image_size = tuple(int(v) for v in getattr(args, "image_size"))
    deadline = _deadline_after(float(getattr(args, "layer_timeout_seconds", 30.0)))
    allow_partial_data = bool(getattr(args, "allow_partial_data", False))

    use_lod6_catalog = float(user_options.vmag_limit) <= 6.0
    star_catalog = catalogs.star_catalog_np
    star_subset_indices = catalogs.star_catalog_lod6_indices if use_lod6_catalog else None
    star_vmag_limit = None if use_lod6_catalog else float(user_options.vmag_limit)
    sky_payload = compute_sky_snapshot(
        lat=float(viewer_data.lat_deg),
        lon=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        view_center=tuple(viewer_data.view_center),
        star_catalog=star_catalog,
        dso_catalog=catalogs.dso_catalog_np,
        star_vmag_limit=star_vmag_limit,
        star_subset_indices=star_subset_indices,
        delta_t=runtime_options.delta_t,
        sky_disc_alpha=float(user_options.sky_disc_alpha),
        sky_disc_base_size=max(image_size),
        content_fov_deg=float(viewer_data.content_fov_deg),
        visual_preset=str(user_options.visual_preset),
        star_catalog_meta=getattr(catalogs, "star_catalog_meta", None),
        render_width_px=int(image_size[0]),
        render_height_px=int(image_size[1]),
        render_generation=0,
    )
    celestial_data = sky_payload["celestial"]
    sky_disc_image = sky_payload["sky_disc"]

    layer_failures: list[str] = []
    cloud_image = None
    cloud_missing_mask = None
    cloud_amount_field = None
    if user_options.cloud_disc_alpha > 0.0:
        try:
            cloud_image, cloud_missing_mask, cloud_amount_field = _fetch_cloud_layer(
                viewer_data=viewer_data,
                user_options=user_options,
                deadline=deadline,
            )
        except Exception as exc:
            logger.warning("Export layer unavailable: cloud (%s)", exc)
            layer_failures.append("cloud")
            if not allow_partial_data:
                raise SystemExit(1)

    terrain_horizon_profile = None
    if user_options.terrain_horizon_opacity > 0.0:
        try:
            terrain_horizon_profile = _fetch_terrain_horizon_layer(
                viewer_data=viewer_data,
                deadline=deadline,
            )
        except Exception as exc:
            logger.warning("Export layer unavailable: terrain (%s)", exc)
            layer_failures.append("terrain")
            if not allow_partial_data:
                raise SystemExit(1)

    urban_outlines = None
    if user_options.urban_outline_opacity > 0.0:
        try:
            urban_outlines = _fetch_urban_outline_layer(
                viewer_data=viewer_data,
                runtime_options=runtime_options,
                deadline=deadline,
            )
        except Exception as exc:
            logger.warning("Export layer unavailable: urban (%s)", exc)
            layer_failures.append("urban")
            if not allow_partial_data:
                raise SystemExit(1)

    aircraft_overlay_points = None
    if user_options.aircraft_opacity > 0.0:
        try:
            aircraft_overlay_points = _fetch_aircraft_overlay_points(
                viewer_data=viewer_data,
                celestial_time_obj=celestial_data.time,
                deadline=deadline,
            )
        except Exception as exc:
            logger.warning("Export layer unavailable: aircraft (%s)", exc)
            layer_failures.append("aircraft")
            if not allow_partial_data:
                raise SystemExit(1)

    satellite_overlay_points = None
    if user_options.satellite_opacity > 0.0:
        try:
            satellite_overlay_points = _fetch_satellite_overlay_points(
                viewer_data=viewer_data,
                celestial_time_obj=celestial_data.time,
                target_time_utc=celestial_data.time.to_datetime(timezone=timezone.utc),
                deadline=deadline,
            )
        except Exception as exc:
            logger.warning("Export layer unavailable: satellites (%s)", exc)
            layer_failures.append("satellites")
            if not allow_partial_data:
                raise SystemExit(1)

    if layer_failures and not allow_partial_data:
        logger.error("Export aborted because partial data is not allowed.")
        raise SystemExit(1)

    style = _build_render_style(
        text_font=text_font,
        status_line_font=status_line_font,
        catalogs=catalogs,
        user_options=user_options,
        runtime_options=runtime_options,
    )
    scene = RenderSceneData(
        viewer=viewer_data,
        celestial_data=celestial_data,
        sky_disc_image=sky_disc_image,
        cloud_image=cloud_image,
        cloud_missing_mask=cloud_missing_mask,
        cloud_amount_field=cloud_amount_field,
        terrain_horizon_profile=terrain_horizon_profile,
        urban_outlines=urban_outlines,
        satellite_overlay_points=satellite_overlay_points,
        aircraft_overlay_points=aircraft_overlay_points,
    )
    image = _render_image(
        image_size=image_size,
        scene=scene,
        style=style,
        compositor=compositor,
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
    )


if __name__ == "__main__":
    main()
