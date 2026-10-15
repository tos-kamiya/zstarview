"""Render cached history as a minute-stepped PNG sequence."""

from __future__ import annotations

import hashlib
import json
import math
import os
import pickle
import sys
import tempfile
from collections import deque
from concurrent.futures import Future, ThreadPoolExecutor, TimeoutError as FutureTimeout
from dataclasses import replace
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, cast
from zoneinfo import ZoneInfo

import numpy as np
from PySide6.QtGui import QColor, QImage, QPainter
from PySide6.QtWidgets import QApplication

from ..__about__ import __version__
from ..aircraft.cache import _cached_set_from_payload as decode_aircraft
from ..catalog import load_dso_catalog, load_star_catalog
from ..clouddisc.altaz_grid import CloudAltAzGrid
from ..clouddisc.types import CloudSourceData
from ..paths import (
    CACHE_PATH,
    DSO_CSV_FILE,
    EPHEMERIS_FILENAME,
    STARS_CSV_FILE,
    TEXT_FONT_PATH,
    THEME_STYLES_BY_PRESET,
)
from ..data.derived_tile_cache import (
    parse_derived_tile_buildings,
    select_derived_tile_envelopes,
)
from ..data.building_source import select_prepared_building_source
from ..data.urban_outline_from_buildings import compute_urban_outlines
from ..data.urban_outline_common import BuildingFootprint
from ..night_light_source import (
    _manifest_path as night_light_manifest_path,
    _required_tile_names,
    _sha256_file,
    _tile_path as night_light_tile_path,
    _validate_geotiff as validate_night_light_geotiff,
    _validate_manifest as validate_night_light_manifest,
)
from ..night_lights import compute_night_light_glow_profile, is_night_light_enabled
from ..night_lights_constants import NIGHT_LIGHTS_MAX_DISTANCE_KM
from ..satellites.cache import _load_cached_set_from_payload as decode_satellites
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
from ..terrain.dem import collect_copernicus_tile_keys
from ..urban_outline_layer import (
    _merge_building_footprints,
    _project_observer_centric_urban_outline_result,
)
from ..types import CelestialData, ViewerData
from ..gui.application_services import ApplicationServices
from ..gui.sky_worker import compute_sky_snapshot
from ..gui.window_inputs import (
    SkyWindowRuntimeOptions,
    prepare_window_catalogs,
    prepare_window_user_options,
    prepare_window_viewer_data,
)
from ..location_resolver.viewpoints import Viewpoint
from ..render import geometry as render_geometry
from ..render import terrain as render_terrain
from ..render.pipeline import compute_star_render_upscale_factor
from ..render.render_types import RenderSceneData
from ..replay.inventory import Inventory, contained
from ..location_resolver.resolve import ResolvedLocation
from ..replay.location import location_payload
from ..replay.models import MINUTE, HistoryRecord, iso
from ..replay.selection import HistoryIndex
from ..cli.export_image_render import (
    _build_compositor,
    _build_render_style,
    _load_fonts,
    _render_image,
)
from ..cli.export_image_support import (
    DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL,
    DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG,
    DEFAULT_EXPORT_IMAGE_TERRAIN_SAMPLE_STEP_M,
    TerrainHorizonPayload,
)


_STATIC_UNSUPPORTED = {
    "water_surface_opacity": "water-overlay",
    "road_light_opacity": "road-lights",
    "diffuse_sky_opacity": "diffuse-sky",
    "akari_ir_bands_opacity": "diffuse-sky",
}
_MAX_VOXEL_CHUNK_FRAMES = 10


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    ).encode("ascii")


def _frame_times(start: datetime, end: datetime) -> list[datetime]:
    times: list[datetime] = []
    value = start
    while value < end:
        times.append(value)
        value += MINUTE
    return times


def _record_key(record: HistoryRecord) -> tuple[str, str]:
    return record.layer, str(record.path)


def _chunk_frame_indices(
    render_indices: list[int],
    plans: list[dict[str, HistoryRecord | None]],
    decoded: dict[tuple[str, str], Any],
    *,
    voxel_enabled: bool,
    max_frames: int = _MAX_VOXEL_CHUNK_FRAMES,
) -> list[list[int]]:
    """Group adjacent frames by selected native-cloud record and chunk size."""
    chunks: list[list[int]] = []

    def source_key(index: int) -> tuple[str, str] | None:
        record = plans[index].get("cloud")
        data = decoded.get(_record_key(record)) if record is not None else None
        if voxel_enabled and record is not None and isinstance(data, CloudSourceData):
            return _record_key(record)
        return None

    for index in render_indices:
        if (
            chunks
            and len(chunks[-1]) < max_frames
            and chunks[-1][-1] + 1 == index
            and source_key(index) == source_key(chunks[-1][-1])
        ):
            chunks[-1].append(index)
        else:
            chunks.append([index])
    return chunks


def _selected_layers(args: Any, location: Any) -> dict[str, str]:
    selected: dict[str, str] = {}
    if args.cloud_opacity > 0:
        selected["cloud"] = "cloud-geo" if args.geo_satellite else "cloud-native"
    if args.aircraft_opacity > 0:
        selected["aircraft"] = "aircraft"
    if args.satellite_opacity > 0:
        selected["satellite-iss"] = "satellite-iss"
        selected["satellite-horizons"] = "satellite-horizons"
    if args.geo_satellite and selected.get("cloud") == "cloud-geo":
        selected["cloud"] = "cloud-geo"
    del location
    return selected


def _unsupported_layers(args: Any) -> list[str]:
    values = []
    for attribute, label in _STATIC_UNSUPPORTED.items():
        if float(getattr(args, attribute, 0.0) or 0.0) > 0.0 and label not in values:
            values.append(label)
    if args.geo_satellite and float(args.cloud_opacity) > 0.0:
        if "geo-satellite-cloud" not in values:
            values.append("geo-satellite-cloud")
    return values


def _read_verified_sources(
    inventory: Inventory,
    selected: dict[str, str],
    times: list[datetime],
) -> tuple[
    list[dict[str, HistoryRecord | None]],
    dict[tuple[str, str], Any],
    dict[str, dict[str, str]],
    dict[str, bytes],
]:
    indexes = {
        name: HistoryIndex(records) for name, records in inventory.records.items()
    }
    plans: list[dict[str, HistoryRecord | None]] = []
    unique: dict[tuple[str, str], HistoryRecord] = {}
    for frame_time in times:
        frame: dict[str, HistoryRecord | None] = {}
        for key, layer in selected.items():
            record = indexes[layer].select(frame_time)
            frame[key] = record
            if record is not None:
                unique[_record_key(record)] = record
        plans.append(frame)

    decoded: dict[tuple[str, str], Any] = {}
    provenance: dict[str, dict[str, str]] = {}
    frozen_bytes: dict[str, bytes] = {}
    for record_key, record in unique.items():
        try:
            path = contained(inventory.root, record.path)
            raw = path.read_bytes()
            metadata = json.loads(raw)
            if not isinstance(metadata, dict):
                raise ValueError("Cache record is not a JSON object")
            source_ref = {"record": str(path.relative_to(inventory.root))}
            frozen_bytes[frozen_bytes_key(path)] = raw
            if record.layer == "cloud-native":
                artifact = contained(
                    inventory.root, inventory.root / metadata["artifact_path"]
                )
                payload = artifact.read_bytes()
                if _digest(payload).lower() != str(metadata["artifact_sha256"]).lower():
                    raise ValueError("Native cloud artifact digest mismatch")
                source = pickle.loads(payload)
                if not isinstance(source, CloudSourceData) or not isinstance(
                    source.altaz_grid, CloudAltAzGrid
                ):
                    raise ValueError("Native cloud artifact has no renderable grid")
                decoded[record_key] = source
                frozen_bytes[frozen_bytes_key(artifact)] = payload
                source_ref["artifact"] = str(artifact.relative_to(inventory.root))
                source_ref["artifact_sha256"] = _digest(payload)
            elif record.layer == "aircraft":
                decoded[record_key] = decode_aircraft(metadata)
            elif record.layer.startswith("satellite-"):
                group = record.layer.removeprefix("satellite-")
                cached = decode_satellites(metadata, group_key=group)
                if cached is None:
                    raise ValueError("Satellite cache record could not be decoded")
                decoded[record_key] = cached
            elif record.layer == "cloud-geo":
                raise ValueError("Geo-satellite cloud rendering is not yet supported")
            else:
                raise ValueError("Unsupported replay layer")
            source_ref["record_sha256"] = _digest(raw)
            source_ref["fetched_at_utc"] = iso(record.fetched_at)
            source_ref["source_time_utc"] = iso(record.source_time)
            provenance[f"{record.layer}:{path.relative_to(inventory.root)}"] = (
                source_ref
            )
        except Exception as exc:
            decoded[record_key] = exc
    return plans, decoded, provenance, frozen_bytes


def frozen_bytes_key(path: Path) -> str:
    return str(path)


def _build_viewer_data(args: Any, location: Any) -> ViewerData:
    return prepare_window_viewer_data(
        location.display_name,
        (location.lat, location.lon, location.tz),
        (args.view_center_alt, args.view_center_az % 360.0),
        edge_fov_deg=args.edge_fov_deg,
        content_fov_deg=args.content_fov_deg,
        ground_elevation_m=location.ground_elevation_m,
        location_height_label=location.location_height_label,
        location_height_m=location.location_height_m,
        height_add_m=location.height_add_m,
    )


def _build_render_inputs(args: Any, location: Any, viewer: ViewerData | None = None):
    star_catalog = load_star_catalog(STARS_CSV_FILE, vmag_threshold=args.vmag_limit)
    try:
        dso_catalog = load_dso_catalog(DSO_CSV_FILE)
    except (OSError, ValueError):
        dso_catalog = None
    brightness_scale = -math.log10(float(args.vmag_brightness_multiplier))
    catalogs = prepare_window_catalogs(
        star_catalog,
        dso_catalog=dso_catalog,
        vmag_brightness_scale=brightness_scale,
    )
    if viewer is None:
        viewer = _build_viewer_data(args, location)
    user_options = prepare_window_user_options(
        sky_disc_alpha=args.sky_opacity,
        sky_disc_style=args.sky_disc_style,
        sky_disc_altaz_rings=args.sky_disc_altaz_rings,
        sky_disc_altaz_rings_hover=args.sky_disc_altaz_rings_hover,
        cloud_disc_alpha=args.cloud_opacity,
        satellite_opacity=args.satellite_opacity,
        aircraft_opacity=args.aircraft_opacity,
        terrain_horizon_opacity=args.terrain_horizon_opacity,
        earth_guide_opacity=args.earth_guide_opacity,
        urban_outline_opacity=args.urban_outline_opacity,
        bright_bodies_mode=args.bright_bodies,
        star_base_radius=args.star_base_radius,
        vmag_limit=args.vmag_limit,
        visual_preset=args.theme,
        star_visibility_boost=THEME_STYLES_BY_PRESET[args.theme].star_visibility_boost,
        visibility_boost=args.visibility_boost,
        show_dso_initial=args.show_dso_initial,
        show_asterisms_initial=args.show_asterisms_initial,
        show_guidelines_initial=args.show_guidelines_initial,
        observation_info_mode=args.observation_info,
        sky_disc_gui_allowed=args.sky_opacity > 0,
        cloud_gui_allowed=args.cloud_opacity > 0,
        satellite_gui_allowed=args.satellite_opacity > 0,
        aircraft_gui_allowed=args.aircraft_opacity > 0,
        tropical_cyclone_gui_allowed=False,
        terrain_horizon_gui_allowed=args.terrain_horizon_opacity > 0,
        earth_guide_gui_allowed=args.earth_guide_opacity > 0,
        geo_satellite=args.geo_satellite,
        night_light_opacity=float(getattr(args, "night_light_opacity", 0.0)),
        road_light_opacity=float(getattr(args, "road_light_opacity", 0.0)),
        diffuse_sky_source=str(getattr(args, "diffuse_sky_source", "gaia")),
        akari_ir_bands_opacity=float(
            getattr(args, "akari_ir_bands_opacity", 0.0) or 0.0
        ),
        ridge_glow_opacity=float(getattr(args, "ridge_glow_opacity", 0.0)),
        water_overlay_opacity=float(getattr(args, "water_surface_opacity", 0.0)),
        overlay_font_size=float(args.overlay_font_size),
        moon_style=str(args.moon_style),
        moon_scale=int(args.moon_scale),
        asterism_opacity=args.asterism_opacity,
        meteor_trails_gui_allowed=False,
        night_light_gui_allowed=False,
        road_light_gui_allowed=False,
        akari_ir_bands_gui_allowed=False,
        urban_outline_gui_allowed=False,
    )
    cloud_stripe_mode, stripe_count, stripe_width = args.cloud_stripe
    runtime = replace(
        SkyWindowRuntimeOptions(),
        cloud_stripe_style=(stripe_count, stripe_width),
        cloud_stripe_mode=cloud_stripe_mode,
        cloud_voxel=bool(args.cloud_voxel),
        cloud_voxel_cutout_count=(
            stripe_count if cloud_stripe_mode == "cutout" else 16
        ),
        cloud_missing_tint_opacity=float(args.cloud_missing_tint_opacity),
        star_render_expected_width=int(args.expected_render_width),
        content_fov_deg=float(args.content_fov_deg),
    )
    return catalogs, viewer, user_options, runtime


def _prepare_terrain_payload(
    *, root: Path, viewer: ViewerData, source_paths: set[Path]
) -> TerrainHorizonPayload:
    bbox = build_download_bbox(
        lat_deg=float(viewer.lat_deg), lon_deg=float(viewer.lon_deg), radius_km=138.0
    )
    dem_root = root / "copernicus-dem"
    paths = tuple(
        dem_root / key
        for key in collect_copernicus_tile_keys(bbox)
        if (dem_root / key).is_file()
    )
    distance_samples = build_distance_samples(
        128.0, DEFAULT_EXPORT_IMAGE_TERRAIN_SAMPLE_STEP_M
    )
    if paths:
        dem = GeoTiffDem(paths, default_elevation_m=0.0)
        try:
            dem_grid = dem.build_grid(bbox)
            ground_m = sample_ground_elevation(
                dem_grid,
                latitude_deg=float(viewer.lat_deg),
                longitude_deg=float(viewer.lon_deg),
                dem_resampling="bilinear",
            )
            layers = compute_horizon_layers(
                dem_grid=dem_grid,
                geod=WGS84_GEOD,
                observer=ObserverLocation(
                    latitude_deg=float(viewer.lat_deg),
                    longitude_deg=float(viewer.lon_deg),
                    observer_ground_m=ground_m,
                    observer_eye_m=float(viewer.observer_height_m),
                ),
                azimuth_step_deg=DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG,
                distance_samples_m=distance_samples,
                dem_resampling="bilinear",
                earth_radius_m=EARTH_MEAN_RADIUS_M,
                refraction_coefficient=0.13,
            )
            source_paths.update(paths)
        finally:
            dem.close()
    else:
        layers = compute_flat_ground_horizon_layers(
            geod=WGS84_GEOD,
            observer=ObserverLocation(
                latitude_deg=float(viewer.lat_deg),
                longitude_deg=float(viewer.lon_deg),
                observer_ground_m=0.0,
                observer_eye_m=float(viewer.observer_height_m),
            ),
            azimuth_step_deg=DEFAULT_EXPORT_IMAGE_TERRAIN_AZIMUTH_STEP_DEG,
            distance_samples_m=distance_samples,
            earth_radius_m=EARTH_MEAN_RADIUS_M,
            refraction_coefficient=0.13,
        )
    sample_distances = layers.sample_distances_m
    sample_elevations = layers.sample_terrain_elevation_m
    if sample_distances is None or sample_elevations is None:
        raise ValueError("Terrain horizon source arrays are unavailable")
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
        "sample_distances_m": sample_distances,
        "sample_terrain_elevation_m": sample_elevations,
    }


def _prepare_urban_outlines(
    *, root: Path, viewer: ViewerData, args: Any, source_paths: set[Path]
) -> list[Any] | None:
    radius = float(args.urban_outline_radius_km)
    plateau_root = root / "plateau_buildings"
    selection = select_prepared_building_source(
        observer_lat_deg=float(viewer.lat_deg),
        observer_lon_deg=float(viewer.lon_deg),
        radius_km=radius,
        plateau_root_dir=plateau_root,
    )
    if selection.source == "plateau":
        derived_dirs = selection.derived_dirs
    else:
        overture_root = root / "overture_buildings"
        derived_dirs = tuple(sorted(overture_root.glob("*/bldg")))
    buildings: list[BuildingFootprint] = []
    for derived_dir in derived_dirs:
        try:
            envelopes = select_derived_tile_envelopes(
                derived_dir,
                observer_lat_deg=float(viewer.lat_deg),
                observer_lon_deg=float(viewer.lon_deg),
                radius_km=radius,
            )
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        index_path = derived_dir / "tile_index.json"
        if index_path.is_file():
            source_paths.add(index_path)
        for envelope in envelopes:
            source_paths.add(envelope.path)
            buildings.extend(
                parse_derived_tile_buildings(
                    envelope.path,
                    min_building_height_m=float(args.urban_outline_min_height_m),
                )
            )
    if not buildings:
        return None
    merged = _merge_building_footprints(tuple(buildings))
    tower = Viewpoint(
        id="replay",
        qid="",
        kind="tower",
        name=viewer.city_name,
        labels={},
        names=(),
        latitude_deg=float(viewer.lat_deg),
        longitude_deg=float(viewer.lon_deg),
        height_m=0.0,
        viewpoint_height_m=float(viewer.observer_height_m),
        meta={},
    )
    result = compute_urban_outlines(
        tower,
        merged,
        radius_km=radius,
        min_distance_km=0.0,
        observer_ground_elevation_m=float(viewer.ground_elevation_m),
        view_center=None,
        edge_fov_deg=float(viewer.edge_fov_deg),
        edge_sample_step_m=10.0,
        max_candidates=int(args.urban_outline_max_candidates),
    )
    return _project_observer_centric_urban_outline_result(
        result,
        view_center=(float(viewer.view_center[0]), float(viewer.view_center[1])),
        edge_fov_deg=float(viewer.edge_fov_deg),
        front_hemisphere_view_center=(
            float(viewer.view_center[0]),
            float(viewer.view_center[1]),
        ),
        front_hemisphere_fov_deg=float(viewer.content_fov_deg),
    )


def _cached_night_lights_ready(
    *, root: Path, viewer: ViewerData, source_paths: set[Path]
) -> bool:
    cache_root = root / "night_lights"
    manifest_path = night_light_manifest_path(cache_root)
    if not manifest_path.is_file():
        return False
    try:
        manifest = validate_night_light_manifest(
            json.loads(manifest_path.read_text(encoding="utf-8"))
        )
        tiles = manifest["tiles"]
        if not isinstance(tiles, dict):
            return False
        required = _required_tile_names(
            observer_lat_deg=float(viewer.lat_deg),
            observer_lon_deg=float(viewer.lon_deg),
            max_distance_km=NIGHT_LIGHTS_MAX_DISTANCE_KM,
        )
        for name in required:
            entry = tiles.get(name)
            if not isinstance(entry, dict):
                return False
            path = night_light_tile_path(name, manifest, cache_root)
            if (
                not path.is_file()
                or _sha256_file(path).lower() != str(entry["sha256"]).lower()
            ):
                return False
            validate_night_light_geotiff(path, name, entry)
            source_paths.add(path)
        source_paths.add(manifest_path)
        return True
    except (OSError, ValueError, KeyError, TypeError, RuntimeError):
        return False


def _draw_selected_hud(
    image: Any, args: Any, viewer: ViewerData, target: datetime
) -> None:
    if not args.hud:
        return
    try:
        local_time = target.astimezone(ZoneInfo(viewer.timezone_name))
    except Exception:
        local_time = target
    values = {
        "datetime": local_time.strftime("%Y-%m-%d %H:%M %Z"),
        "location": viewer.city_name,
        "coordinates": f"{viewer.lat_deg:.5f}, {viewer.lon_deg:.5f}",
        "height": f"Observer height: {viewer.observer_height_m:.1f} m",
        "view-direction": f"Alt {viewer.view_alt_deg:.1f} deg / Az {viewer.view_az_deg:.1f} deg",
    }
    hud_items = tuple(item for item in args.hud if item != "time-of-day")
    if not hud_items:
        return
    painter = QPainter(image)
    try:
        painter.setPen(QColor(255, 255, 255, 230))
        font = painter.font()
        font.setPointSize(10)
        painter.setFont(font)
        x = 16
        y = (
            24
            if args.observation_info == "top"
            else image.height() - 16 * len(hud_items)
        )
        for key in ("datetime", "location", "coordinates", "height", "view-direction"):
            if key in hud_items:
                painter.drawText(x, y, values[key])
                y += 18
    finally:
        painter.end()


def export_replay_sequence(args: Any) -> int:
    from ..replay.location import resolve_location

    root = Path(args.cache_dir or CACHE_PATH).expanduser().resolve()
    output = Path(args.output_dir).expanduser().resolve() if args.output_dir else None
    resume_manifest: dict[str, Any] | None = None
    if args.resume:
        if output is None or not (output / "manifest.json").is_file():
            raise ValueError("--resume requires an existing replay manifest")
        loaded_manifest = json.loads(
            (output / "manifest.json").read_text(encoding="utf-8")
        )
        if not isinstance(loaded_manifest, dict):
            raise ValueError("Existing replay manifest is invalid")
        resume_manifest = loaded_manifest
        location_data = loaded_manifest.get("resolved_location")
        if not isinstance(location_data, dict):
            raise ValueError("Existing replay manifest has no resolved location")
        location = ResolvedLocation(**location_data)
        location_issues = list(loaded_manifest.get("location_issues", []))
    else:
        location_source_paths: set[Path] = set()
        location, location_issues = resolve_location(
            args,
            root,
            inspect_local_sources=True,
            source_paths=location_source_paths,
        )
    if args.resume:
        location_source_paths = set()
    ephemeris_path = root / EPHEMERIS_FILENAME
    if not ephemeris_path.is_file():
        raise ValueError(f"Ephemeris file is missing: {ephemeris_path}")
    for asset_path, label in (
        (Path(STARS_CSV_FILE), "Star catalog"),
        (Path(TEXT_FONT_PATH), "Overlay font"),
    ):
        if not asset_path.is_file():
            raise ValueError(f"{label} is missing: {asset_path}")
    height = (
        None
        if any(issue.startswith("building-") for issue in location_issues)
        else location.observer_height_m
    )
    inventory = Inventory(root, location.lat, location.lon, height)
    inventory.scan()
    times = _frame_times(args.start, args.end)
    selected = _selected_layers(args, location)
    plans, decoded, provenance, frozen_bytes = _read_verified_sources(
        inventory, selected, times
    )
    static_source_paths = set(location_source_paths)
    viewer = _build_viewer_data(args, location)
    terrain_payload = None
    if (
        float(args.terrain_horizon_opacity) > 0.0
        or float(args.night_light_opacity) > 0.0
        or float(args.ridge_glow_opacity) > 0.0
    ):
        terrain_payload = _prepare_terrain_payload(
            root=root, viewer=viewer, source_paths=static_source_paths
        )
    urban_outlines = None
    static_unsupported: list[str] = []
    if float(args.urban_outline_opacity) > 0.0:
        urban_outlines = _prepare_urban_outlines(
            root=root, viewer=viewer, args=args, source_paths=static_source_paths
        )
        if urban_outlines is None:
            static_unsupported.append("urban-outline-cache-unavailable")
    night_lights_ready = True
    if float(args.night_light_opacity) > 0.0:
        night_lights_ready = _cached_night_lights_ready(
            root=root, viewer=viewer, source_paths=static_source_paths
        )
        if not night_lights_ready:
            static_unsupported.append("night-lights-cache-unavailable")
    static_paths = [ephemeris_path, Path(STARS_CSV_FILE), Path(TEXT_FONT_PATH)]
    if Path(DSO_CSV_FILE).is_file():
        static_paths.append(Path(DSO_CSV_FILE))
    static_paths.extend(sorted(static_source_paths))
    for asset_path in static_paths:
        if not asset_path.is_file():
            raise ValueError(f"Could not read required replay input: {asset_path}")
        try:
            source_digest = _file_digest(asset_path)
        except OSError as exc:
            raise ValueError(
                f"Could not read required replay input: {asset_path}"
            ) from exc
        provenance_key = (
            f"static:{asset_path.relative_to(root)}"
            if asset_path.is_relative_to(root)
            else f"static:{asset_path.name}"
        )
        provenance[provenance_key] = {
            "path": str(asset_path),
            "sha256": source_digest,
        }
        # Very large rasters are fingerprinted in place instead of copied into
        # each export job. Resume still verifies their content digest.
        if asset_path.stat().st_size <= 16 * 1024 * 1024:
            raw = asset_path.read_bytes()
            frozen_bytes[frozen_bytes_key(asset_path)] = raw
    if resume_manifest is not None:
        old_inputs = resume_manifest.get("inputs", {})
        if isinstance(old_inputs, dict):
            for key, value in old_inputs.items():
                if str(key).startswith("static:") and isinstance(value, dict):
                    provenance.setdefault(str(key), value)
    unsupported = _unsupported_layers(args) + static_unsupported
    missing: list[str] = []
    for index, frame in enumerate(plans):
        for key, record in frame.items():
            if record is None:
                missing.append(
                    f"{iso(times[index])}: {key} has no matching cache record"
                )
            elif isinstance(decoded.get(_record_key(record)), Exception):
                missing.append(
                    f"{iso(times[index])}: {key} cache record cannot be rendered"
                )
    if location_issues:
        unsupported.extend(
            issue for issue in location_issues if issue not in unsupported
        )
    if (missing or unsupported) and not args.allow_partial_data and not args.plan_only:
        details = missing[:8] + [
            f"Unsupported or unverified input: {item}" for item in unsupported
        ]
        raise ValueError("Replay preflight failed: " + "; ".join(details))
    if args.plan_only:
        result = {
            "schema_version": 1,
            "start_utc": iso(args.start),
            "end_utc_exclusive": iso(args.end),
            "frame_count": len(times),
            "location": {
                "name": location.display_name,
                "latitude_deg": location.lat,
                "longitude_deg": location.lon,
            },
            "missing": missing,
            "unsupported": unsupported,
            "image_export_ready": not missing and not unsupported,
        }
        report = (
            json.dumps(result, ensure_ascii=True, indent=2)
            if args.json
            else _format_plan(result)
        )
        print(report.encode("ascii", "backslashreplace").decode("ascii"))
        return 0 if not missing and not unsupported else 2

    if output is None:
        raise ValueError("Output directory is required for image export")
    if output == root or root in output.parents or output in root.parents:
        raise ValueError("Output directory and cache root must not overlap")
    output_exists = output.exists()
    if output_exists and not output.is_dir():
        raise ValueError("Output path is not a directory")
    manifest_path = output / "manifest.json"
    operational_args = {
        "city",
        "place",
        "place_countrycode",
        "place_lang",
        "cache_dir",
        "start",
        "end",
        "output_dir",
        "plan_only",
        "allow_partial_data",
        "resume",
        "list_cache_ranges",
        "json",
        "time_reference_utc",
        "jobs",
    }
    settings = {
        "location": {
            "name": location.display_name,
            "latitude_deg": location.lat,
            "longitude_deg": location.lon,
            "timezone": location.tz,
            "observer_height_m": location.observer_height_m,
        },
        "start_utc": iso(args.start),
        "end_utc": iso(args.end),
        "render_options": {
            name: value
            for name, value in vars(args).items()
            if name not in operational_args and not name.startswith("_")
        },
    }
    settings_digest = _digest(_json_bytes(settings))
    plan_digest = _digest(_json_bytes(provenance))
    if output_exists and any(output.iterdir()):
        if not args.resume or not manifest_path.is_file():
            raise ValueError(
                "Output directory is non-empty; pass --resume for a matching job"
            )
        old_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            old_manifest.get("settings_sha256") != settings_digest
            or old_manifest.get("plan_sha256") != plan_digest
        ):
            raise ValueError(
                "Existing replay job does not match the requested settings or inputs"
            )
    else:
        output.mkdir(parents=True, exist_ok=True)
        input_dir = output / "inputs"
        for source_name, raw in frozen_bytes.items():
            source_digest = _digest(raw)
            destination = input_dir / source_digest
            if not destination.exists():
                _atomic_write(destination, raw)
        manifest = {
            "schema_version": 1,
            "application_version": __version__,
            "renderer_version": __version__,
            "resolved_location": location_payload(location),
            "location_issues": location_issues,
            "settings": settings,
            "settings_sha256": settings_digest,
            "plan_sha256": plan_digest,
            "time_reference_utc": iso(args.time_reference_utc),
            "start_utc": iso(args.start),
            "end_utc_exclusive": iso(args.end),
            "frame_count": len(times),
            "inputs": provenance,
            "missing": missing,
            "unsupported": unsupported,
            "frames": {},
        }
        _atomic_write(manifest_path, _json_bytes(manifest))

    if not Path(TEXT_FONT_PATH).is_file():
        raise ValueError(f"Overlay font is missing: {TEXT_FONT_PATH}")
    if not os.environ.get("DISPLAY") and not os.environ.get("WAYLAND_DISPLAY"):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QApplication.instance() or QApplication([])
    app.setQuitOnLastWindowClosed(False)
    catalogs, viewer, user_options, runtime = _build_render_inputs(
        args, location, viewer
    )
    text_font, status_font = _load_fonts(user_options.overlay_font_size)
    theme = THEME_STYLES_BY_PRESET[args.theme]
    style = _build_render_style(
        text_font=text_font,
        status_line_font=status_font,
        catalogs=catalogs,
        user_options=user_options,
        runtime_options=runtime,
        theme=theme,
    )
    compositor = _build_compositor(runtime)
    # Load only the ephemeris inside the selected cache root. The global
    # application loader may create or populate the default cache.
    from skyfield.api import Loader

    ephemeris = Loader(str(root))(EPHEMERIS_FILENAME)
    geometry = render_geometry.get_screen_geometry(
        int(args.image_size[0]),
        int(args.image_size[1]),
        viewer.view_alt_deg,
        edge_fov_deg=viewer.edge_fov_deg,
        content_fov_deg=viewer.content_fov_deg,
    )
    ridge_overlay_image = None
    ridge_layers = (
        terrain_payload["secondary_ridges_altaz_layers"]
        if terrain_payload is not None
        else None
    )
    ridge_distances = (
        terrain_payload["secondary_ridges_distances_m_layers"]
        if terrain_payload is not None
        else None
    )
    if ridge_layers and float(style.terrain_horizon_opacity) > 0.0:
        ridge_overlay_image = QImage(
            int(args.image_size[0]),
            int(args.image_size[1]),
            QImage.Format.Format_ARGB32_Premultiplied,
        )
        ridge_overlay_image.fill(0)
        ridge_painter = QPainter(ridge_overlay_image)
        ridge_painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        try:
            render_terrain.draw_terrain_secondary_ridges(
                ridge_painter,
                geometry,
                viewer,
                ridge_layers,
                ridge_distances,
                opacity=max(0.0, float(style.terrain_horizon_opacity) * 0.72),
                line_width_scale=compute_star_render_upscale_factor(
                    geometry.radius * 2,
                    style.star_render_expected_width,
                ),
                layer_style=style.theme.overlays.terrain_horizon,
            )
        finally:
            ridge_painter.end()
    services = ApplicationServices(max_workers=1)
    failures = 0
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    render_indices: list[int] = []
    for index, target in enumerate(times):
        frame_name = target.strftime("%Y%m%dT%H%MZ")
        png_path = output / "frames" / f"{frame_name}.png"
        sidecar_path = output / "frames" / f"{frame_name}.json"
        if args.resume and _valid_completed_frame(png_path, sidecar_path, plan_digest):
            continue
        if args.resume and (png_path.exists() or sidecar_path.exists()):
            raise ValueError(
                f"Incomplete or corrupt frame artifacts require manual review: {frame_name}"
            )
        render_indices.append(index)

    # A chunk ends at every selected cloud-record change, with a ten-frame cap
    # to bound retained sky snapshots for unusually long acquisition holds.
    voxel_executor = ThreadPoolExecutor(
        max_workers=int(args.jobs), thread_name_prefix="replay-cloud"
    )
    PreparedFrame = tuple[int, dict[str, Any], int]
    PreparedChunk = tuple[
        list[PreparedFrame],
        Future[None] | None,
        tempfile.TemporaryDirectory[str] | None,
    ]
    prepared_chunks: deque[PreparedChunk] = deque()
    chunk_workdirs: list[tempfile.TemporaryDirectory[str]] = []
    render_chunks = _chunk_frame_indices(
        render_indices,
        plans,
        decoded,
        voxel_enabled=bool(runtime.cloud_voxel),
    )

    next_chunk_position = 0

    def prepare_chunk(indices: list[int]) -> PreparedChunk:
        skies: list[PreparedFrame] = []
        frame_inputs: list[dict[str, object]] = []
        cloud_record = plans[indices[0]].get("cloud")
        cloud_data = (
            decoded.get(_record_key(cloud_record)) if cloud_record is not None else None
        )
        for index in indices:
            frame_name = times[index].strftime("%Y%m%dT%H%MZ")
            print(f"Rendering {index + 1}/{len(times)}: {frame_name}", file=sys.stderr)
            with services.native_work_lock:
                sky = compute_sky_snapshot(
                    ephemeris=ephemeris,
                    viewer_data=viewer,
                    geometry=geometry,
                    star_catalog=catalogs.star_catalog_np,
                    dso_catalog=catalogs.dso_catalog_np,
                    star_vmag_limit=None if args.vmag_limit <= 6.0 else args.vmag_limit,
                    star_subset_indices=(
                        catalogs.star_catalog_lod6_indices
                        if args.vmag_limit <= 6.0
                        else None
                    ),
                    delta_t=timedelta(0),
                    sky_update_interval=DEFAULT_EXPORT_IMAGE_SKY_UPDATE_INTERVAL,
                    sky_disc_alpha=user_options.sky_disc_alpha,
                    theme=theme,
                    star_catalog_meta=catalogs.star_catalog_meta,
                    image_size=(int(args.image_size[0]), int(args.image_size[1])),
                    target_time_utc=times[index],
                )
            skies.append((index, sky, len(skies)))
            celestial = cast(CelestialData, sky["celestial"])
            sun_body = next(
                (body for body in celestial.planets if body.name == "sun"), None
            )
            if sun_body is not None:
                from ..night_lights import night_light_strength_factor

                frame_inputs.append(
                    {
                        "request_id": index,
                        "sun_alt_deg": float(sun_body.alt),
                        "sun_az_deg": float(sun_body.az),
                        "sunlight_mix": 1.0
                        - night_light_strength_factor(float(sun_body.alt)),
                    }
                )
        voxel_future = None
        work_manager = None
        if runtime.cloud_voxel and isinstance(cloud_data, CloudSourceData):
            from ..cloud_voxel_style import CLOUD_VOXEL_NIGHT_COLOR_RGB
            from ..clouddisc.workers.cloud_voxel_worker import (
                render_cloud_voxel_chunk_in_subprocess,
            )

            if len(frame_inputs) != len(indices):
                raise ValueError("Could not determine the Sun position for every frame")
            work_manager = tempfile.TemporaryDirectory(
                prefix="zstarview-replay-cloud-chunk-"
            )
            chunk_workdirs.append(work_manager)
            voxel_future = voxel_executor.submit(
                render_cloud_voxel_chunk_in_subprocess,
                cloud_data,
                lat=viewer.lat_deg,
                lon=viewer.lon_deg,
                view_center=viewer.view_center,
                edge_fov_deg=viewer.edge_fov_deg,
                content_fov_deg=viewer.content_fov_deg,
                radius_px=max(1, int(geometry.radius)),
                frames=frame_inputs,
                night_color_rgb=CLOUD_VOXEL_NIGHT_COLOR_RGB,
                work_dir=Path(work_manager.name),
                timeout_s=120.0 * len(frame_inputs),
                native_work_lock=services.native_work_lock,
            )
        return skies, voxel_future, work_manager

    def prepare_available_chunks() -> None:
        nonlocal next_chunk_position
        while next_chunk_position < len(render_chunks) and len(prepared_chunks) < int(
            args.jobs
        ):
            indices = render_chunks[next_chunk_position]
            next_chunk_position += 1
            prepared_chunks.append(prepare_chunk(indices))

    try:
        active_frames: deque[PreparedFrame] = deque()
        active_future: Future[None] | None = None
        active_work_manager: tempfile.TemporaryDirectory[str] | None = None
        for index in render_indices:
            target = times[index]
            frame_name = target.strftime("%Y%m%dT%H%MZ")
            if not active_frames:
                if active_future is not None:
                    active_future.result()
                if active_work_manager is not None:
                    active_work_manager.cleanup()
                    chunk_workdirs.remove(active_work_manager)
                active_future = None
                active_work_manager = None
                prepare_available_chunks()
                if not prepared_chunks:
                    raise RuntimeError("Replay frame chunk was not prepared")
                frames, voxel_future, work_manager = prepared_chunks.popleft()
                active_frames = deque(frames)
                active_future = voxel_future
                active_work_manager = work_manager
            prepared_index, sky, chunk_offset = active_frames.popleft()
            if prepared_index != index:
                raise RuntimeError("Replay frame preparation order changed")
            png_path = output / "frames" / f"{frame_name}.png"
            sidecar_path = output / "frames" / f"{frame_name}.json"
            frame = plans[index]
            frame_records: dict[str, Any] = {}
            frame_status: dict[str, str] = {}
            for reason in unsupported:
                frame_status[f"unsupported:{reason}"] = "unsupported"
            if (
                terrain_payload is not None
                and float(args.terrain_horizon_opacity) > 0.0
            ):
                frame_status["terrain-horizon"] = "ok"
            if urban_outlines and float(args.urban_outline_opacity) > 0.0:
                frame_status["urban-outline"] = "ok"
            cloud_grid = None
            aircraft = None
            satellites: dict[str, list[dict[str, Any]]] = {}
            for key, record in frame.items():
                if record is None or isinstance(
                    decoded.get(_record_key(record)), Exception
                ):
                    frame_status[key] = "missing"
                    continue
                value = decoded[_record_key(record)]
                frame_records[key] = value
                frame_status[key] = "ok"
                if key == "cloud" and isinstance(value, CloudSourceData):
                    cloud_grid = value.altaz_grid
                elif key == "aircraft":
                    aircraft = value.snapshots
                    if not aircraft:
                        frame_status[key] = "empty"
                elif key.startswith("satellite-"):
                    satellites[key.removeprefix("satellite-")] = value.records
            target_viewer = viewer
            celestial = cast(CelestialData, sky["celestial"])
            sun_body = next(
                (body for body in celestial.planets if body.name == "sun"), None
            )
            night_light_profile = None
            night_light_active = bool(
                night_lights_ready
                and terrain_payload is not None
                and sun_body is not None
                and is_night_light_enabled(float(sun_body.alt))
                and (
                    float(args.night_light_opacity) > 0.0
                    or float(args.ridge_glow_opacity) > 0.0
                )
            )
            if night_light_active:
                assert sun_body is not None and terrain_payload is not None
                night_light_profile = compute_night_light_glow_profile(
                    observer_lat_deg=float(viewer.lat_deg),
                    observer_lon_deg=float(viewer.lon_deg),
                    observer_elevation_m=max(
                        0.0,
                        float(viewer.ground_elevation_m)
                        + float(viewer.observer_height_m),
                    ),
                    sun_alt_deg=float(sun_body.alt),
                    terrain_profile_altaz=terrain_payload["profile_altaz"],
                    terrain_profile_distances_m=terrain_payload["profile_distances_m"],
                    terrain_secondary_ridges_altaz_layers=terrain_payload[
                        "secondary_ridges_altaz_layers"
                    ],
                    terrain_secondary_ridges_distances_m_layers=terrain_payload[
                        "secondary_ridges_distances_m_layers"
                    ],
                    terrain_sample_distances_m=terrain_payload["sample_distances_m"],
                    terrain_sample_terrain_elevation_m=terrain_payload[
                        "sample_terrain_elevation_m"
                    ],
                    include_night_light_tiles=True,
                    cache_root=root / "night_lights",
                )
            if float(args.night_light_opacity) > 0.0:
                frame_status["night-lights"] = (
                    "ok"
                    if night_light_profile is not None
                    else "empty"
                    if night_light_active
                    else "inactive-daylight"
                )
            if float(args.ridge_glow_opacity) > 0.0:
                frame_status["ridge-glow"] = (
                    "ok"
                    if night_light_profile is not None
                    else "empty"
                    if night_light_active
                    else "inactive-daylight"
                )
            voxel_image = None
            if active_future is not None:
                assert active_work_manager is not None
                cloud_frame_path = (
                    Path(active_work_manager.name)
                    / "frames"
                    / f"frame-{chunk_offset:04d}.npy"
                )
                while not cloud_frame_path.is_file():
                    try:
                        active_future.result(timeout=0.1)
                    except FutureTimeout:
                        continue
                    raise RuntimeError("Cloud voxel chunk completed without its frame")
                with services.native_work_lock:
                    with cloud_frame_path.open("rb") as cloud_stream:
                        voxel_image = np.load(cloud_stream, allow_pickle=False)
                expected_cloud_shape = (
                    max(1, int(geometry.radius)) * 2 + 1,
                    max(1, int(geometry.radius)) * 2 + 1,
                    4,
                )
                if (
                    voxel_image.shape != expected_cloud_shape
                    or voxel_image.dtype != np.uint8
                ):
                    raise RuntimeError(
                        "Cloud voxel chunk returned an invalid frame image"
                    )
                from ..cloud_voxel_style import apply_cloud_voxel_cutout

                voxel_image = apply_cloud_voxel_cutout(
                    voxel_image, count=runtime.cloud_voxel_cutout_count
                )
            scene = RenderSceneData(
                celestial_data=celestial,
                sky_disc_image=cast(QImage | None, sky.get("sky_disc")),
                cloud_missing_mask=None,
                cloud_altaz_grid=None if voxel_image is not None else cloud_grid,
                terrain_horizon_profile=(
                    terrain_payload["profile_altaz"]
                    if terrain_payload is not None
                    else None
                ),
                terrain_horizon_profile_distances_m=(
                    terrain_payload["profile_distances_m"]
                    if terrain_payload is not None
                    else None
                ),
                terrain_secondary_ridges_altaz_layers=(
                    terrain_payload["secondary_ridges_altaz_layers"]
                    if terrain_payload is not None
                    else None
                ),
                terrain_secondary_ridges_distances_m_layers=(
                    terrain_payload["secondary_ridges_distances_m_layers"]
                    if terrain_payload is not None
                    else None
                ),
                urban_outlines=urban_outlines,
                terrain_secondary_ridges_image=ridge_overlay_image,
                cloud_voxel_image=voxel_image,
                satellite_records_by_group=satellites,
                aircraft_snapshots=aircraft,
                night_light_glow_profile=night_light_profile,
            )
            frame_style = style
            if "cloud" not in frame_records:
                frame_style = replace(frame_style, cloud_disc_alpha=0.0)
            if "aircraft" not in frame_records:
                frame_style = replace(frame_style, aircraft_opacity=0.0)
            if not satellites:
                frame_style = replace(frame_style, satellite_opacity=0.0)
            image = _render_image(
                image_size=(int(args.image_size[0]), int(args.image_size[1])),
                viewer=target_viewer,
                scene=scene,
                style=frame_style,
                compositor=compositor,
                draw_direction_grid=bool(args.include_direction_grid),
                time_of_day_marker_bottom_left=(
                    args.observation_info == "top"
                    if "time-of-day" in args.hud
                    else None
                ),
            )
            _draw_selected_hud(image, args, target_viewer, target)
            png_path.parent.mkdir(parents=True, exist_ok=True)
            temporary_png = png_path.with_name(f".{png_path.name}.tmp")
            if not image.save(str(temporary_png), "PNG"):
                failures += 1
                frame_status["render"] = "failed"
                _atomic_write(
                    sidecar_path,
                    _json_bytes({"time_utc": iso(target), "status": frame_status}),
                )
                manifest.setdefault("frames", {})[frame_name] = {
                    "png": str(png_path.relative_to(output)),
                    "sidecar": str(sidecar_path.relative_to(output)),
                    "status": {**frame_status, "render": "failed"},
                }
                _atomic_write(manifest_path, _json_bytes(manifest))
                continue
            os.replace(temporary_png, png_path)
            image_hash = _digest(png_path.read_bytes())
            sidecar = {
                "schema_version": 1,
                "time_utc": iso(target),
                "plan_sha256": plan_digest,
                "image_sha256": image_hash,
                "layers": frame_status,
            }
            _atomic_write(sidecar_path, _json_bytes(sidecar))
            manifest.setdefault("frames", {})[frame_name] = {
                "png": str(png_path.relative_to(output)),
                "sidecar": str(sidecar_path.relative_to(output)),
                "image_sha256": image_hash,
                "status": frame_status,
            }
            _atomic_write(manifest_path, _json_bytes(manifest))
            print(
                f"Rendered {index + 1}/{len(times)}: {frame_name}",
                file=sys.stderr,
            )
        if active_future is not None:
            active_future.result()
        if active_work_manager is not None:
            active_work_manager.cleanup()
            chunk_workdirs.remove(active_work_manager)
    except Exception as exc:
        manifest["failure"] = {
            "type": type(exc).__name__,
            "message": str(exc).encode("ascii", "backslashreplace").decode("ascii"),
        }
        _atomic_write(manifest_path, _json_bytes(manifest))
        raise
    finally:
        voxel_executor.shutdown(wait=True, cancel_futures=True)
        for work_manager in chunk_workdirs:
            work_manager.cleanup()
        services.shutdown()
    return 1 if failures else (2 if missing or unsupported else 0)


def _valid_completed_frame(
    png_path: Path, sidecar_path: Path, plan_digest: str
) -> bool:
    try:
        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
        if sidecar.get("plan_sha256") != plan_digest:
            return False
        return png_path.is_file() and _digest(png_path.read_bytes()) == sidecar.get(
            "image_sha256"
        )
    except (OSError, ValueError, TypeError):
        return False


def _format_plan(plan: dict[str, Any]) -> str:
    return "\n".join(
        (
            f"Location: {plan['location']['name']}",
            f"Frames: {plan['frame_count']}",
            f"Range: [{plan['start_utc']}, {plan['end_utc_exclusive']})",
            f"Image export ready: {str(plan['image_export_ready']).lower()}",
            "Missing: " + ("; ".join(plan["missing"]) or "none"),
            "Unsupported: " + ("; ".join(plan["unsupported"]) or "none"),
        )
    )
