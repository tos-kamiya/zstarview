from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace
import math

import numpy as np

from .data.derived_tile_cache import parse_derived_tile_buildings, select_derived_tile_envelopes
from .data.import_overture_buildings import DEFAULT_FETCH_RADIUS_KM
from .data.urban_outline_common import BuildingFootprint
from .data.urban_outline_from_buildings import compute_urban_outlines
from .paths import COPERNICUS_DEM_CACHE_DIR, OVERTURE_DERIVED_ROOT_DIR
from .terrain import GeoTiffDem, build_download_bbox, fetch_copernicus_dem
from .types import UrbanOutlinePolyline, ViewerData

DEM_DOWNLOAD_MARGIN_KM = 5.0


def resolve_urban_outline_layer_for_viewer(
    viewer_data: ViewerData,
    *,
    derived_root_dir: str | Path = OVERTURE_DERIVED_ROOT_DIR,
    dem_cache_dir: str | Path = COPERNICUS_DEM_CACHE_DIR,
    derived_dir: str | Path | None = None,
    derived_dirs: tuple[str | Path, ...] | None = None,
    radius_km: float = DEFAULT_FETCH_RADIUS_KM,
    min_distance_km: float = 0.0,
    min_height_m: float = 0.0,
) -> list[UrbanOutlinePolyline] | None:
    return _build_dynamic_urban_outline_layer(
        lat_deg=float(viewer_data.lat_deg),
        lon_deg=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        edge_fov_deg=float(viewer_data.edge_fov_deg),
        derived_root_dir=Path(derived_root_dir),
        dem_cache_dir=Path(dem_cache_dir),
        derived_dir=None if derived_dir is None else Path(derived_dir),
        derived_dirs=None if derived_dirs is None else tuple(Path(path) for path in derived_dirs),
        radius_km=float(radius_km),
        min_distance_km=float(min_distance_km),
        min_height_m=float(min_height_m),
        view_center=tuple(float(v) for v in viewer_data.view_center),
    )


@lru_cache(maxsize=64)
def _build_dynamic_urban_outline_layer(
    *,
    lat_deg: float,
    lon_deg: float,
    observer_height_m: float,
    edge_fov_deg: float,
    derived_root_dir: Path,
    dem_cache_dir: Path,
    derived_dir: Path | None = None,
    derived_dirs: tuple[Path, ...] | None = None,
    radius_km: float = DEFAULT_FETCH_RADIUS_KM,
    min_distance_km: float = 0.0,
    min_height_m: float = 0.0,
    view_center: tuple[float, float] | None = None,
) -> list[UrbanOutlinePolyline] | None:
    if derived_dirs is not None:
        candidate_dirs = tuple(path for path in derived_dirs if path.exists())
    elif derived_dir is not None:
        candidate_dirs = (derived_dir,) if derived_dir.exists() else ()
    elif derived_root_dir.exists():
        candidate_dirs = _list_derived_dirs(derived_root_dir)
    else:
        candidate_dirs = ()
    if not candidate_dirs:
        return None

    buildings = []
    for derived_dir in candidate_dirs:
        try:
            envelopes = select_derived_tile_envelopes(
                derived_dir,
                observer_lat_deg=lat_deg,
                observer_lon_deg=lon_deg,
                radius_km=radius_km,
            )
        except ValueError:
            continue
        for envelope in envelopes:
            buildings.extend(parse_derived_tile_buildings(envelope.path))
    buildings = _merge_building_footprints(tuple(buildings))
    if min_height_m > 0.0:
        buildings = tuple(
            building for building in buildings if float(building.height_m) >= float(min_height_m)
        )
    if not buildings:
        return None
    observer_ground_elevation_m, buildings = _resolve_building_ground_elevations(
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        buildings=tuple(buildings),
        radius_km=radius_km,
        dem_cache_dir=dem_cache_dir,
    )

    result = compute_urban_outlines(
        SimpleNamespace(
            id="coords",
            name="coords",
            latitude_deg=lat_deg,
            longitude_deg=lon_deg,
            viewpoint_height_m=observer_height_m,
            observer_height_m=observer_height_m,
        ),
        buildings,
        radius_km=radius_km,
        min_distance_km=min_distance_km,
        observer_ground_elevation_m=observer_ground_elevation_m,
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
        edge_sample_step_m=10.0,
    )
    outlines = [
        UrbanOutlinePolyline(
            points=[(point.altitude_deg, point.azimuth_deg) for point in outline.points],
            height_m=float(outline.height_m),
            distance_km=float(outline.distance_km),
            source="base",
        )
        for outline in result.outlines
    ]
    return outlines or None


def _merge_building_footprints(buildings: tuple[BuildingFootprint, ...]) -> tuple[BuildingFootprint, ...]:
    parent_ids_with_parts = {
        building.parent_building_id
        for building in buildings
        if building.parent_building_id
    }
    if not parent_ids_with_parts:
        return buildings
    return tuple(
        building
        for building in buildings
        if building.parent_building_id is not None or building.building_id not in parent_ids_with_parts
    )


@lru_cache(maxsize=8)
def _list_derived_dirs(derived_root_dir: Path) -> tuple[Path, ...]:
    return tuple(
        path
        for path in sorted(derived_root_dir.glob("*/bldg"))
        if path.is_dir()
    )


def _resolve_building_ground_elevations(
    *,
    lat_deg: float,
    lon_deg: float,
    buildings: tuple[BuildingFootprint, ...],
    radius_km: float,
    dem_cache_dir: Path,
) -> tuple[float, tuple[BuildingFootprint, ...]]:
    if not buildings:
        return 0.0, ()

    try:
        download = fetch_copernicus_dem(
            observer_lat_deg=lat_deg,
            observer_lon_deg=lon_deg,
            max_distance_km=radius_km,
            margin_km=DEM_DOWNLOAD_MARGIN_KM,
            cache_dir=dem_cache_dir,
        )
    except RuntimeError as exc:
        if str(exc) == "No Copernicus DEM tiles were downloaded for the requested area.":
            return 0.0, buildings
        raise

    dem = GeoTiffDem(download.paths, default_elevation_m=0.0)
    try:
        bbox = build_download_bbox(
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            radius_km=radius_km + DEM_DOWNLOAD_MARGIN_KM,
        )
        dem_grid = dem.build_grid(bbox)
        observer_ground_elevation_m = float(
            dem_grid.sample_lonlat(
                np.array([lon_deg], dtype=np.float64),
                np.array([lat_deg], dtype=np.float64),
                method="bilinear",
            )[0]
        )
        building_lonlat = tuple(_representative_point_lonlat(building) for building in buildings)
        building_lon = np.array([point[0] for point in building_lonlat], dtype=np.float64)
        building_lat = np.array([point[1] for point in building_lonlat], dtype=np.float64)
        sampled_ground = dem_grid.sample_lonlat(building_lon, building_lat, method="bilinear")
        resolved_buildings = tuple(
            BuildingFootprint(
                building_id=building.building_id,
                height_m=building.height_m,
                rings_lonlat=building.rings_lonlat,
                parent_building_id=building.parent_building_id,
                ground_elevation_m=(
                    float(sampled_ground[index])
                    if math.isfinite(float(sampled_ground[index]))
                    else 0.0
                ),
                min_height_m=float(building.min_height_m),
            )
            for index, building in enumerate(buildings)
        )
        if not math.isfinite(observer_ground_elevation_m):
            observer_ground_elevation_m = 0.0
        return observer_ground_elevation_m, resolved_buildings
    finally:
        dem.close()


def _representative_point_lonlat(building: BuildingFootprint) -> tuple[float, float]:
    ring = _select_representative_ring(building.rings_lonlat)
    centroid = _polygon_centroid_lonlat(ring)
    if centroid is not None and _point_in_ring(centroid, ring):
        return centroid
    bbox_center = _bbox_center_lonlat(ring)
    if _point_in_ring(bbox_center, ring):
        return bbox_center
    return _mean_point_lonlat(ring)


def _select_representative_ring(
    rings_lonlat: tuple[tuple[tuple[float, float], ...], ...],
) -> tuple[tuple[float, float], ...]:
    best_ring = rings_lonlat[0]
    best_area = abs(_signed_area(best_ring))
    for ring in rings_lonlat[1:]:
        area = abs(_signed_area(ring))
        if area > best_area:
            best_ring = ring
            best_area = area
    return best_ring


def _polygon_centroid_lonlat(
    ring_lonlat: tuple[tuple[float, float], ...],
) -> tuple[float, float] | None:
    signed_area = _signed_area(ring_lonlat)
    if abs(signed_area) < 1e-12:
        return None
    factor = 1.0 / (6.0 * signed_area)
    cx = 0.0
    cy = 0.0
    for (x0, y0), (x1, y1) in zip(ring_lonlat[:-1], ring_lonlat[1:]):
        cross = x0 * y1 - x1 * y0
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    return (cx * factor, cy * factor)


def _signed_area(ring_lonlat: tuple[tuple[float, float], ...]) -> float:
    if len(ring_lonlat) < 3:
        return 0.0
    area = 0.0
    pairs = zip(ring_lonlat, ring_lonlat[1:]) if ring_lonlat[0] == ring_lonlat[-1] else zip(ring_lonlat, ring_lonlat[1:] + ring_lonlat[:1])
    for (x0, y0), (x1, y1) in pairs:
        area += x0 * y1 - x1 * y0
    return 0.5 * area


def _point_in_ring(point_lonlat: tuple[float, float], ring_lonlat: tuple[tuple[float, float], ...]) -> bool:
    px, py = point_lonlat
    if len(ring_lonlat) < 3:
        return False
    inside = False
    points = ring_lonlat if ring_lonlat[0] == ring_lonlat[-1] else ring_lonlat + (ring_lonlat[0],)
    for (x0, y0), (x1, y1) in zip(points[:-1], points[1:]):
        intersects = ((y0 > py) != (y1 > py)) and (px < (x1 - x0) * (py - y0) / ((y1 - y0) or 1e-12) + x0)
        if intersects:
            inside = not inside
    return inside


def _bbox_center_lonlat(ring_lonlat: tuple[tuple[float, float], ...]) -> tuple[float, float]:
    lon = [point[0] for point in ring_lonlat]
    lat = [point[1] for point in ring_lonlat]
    return ((min(lon) + max(lon)) * 0.5, (min(lat) + max(lat)) * 0.5)


def _mean_point_lonlat(ring_lonlat: tuple[tuple[float, float], ...]) -> tuple[float, float]:
    points = ring_lonlat[:-1] if len(ring_lonlat) > 1 and ring_lonlat[0] == ring_lonlat[-1] else ring_lonlat
    lon = sum(point[0] for point in points) / max(1, len(points))
    lat = sum(point[1] for point in points) / max(1, len(points))
    return (lon, lat)
