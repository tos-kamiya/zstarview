from __future__ import annotations

import json
import logging
import math
from functools import lru_cache
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np

from .paths import CACHE_PATH
from .terrain import GeoTiffDem, build_download_bbox, fetch_copernicus_dem
from .types import ViewerData
from .water_overlay import (
    WaterPolygonFootprint,
    WaterSurfacePatch,
    bbox_from_point,
    build_overpass_query,
    extract_water_polygons,
)
from .water_surface_mesh import WaterSurfaceMesh, build_water_surface_mesh

logger = logging.getLogger(__name__)

DEFAULT_FETCH_RADIUS_KM = 2.5
DEM_DOWNLOAD_MARGIN_KM = 5.0
OVERPASS_API_URL = "https://overpass-api.de/api/interpreter"
OVERPASS_TIMEOUT_SECONDS = 60.0


def resolve_water_overlay_for_viewer(
    viewer_data: ViewerData,
    *,
    radius_km: float = DEFAULT_FETCH_RADIUS_KM,
    dem_cache_dir: str | Path = Path(CACHE_PATH) / "copernicus-dem",
    overpass_url: str = OVERPASS_API_URL,
    timeout_seconds: float = OVERPASS_TIMEOUT_SECONDS,
) -> list[WaterSurfaceMesh] | None:
    return _build_dynamic_water_overlay_layer(
        lat_deg=float(viewer_data.lat_deg),
        lon_deg=float(viewer_data.lon_deg),
        radius_km=float(radius_km),
        dem_cache_dir=Path(dem_cache_dir),
        overpass_url=str(overpass_url),
        timeout_seconds=float(timeout_seconds),
    )


@lru_cache(maxsize=64)
def _build_dynamic_water_overlay_layer(
    *,
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    dem_cache_dir: Path,
    overpass_url: str,
    timeout_seconds: float,
) -> list[WaterSurfaceMesh] | None:
    bbox = bbox_from_point(lat_deg, lon_deg, radius_km)
    elements = _fetch_overpass_elements(
        bbox,
        overpass_url=overpass_url,
        timeout_seconds=timeout_seconds,
    )
    footprints = extract_water_polygons(elements)
    if not footprints:
        return None

    dem_grid = _resolve_dem_grid(
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        radius_km=radius_km,
        dem_cache_dir=dem_cache_dir,
    )
    meshes: list[WaterSurfaceMesh] = []
    for footprint in footprints:
        patch = _build_water_surface_patch(footprint, dem_grid=dem_grid)
        mesh = build_water_surface_mesh(
            footprint,
            center_lat_deg=lat_deg,
            center_lon_deg=lon_deg,
            patch=patch,
            grid_m=1.0,
            simplify_tolerance_m=20.0,
        )
        if mesh is not None:
            meshes.append(mesh)
    return meshes or None


def _fetch_overpass_elements(
    bbox: tuple[float, float, float, float],
    *,
    overpass_url: str,
    timeout_seconds: float,
) -> list[dict[str, object]]:
    query = build_overpass_query(bbox)
    request = Request(
        overpass_url,
        data=query.encode("utf-8"),
        headers={
            "User-Agent": "zstarview/1.0",
            "Content-Type": "application/x-www-form-urlencoded; charset=utf-8",
            "Accept": "application/json",
        },
    )
    with urlopen(request, timeout=float(timeout_seconds)) as response:  # nosec: B310
        payload = json.loads(response.read().decode("utf-8"))
    if not isinstance(payload, dict):
        return []
    elements = payload.get("elements", [])
    if not isinstance(elements, list):
        return []
    return [element for element in elements if isinstance(element, dict)]


def _resolve_dem_grid(
    *,
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    dem_cache_dir: Path,
) -> object | None:
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
            return None
        raise

    dem = GeoTiffDem(download.paths, default_elevation_m=0.0)
    try:
        bbox = build_download_bbox(
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            radius_km=radius_km + DEM_DOWNLOAD_MARGIN_KM,
        )
        return dem.build_grid(bbox)
    finally:
        dem.close()


def _build_water_surface_patch(
    footprint: WaterPolygonFootprint,
    *,
    dem_grid: object | None,
) -> WaterSurfacePatch | None:
    ring = _select_representative_ring(footprint.outer_rings_lonlat)
    anchor_points = _select_anchor_points_lonlat(ring)
    if anchor_points is None:
        return None
    anchor_elevations = _sample_anchor_elevations(anchor_points, dem_grid=dem_grid)
    if anchor_elevations is None:
        return None
    return WaterSurfacePatch(
        patch_id=footprint.water_id,
        polygon_id=footprint.water_id,
        anchor_points_lonlat=anchor_points,
        anchor_elevations_m=anchor_elevations,
        flat_threshold_m=1.0,
    )


def _select_representative_ring(
    rings_lonlat: tuple[tuple[tuple[float, float], ...], ...],
) -> tuple[tuple[float, float], ...]:
    best_ring = rings_lonlat[0]
    best_area = abs(_signed_area_lonlat(best_ring))
    for ring in rings_lonlat[1:]:
        area = abs(_signed_area_lonlat(ring))
        if area > best_area:
            best_ring = ring
            best_area = area
    return best_ring


def _select_anchor_points_lonlat(
    ring_lonlat: tuple[tuple[float, float], ...],
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None:
    body = ring_lonlat[:-1] if len(ring_lonlat) >= 2 and ring_lonlat[0] == ring_lonlat[-1] else ring_lonlat
    if len(body) < 3:
        return None
    step = max(1, len(body) // 3)
    anchor_0 = body[0]
    anchor_1 = body[step % len(body)]
    anchor_2 = body[(2 * step) % len(body)]
    if len({anchor_0, anchor_1, anchor_2}) < 3:
        anchor_1 = body[len(body) // 2]
        anchor_2 = body[-1]
    if len({anchor_0, anchor_1, anchor_2}) < 3:
        return None
    return anchor_0, anchor_1, anchor_2


def _sample_anchor_elevations(
    anchor_points_lonlat: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    *,
    dem_grid: object | None,
) -> tuple[float, float, float] | None:
    if dem_grid is None:
        return (0.0, 0.0, 0.0)
    lon = np.array([point[0] for point in anchor_points_lonlat], dtype=np.float64)
    lat = np.array([point[1] for point in anchor_points_lonlat], dtype=np.float64)
    try:
        sampled = dem_grid.sample_lonlat(lon, lat, method="bilinear")
    except Exception:
        return None
    values = tuple(
        float(value) if math.isfinite(float(value)) else 0.0
        for value in np.asarray(sampled, dtype=np.float64)
    )
    if len(values) != 3:
        return None
    return values


def _signed_area_lonlat(ring_lonlat: tuple[tuple[float, float], ...]) -> float:
    if len(ring_lonlat) < 3:
        return 0.0
    points = ring_lonlat
    if points[0] != points[-1]:
        points = points + (points[0],)
    area = 0.0
    for (x0, y0), (x1, y1) in zip(points[:-1], points[1:]):
        area += x0 * y1 - x1 * y0
    return 0.5 * area
