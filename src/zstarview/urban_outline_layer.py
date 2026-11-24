from __future__ import annotations

import math
from dataclasses import replace
from functools import lru_cache
from pathlib import Path

import numpy as np

from .astro import is_in_fov_vectorized
from .data.derived_tile_cache import (
    parse_derived_tile_buildings,
    select_derived_tile_envelopes,
)
from .data.import_overture_buildings import DEFAULT_FETCH_RADIUS_KM
from .data.urban_outline_common import (
    BuildingFootprint,
    make_local_transformer,
    project_ring_xy,
)
from .data.urban_outline_from_buildings import (
    UrbanOutlinePolyline as ObserverOutline,
)
from .data.urban_outline_from_buildings import (
    UrbanOutlineResult,
    UrbanPolylinePoint,
    _maybe_linearize_run_points,
    compute_urban_outlines,
    ring_contains_origin_xy,
)
from .location_resolver.viewpoints import Viewpoint
from .nearby_walls import (
    NEARBY_WALL_MAX_DISTANCE_M,
    NEARBY_WALL_MIN_DISTANCE_M,
    WallSegment,
    roof_floor_prisms,
    visible_wall_segments,
    wall_opacity_level,
)
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
    max_candidates: int = 5000,
    front_hemisphere_view_center: tuple[float, float] | None = None,
    front_hemisphere_fov_deg: float = 90.0,
    include_nearby_walls: bool = False,
    wall_radius_km: float = NEARBY_WALL_MAX_DISTANCE_M / 1000.0,
) -> list[UrbanOutlinePolyline] | None:
    observer_centric_result = _build_observer_centric_urban_outline_result(
        lat_deg=float(viewer_data.lat_deg),
        lon_deg=float(viewer_data.lon_deg),
        observer_height_m=float(viewer_data.observer_height_m),
        observer_ground_elevation_m=float(viewer_data.ground_elevation_m),
        derived_root_dir=Path(derived_root_dir),
        dem_cache_dir=Path(dem_cache_dir),
        derived_dir=None if derived_dir is None else Path(derived_dir),
        derived_dirs=None
        if derived_dirs is None
        else tuple(Path(path) for path in derived_dirs),
        radius_km=float(radius_km),
        min_distance_km=float(min_distance_km),
        min_height_m=float(min_height_m),
        max_candidates=max(0, int(max_candidates)),
        include_nearby_walls=include_nearby_walls,
        wall_radius_km=min(float(wall_radius_km), float(radius_km)),
    )
    if observer_centric_result is None:
        return None
    return _project_observer_centric_urban_outline_result(
        observer_centric_result,
        view_center=(
            float(viewer_data.view_center[0]),
            float(viewer_data.view_center[1]),
        ),
        edge_fov_deg=float(viewer_data.edge_fov_deg),
        front_hemisphere_view_center=front_hemisphere_view_center,
        front_hemisphere_fov_deg=front_hemisphere_fov_deg,
    )


@lru_cache(maxsize=64)
def _build_observer_centric_urban_outline_result(
    *,
    lat_deg: float,
    lon_deg: float,
    observer_height_m: float,
    derived_root_dir: Path,
    dem_cache_dir: Path,
    derived_dir: Path | None = None,
    derived_dirs: tuple[Path, ...] | None = None,
    radius_km: float = DEFAULT_FETCH_RADIUS_KM,
    min_distance_km: float = 0.0,
    min_height_m: float = 0.0,
    max_candidates: int = 5000,
    include_nearby_walls: bool = False,
    wall_radius_km: float = NEARBY_WALL_MAX_DISTANCE_M / 1000.0,
    observer_ground_elevation_m: float | None = None,
) -> UrbanOutlineResult | None:
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

    collected_buildings: list[BuildingFootprint] = []
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
            collected_buildings.extend(parse_derived_tile_buildings(envelope.path))
    buildings = _merge_building_footprints(tuple(collected_buildings))
    if min_height_m > 0.0:
        buildings = tuple(
            building
            for building in buildings
            if float(building.height_m) >= float(min_height_m)
        )
    if not buildings:
        return None
    sampled_observer_ground_m, buildings = _resolve_building_ground_elevations(
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        buildings=tuple(buildings),
        radius_km=radius_km,
        dem_cache_dir=dem_cache_dir,
        fallback_ground_elevation_m=(
            observer_ground_elevation_m
            if observer_ground_elevation_m is not None
            else 0.0
        ),
    )
    if observer_ground_elevation_m is None:
        observer_ground_elevation_m = sampled_observer_ground_m

    result = compute_urban_outlines(
        Viewpoint(
            id="coords",
            qid="",
            kind="coords",
            name="coords",
            labels={},
            names=(),
            latitude_deg=lat_deg,
            longitude_deg=lon_deg,
            height_m=observer_height_m,
            viewpoint_height_m=observer_height_m,
            meta={},
        ),
        buildings,
        radius_km=radius_km,
        min_distance_km=min_distance_km,
        observer_ground_elevation_m=observer_ground_elevation_m,
        view_center=None,
        edge_fov_deg=90.0,
        edge_sample_step_m=10.0,
        max_candidates=max_candidates,
    )
    if include_nearby_walls and wall_radius_km > 0.0:
        result = replace(
            result,
            walls=_nearby_wall_faces(
                buildings,
                lat_deg=lat_deg,
                lon_deg=lon_deg,
                observer_elevation_m=observer_ground_elevation_m + observer_height_m,
                max_distance_m=wall_radius_km * 1000.0,
            ),
        )
    return result


def _project_observer_centric_urban_outline_result(
    result: UrbanOutlineResult,
    *,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    front_hemisphere_view_center: tuple[float, float] | None = None,
    front_hemisphere_fov_deg: float = 90.0,
) -> list[UrbanOutlinePolyline] | None:
    if result is None:
        return None

    outlines = [
        _project_observer_centric_outline(
            outline,
            view_center=view_center,
            edge_fov_deg=edge_fov_deg,
            front_hemisphere_view_center=front_hemisphere_view_center,
            front_hemisphere_fov_deg=front_hemisphere_fov_deg,
        )
        for outline in result.outlines
    ]
    projected = [outline for outline in outlines if outline is not None]
    projected.extend(
        UrbanOutlinePolyline(
            points=[(point.altitude_deg, point.azimuth_deg) for point in wall.points],
            height_m=wall.height_m,
            distance_km=wall.distance_km,
            source="nearby_wall",
            wall_opacity_level=wall.wall_opacity_level,
        )
        for wall in getattr(result, "walls", ())
    )
    return projected or None


def _building_wall_prisms(
    building: BuildingFootprint,
    footprint_rings_xy: list[np.ndarray],
    transformer,
) -> list[tuple[list[np.ndarray], float, float]]:
    """Return (rings_xy, bottom_m, top_m) extrusions for one building's walls.

    PLATEAU buildings with LOD2 roof surfaces are extruded per roof polygon down
    to the next lower floor (see ``roof_floor_prisms``).  Other buildings,
    including Overture buildings and building parts, extrude their footprint
    from ``min_height_m`` to ``height_m``.
    """

    if building.geometry_lod >= 2 and building.roof_surfaces_lonlat:
        prisms = []
        for prism in roof_floor_prisms(
            building.roof_surfaces_lonlat, building.ground_elevation_m
        ):
            ring_xy = project_ring_xy(prism.ring_lonlat, transformer)
            if ring_contains_origin_xy(ring_xy):
                continue
            prisms.append(([ring_xy], prism.bottom_m, prism.top_m))
        if prisms:
            return prisms
    return [(footprint_rings_xy, building.min_height_m, building.height_m)]


def _nearby_wall_faces(
    buildings: tuple[BuildingFootprint, ...],
    *,
    lat_deg: float,
    lon_deg: float,
    observer_elevation_m: float,
    max_distance_m: float = NEARBY_WALL_MAX_DISTANCE_M,
) -> tuple[ObserverOutline, ...]:
    """Build small observer-centric wall patches from unmodified footprints."""
    from .location_resolver.resolve import _ring_min_distance_to_lonlat_m

    viewpoint = type(
        "WallViewpoint", (), {"latitude_deg": lat_deg, "longitude_deg": lon_deg}
    )()
    transformer = make_local_transformer(viewpoint)
    faces: list[ObserverOutline] = []
    segments: list[WallSegment] = []
    lat_margin = max_distance_m / 110000.0
    lon_margin = max_distance_m / max(
        1.0, 111000.0 * abs(math.cos(math.radians(lat_deg)))
    )
    for building in buildings:
        points = [point for ring in building.rings_lonlat for point in ring]
        if not points:
            continue
        if (
            min(p[0] for p in points) > lon_deg + lon_margin
            or max(p[0] for p in points) < lon_deg - lon_margin
            or min(p[1] for p in points) > lat_deg + lat_margin
            or max(p[1] for p in points) < lat_deg - lat_margin
        ):
            continue
        rings = [project_ring_xy(ring, transformer) for ring in building.rings_lonlat]
        if not rings or any(ring_contains_origin_xy(ring) for ring in rings):
            continue
        distance = min(
            _ring_min_distance_to_lonlat_m(ring, lon_deg=lon_deg, lat_deg=lat_deg)
            for ring in building.rings_lonlat
        )
        if not NEARBY_WALL_MIN_DISTANCE_M <= distance <= max_distance_m:
            continue
        if wall_opacity_level(distance, max_distance_m) <= 0:
            continue
        for prism_rings, prism_bottom_m, prism_top_m in _building_wall_prisms(
            building, rings, transformer
        ):
            bottom = building.ground_elevation_m + prism_bottom_m - observer_elevation_m
            top = building.ground_elevation_m + prism_top_m - observer_elevation_m
            if top <= bottom:
                continue
            perimeter = sum(
                float(np.linalg.norm(end - start))
                for ring in prism_rings
                for start, end in zip(ring[:-1], ring[1:])
            )
            step_m = max(2.0, math.sqrt(perimeter * (top - bottom) / 10000.0))
            for ring in prism_rings:
                for start, end in zip(ring[:-1], ring[1:]):
                    segments.append(
                        WallSegment(
                            start=(float(start[0]), float(start[1])),
                            end=(float(end[0]), float(end[1])),
                            bottom=bottom,
                            top=top,
                            distance_m=distance,
                            height_m=building.height_m,
                            step_m=step_m,
                        )
                    )
    for segment in visible_wall_segments(segments, max_distance_m):
        start = np.asarray(segment.start)
        end = np.asarray(segment.end)
        bottom, top = segment.bottom, segment.top
        step_m = segment.step_m
        length = float(np.linalg.norm(end - start))
        horizontal_steps = max(1, math.ceil(length / step_m))
        vertical_steps = max(1, math.ceil((top - bottom) / step_m))
        for i in range(horizontal_steps):
            a = start + (end - start) * (i / horizontal_steps)
            b = start + (end - start) * ((i + 1) / horizontal_steps)
            for j in range(vertical_steps):
                lo = bottom + (top - bottom) * j / vertical_steps
                hi = bottom + (top - bottom) * (j + 1) / vertical_steps
                wall_points: list[UrbanPolylinePoint] = []
                for xy, z in ((a, lo), (b, lo), (b, hi), (a, hi)):
                    d = float(np.hypot(xy[0], xy[1]))
                    wall_points.append(
                        UrbanPolylinePoint(
                            azimuth_deg=math.degrees(math.atan2(xy[0], xy[1])) % 360,
                            altitude_deg=math.degrees(math.atan2(z, d)),
                        )
                    )
                faces.append(
                    ObserverOutline(
                        height_m=segment.height_m,
                        distance_km=segment.distance_m / 1000.0,
                        points=tuple(wall_points),
                        wall_opacity_level=wall_opacity_level(segment.distance_m, max_distance_m),
                    )
                )
    return tuple(faces)


def _project_observer_centric_outline(
    outline,
    *,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    front_hemisphere_view_center: tuple[float, float] | None = None,
    front_hemisphere_fov_deg: float = 90.0,
) -> UrbanOutlinePolyline | None:
    run_points = list(outline.points)
    if front_hemisphere_view_center is not None and run_points:
        alt = np.array([point.altitude_deg for point in run_points], dtype=np.float64)
        az = np.array([point.azimuth_deg for point in run_points], dtype=np.float64)
        visible = is_in_fov_vectorized(
            alt,
            az,
            (
                float(front_hemisphere_view_center[0]),
                float(front_hemisphere_view_center[1]),
            ),
            fov_deg=float(front_hemisphere_fov_deg),
        )
        run_points = [
            point
            for point, keep in zip(run_points, visible, strict=False)
            if bool(keep)
        ]
    run_points = _maybe_linearize_run_points(
        run_points,
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
    )
    if len(run_points) == 0:
        return None
    return UrbanOutlinePolyline(
        points=[(point.altitude_deg, point.azimuth_deg) for point in run_points],
        height_m=float(outline.height_m),
        distance_km=float(outline.distance_km),
        source="base",
    )


def _merge_building_footprints(
    buildings: tuple[BuildingFootprint, ...],
) -> tuple[BuildingFootprint, ...]:
    seen: set[BuildingFootprint] = set()
    unique: list[BuildingFootprint] = []
    for building in buildings:
        if not building.id_is_synthetic:
            if building in seen:
                continue
            seen.add(building)
        unique.append(building)
    parent_ids_with_parts = {
        building.parent_building_id
        for building in unique
        if building.parent_building_id
    }
    if not parent_ids_with_parts:
        return tuple(unique)
    return tuple(
        building
        for building in unique
        if building.id_is_synthetic
        or building.parent_building_id is not None
        or building.building_id not in parent_ids_with_parts
    )


@lru_cache(maxsize=8)
def _list_derived_dirs(derived_root_dir: Path) -> tuple[Path, ...]:
    return tuple(
        path for path in sorted(derived_root_dir.glob("*/bldg")) if path.is_dir()
    )


def _resolve_building_ground_elevations(
    *,
    lat_deg: float,
    lon_deg: float,
    buildings: tuple[BuildingFootprint, ...],
    radius_km: float,
    dem_cache_dir: Path,
    fallback_ground_elevation_m: float = 0.0,
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
        if (
            str(exc)
            == "No Copernicus DEM tiles were downloaded for the requested area."
        ):
            return fallback_ground_elevation_m, tuple(
                replace(building, ground_elevation_m=fallback_ground_elevation_m)
                for building in buildings
            )
        raise

    dem = GeoTiffDem(download.paths, default_elevation_m=fallback_ground_elevation_m)
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
        building_lonlat = tuple(
            _representative_point_lonlat(building) for building in buildings
        )
        building_lon = np.array(
            [point[0] for point in building_lonlat], dtype=np.float64
        )
        building_lat = np.array(
            [point[1] for point in building_lonlat], dtype=np.float64
        )
        sampled_ground = dem_grid.sample_lonlat(
            building_lon, building_lat, method="bilinear"
        )
        resolved_buildings = tuple(
            BuildingFootprint(
                building_id=building.building_id,
                height_m=building.height_m,
                rings_lonlat=building.rings_lonlat,
                geometry_lod=building.geometry_lod,
                roof_surfaces_lonlat=building.roof_surfaces_lonlat,
                parent_building_id=building.parent_building_id,
                ground_elevation_m=(
                    float(sampled_ground[index])
                    if math.isfinite(float(sampled_ground[index]))
                    else fallback_ground_elevation_m
                ),
                min_height_m=float(building.min_height_m),
                id_is_synthetic=building.id_is_synthetic,
            )
            for index, building in enumerate(buildings)
        )
        if not math.isfinite(observer_ground_elevation_m):
            observer_ground_elevation_m = fallback_ground_elevation_m
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
    pairs = (
        zip(ring_lonlat, ring_lonlat[1:])
        if ring_lonlat[0] == ring_lonlat[-1]
        else zip(ring_lonlat, ring_lonlat[1:] + ring_lonlat[:1])
    )
    for (x0, y0), (x1, y1) in pairs:
        area += x0 * y1 - x1 * y0
    return 0.5 * area


def _point_in_ring(
    point_lonlat: tuple[float, float], ring_lonlat: tuple[tuple[float, float], ...]
) -> bool:
    px, py = point_lonlat
    if len(ring_lonlat) < 3:
        return False
    inside = False
    points = (
        ring_lonlat
        if ring_lonlat[0] == ring_lonlat[-1]
        else ring_lonlat + (ring_lonlat[0],)
    )
    for (x0, y0), (x1, y1) in zip(points[:-1], points[1:]):
        intersects = ((y0 > py) != (y1 > py)) and (
            px < (x1 - x0) * (py - y0) / ((y1 - y0) or 1e-12) + x0
        )
        if intersects:
            inside = not inside
    return inside


def _bbox_center_lonlat(
    ring_lonlat: tuple[tuple[float, float], ...],
) -> tuple[float, float]:
    lon = [point[0] for point in ring_lonlat]
    lat = [point[1] for point in ring_lonlat]
    return ((min(lon) + max(lon)) * 0.5, (min(lat) + max(lat)) * 0.5)


def _mean_point_lonlat(
    ring_lonlat: tuple[tuple[float, float], ...],
) -> tuple[float, float]:
    points = (
        ring_lonlat[:-1]
        if len(ring_lonlat) > 1 and ring_lonlat[0] == ring_lonlat[-1]
        else ring_lonlat
    )
    lon = sum(point[0] for point in points) / max(1, len(points))
    lat = sum(point[1] for point in points) / max(1, len(points))
    return (lon, lat)
