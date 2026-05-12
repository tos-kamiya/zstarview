from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from pyproj import CRS, Transformer

try:
    from shapely import Point, Polygon, constrained_delaunay_triangles, set_precision, simplify
except ImportError:  # pragma: no cover - optional dependency in dev env
    Point = None  # type: ignore[assignment]
    Polygon = None  # type: ignore[assignment]
    constrained_delaunay_triangles = None  # type: ignore[assignment]
    set_precision = None  # type: ignore[assignment]
    simplify = None  # type: ignore[assignment]

from .water_overlay import WaterPolygonFootprint, WaterSurfacePatch


EARTH_RADIUS_M = 6371.0088 * 1000.0


@dataclass(frozen=True, slots=True)
class WaterSurfaceMesh:
    water_id: str
    kind: str
    surface_mode: str
    surface_elevation_m: float
    triangles_xy_m: tuple[tuple[tuple[float, float], tuple[float, float], tuple[float, float]], ...]
    source: str
    simplified_grid_m: float
    simplified_tolerance_m: float


def make_local_transformer(lat_deg: float, lon_deg: float) -> Transformer:
    local_crs = CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat_deg} +lon_0={lon_deg} +datum=WGS84 +units=m +no_defs"
    )
    return Transformer.from_crs("EPSG:4326", local_crs, always_xy=True)


def project_ring_xy(
    ring_lonlat: tuple[tuple[float, float], ...],
    transformer: Transformer,
) -> list[tuple[float, float]]:
    lon = [point[0] for point in ring_lonlat]
    lat = [point[1] for point in ring_lonlat]
    x, y = transformer.transform(lon, lat)
    return [(float(px), float(py)) for px, py in zip(x, y)]


def _dedupe_closed_ring(points_xy: list[tuple[float, float]]) -> list[tuple[float, float]]:
    if not points_xy:
        return []
    deduped: list[tuple[float, float]] = []
    for point in points_xy:
        if not deduped or point != deduped[-1]:
            deduped.append(point)
    if len(deduped) >= 2 and deduped[0] == deduped[-1]:
        return deduped
    return deduped + [deduped[0]]


def _build_polygon_geometry(
    shell_xy: list[tuple[float, float]],
    hole_xy_list: list[list[tuple[float, float]]],
    *,
    grid_m: float,
    simplify_m: float,
) -> Any:
    if Polygon is None:
        raise RuntimeError("shapely is required for water surface mesh generation")
    polygon = Polygon(shell_xy, holes=hole_xy_list or None)
    if polygon.is_empty or not polygon.is_valid or polygon.area <= 0.0:
        return None
    if set_precision is not None and grid_m > 0.0:
        polygon = set_precision(polygon, grid_m)
    if simplify is not None and simplify_m > 0.0:
        polygon = simplify(polygon, simplify_m, preserve_topology=True)
    if polygon.is_empty or not polygon.is_valid or polygon.area <= 0.0:
        return None
    return polygon


def _triangles_from_geometry(geometry: Any) -> tuple[tuple[tuple[float, float], tuple[float, float], tuple[float, float]], ...]:
    if constrained_delaunay_triangles is None:
        raise RuntimeError("shapely is required for water surface mesh generation")
    triangles = constrained_delaunay_triangles(geometry)
    result: list[tuple[tuple[float, float], tuple[float, float], tuple[float, float]]] = []
    for triangle in triangles.geoms:
        if triangle.is_empty or triangle.geom_type != "Polygon":
            continue
        coords = list(triangle.exterior.coords)
        if len(coords) < 4:
            continue
        if coords[0] == coords[-1]:
            coords = coords[:-1]
        result.append((coords[0], coords[1], coords[2]))
    return tuple(result)


def build_water_surface_mesh(
    footprint: WaterPolygonFootprint,
    *,
    center_lat_deg: float,
    center_lon_deg: float,
    patch: WaterSurfacePatch | None = None,
    grid_m: float = 1.0,
    simplify_tolerance_m: float = 20.0,
) -> WaterSurfaceMesh | None:
    if not footprint.outer_rings_lonlat:
        return None
    transformer = make_local_transformer(center_lat_deg, center_lon_deg)
    triangles: list[tuple[tuple[float, float], tuple[float, float], tuple[float, float]]] = []

    for outer_ring in footprint.outer_rings_lonlat:
        shell_xy = _dedupe_closed_ring(project_ring_xy(outer_ring, transformer))
        if len(shell_xy) < 4:
            continue
        hole_xy_list: list[list[tuple[float, float]]] = []
        for inner_ring in footprint.inner_rings_lonlat:
            inner_xy = _dedupe_closed_ring(project_ring_xy(inner_ring, transformer))
            if len(inner_xy) < 4:
                continue
            if Polygon is not None:
                try:
                    anchor = Point(inner_xy[0])
                except Exception:
                    anchor = None
                if anchor is None:
                    continue
                shell_geom = Polygon(shell_xy)
                if not shell_geom.covers(anchor):
                    continue
            hole_xy_list.append(inner_xy)
        geometry = _build_polygon_geometry(
            shell_xy,
            hole_xy_list,
            grid_m=grid_m,
            simplify_m=simplify_tolerance_m,
        )
        if geometry is None:
            continue
        triangles.extend(_triangles_from_geometry(geometry))

    if not triangles:
        return None

    surface_mode = patch.surface_mode if patch is not None else "flat"
    if patch is not None:
        surface_elevation_m = float(sum(float(value) for value in patch.anchor_elevations_m) / 3.0)
    else:
        surface_elevation_m = 0.0
    return WaterSurfaceMesh(
        water_id=footprint.water_id,
        kind=footprint.kind,
        surface_mode=surface_mode,
        surface_elevation_m=surface_elevation_m,
        triangles_xy_m=tuple(triangles),
        source=footprint.source,
        simplified_grid_m=float(grid_m),
        simplified_tolerance_m=float(simplify_tolerance_m),
    )
