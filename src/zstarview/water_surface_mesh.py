from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import mapbox_earcut as earcut
import numpy as np
from pyproj import CRS, Transformer

from .water_overlay import WaterPolygonFootprint, WaterSurfacePatch


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


def _dedupe_closed_ring(points_xy: Sequence[tuple[float, float]]) -> list[tuple[float, float]]:
    if not points_xy:
        return []
    deduped: list[tuple[float, float]] = []
    for point in points_xy:
        if not deduped or point != deduped[-1]:
            deduped.append(point)
    if len(deduped) >= 2 and deduped[0] == deduped[-1]:
        return deduped
    return deduped + [deduped[0]]


def _ring_body(points_xy: Sequence[tuple[float, float]]) -> list[tuple[float, float]]:
    ring = _dedupe_closed_ring(points_xy)
    if len(ring) < 4:
        return []
    return ring[:-1]


def _ring_signed_area(points_xy: Sequence[tuple[float, float]]) -> float:
    body = _ring_body(points_xy)
    if len(body) < 3:
        return 0.0
    area = 0.0
    for index in range(len(body)):
        x0, y0 = body[index]
        x1, y1 = body[(index + 1) % len(body)]
        area += x0 * y1 - x1 * y0
    return 0.5 * area


def _round_ring(
    points_xy: Sequence[tuple[float, float]],
    *,
    grid_m: float,
) -> list[tuple[float, float]]:
    ring = _dedupe_closed_ring(points_xy)
    if grid_m <= 0.0 or len(ring) < 4:
        return ring
    rounded: list[tuple[float, float]] = []
    for x, y in ring[:-1]:
        rounded_point = (round(float(x) / grid_m) * grid_m, round(float(y) / grid_m) * grid_m)
        if not rounded or rounded_point != rounded[-1]:
            rounded.append(rounded_point)
    if len(rounded) < 3:
        return ring
    if rounded[0] != rounded[-1]:
        rounded.append(rounded[0])
    return rounded


def _simplify_body(
    body: Sequence[tuple[float, float]],
    *,
    tolerance_m: float,
) -> list[tuple[float, float]]:
    if tolerance_m <= 0.0 or len(body) <= 4:
        return list(body)

    simplified = _rdp(list(body), tolerance_m)
    if len(simplified) < 3:
        return list(body)
    return simplified


def _rdp(
    points: list[tuple[float, float]],
    tolerance_m: float,
) -> list[tuple[float, float]]:
    if len(points) <= 2:
        return list(points)

    keep = [False] * len(points)
    keep[0] = True
    keep[-1] = True
    stack: list[tuple[int, int]] = [(0, len(points) - 1)]
    while stack:
        start_index, end_index = stack.pop()
        max_distance = -1.0
        split_index = -1
        ax, ay = points[start_index]
        bx, by = points[end_index]
        dx = bx - ax
        dy = by - ay
        segment_length_sq = dx * dx + dy * dy
        for index in range(start_index + 1, end_index):
            px, py = points[index]
            if segment_length_sq == 0.0:
                distance = math.hypot(px - ax, py - ay)
            else:
                t = ((px - ax) * dx + (py - ay) * dy) / segment_length_sq
                t = max(0.0, min(1.0, t))
                cx = ax + t * dx
                cy = ay + t * dy
                distance = math.hypot(px - cx, py - cy)
            if distance > max_distance:
                max_distance = distance
                split_index = index
        if split_index >= 0 and max_distance > tolerance_m:
            keep[split_index] = True
            stack.append((start_index, split_index))
            stack.append((split_index, end_index))

    simplified = [point for point, selected in zip(points, keep) if selected]
    return simplified if len(simplified) >= 3 else list(points)


def _point_in_ring(point: tuple[float, float], ring: Sequence[tuple[float, float]]) -> bool:
    body = _ring_body(ring)
    if len(body) < 3:
        return False
    x, y = point
    inside = False
    for index in range(len(body)):
        x0, y0 = body[index]
        x1, y1 = body[(index + 1) % len(body)]
        if ((y0 > y) != (y1 > y)) and (y1 != y0):
            x_cross = ((x1 - x0) * (y - y0) / (y1 - y0)) + x0
            if x < x_cross:
                inside = not inside
    return inside


def _assign_holes_to_shells(
    shells: Sequence[list[tuple[float, float]]],
    holes: Sequence[list[tuple[float, float]]],
) -> list[list[list[tuple[float, float]]]]:
    assigned: list[list[list[tuple[float, float]]]] = [[] for _ in shells]
    shell_areas = [abs(_ring_signed_area(shell)) for shell in shells]
    for hole in holes:
        if len(hole) < 4:
            continue
        anchor = hole[0]
        matches: list[int] = []
        for index, shell in enumerate(shells):
            if _point_in_ring(anchor, shell):
                matches.append(index)
        if not matches:
            continue
        target_index = min(matches, key=lambda index: shell_areas[index])
        assigned[target_index].append(hole)
    return assigned


def _triangulate_shell(
    shell: list[tuple[float, float]],
    holes: Sequence[list[tuple[float, float]]],
) -> tuple[tuple[tuple[float, float], tuple[float, float], tuple[float, float]], ...]:
    vertices: list[tuple[float, float]] = []
    ring_ends: list[int] = []

    def add_ring(ring: Sequence[tuple[float, float]]) -> None:
        body = _ring_body(ring)
        if len(body) < 3:
            return
        vertices.extend(body)
        ring_ends.append(len(vertices))

    add_ring(shell)
    for hole in holes:
        add_ring(hole)

    if not ring_ends or ring_ends[-1] != len(vertices):
        return ()

    vertices_array = np.asarray(vertices, dtype=np.float64)
    ring_ends_array = np.asarray(ring_ends, dtype=np.uint32)
    if vertices_array.shape[0] < 3:
        return ()

    indices = earcut.triangulate_float64(vertices_array, ring_ends_array)
    if indices.size == 0:
        return ()

    result: list[tuple[tuple[float, float], tuple[float, float], tuple[float, float]]] = []
    for tri_start in range(0, int(indices.size), 3):
        tri_indices = indices[tri_start : tri_start + 3]
        if len(tri_indices) < 3:
            continue
        triangle = tuple(
            tuple(vertices_array[int(index)])
            for index in tri_indices
        )
        result.append(
            tuple((float(x), float(y)) for x, y in triangle)
        )
    return tuple(result)


def _mean_surface_elevation_m(patch: WaterSurfacePatch | None) -> float:
    if patch is None:
        return 0.0
    values = [float(value) for value in patch.anchor_elevations_m]
    return sum(values) / float(len(values))


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

    shells: list[list[tuple[float, float]]] = []
    for outer_ring in footprint.outer_rings_lonlat:
        projected = project_ring_xy(outer_ring, transformer)
        rounded = _round_ring(projected, grid_m=grid_m)
        body = _simplify_body(_ring_body(rounded), tolerance_m=simplify_tolerance_m)
        if len(body) < 3:
            continue
        if _ring_signed_area(body) < 0.0:
            body = list(reversed(body))
        shells.append(body + [body[0]])

    if not shells:
        return None

    holes: list[list[tuple[float, float]]] = []
    for inner_ring in footprint.inner_rings_lonlat:
        projected = project_ring_xy(inner_ring, transformer)
        rounded = _round_ring(projected, grid_m=grid_m)
        body = _simplify_body(_ring_body(rounded), tolerance_m=simplify_tolerance_m)
        if len(body) < 3:
            continue
        if _ring_signed_area(body) > 0.0:
            body = list(reversed(body))
        holes.append(body + [body[0]])

    assigned_holes = _assign_holes_to_shells(shells, holes)
    triangles: list[tuple[tuple[float, float], tuple[float, float], tuple[float, float]]] = []
    for shell, shell_holes in zip(shells, assigned_holes):
        triangles.extend(_triangulate_shell(shell, shell_holes))

    if not triangles:
        return None

    surface_mode = patch.surface_mode if patch is not None else "flat"
    return WaterSurfaceMesh(
        water_id=footprint.water_id,
        kind=footprint.kind,
        surface_mode=surface_mode,
        surface_elevation_m=_mean_surface_elevation_m(patch),
        triangles_xy_m=tuple(triangles),
        source=footprint.source,
        simplified_grid_m=float(grid_m),
        simplified_tolerance_m=float(simplify_tolerance_m),
    )
