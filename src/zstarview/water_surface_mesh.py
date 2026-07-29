from __future__ import annotations

import math
from collections.abc import Callable, Sequence

from pyproj import CRS, Transformer

DEFAULT_SPLIT_CELL_M = 300.0


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


def _clip_ring_to_edge(
    points: Sequence[tuple[float, float]],
    *,
    edge_value: float,
    keep_inside: Callable[[float, float, float], bool],
    intersect: Callable[[tuple[float, float], tuple[float, float], float], tuple[float, float]],
) -> list[tuple[float, float]]:
    if not points:
        return []
    clipped: list[tuple[float, float]] = []
    previous = points[-1]
    previous_inside = keep_inside(previous[0], previous[1], edge_value)
    for current in points:
        current_inside = keep_inside(current[0], current[1], edge_value)
        if current_inside:
            if not previous_inside:
                clipped.append(intersect(previous, current, edge_value))
            clipped.append((float(current[0]), float(current[1])))
        elif previous_inside:
            clipped.append(intersect(previous, current, edge_value))
        previous = current
        previous_inside = current_inside
    if len(clipped) >= 2 and clipped[0] == clipped[-1]:
        clipped = clipped[:-1]
    return clipped


def _clip_ring_to_rect(
    ring: Sequence[tuple[float, float]],
    *,
    min_x: float,
    min_y: float,
    max_x: float,
    max_y: float,
) -> list[tuple[float, float]]:
    clipped = list(ring)
    if len(clipped) < 3:
        return []

    def keep_left(x: float, _: float, edge_value: float) -> bool:
        return x >= edge_value

    def keep_right(x: float, _: float, edge_value: float) -> bool:
        return x <= edge_value

    def keep_bottom(_: float, y: float, edge_value: float) -> bool:
        return y >= edge_value

    def keep_top(_: float, y: float, edge_value: float) -> bool:
        return y <= edge_value

    def intersect_vertical(
        p0: tuple[float, float],
        p1: tuple[float, float],
        x_value: float,
    ) -> tuple[float, float]:
        x0, y0 = p0
        x1, y1 = p1
        dx = x1 - x0
        if dx == 0.0:
            return (float(x_value), float(y0))
        t = (x_value - x0) / dx
        return (float(x_value), float(y0 + (t * (y1 - y0))))

    def intersect_horizontal(
        p0: tuple[float, float],
        p1: tuple[float, float],
        y_value: float,
    ) -> tuple[float, float]:
        x0, y0 = p0
        x1, y1 = p1
        dy = y1 - y0
        if dy == 0.0:
            return (float(x0), float(y_value))
        t = (y_value - y0) / dy
        return (float(x0 + (t * (x1 - x0))), float(y_value))

    clipped = _clip_ring_to_edge(
        clipped,
        edge_value=min_x,
        keep_inside=keep_left,
        intersect=intersect_vertical,
    )
    clipped = _clip_ring_to_edge(
        clipped,
        edge_value=max_x,
        keep_inside=keep_right,
        intersect=intersect_vertical,
    )
    clipped = _clip_ring_to_edge(
        clipped,
        edge_value=min_y,
        keep_inside=keep_bottom,
        intersect=intersect_horizontal,
    )
    clipped = _clip_ring_to_edge(
        clipped,
        edge_value=max_y,
        keep_inside=keep_top,
        intersect=intersect_horizontal,
    )
    if len(clipped) < 3:
        return []
    if clipped[0] != clipped[-1]:
        clipped.append(clipped[0])
    return clipped


def _ring_bounds(points_xy: Sequence[tuple[float, float]]) -> tuple[float, float, float, float]:
    xs = [point[0] for point in points_xy]
    ys = [point[1] for point in points_xy]
    if not xs or not ys:
        return 0.0, 0.0, 0.0, 0.0
    return min(xs), min(ys), max(xs), max(ys)


def _ring_centroid(points_xy: Sequence[tuple[float, float]]) -> tuple[float, float]:
    body = _ring_body(points_xy)
    if not body:
        return 0.0, 0.0
    sum_x = sum(point[0] for point in body)
    sum_y = sum(point[1] for point in body)
    count = float(len(body))
    return sum_x / count, sum_y / count


def split_local_polygon_by_grid(
    shell: Sequence[tuple[float, float]],
    holes: Sequence[Sequence[tuple[float, float]]],
    *,
    cell_m: float = DEFAULT_SPLIT_CELL_M,
) -> list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]]:
    if cell_m <= 0.0:
        return [(list(shell), [list(hole) for hole in holes])]
    if len(shell) < 4:
        return []

    min_x, min_y, max_x, max_y = _ring_bounds(shell)
    if min_x >= max_x or min_y >= max_y:
        return []

    pieces: list[tuple[list[tuple[float, float]], list[list[tuple[float, float]]]]] = []
    epsilon = max(1.0e-9, abs(cell_m) * 1.0e-12)
    start_x = int(math.floor(min_x / cell_m))
    end_x = int(math.floor((max_x - epsilon) / cell_m))
    start_y = int(math.floor(min_y / cell_m))
    end_y = int(math.floor((max_y - epsilon) / cell_m))
    for cell_x in range(start_x, end_x + 1):
        cell_min_x = float(cell_x) * cell_m
        cell_max_x = cell_min_x + cell_m
        for cell_y in range(start_y, end_y + 1):
            cell_min_y = float(cell_y) * cell_m
            cell_max_y = cell_min_y + cell_m
            clipped_shell = _clip_ring_to_rect(
                shell,
                min_x=cell_min_x,
                min_y=cell_min_y,
                max_x=cell_max_x,
                max_y=cell_max_y,
            )
            if len(clipped_shell) < 4:
                continue
            clipped_holes: list[list[tuple[float, float]]] = []
            for hole in holes:
                clipped_hole = _clip_ring_to_rect(
                    hole,
                    min_x=cell_min_x,
                    min_y=cell_min_y,
                    max_x=cell_max_x,
                    max_y=cell_max_y,
                )
                if len(clipped_hole) < 4:
                    continue
                if _point_in_ring(_ring_centroid(clipped_hole), clipped_shell):
                    clipped_holes.append(clipped_hole)
            pieces.append((clipped_shell, clipped_holes))
    return pieces
