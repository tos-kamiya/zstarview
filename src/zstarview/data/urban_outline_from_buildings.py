#!/usr/bin/env python3
"""Build urban outline polylines from derived building tiles."""

from __future__ import annotations

import argparse
import heapq
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.astro import altaz_to_normalized_xy  # noqa: E402
from zstarview.data.derived_tile_cache import (  # noqa: E402
    parse_derived_tile_buildings,
    select_derived_tile_envelopes,
)
from zstarview.data.urban_outline_common import (  # noqa: E402
    BuildingFootprint,
    bbox_min_distance_m,
    iter_true_runs,
    make_local_transformer,
    project_ring_xy,
    sample_ring_points_xy,
)
from zstarview.location_resolver import (  # noqa: E402
    TowerViewpoint,
    load_tower_viewpoints,
    resolve_tower_viewpoint,
)
from zstarview.location_resolver.viewpoints import normalize_viewpoint_name  # noqa: E402

DEFAULT_RADIUS_KM = 3.0
DEFAULT_MIN_BUILDING_HEIGHT_M = 40.0
DEFAULT_EDGE_SAMPLE_STEP_M = 10.0
MAX_URBAN_OUTLINE_CANDIDATES = 5000
HOLE_RING_SUPPRESSION_MIN_DISTANCE_M = 1000.0
HOLE_RING_SUPPRESSION_MAX_SPAN_M = 250.0
VERTICAL_THIN_RUN_MAX_NORMALIZED_HEIGHT = 0.01
VERTICAL_THIN_RUN_MAX_ABS_ALTITUDE_DEG = 5.0
DETAILED_ROOF_MIN_BUILDING_HEIGHT_M = 40.0
ROOF_SURFACE_ELEVATION_TOLERANCE_M = -1.0
ROOF_SURFACE_FAR_ELEVATION_TOLERANCE_M = 3.0
ROOF_SURFACE_TOLERANCE_FAR_DISTANCE_M = 4000.0
ROOF_SURFACE_EDGE_SNAP_M = 0.1


@dataclass(frozen=True)
class UrbanPolylinePoint:
    azimuth_deg: float
    altitude_deg: float


@dataclass(frozen=True)
class UrbanOutlinePolyline:
    height_m: float
    distance_km: float
    points: tuple[UrbanPolylinePoint, ...]


@dataclass(frozen=True)
class UrbanOutlineResult:
    tower: TowerViewpoint
    outlines: tuple[UrbanOutlinePolyline, ...]
    buildings_considered: int
    outlines_emitted: int


@dataclass(frozen=True)
class _UrbanOutlineRingCandidate:
    score: float
    order: int
    building: BuildingFootprint
    ring_xy: np.ndarray
    building_distance_m: float
    min_distance_m: float


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate urban outline overlays from derived building tiles."
    )
    parser.add_argument(
        "--derived-dir",
        type=Path,
        required=True,
        help="Directory containing derived building tile JSON files.",
    )
    parser.add_argument(
        "--tower",
        action="append",
        default=[],
        help="Bundled tower viewpoint name or wikidata:Q... identifier. May be repeated.",
    )
    parser.add_argument(
        "--all-covered-towers",
        action="store_true",
        help="Process every bundled Japan viewpoint whose search circle intersects the provided derived tiles.",
    )
    parser.add_argument(
        "--radius-km",
        type=float,
        default=DEFAULT_RADIUS_KM,
        help="Maximum building search radius around each viewpoint (default: 3.0).",
    )
    parser.add_argument(
        "--min-building-height-m",
        type=float,
        default=DEFAULT_MIN_BUILDING_HEIGHT_M,
        help="Ignore buildings lower than this height in meters (default: 40.0).",
    )
    parser.add_argument(
        "--edge-sample-step-m",
        type=float,
        default=DEFAULT_EDGE_SAMPLE_STEP_M,
        help="Approximate spacing for sampling building edges in meters (default: 10.0).",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        required=True,
        help="Output JSON path for one-off outline export.",
    )
    return parser


def select_towers(
    *,
    tower_queries: Sequence[str],
    all_covered_towers: bool,
    derived_dir: Path,
    radius_km: float,
) -> tuple[TowerViewpoint, ...]:
    if tower_queries:
        selected: list[TowerViewpoint] = []
        for query in tower_queries:
            tower = resolve_tower_viewpoint(query)
            if tower is None:
                raise ValueError(f"Viewpoint not found: {query}")
            if not _is_japan_tower(tower):
                continue
            selected.append(tower)
        return tuple(selected)

    if not all_covered_towers:
        raise ValueError("Specify --tower or --all-covered-towers.")

    selected: list[TowerViewpoint] = []
    for tower in load_tower_viewpoints():
        if not _is_japan_tower(tower):
            continue
        try:
            select_derived_tile_envelopes(
                derived_dir,
                observer_lat_deg=tower.latitude_deg,
                observer_lon_deg=tower.longitude_deg,
                radius_km=radius_km,
            )
        except ValueError:
            continue
        selected.append(tower)
    return tuple(selected)


def _is_japan_tower(tower: TowerViewpoint) -> bool:
    country = str(tower.meta.get("country", "")).strip().casefold()
    if country in {"jp", "jpn", "japan", "\u65e5\u672c"}:
        return True
    labels = " ".join(tower.labels.values())
    names = " ".join(tower.names)
    haystack = " ".join((tower.name, labels, names))
    return "\u65e5\u672c" in haystack or normalize_viewpoint_name(haystack).endswith(" japan")


def compute_urban_outlines(
    tower: TowerViewpoint,
    buildings: Sequence[BuildingFootprint],
    *,
    radius_km: float,
    min_distance_km: float = 0.0,
    observer_ground_elevation_m: float = 0.0,
    view_center: tuple[float, float] | None = None,
    edge_fov_deg: float = 90.0,
    edge_sample_step_m: float,
    max_candidates: int = MAX_URBAN_OUTLINE_CANDIDATES,
) -> UrbanOutlineResult:
    if radius_km <= 0.0:
        raise ValueError("--radius-km must be positive.")
    if min_distance_km < 0.0:
        raise ValueError("--min-distance-km must be non-negative.")
    if min_distance_km >= radius_km:
        raise ValueError("--min-distance-km must be smaller than --radius-km.")
    if edge_sample_step_m <= 0.0:
        raise ValueError("--edge-sample-step-m must be positive.")

    transformer = make_local_transformer(tower)
    radius_m = radius_km * 1000.0
    min_distance_m = min_distance_km * 1000.0
    observer_height_m = float(getattr(tower, "viewpoint_height_m", getattr(tower, "observer_height_m", 0.0)) or 0.0)
    observer_elevation_m = float(observer_ground_elevation_m) + observer_height_m
    buildings_considered = 0
    candidate_heap: list[tuple[float, int, _UrbanOutlineRingCandidate]] = []
    candidate_order = 0
    candidate_capacity = max(0, int(max_candidates))

    for building in buildings:
        projected_rings = tuple(project_ring_xy(ring, transformer) for ring in building.rings_lonlat)
        if not projected_rings:
            continue
        if any(ring_contains_origin_xy(ring_xy) for ring_xy in projected_rings if ring_xy.size > 0):
            continue
        min_distance = min(bbox_min_distance_m(ring_xy) for ring_xy in projected_rings if ring_xy.size > 0)
        if min_distance > radius_m:
            continue
        max_distance = max(float(np_hypot_xy(ring_xy).max()) for ring_xy in projected_rings if ring_xy.size > 0)
        if max_distance <= min_distance_m:
            continue
        buildings_considered += 1
        for ring_index, ring_xy in enumerate(projected_rings):
            if ring_index > 0 and _should_skip_hole_ring(ring_xy, building_distance_m=min_distance):
                continue
            score = _urban_outline_ring_candidate_score(
                building=building,
                ring_xy=ring_xy,
                building_distance_m=min_distance,
            )
            if score <= 0.0:
                continue
            candidate = _UrbanOutlineRingCandidate(
                score=score,
                order=candidate_order,
                building=building,
                ring_xy=ring_xy,
                building_distance_m=min_distance,
                min_distance_m=min_distance_m,
            )
            candidate_order += 1
            if candidate_capacity == 0:
                continue
            if len(candidate_heap) < candidate_capacity:
                heapq.heappush(candidate_heap, (candidate.score, candidate.order, candidate))
            elif candidate.score > candidate_heap[0][0]:
                heapq.heapreplace(candidate_heap, (candidate.score, candidate.order, candidate))

    selected_candidates = [entry[2] for entry in candidate_heap]
    selected_candidates.sort(key=lambda candidate: candidate.order)
    outlines: list[UrbanOutlinePolyline] = []
    for candidate in selected_candidates:
        if (
            candidate.building.geometry_lod >= 2
            and candidate.building.roof_surfaces_lonlat
            and candidate.building.height_m >= DETAILED_ROOF_MIN_BUILDING_HEIGHT_M
        ):
            _emit_lod2_roof_surface_outlines(
                building=candidate.building,
                transformer=transformer,
                observer_elevation_m=observer_elevation_m,
                building_distance_m=candidate.building_distance_m,
                min_distance_m=candidate.min_distance_m,
                radius_m=radius_m,
                edge_sample_step_m=edge_sample_step_m,
                view_center=view_center,
                edge_fov_deg=edge_fov_deg,
                outlines=outlines,
            )
        else:
            emitter = (
                _emit_lod1_roof_polygon_outline
                if candidate.building.geometry_lod >= 1
                else _emit_ring_outlines
            )
            emitter(
                building=candidate.building,
                ring_xy=candidate.ring_xy,
                observer_elevation_m=observer_elevation_m,
                building_distance_m=candidate.building_distance_m,
                min_distance_m=candidate.min_distance_m,
                radius_m=radius_m,
                edge_sample_step_m=edge_sample_step_m,
                view_center=view_center,
                edge_fov_deg=edge_fov_deg,
                outlines=outlines,
            )
    return UrbanOutlineResult(
        tower=tower,
        outlines=tuple(outlines),
        buildings_considered=buildings_considered,
        outlines_emitted=len(outlines),
    )


def _emit_ring_outlines(
    *,
    building: BuildingFootprint,
    ring_xy: np.ndarray,
    observer_elevation_m: float,
    building_distance_m: float,
    min_distance_m: float,
    radius_m: float,
    edge_sample_step_m: float,
    view_center: tuple[float, float] | None,
    edge_fov_deg: float,
    outlines: list[UrbanOutlinePolyline],
) -> None:
    _emit_roof_polygon_outline(
        building=building,
        ring_xy=ring_xy,
        observer_elevation_m=observer_elevation_m,
        building_distance_m=building_distance_m,
        min_distance_m=min_distance_m,
        radius_m=radius_m,
        edge_sample_step_m=edge_sample_step_m,
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
        outlines=outlines,
    )


def _emit_lod1_roof_polygon_outline(
    *,
    building: BuildingFootprint,
    ring_xy: np.ndarray,
    observer_elevation_m: float,
    building_distance_m: float,
    min_distance_m: float,
    radius_m: float,
    edge_sample_step_m: float,
    view_center: tuple[float, float] | None,
    edge_fov_deg: float,
    outlines: list[UrbanOutlinePolyline],
) -> None:
    """Emit only the projected boundary of an LOD1 roof polygon."""
    _emit_roof_polygon_outline(
        building=building,
        ring_xy=ring_xy,
        observer_elevation_m=observer_elevation_m,
        building_distance_m=building_distance_m,
        min_distance_m=min_distance_m,
        radius_m=radius_m,
        edge_sample_step_m=edge_sample_step_m,
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
        outlines=outlines,
    )


def _emit_lod2_roof_surface_outlines(
    *,
    building: BuildingFootprint,
    transformer,
    observer_elevation_m: float,
    building_distance_m: float,
    min_distance_m: float,
    radius_m: float,
    edge_sample_step_m: float,
    view_center: tuple[float, float] | None,
    edge_fov_deg: float,
    outlines: list[UrbanOutlinePolyline],
) -> None:
    roof_surfaces = _select_outer_roof_surfaces(
        building.roof_surfaces_lonlat,
        transformer=transformer,
    )
    elevation_tolerance_m = _roof_surface_elevation_tolerance(
        building_distance_m
    )
    projected_surfaces = []
    for surface in roof_surfaces:
        if len(surface) < 4:
            continue
        ring_lonlat = tuple((lon, lat) for lon, lat, _elevation in surface)
        ring_xy = project_ring_xy(ring_lonlat, transformer)
        elevations = np.asarray(
            [elevation for _lon, _lat, elevation in surface], dtype=np.float64
        )
        projected_surfaces.append((ring_xy, elevations))

    for ring_xy, elevations in _merge_projected_roof_surface_groups(
        projected_surfaces,
        elevation_tolerance_m=elevation_tolerance_m,
    ):
        sampled_points = _sample_surface_points_xy(
            ring_xy,
            elevations,
            sample_step_m=edge_sample_step_m,
        )
        if sampled_points.size == 0:
            continue
        distances = np.hypot(sampled_points[:, 0], sampled_points[:, 1])
        valid = (distances > max(0.1, min_distance_m)) & (distances <= radius_m)
        if not valid.any():
            continue
        azimuth_deg = (
            np.degrees(np.arctan2(sampled_points[:, 0], sampled_points[:, 1])) + 360.0
        ) % 360.0
        altitude_deg = np.degrees(
            np.arctan2(sampled_points[:, 2] - observer_elevation_m, distances)
        )
        for run in iter_true_runs(valid):
            run_points = [
                UrbanPolylinePoint(azimuth_deg=float(az), altitude_deg=float(alt))
                for az, alt in zip(azimuth_deg[run], altitude_deg[run])
            ]
            run_points = _maybe_linearize_run_points(
                run_points,
                view_center=view_center,
                edge_fov_deg=edge_fov_deg,
            )
            if len(run_points) >= 2:
                outlines.append(
                    UrbanOutlinePolyline(
                        height_m=float(building.height_m),
                        distance_km=float(building_distance_m / 1000.0),
                        points=tuple(run_points),
                    )
                )


def _merge_projected_roof_surface_groups(
    surfaces: Sequence[tuple[np.ndarray, np.ndarray]],
    *,
    elevation_tolerance_m: float = ROOF_SURFACE_ELEVATION_TOLERANCE_M,
) -> tuple[tuple[np.ndarray, np.ndarray], ...]:
    """Merge shared edges of nearby-height roof surfaces into boundary rings."""
    if elevation_tolerance_m < 0.0:
        return tuple(surfaces)
    groups: list[list[tuple[np.ndarray, np.ndarray]]] = []
    group_elevations: list[float] = []
    ordered = sorted(
        surfaces,
        key=lambda item: float(np.mean(item[1])),
    )
    for ring_xy, elevations in ordered:
        mean_elevation = float(np.mean(elevations))
        matching_group = next(
            (
                index
                for index, group_elevation in enumerate(group_elevations)
                if abs(mean_elevation - group_elevation)
                <= elevation_tolerance_m
            ),
            None,
        )
        if matching_group is None:
            groups.append([(ring_xy, elevations)])
            group_elevations.append(mean_elevation)
        else:
            groups[matching_group].append((ring_xy, elevations))

    merged: list[tuple[np.ndarray, np.ndarray]] = []
    for group in groups:
        merged.extend(_merge_roof_surface_group_edges(group))
    return tuple(merged)


def _merge_roof_surface_group_edges(
    surfaces: Sequence[tuple[np.ndarray, np.ndarray]],
) -> tuple[tuple[np.ndarray, np.ndarray], ...]:
    if len(surfaces) <= 1:
        return tuple(surfaces)

    nodes: dict[tuple[int, int], int] = {}
    node_points: list[np.ndarray] = []
    edge_records: dict[tuple[int, int], list[tuple[int, int]]] = {}

    def node_id(point: np.ndarray, elevation: float) -> int:
        key = (
            int(round(float(point[0]) / ROOF_SURFACE_EDGE_SNAP_M)),
            int(round(float(point[1]) / ROOF_SURFACE_EDGE_SNAP_M)),
        )
        existing = nodes.get(key)
        if existing is not None:
            return existing
        nodes[key] = len(node_points)
        node_points.append(
            np.array((point[0], point[1], elevation), dtype=np.float64)
        )
        return len(node_points) - 1

    for ring_xy, elevations in surfaces:
        for index, (start, end) in enumerate(zip(ring_xy[:-1], ring_xy[1:])):
            start_id = node_id(start, float(elevations[index]))
            end_id = node_id(end, float(elevations[index + 1]))
            if start_id == end_id:
                continue
            edge_key = tuple(sorted((start_id, end_id)))
            edge_records.setdefault(edge_key, []).append((start_id, end_id))

    boundary_edges = [
        record[0]
        for record in edge_records.values()
        if len(record) % 2 == 1
    ]
    if not boundary_edges:
        return ()

    adjacency: dict[int, list[int]] = {}
    for edge_index, (start_id, end_id) in enumerate(boundary_edges):
        adjacency.setdefault(start_id, []).append(edge_index)
        adjacency.setdefault(end_id, []).append(edge_index)

    unused = set(range(len(boundary_edges)))
    rings: list[tuple[np.ndarray, np.ndarray]] = []
    while unused:
        edge_index = min(unused)
        start_id, current_id = boundary_edges[edge_index]
        path = [start_id, current_id]
        unused.remove(edge_index)
        while current_id != start_id:
            next_edges = [
                candidate
                for candidate in adjacency.get(current_id, [])
                if candidate in unused
            ]
            if not next_edges:
                break
            next_edge = next_edges[0]
            unused.remove(next_edge)
            left, right = boundary_edges[next_edge]
            current_id = right if left == current_id else left
            path.append(current_id)
        if path[-1] != start_id or len(path) < 4:
            continue
        points = np.vstack([node_points[index] for index in path])
        rings.append((points[:, :2], points[:, 2]))

    if len(rings) <= 1:
        return tuple(rings)
    outer_rings = []
    for index, (ring_xy, elevations) in enumerate(rings):
        if any(
            index != other_index
            and abs(_ring_area_xy(other_ring_xy)) > abs(_ring_area_xy(ring_xy))
            and _ring_contains_ring_xy(ring_xy, other_ring_xy)
            for other_index, (other_ring_xy, _other_elevations) in enumerate(rings)
        ):
            continue
        outer_rings.append((ring_xy, elevations))
    return tuple(outer_rings)


def _roof_surface_elevation_tolerance(distance_m: float) -> float:
    distance_m = max(0.0, float(distance_m))
    if distance_m >= ROOF_SURFACE_TOLERANCE_FAR_DISTANCE_M:
        return ROOF_SURFACE_FAR_ELEVATION_TOLERANCE_M
    return distance_m / 1000.0 - 1.0


def _select_outer_roof_surfaces(
    surfaces: Sequence[tuple[tuple[float, float, float], ...]],
    *,
    transformer,
) -> tuple[tuple[tuple[float, float, float], ...], ...]:
    """Drop same-height roof rings fully contained by a larger ring."""
    projected: list[tuple[int, np.ndarray, float, float]] = []
    for index, surface in enumerate(surfaces):
        if len(surface) < 4:
            continue
        ring_xy = project_ring_xy(
            tuple((lon, lat) for lon, lat, _elevation in surface),
            transformer,
        )
        area_m2 = abs(_ring_area_xy(ring_xy))
        if area_m2 <= 1e-6:
            continue
        mean_elevation_m = float(
            np.mean([elevation for _lon, _lat, elevation in surface])
        )
        projected.append((index, ring_xy, area_m2, mean_elevation_m))

    kept: list[tuple[int, np.ndarray, float, float]] = []
    for candidate in sorted(projected, key=lambda item: item[2], reverse=True):
        _index, ring_xy, area_m2, elevation_m = candidate
        if any(
            abs(elevation_m - outer_elevation_m)
            <= ROOF_SURFACE_ELEVATION_TOLERANCE_M
            and area_m2 < outer_area_m2
            and _ring_contains_ring_xy(ring_xy, outer_ring_xy)
            for _outer_index, outer_ring_xy, outer_area_m2, outer_elevation_m in kept
        ):
            continue
        kept.append(candidate)

    kept.sort(key=lambda item: item[0])
    return tuple(surfaces[index] for index, _ring, _area, _elevation in kept)


def _ring_area_xy(ring_xy: np.ndarray) -> float:
    if ring_xy.ndim != 2 or ring_xy.shape[0] < 4:
        return 0.0
    return 0.5 * float(
        np.sum(ring_xy[:-1, 0] * ring_xy[1:, 1])
        - np.sum(ring_xy[1:, 0] * ring_xy[:-1, 1])
    )


def _ring_contains_ring_xy(inner_ring_xy: np.ndarray, outer_ring_xy: np.ndarray) -> bool:
    if inner_ring_xy.ndim != 2 or outer_ring_xy.ndim != 2:
        return False
    if inner_ring_xy.shape[0] < 4 or outer_ring_xy.shape[0] < 4:
        return False
    outer_min = np.min(outer_ring_xy, axis=0)
    outer_max = np.max(outer_ring_xy, axis=0)
    inner_points = inner_ring_xy[:-1]
    if np.any(inner_points < outer_min) or np.any(inner_points > outer_max):
        return False
    return all(
        _point_in_or_on_ring_xy(float(point[0]), float(point[1]), outer_ring_xy)
        for point in inner_points
    )


def _point_in_or_on_ring_xy(x: float, y: float, ring_xy: np.ndarray) -> bool:
    inside = False
    for start, end in zip(ring_xy[:-1], ring_xy[1:]):
        x0, y0 = float(start[0]), float(start[1])
        x1, y1 = float(end[0]), float(end[1])
        cross = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0)
        if abs(cross) <= 1e-7 and min(x0, x1) - 1e-7 <= x <= max(x0, x1) + 1e-7 and min(y0, y1) - 1e-7 <= y <= max(y0, y1) + 1e-7:
            return True
        if (y0 > y) != (y1 > y):
            crossing_x = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if x < crossing_x:
                inside = not inside
    return inside


def _sample_surface_points_xy(
    ring_xy: np.ndarray,
    elevations: np.ndarray,
    *,
    sample_step_m: float,
) -> np.ndarray:
    if ring_xy.ndim != 2 or ring_xy.shape[0] < 2:
        return np.empty((0, 3), dtype=np.float64)
    samples: list[np.ndarray] = []
    for index, (start_xy, end_xy) in enumerate(zip(ring_xy[:-1], ring_xy[1:])):
        delta = end_xy - start_xy
        length = float(np.hypot(delta[0], delta[1]))
        count = max(2, int(np.ceil(length / max(sample_step_m, 0.1))) + 1)
        t = np.linspace(0.0, 1.0, num=count, dtype=np.float64)
        segment = np.column_stack(
            (
                start_xy[0] + (end_xy[0] - start_xy[0]) * t,
                start_xy[1] + (end_xy[1] - start_xy[1]) * t,
                elevations[index] + (elevations[index + 1] - elevations[index]) * t,
            )
        )
        if samples:
            segment = segment[1:]
        samples.append(segment)
    return np.vstack(samples) if samples else np.empty((0, 3), dtype=np.float64)


def _emit_roof_polygon_outline(
    *,
    building: BuildingFootprint,
    ring_xy: np.ndarray,
    observer_elevation_m: float,
    building_distance_m: float,
    min_distance_m: float,
    radius_m: float,
    edge_sample_step_m: float,
    view_center: tuple[float, float] | None,
    edge_fov_deg: float,
    outlines: list[UrbanOutlinePolyline],
) -> None:
    sampled_points = sample_ring_points_xy(ring_xy, sample_step_m=edge_sample_step_m)
    if sampled_points.size == 0:
        return
    distances = np_hypot_xy(sampled_points)
    valid = (distances > max(0.1, min_distance_m)) & (distances <= radius_m)
    if not valid.any():
        return
    azimuth_deg = (np_degrees_arctan2(sampled_points[:, 0], sampled_points[:, 1]) + 360.0) % 360.0
    building_top_elevation_m = float(building.ground_elevation_m) + float(building.height_m)
    altitude_deg = np_degrees_arctan2_scalar(
        building_top_elevation_m - observer_elevation_m,
        distances,
    )
    for run in iter_true_runs(valid):
        run_points: list[UrbanPolylinePoint] = []
        for az, alt in zip(azimuth_deg[run], altitude_deg[run]):
            run_points.append(UrbanPolylinePoint(azimuth_deg=float(az), altitude_deg=float(alt)))
        run_points = _maybe_linearize_run_points(
            run_points,
            view_center=view_center,
            edge_fov_deg=edge_fov_deg,
        )
        if len(run_points) >= 2:
            outlines.append(
                UrbanOutlinePolyline(
                    height_m=float(building.height_m),
                    distance_km=float(building_distance_m / 1000.0),
                    points=tuple(run_points),
                )
            )


def _should_skip_hole_ring(ring_xy: np.ndarray, *, building_distance_m: float) -> bool:
    if building_distance_m < HOLE_RING_SUPPRESSION_MIN_DISTANCE_M:
        return False
    if ring_xy.ndim != 2 or ring_xy.shape[0] == 0:
        return False
    span_x_m = float(np.ptp(ring_xy[:, 0]))
    span_y_m = float(np.ptp(ring_xy[:, 1]))
    return max(span_x_m, span_y_m) <= HOLE_RING_SUPPRESSION_MAX_SPAN_M


def _maybe_linearize_run_points(
    run_points: list[UrbanPolylinePoint],
    *,
    view_center: tuple[float, float] | None,
    edge_fov_deg: float,
) -> list[UrbanPolylinePoint]:
    if view_center is None or len(run_points) < 2:
        return run_points

    alt = np.array([point.altitude_deg for point in run_points], dtype=np.float64)
    if float(np.max(np.abs(alt))) > VERTICAL_THIN_RUN_MAX_ABS_ALTITUDE_DEG:
        return run_points
    az = np.array([point.azimuth_deg for point in run_points], dtype=np.float64)
    ny = np.array(
        [
            altaz_to_normalized_xy(float(alt_deg), float(az_deg), view_center, edge_fov_deg=edge_fov_deg)[1]
            for alt_deg, az_deg in zip(alt, az, strict=True)
        ],
        dtype=np.float64,
    )
    if float(np.ptp(ny)) > VERTICAL_THIN_RUN_MAX_NORMALIZED_HEIGHT:
        return run_points

    start_index = int(np.argmin(ny))
    end_index = int(np.argmax(ny))
    if start_index == end_index:
        start_index = 0
        end_index = len(run_points) - 1

    if start_index == end_index:
        return run_points

    return [
        run_points[start_index],
        run_points[end_index],
    ]


def _urban_outline_ring_candidate_score(
    *,
    building: BuildingFootprint,
    ring_xy: np.ndarray,
    building_distance_m: float,
) -> float:
    width_deg = _urban_outline_ring_width_angle_deg(
        ring_xy,
        building_distance_m=building_distance_m,
    )
    height_deg = _urban_outline_height_angle_deg(building.height_m, building_distance_m)
    return max(width_deg, height_deg)


def _urban_outline_ring_width_angle_deg(
    ring_xy: np.ndarray,
    *,
    building_distance_m: float,
) -> float:
    if ring_xy.ndim != 2 or ring_xy.shape[0] < 2:
        return 0.0
    if np.array_equal(ring_xy[0], ring_xy[-1]):
        ring_xy = ring_xy[:-1]
    if ring_xy.shape[0] < 2:
        return 0.0
    center_xy = np.mean(ring_xy, axis=0)
    center_distance_m = float(np.hypot(center_xy[0], center_xy[1]))
    if center_distance_m <= 1e-6:
        width_m = float(max(np.ptp(ring_xy[:, 0]), np.ptp(ring_xy[:, 1])))
        return float(np.degrees(np.arctan2(width_m, max(1.0, float(building_distance_m)))))

    forward_xy = center_xy / center_distance_m
    side_xy = np.array((-forward_xy[1], forward_xy[0]), dtype=np.float64)
    projected = ring_xy @ side_xy
    width_m = float(np.ptp(projected))
    return float(np.degrees(np.arctan2(width_m, max(1.0, float(building_distance_m)))))


def _urban_outline_height_angle_deg(height_m: float, distance_km: float) -> float:
    if height_m <= 0.0:
        return 0.0
    distance_m = max(1.0, float(distance_km) * 1000.0)
    return float(np.degrees(np.arctan2(float(height_m), distance_m)))


def np_hypot_xy(points_xy):
    return np.hypot(points_xy[:, 0], points_xy[:, 1])


def np_degrees_arctan2(x, y):
    return np.degrees(np.arctan2(x, y))


def np_degrees_arctan2_scalar(delta_height_m, distances_m):
    return np.degrees(np.arctan2(delta_height_m, distances_m))


def ring_contains_origin_xy(ring_xy: np.ndarray) -> bool:
    if ring_xy.ndim != 2 or ring_xy.shape[0] < 4:
        return False
    x = ring_xy[:, 0]
    y = ring_xy[:, 1]

    inside = False
    for x0, y0, x1, y1 in zip(x[:-1], y[:-1], x[1:], y[1:]):
        intersects = ((y0 > 0.0) != (y1 > 0.0)) and (0.0 < (x1 - x0) * (-y0) / (y1 - y0 + 1e-12) + x0)
        if intersects:
            inside = not inside
    return inside


def collect_buildings_for_tower(
    derived_dir: Path,
    tower: TowerViewpoint,
    *,
    radius_km: float,
    min_building_height_m: float,
) -> tuple[BuildingFootprint, ...]:
    envelopes = select_derived_tile_envelopes(
        derived_dir,
        observer_lat_deg=tower.latitude_deg,
        observer_lon_deg=tower.longitude_deg,
        radius_km=radius_km,
    )
    buildings: list[BuildingFootprint] = []
    for envelope in envelopes:
        buildings.extend(
            parse_derived_tile_buildings(
                envelope.path,
                min_building_height_m=min_building_height_m,
            )
        )
    return tuple(buildings)


def build_output_payload(results: Sequence[UrbanOutlineResult]) -> dict[str, object]:
    towers_payload: list[dict[str, object]] = []
    for result in results:
        towers_payload.append(
            {
                "tower": {
                    "id": result.tower.id,
                    "name": result.tower.name,
                    "latitude_deg": result.tower.latitude_deg,
                    "longitude_deg": result.tower.longitude_deg,
                    "viewpoint_height_m": result.tower.viewpoint_height_m,
                },
                "buildings_considered": result.buildings_considered,
                "outlines_emitted": result.outlines_emitted,
                "outlines": [
                    [
                        {
                            "azimuth_deg": point.azimuth_deg,
                            "altitude_deg": point.altitude_deg,
                        }
                        for point in outline
                    ]
                    for outline in result.outlines
                ],
            }
        )
    return {
        "schema_version": 1,
        "towers": towers_payload,
    }


def write_output_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    towers = select_towers(
        tower_queries=args.tower,
        all_covered_towers=bool(args.all_covered_towers),
        derived_dir=args.derived_dir,
        radius_km=float(args.radius_km),
    )
    results: list[UrbanOutlineResult] = []
    for tower in towers:
        buildings = collect_buildings_for_tower(
            args.derived_dir,
            tower,
            radius_km=float(args.radius_km),
            min_building_height_m=float(args.min_building_height_m),
        )
        results.append(
            compute_urban_outlines(
                tower,
                buildings,
                radius_km=float(args.radius_km),
                edge_sample_step_m=float(args.edge_sample_step_m),
            )
        )

    payload = build_output_payload(results)
    write_output_json(args.output_json, payload)
    print(f"[ok] wrote urban outline json: {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
