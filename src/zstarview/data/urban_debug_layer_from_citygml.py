#!/usr/bin/env python3
"""Build a Tokyo23-only debug urban layer from PLATEAU CityGML tiles.

This tool samples the tops of building footprints within 3 km of a bundled
viewpoint and stores those projected roof-edge polylines as a debug overlay.
Unlike the urban skyline generator, it does not collapse buildings into a
single max-altitude profile.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
DATA_ROOT = SRC_ROOT / "zstarview" / "data"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.data.urban_skyline_demo import (  # noqa: E402
    BuildingFootprint,
    bbox_min_distance_m,
    is_japan_tower,
    iter_true_runs,
    make_local_transformer,
    project_ring_xy,
    sample_ring_points_xy,
)
from zstarview.data.urban_skyline_from_citygml import (  # noqa: E402
    parse_citygml_buildings,
    select_tile_envelopes,
)
from zstarview.paths import URBAN_DEBUG_LAYER_FILE  # noqa: E402
from zstarview.tower_viewpoints import (  # noqa: E402
    TowerViewpoint,
    load_tower_viewpoints,
    resolve_tower_viewpoint,
)

DEFAULT_RADIUS_KM = 3.0
DEFAULT_MIN_BUILDING_HEIGHT_M = 40.0
DEFAULT_EDGE_SAMPLE_STEP_M = 10.0


@dataclass(frozen=True)
class DebugPolylinePoint:
    azimuth_deg: float
    altitude_deg: float


@dataclass(frozen=True)
class DebugOutlineResult:
    tower: TowerViewpoint
    outlines: tuple[tuple[DebugPolylinePoint, ...], ...]
    buildings_considered: int
    outlines_emitted: int


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate Tokyo23 debug urban outline overlays from PLATEAU CityGML building tiles."
    )
    parser.add_argument(
        "--citygml-dir",
        type=Path,
        required=True,
        help="Directory containing Tokyo23 PLATEAU building GML files such as raw-data/.../udx/bldg.",
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
        help="Process every bundled Japan viewpoint whose 3km search circle intersects the provided Tokyo23 tiles.",
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
        default=Path(URBAN_DEBUG_LAYER_FILE),
        help="Output JSON path.",
    )
    return parser


def select_towers(
    *,
    tower_queries: Sequence[str],
    all_covered_towers: bool,
    citygml_dir: Path,
    radius_km: float,
) -> tuple[TowerViewpoint, ...]:
    if tower_queries:
        selected: list[TowerViewpoint] = []
        for query in tower_queries:
            tower = resolve_tower_viewpoint(query)
            if tower is None:
                raise ValueError(f"Viewpoint not found: {query}")
            if not is_japan_tower(tower):
                continue
            selected.append(tower)
        return tuple(selected)

    if not all_covered_towers:
        raise ValueError("Specify --tower or --all-covered-towers.")

    selected: list[TowerViewpoint] = []
    for tower in load_tower_viewpoints():
        if not is_japan_tower(tower):
            continue
        try:
            select_tile_envelopes(
                citygml_dir,
                observer_lat_deg=tower.latitude_deg,
                observer_lon_deg=tower.longitude_deg,
                radius_km=radius_km,
            )
        except ValueError:
            continue
        selected.append(tower)
    return tuple(selected)


def compute_debug_outlines(
    tower: TowerViewpoint,
    buildings: Sequence[BuildingFootprint],
    *,
    radius_km: float,
    edge_sample_step_m: float,
) -> DebugOutlineResult:
    if radius_km <= 0.0:
        raise ValueError("--radius-km must be positive.")
    if edge_sample_step_m <= 0.0:
        raise ValueError("--edge-sample-step-m must be positive.")

    transformer = make_local_transformer(tower)
    radius_m = radius_km * 1000.0
    observer_height_m = tower.observer_height_m
    buildings_considered = 0
    outlines: list[tuple[DebugPolylinePoint, ...]] = []

    for building in buildings:
        projected_rings = tuple(project_ring_xy(ring, transformer) for ring in building.rings_lonlat)
        if not projected_rings:
            continue
        if any(ring_contains_origin_xy(ring_xy) for ring_xy in projected_rings if ring_xy.size > 0):
            continue
        min_distance = min(bbox_min_distance_m(ring_xy) for ring_xy in projected_rings if ring_xy.size > 0)
        if min_distance > radius_m:
            continue
        buildings_considered += 1
        for ring_xy in projected_rings:
            sampled_points = sample_ring_points_xy(ring_xy, sample_step_m=edge_sample_step_m)
            if sampled_points.size == 0:
                continue
            distances = np_hypot_xy(sampled_points)
            valid = (distances > 0.1) & (distances <= radius_m)
            if not valid.any():
                continue
            azimuth_deg = (np_degrees_arctan2(sampled_points[:, 0], sampled_points[:, 1]) + 360.0) % 360.0
            altitude_deg = np_degrees_arctan2_scalar(building.height_m - observer_height_m, distances)
            for run in iter_true_runs(valid):
                run_points: list[DebugPolylinePoint] = []
                for az, alt in zip(azimuth_deg[run], altitude_deg[run]):
                    run_points.append(
                        DebugPolylinePoint(
                            azimuth_deg=float(az),
                            altitude_deg=float(alt),
                        )
                    )
                if len(run_points) >= 2:
                    outlines.append(tuple(run_points))

    return DebugOutlineResult(
        tower=tower,
        outlines=tuple(outlines),
        buildings_considered=buildings_considered,
        outlines_emitted=len(outlines),
    )


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
    if np.any(np.isclose(x, 0.0) & np.isclose(y, 0.0)):
        return True

    inside = False
    px = 0.0
    py = 0.0
    for i in range(ring_xy.shape[0] - 1):
        x1 = float(x[i])
        y1 = float(y[i])
        x2 = float(x[i + 1])
        y2 = float(y[i + 1])
        intersects = ((y1 > py) != (y2 > py)) and (
            px < ((x2 - x1) * (py - y1) / ((y2 - y1) or 1.0e-12)) + x1
        )
        if intersects:
            inside = not inside
    return inside


def load_buildings_for_tower(
    citygml_dir: Path,
    *,
    tower: TowerViewpoint,
    radius_km: float,
    min_building_height_m: float,
) -> tuple[BuildingFootprint, ...]:
    envelopes = select_tile_envelopes(
        citygml_dir,
        observer_lat_deg=tower.latitude_deg,
        observer_lon_deg=tower.longitude_deg,
        radius_km=radius_km,
    )
    buildings: list[BuildingFootprint] = []
    for envelope in envelopes:
        buildings.extend(
            parse_citygml_buildings(
                envelope.path,
                min_building_height_m=min_building_height_m,
            )
        )
    return tuple(buildings)


def write_debug_layer_json(path: Path, results: Sequence[DebugOutlineResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, object] = {}
    for result in results:
        payload[result.tower.id] = {
            "name": result.tower.name,
            "outlines": [
                [{"az": point.azimuth_deg, "alt": point.altitude_deg} for point in outline]
                for outline in result.outlines
            ],
        }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    towers = select_towers(
        tower_queries=args.tower,
        all_covered_towers=bool(args.all_covered_towers),
        citygml_dir=args.citygml_dir,
        radius_km=float(args.radius_km),
    )
    if not towers:
        raise ValueError("No bundled Tokyo23 viewpoints matched the provided inputs.")

    results: list[DebugOutlineResult] = []
    for tower in towers:
        buildings = load_buildings_for_tower(
            args.citygml_dir,
            tower=tower,
            radius_km=float(args.radius_km),
            min_building_height_m=float(args.min_building_height_m),
        )
        result = compute_debug_outlines(
            tower,
            buildings,
            radius_km=float(args.radius_km),
            edge_sample_step_m=float(args.edge_sample_step_m),
        )
        print(
            f"[ok] {tower.name}: buildings={result.buildings_considered} outlines={result.outlines_emitted}"
        )
        results.append(result)

    write_debug_layer_json(args.output_json, results)
    print(f"[ok] debug-json: {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
