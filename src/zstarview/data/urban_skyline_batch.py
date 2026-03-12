#!/usr/bin/env python3
"""Batch-build urban skyline profiles for bundled viewpoints covered by PLATEAU."""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.data.urban_skyline_demo import (  # noqa: E402
    DEFAULT_CUMULATIVE_RADII_KM,
    DEFAULT_MIN_BUILDING_HEIGHT_M,
    DEFAULT_RADIUS_BAND_WIDTH_M,
    is_japan_tower,
    sanitize_slug,
    write_preview_png,
    write_profiles_json,
)
from zstarview.data.urban_skyline_from_citygml import (  # noqa: E402
    TileEnvelope,
    combine_tile_results,
    compute_tile_skylines,
    envelope_min_distance_km,
    load_tile_envelope,
    normalize_cumulative_radii_km,
    select_tile_envelopes,
)
from zstarview.tower_viewpoints import (  # noqa: E402
    TowerViewpoint,
    load_tower_viewpoints,
    resolve_tower_viewpoint,
)


@dataclass(frozen=True)
class CityGmlCoverage:
    bldg_dir: Path
    min_lat_deg: float
    min_lon_deg: float
    max_lat_deg: float
    max_lon_deg: float

    def as_tile_envelope(self) -> TileEnvelope:
        return TileEnvelope(
            path=self.bldg_dir,
            min_lat_deg=self.min_lat_deg,
            min_lon_deg=self.min_lon_deg,
            max_lat_deg=self.max_lat_deg,
            max_lon_deg=self.max_lon_deg,
        )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate urban skyline profiles for every bundled viewpoint covered "
            "by one or more PLATEAU CityGML building directories."
        )
    )
    parser.add_argument(
        "--citygml-root",
        action="append",
        default=[],
        help="Root directory to scan for */udx/bldg directories. May be repeated.",
    )
    parser.add_argument(
        "--citygml-dir",
        action="append",
        default=[],
        help="Explicit PLATEAU building directory such as raw-data/.../udx/bldg. May be repeated.",
    )
    parser.add_argument(
        "--tower",
        action="append",
        default=[],
        help="Specific bundled viewpoint name or wikidata:Q... identifier. May be repeated.",
    )
    parser.add_argument(
        "--all-covered-towers",
        action="store_true",
        help="Process every bundled Japan viewpoint covered by the discovered PLATEAU directories.",
    )
    parser.add_argument(
        "--radius-km",
        type=float,
        default=max(DEFAULT_CUMULATIVE_RADII_KM),
        help="Maximum building search radius around each viewpoint (default: 5.76650390625).",
    )
    parser.add_argument(
        "--cumulative-radius-km",
        action="append",
        default=[],
        help=(
            "Skyline scan radius in km. May be repeated. Defaults to: "
            + ", ".join(str(value) for value in DEFAULT_CUMULATIVE_RADII_KM)
        ),
    )
    parser.add_argument(
        "--radius-band-width-m",
        type=float,
        default=DEFAULT_RADIUS_BAND_WIDTH_M,
        help="Radial scan band width in meters for each skyline layer (default: 90.0).",
    )
    parser.add_argument(
        "--azimuth-step",
        type=float,
        default=0.1,
        help="Azimuth sampling step in degrees (default: 0.1).",
    )
    parser.add_argument(
        "--edge-sample-step-m",
        type=float,
        default=5.0,
        help="Approximate spacing for sampling building edges in meters (default: 5.0).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, ((__import__("os").cpu_count() or 1) // 2)),
        help="Number of worker processes for per-tile processing (default: half of CPU cores).",
    )
    parser.add_argument(
        "--min-building-height-m",
        type=float,
        default=DEFAULT_MIN_BUILDING_HEIGHT_M,
        help="Ignore buildings lower than this height in meters (default: 50.0).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "src" / "zstarview" / "data" / "viewpoints" / "urban_skyline",
        help="Directory where PNG previews and JSON outputs will be written.",
    )
    parser.add_argument(
        "--write-json",
        action="store_true",
        help="Write a consolidated skyline JSON keyed by viewpoint id.",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=None,
        help="Override consolidated JSON output path.",
    )
    parser.add_argument(
        "--list-covered-towers",
        action="store_true",
        help="List discovered covered bundled viewpoints and exit.",
    )
    return parser


def discover_citygml_dirs(*, roots: Sequence[Path], explicit_dirs: Sequence[Path]) -> tuple[Path, ...]:
    found: set[Path] = set()
    for path in explicit_dirs:
        if path.is_dir():
            found.add(path.resolve())
    for root in roots:
        if not root.is_dir():
            continue
        for path in root.rglob("udx/bldg"):
            if path.is_dir():
                found.add(path.resolve())
    return tuple(sorted(found))


def load_citygml_coverage(bldg_dir: Path) -> CityGmlCoverage | None:
    min_lat_deg = math.inf
    min_lon_deg = math.inf
    max_lat_deg = -math.inf
    max_lon_deg = -math.inf
    seen = 0
    for path in bldg_dir.glob("*.gml"):
        envelope = load_tile_envelope(path)
        if envelope is None:
            continue
        seen += 1
        min_lat_deg = min(min_lat_deg, envelope.min_lat_deg)
        min_lon_deg = min(min_lon_deg, envelope.min_lon_deg)
        max_lat_deg = max(max_lat_deg, envelope.max_lat_deg)
        max_lon_deg = max(max_lon_deg, envelope.max_lon_deg)
    if seen == 0:
        return None
    return CityGmlCoverage(
        bldg_dir=bldg_dir,
        min_lat_deg=min_lat_deg,
        min_lon_deg=min_lon_deg,
        max_lat_deg=max_lat_deg,
        max_lon_deg=max_lon_deg,
    )


def load_coverages(bldg_dirs: Sequence[Path]) -> tuple[CityGmlCoverage, ...]:
    coverages = [coverage for coverage in (load_citygml_coverage(path) for path in bldg_dirs) if coverage is not None]
    return tuple(sorted(coverages, key=lambda coverage: str(coverage.bldg_dir)))


def select_candidate_coverages(
    tower: TowerViewpoint,
    coverages: Sequence[CityGmlCoverage],
    *,
    radius_km: float,
) -> tuple[CityGmlCoverage, ...]:
    matched = []
    for coverage in coverages:
        if envelope_min_distance_km(
            coverage.as_tile_envelope(),
            observer_lat_deg=tower.latitude_deg,
            observer_lon_deg=tower.longitude_deg,
        ) <= radius_km:
            matched.append(coverage)
    return tuple(matched)


def select_towers(
    *,
    tower_queries: Sequence[str],
    all_covered_towers: bool,
    coverages: Sequence[CityGmlCoverage],
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

    selected = []
    search_radius_km = float(radius_km)
    for tower in load_tower_viewpoints():
        if not is_japan_tower(tower):
            continue
        if select_candidate_coverages(tower, coverages, radius_km=search_radius_km):
            selected.append(tower)
    return tuple(selected)


def compute_tower_skyline(
    tower: TowerViewpoint,
    *,
    coverages: Sequence[CityGmlCoverage],
    cumulative_radii_km: Sequence[float],
    radius_km: float,
    radius_band_width_m: float,
    azimuth_step_deg: float,
    edge_sample_step_m: float,
    workers: int,
    min_building_height_m: float,
):
    selected_coverages = select_candidate_coverages(tower, coverages, radius_km=radius_km)
    envelopes: list[TileEnvelope] = []
    for coverage in selected_coverages:
        try:
            matched = select_tile_envelopes(
                coverage.bldg_dir,
                observer_lat_deg=tower.latitude_deg,
                observer_lon_deg=tower.longitude_deg,
                radius_km=radius_km,
            )
        except ValueError:
            continue
        envelopes.extend(matched)
    if not envelopes:
        return None
    tile_results = compute_tile_skylines(
        tuple(envelopes),
        tower=tower,
        cumulative_radii_km=cumulative_radii_km,
        radius_band_width_m=radius_band_width_m,
        azimuth_step_deg=azimuth_step_deg,
        edge_sample_step_m=edge_sample_step_m,
        workers=workers,
        min_building_height_m=min_building_height_m,
    )
    if not tile_results:
        return None
    return combine_tile_results(
        tower,
        tile_results,
        radii_km=cumulative_radii_km,
        azimuth_step_deg=azimuth_step_deg,
    )


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        bldg_dirs = discover_citygml_dirs(
            roots=tuple(Path(text) for text in args.citygml_root),
            explicit_dirs=tuple(Path(text) for text in args.citygml_dir),
        )
        if not bldg_dirs:
            raise ValueError("No PLATEAU building directories were found.")
        print(f"[batch] discovered-bldg-dirs={len(bldg_dirs)}")
        coverages = load_coverages(bldg_dirs)
        if not coverages:
            raise ValueError("No usable CityGML building coverages were found.")
        print(f"[batch] usable-coverages={len(coverages)}")

        search_radius_km = float(args.radius_km) + (float(args.radius_band_width_m) / 1000.0)
        towers = select_towers(
            tower_queries=tuple(args.tower),
            all_covered_towers=bool(args.all_covered_towers),
            coverages=coverages,
            radius_km=search_radius_km,
        )
        print(f"[batch] selected-viewpoints={len(towers)}")
        if args.list_covered_towers:
            for tower in towers:
                print(f"{tower.name}\t{tower.id}")
            return 0

        cumulative_radii_km = normalize_cumulative_radii_km(
            [float(value) for value in args.cumulative_radius_km],
            max_radius_km=float(args.radius_km),
        )
        results = []
        for index, tower in enumerate(towers, start=1):
            selected_coverages = select_candidate_coverages(
                tower,
                coverages,
                radius_km=search_radius_km,
            )
            print(
                f"[batch] start {index}/{len(towers)}  "
                f"{tower.name}  coverages={len(selected_coverages)}  "
                f"observer_height={tower.observer_height_m:.1f}m"
            )
            radius_results = compute_tower_skyline(
                tower,
                coverages=coverages,
                cumulative_radii_km=cumulative_radii_km,
                radius_km=search_radius_km,
                radius_band_width_m=float(args.radius_band_width_m),
                azimuth_step_deg=float(args.azimuth_step),
                edge_sample_step_m=float(args.edge_sample_step_m),
                workers=int(args.workers),
                min_building_height_m=float(args.min_building_height_m),
            )
            if radius_results is None:
                print(f"[batch] skip {index}/{len(towers)}  {tower.name}: no usable building tiles")
                continue
            results.append((tower, radius_results))
            result = radius_results[-1].result
            stem = sanitize_slug(str(tower.meta.get("slug") or tower.name))
            png_path = args.output_dir / f"{stem}_urban.png"
            write_preview_png(png_path, result)
            print(
                f"[batch] done {index}/{len(towers)}  {tower.name}: {png_path}  "
                f"buildings={result.buildings_considered}/{result.buildings_contributing}  "
                f"peak={result.peak_altitude_deg:.2f}deg@{result.peak_azimuth_deg:.1f}"
            )

        if args.write_json and results:
            json_output = args.json_output or (args.output_dir / "urban_skyline_profiles.json")
            write_profiles_json(json_output, tuple(results))
            print(f"[batch] wrote-json: {json_output}  viewpoints={len(results)}")
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
