#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Show how sea-surface tiles look around an observer location.

This script is intended as a general inspection tool for the water-mask tile
sets. It prints:
- per-resolution tile grid information
- tile suffix counts (.tif/.0/.1)
- the tile hit at the observer location
- a small table of nearby probe points
- an optional ray-scan summary that mirrors the water-mask sampling path

Examples:
  python dev-samples/water_tile_overview.py --lat 0 --lon 135
  python dev-samples/water_tile_overview.py --lat 0 --lon 135 --ray-scan
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import rasterio

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.terrain import WGS84_GEOD, build_ray_scan_grid  # noqa: E402
from zstarview.water_mask_interface import (  # noqa: E402
    DEFAULT_WATER_TILES_ROOT_125M,
    DEFAULT_WATER_TILES_ROOT_250M,
    DEFAULT_WATER_TILES_ROOT_500M,
    _sample_water_mask_for_lonlat_points,
    _tile_key_from_path,
    _tile_marker_value,
    _tile_paths,
)


ROOTS = (
    ("125m", DEFAULT_WATER_TILES_ROOT_125M, 16, 32),
    ("250m", DEFAULT_WATER_TILES_ROOT_250M, 8, 16),
    ("500m", DEFAULT_WATER_TILES_ROOT_500M, 4, 8),
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Inspect how sea-surface tiles appear around an observer point."
    )
    parser.add_argument("--lat", type=float, required=True, help="Observer latitude in degrees.")
    parser.add_argument("--lon", type=float, required=True, help="Observer longitude in degrees.")
    parser.add_argument(
        "--probe-distances-km",
        default="0.25,0.5,1,2,4,8,12",
        help="Comma-separated probe distances in kilometers (default: 0.25,0.5,1,2,4,8,12).",
    )
    parser.add_argument(
        "--probe-azimuths-deg",
        default="0,45,90,135,180,225,270,315",
        help="Comma-separated probe azimuths in degrees (default: cardinal and diagonal directions).",
    )
    parser.add_argument(
        "--max-distance-km",
        type=float,
        default=12.0,
        help="Maximum distance for the optional ray-scan summary (default: 12.0).",
    )
    parser.add_argument(
        "--azimuth-step-deg",
        type=float,
        default=15.0,
        help="Azimuth step for the optional ray-scan summary (default: 15.0).",
    )
    parser.add_argument(
        "--sample-step-m",
        type=float,
        default=125.0,
        help="Sample step for the optional ray-scan summary (default: 125.0).",
    )
    parser.add_argument(
        "--ray-scan",
        action="store_true",
        help="Run the ray-scan summary path in addition to point probes.",
    )
    return parser


def _parse_float_list(text: str) -> list[float]:
    values: list[float] = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        values.append(float(token))
    return values


def _load_tile_map(tile_root: Path) -> dict[tuple[int, int], Path]:
    tile_map: dict[tuple[int, int], Path] = {}
    for tile_path in _tile_paths(tile_root):
        tile_key = _tile_key_from_path(tile_path)
        if tile_key is not None:
            tile_map[tile_key] = tile_path
    return tile_map


def _tile_key_for_lonlat(lon_deg: float, lat_deg: float, *, rows: int, cols: int) -> tuple[int, int]:
    row = int((90.0 - float(lat_deg)) // (180.0 / float(rows)))
    col = int((float(lon_deg) + 180.0) // (360.0 / float(cols)))
    row = max(0, min(rows - 1, row))
    col = max(0, min(cols - 1, col))
    return row, col


def _sample_point(
    tile_map: dict[tuple[int, int], Path],
    lon_deg: float,
    lat_deg: float,
    *,
    rows: int,
    cols: int,
) -> tuple[str, str, float | None]:
    tile_key = _tile_key_for_lonlat(lon_deg, lat_deg, rows=rows, cols=cols)
    tile_path = tile_map.get(tile_key)
    if tile_path is None:
        return (f"{tile_key}", "missing", None)

    marker_value = _tile_marker_value(tile_path)
    if marker_value is not None:
        return (f"{tile_key}", tile_path.name, 1.0 if marker_value else 0.0)

    with rasterio.open(tile_path) as dataset:
        sample = list(dataset.sample([(float(lon_deg), float(lat_deg))]))
    if not sample:
        return (f"{tile_key}", tile_path.name, None)
    value = sample[0]
    if value.size == 0:
        return (f"{tile_key}", tile_path.name, None)
    return (f"{tile_key}", tile_path.name, float(value[0]))


def _suffix_counts(tile_root: Path) -> dict[str, int]:
    counts = {".tif": 0, ".0": 0, ".1": 0}
    for tile_path in _tile_paths(tile_root):
        if tile_path.suffix in counts:
            counts[tile_path.suffix] += 1
    return counts


def _print_root_summary(
    *,
    tile_root_label: str,
    tile_root: Path,
    rows: int,
    cols: int,
    tile_map: dict[tuple[int, int], Path],
    observer_lat_deg: float,
    observer_lon_deg: float,
    probe_azimuths_deg: list[float],
    probe_distances_km: list[float],
) -> None:
    counts = _suffix_counts(tile_root)
    print(
        f"[{tile_root_label}] grid={rows}x{cols} tiles={len(tile_map)} "
        f"tif={counts['.tif']} zero={counts['.0']} one={counts['.1']} root={tile_root}"
    )

    observer_key_text, observer_tile_name, observer_value = _sample_point(
        tile_map,
        float(observer_lon_deg),
        float(observer_lat_deg),
        rows=rows,
        cols=cols,
    )
    print(
        f"  observer lon={float(observer_lon_deg):11.6f} lat={float(observer_lat_deg):11.6f} "
        f"tile={observer_tile_name} key={observer_key_text} value={observer_value}"
    )

    print(f"  [{tile_root_label}] point probes")
    for az_deg in probe_azimuths_deg:
        for distance_km in probe_distances_km:
            lon_deg, lat_deg, _ = WGS84_GEOD.fwd(
                float(observer_lon_deg),
                float(observer_lat_deg),
                float(az_deg),
                float(distance_km) * 1000.0,
            )
            tile_key_text, tile_name, value = _sample_point(
                tile_map,
                lon_deg,
                lat_deg,
                rows=rows,
                cols=cols,
            )
            print(
                f"    az={az_deg:6.1f} dist={distance_km:6.2f} km "
                f"-> lon={lon_deg:11.6f} lat={lat_deg:11.6f} "
                f"tile={tile_name} key={tile_key_text} value={value}"
            )


def _print_ray_scan_summary(
    *,
    tile_root_label: str,
    tile_root: Path,
    tile_map: dict[tuple[int, int], Path],
    rows: int,
    cols: int,
    observer_lat_deg: float,
    observer_lon_deg: float,
    max_distance_km: float,
    azimuth_step_deg: float,
    sample_step_m: float,
) -> None:
    print(f"  [{tile_root_label}] ray-scan summary")
    distance_samples_m = []
    distance_m = float(sample_step_m)
    while distance_m <= float(max_distance_km) * 1000.0:
        distance_samples_m.append(distance_m)
        distance_m *= 1.15
    if not distance_samples_m:
        distance_samples_m = [float(sample_step_m)]

    ray_scan = build_ray_scan_grid(
        geod=WGS84_GEOD,
        observer_latitude_deg=float(observer_lat_deg),
        observer_longitude_deg=float(observer_lon_deg),
        azimuth_step_deg=float(azimuth_step_deg),
        distance_samples_m=np.asarray(distance_samples_m, dtype=float),
    )

    total_hits = 0
    for row_index, azimuth_deg in enumerate(ray_scan.azimuths_deg):
        row_lonlat_points = [
            (float(lon_deg), float(lat_deg))
            for lon_deg, lat_deg in zip(
                ray_scan.ray_lon_deg[row_index],
                ray_scan.ray_lat_deg[row_index],
                strict=False,
            )
        ]
        row_flags = _sample_water_mask_for_lonlat_points(
            row_lonlat_points,
            tile_root=tile_root,
        )
        row_hits = [index for index, flag in enumerate(row_flags) if flag]
        total_hits += len(row_hits)
        if not row_hits:
            continue
        first_index = row_hits[0]
        lon_deg, lat_deg = row_lonlat_points[first_index]
        tile_key_text, tile_name, value = _sample_point(
            tile_map,
            lon_deg,
            lat_deg,
            rows=rows,
            cols=cols,
        )
        print(
            f"    az={float(azimuth_deg):6.1f} hits={len(row_hits):4d}/{len(row_flags):4d} "
            f"first_hit_dist={float(ray_scan.distance_grid_m[row_index][first_index]) / 1000.0:6.2f} km "
            f"tile={tile_name} key={tile_key_text} value={value}"
        )
    print(f"    total_hits={total_hits}")


def main() -> int:
    args = _build_parser().parse_args()
    probe_distances_km = _parse_float_list(args.probe_distances_km)
    probe_azimuths_deg = _parse_float_list(args.probe_azimuths_deg)

    for tile_root_label, tile_root, rows, cols in ROOTS:
        tile_map = _load_tile_map(tile_root)
        _print_root_summary(
            tile_root_label=tile_root_label,
            tile_root=tile_root,
            rows=rows,
            cols=cols,
            tile_map=tile_map,
            observer_lat_deg=float(args.lat),
            observer_lon_deg=float(args.lon),
            probe_azimuths_deg=probe_azimuths_deg,
            probe_distances_km=probe_distances_km,
        )
        if args.ray_scan:
            _print_ray_scan_summary(
                tile_root_label=tile_root_label,
                tile_root=tile_root,
                tile_map=tile_map,
                rows=rows,
                cols=cols,
                observer_lat_deg=float(args.lat),
                observer_lon_deg=float(args.lon),
                max_distance_km=float(args.max_distance_km),
                azimuth_step_deg=float(args.azimuth_step_deg),
                sample_step_m=float(args.sample_step_m),
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
