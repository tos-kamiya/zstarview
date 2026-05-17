#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Find coarse-grid intersection points whose surrounding water tiles are all .tif.

The script checks the 500m coarse boundary grid, then verifies that the four
neighboring tiles around each intersection are full .tif files in the 125m,
250m, and 500m water tile sets.

Examples:
  python dev-samples/find_water_tile_intersections.py
  python dev-samples/find_water_tile_intersections.py --json
  python dev-samples/find_water_tile_intersections.py --tile-root-500m /path/to/water_tiles_500m
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.water_tile_intersections import (  # noqa: E402
    DEFAULT_WATER_TILES_ROOT_125M,
    DEFAULT_WATER_TILES_ROOT_250M,
    DEFAULT_WATER_TILES_ROOT_500M,
    find_common_boundary_intersections,
    load_tile_grid_spec,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Find coarse boundary intersections whose four neighboring tiles are all .tif "
            "in the 125m, 250m, and 500m water tile sets."
        )
    )
    parser.add_argument(
        "--tile-root-125m",
        type=Path,
        default=DEFAULT_WATER_TILES_ROOT_125M,
        help=f"125m tile root (default: {DEFAULT_WATER_TILES_ROOT_125M})",
    )
    parser.add_argument(
        "--tile-root-250m",
        type=Path,
        default=DEFAULT_WATER_TILES_ROOT_250M,
        help=f"250m tile root (default: {DEFAULT_WATER_TILES_ROOT_250M})",
    )
    parser.add_argument(
        "--tile-root-500m",
        type=Path,
        default=DEFAULT_WATER_TILES_ROOT_500M,
        help=f"500m tile root (default: {DEFAULT_WATER_TILES_ROOT_500M})",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit JSON instead of a human-readable table.",
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()

    specs = [
        load_tile_grid_spec("125m", Path(args.tile_root_125m)),
        load_tile_grid_spec("250m", Path(args.tile_root_250m)),
        load_tile_grid_spec("500m", Path(args.tile_root_500m)),
    ]
    matches = find_common_boundary_intersections(specs)

    if args.json:
        payload = [
            {
                "coarse_row_index": match.coarse_row_index,
                "coarse_col_index": match.coarse_col_index,
                "latitude_deg": match.latitude_deg,
                "longitude_deg": match.longitude_deg,
                "tile_paths_by_root": {
                    label: [str(path) for path in quad_paths]
                    for label, quad_paths in match.tile_paths_by_root.items()
                },
            }
            for match in matches
        ]
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0

    print(f"Found {len(matches)} shared boundary intersection(s).")
    for match in matches:
        print(
            f"- lat={match.latitude_deg:.1f}, lon={match.longitude_deg:.1f} "
            f"(coarse_row={match.coarse_row_index}, coarse_col={match.coarse_col_index})"
        )
        for label, quad_paths in match.tile_paths_by_root.items():
            names = ", ".join(path.name for path in quad_paths)
            print(f"  {label}: {names}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
