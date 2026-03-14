#!/usr/bin/env python3
"""Build tile_index.json for a directory of derived building tiles."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from math import inf
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build tile_index.json for a directory of derived building tiles."
    )
    parser.add_argument(
        "--derived-dir",
        type=Path,
        required=True,
        help="Directory containing derived tile JSON files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path for tile_index.json. Defaults to <derived-dir>/tile_index.json.",
    )
    return parser


def detect_city_metadata(derived_dir: Path) -> tuple[str | None, str | None]:
    parent_name = derived_dir.parent.name
    city_code = parent_name.split("_", 1)[0] if "_" in parent_name else None
    city_name = parent_name.split("_", 1)[1] if "_" in parent_name else parent_name
    return city_code, city_name or None


def build_tile_index_payload(derived_dir: Path) -> dict[str, object]:
    city_code, city_name = detect_city_metadata(derived_dir)
    tiles: list[dict[str, object]] = []
    min_lat = inf
    min_lon = inf
    max_lat = -inf
    max_lon = -inf
    min_height_m: float | None = None

    tile_paths = sorted(
        path for path in derived_dir.glob("*.json") if path.name != "tile_index.json"
    )
    if not tile_paths:
        raise ValueError(f"No derived tile JSON files found in {derived_dir}")

    for path in tile_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        bbox = _parse_tile_bbox(payload)
        if bbox is None:
            raise ValueError(f"Derived tile is missing tile bbox: {path}")
        filters = payload.get("filters")
        if isinstance(filters, dict) and min_height_m is None:
            raw_min_height_m = filters.get("min_height_m")
            if isinstance(raw_min_height_m, (int, float)):
                min_height_m = float(raw_min_height_m)
        buildings = payload.get("buildings")
        building_count = len(buildings) if isinstance(buildings, list) else 0

        min_lat = min(min_lat, bbox["min_lat"])
        min_lon = min(min_lon, bbox["min_lon"])
        max_lat = max(max_lat, bbox["max_lat"])
        max_lon = max(max_lon, bbox["max_lon"])

        tiles.append(
            {
                "id": path.stem,
                "path": path.name,
                "bbox": bbox,
                "building_count": building_count,
            }
        )

    return {
        "schema_version": 1,
        "city_code": city_code,
        "city_name": city_name,
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "bbox": {
            "min_lat": min_lat,
            "min_lon": min_lon,
            "max_lat": max_lat,
            "max_lon": max_lon,
        },
        "tile_count": len(tiles),
        "min_height_m": min_height_m,
        "tiles": tiles,
    }


def write_tile_index(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _parse_tile_bbox(payload: object) -> dict[str, float] | None:
    if not isinstance(payload, dict):
        return None
    tile = payload.get("tile")
    if not isinstance(tile, dict):
        return None
    bbox = tile.get("bbox")
    if not isinstance(bbox, dict):
        return None
    try:
        return {
            "min_lat": float(bbox["min_lat"]),
            "min_lon": float(bbox["min_lon"]),
            "max_lat": float(bbox["max_lat"]),
            "max_lon": float(bbox["max_lon"]),
        }
    except (KeyError, TypeError, ValueError):
        return None


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    derived_dir = args.derived_dir
    output = args.output or (derived_dir / "tile_index.json")
    payload = build_tile_index_payload(derived_dir)
    write_tile_index(output, payload)
    print(f"[ok] tile-index: {output}  tiles={payload['tile_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
