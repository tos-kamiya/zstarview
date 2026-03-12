#!/usr/bin/env python3
"""Build tile_index.json for derived PLATEAU building tiles."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, UTC
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.data.plateau_derived_tiles import load_derived_tile_envelope  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build tile_index.json for a directory of derived PLATEAU building tiles."
    )
    parser.add_argument(
        "--derived-dir",
        type=Path,
        required=True,
        help="Directory containing derived PLATEAU tile JSON files.",
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
    min_lat = float("inf")
    min_lon = float("inf")
    max_lat = float("-inf")
    max_lon = float("-inf")
    min_height_m: float | None = None

    tile_paths = sorted(
        path for path in derived_dir.glob("*.json")
        if path.name != "tile_index.json"
    )
    if not tile_paths:
        raise ValueError(f"No derived tile JSON files found in {derived_dir}")

    for path in tile_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        envelope = load_derived_tile_envelope(path)
        if envelope is None:
            raise ValueError(f"Derived tile is missing tile bbox: {path}")
        filters = payload.get("filters")
        if isinstance(filters, dict) and min_height_m is None:
            raw_min_height_m = filters.get("min_height_m")
            if isinstance(raw_min_height_m, (int, float)):
                min_height_m = float(raw_min_height_m)
        buildings = payload.get("buildings")
        building_count = len(buildings) if isinstance(buildings, list) else 0

        min_lat = min(min_lat, envelope.min_lat_deg)
        min_lon = min(min_lon, envelope.min_lon_deg)
        max_lat = max(max_lat, envelope.max_lat_deg)
        max_lon = max(max_lon, envelope.max_lon_deg)

        tiles.append(
            {
                "id": path.stem,
                "path": path.name,
                "bbox": {
                    "min_lat": envelope.min_lat_deg,
                    "min_lon": envelope.min_lon_deg,
                    "max_lat": envelope.max_lat_deg,
                    "max_lon": envelope.max_lon_deg,
                },
                "building_count": building_count,
            }
        )

    return {
        "schema_version": 1,
        "city_code": city_code,
        "city_name": city_name,
        "generated_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
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
