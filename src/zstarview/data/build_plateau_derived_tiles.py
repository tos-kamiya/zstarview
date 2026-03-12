#!/usr/bin/env python3
"""Convert PLATEAU CityGML building tiles into lightweight derived JSON tiles."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.data.plateau_citygml import (  # noqa: E402
    DEFAULT_STOREY_HEIGHT_M,
    load_tile_envelope,
    parse_citygml_buildings,
    parse_citygml_numeric,
    resolve_building_height_m,
)

SCHEMA_VERSION = 1
DEFAULT_MIN_BUILDING_HEIGHT_M = 40.0
_CITY_CODE_RE = re.compile(r"^(?P<code>\d{5})")
DEFAULT_WORKERS = max(1, (os.cpu_count() or 1) // 2)


@dataclass(frozen=True)
class DerivedBuilding:
    building_id: str
    height_m: float
    height_source: str
    bbox: tuple[float, float, float, float]
    rings_lonlat: tuple[tuple[tuple[float, float], ...], ...]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert PLATEAU CityGML building tiles to lightweight derived JSON tiles."
    )
    parser.add_argument(
        "--citygml-dir",
        type=Path,
        required=True,
        help="Directory containing PLATEAU building GML files such as raw-data/.../udx/bldg.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where derived JSON tiles will be written.",
    )
    parser.add_argument(
        "--min-building-height-m",
        type=float,
        default=DEFAULT_MIN_BUILDING_HEIGHT_M,
        help="Ignore buildings lower than this height in meters (default: 40.0).",
    )
    parser.add_argument(
        "--tile",
        action="append",
        default=[],
        help="Optional specific tile basename or filename to convert. May be repeated.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help="Number of worker processes to use for tile conversion (default: half of CPU cores).",
    )
    return parser


def detect_city_code(citygml_dir: Path) -> str | None:
    for part in reversed(citygml_dir.parts):
        match = _CITY_CODE_RE.match(part)
        if match:
            return match.group("code")
    return None


def tile_id_from_path(path: Path) -> str:
    return path.stem


def compute_ring_bbox(ring_lonlat: Sequence[tuple[float, float]]) -> tuple[float, float, float, float]:
    lon = [point[0] for point in ring_lonlat]
    lat = [point[1] for point in ring_lonlat]
    return (min(lat), min(lon), max(lat), max(lon))


def compute_building_bbox(
    rings_lonlat: Sequence[Sequence[tuple[float, float]]],
) -> tuple[float, float, float, float]:
    min_lat = float("inf")
    min_lon = float("inf")
    max_lat = float("-inf")
    max_lon = float("-inf")
    for ring in rings_lonlat:
        ring_min_lat, ring_min_lon, ring_max_lat, ring_max_lon = compute_ring_bbox(ring)
        min_lat = min(min_lat, ring_min_lat)
        min_lon = min(min_lon, ring_min_lon)
        max_lat = max(max_lat, ring_max_lat)
        max_lon = max(max_lon, ring_max_lon)
    return (min_lat, min_lon, max_lat, max_lon)


def detect_height_source(path: Path, building_id: str) -> str:
    import xml.etree.ElementTree as ET

    ns = {
        "bldg": "http://www.opengis.net/citygml/building/2.0",
        "gml": "http://www.opengis.net/gml",
    }
    root = ET.parse(path).getroot()
    for building in root.findall(".//bldg:Building", ns):
        current_id = building.attrib.get("{http://www.opengis.net/gml}id")
        if current_id != building_id:
            continue
        measured_height_m = parse_citygml_numeric(
            building.findtext("bldg:measuredHeight", default="", namespaces=ns)
        )
        if measured_height_m is not None:
            return "measuredHeight"
        estimated_height_m = resolve_building_height_m(building)
        if estimated_height_m is not None:
            return f"storeysAboveGround*{DEFAULT_STOREY_HEIGHT_M}"
        break
    return "unknown"


def convert_tile(
    path: Path,
    *,
    min_building_height_m: float,
    city_code: str | None,
    source_root: Path,
) -> dict[str, object] | None:
    envelope = load_tile_envelope(path)
    if envelope is None:
        return None
    buildings = parse_citygml_buildings(path, min_building_height_m=min_building_height_m)
    if not buildings:
        return None

    derived_buildings: list[dict[str, object]] = []
    for building in buildings:
        min_lat, min_lon, max_lat, max_lon = compute_building_bbox(building.rings_lonlat)
        derived_buildings.append(
            {
                "id": building.building_id,
                "height_m": building.height_m,
                "height_source": detect_height_source(path, building.building_id),
                "bbox": {
                    "min_lat": min_lat,
                    "min_lon": min_lon,
                    "max_lat": max_lat,
                    "max_lon": max_lon,
                },
                "rings": [
                    [[lon, lat] for lon, lat in ring]
                    for ring in building.rings_lonlat
                ],
            }
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "source": {
            "format": "PLATEAU-CityGML",
            "city_code": city_code,
            "path": str(path.relative_to(source_root.parent.parent.parent.parent)),
        },
        "tile": {
            "id": tile_id_from_path(path),
            "bbox": {
                "min_lat": envelope.min_lat_deg,
                "min_lon": envelope.min_lon_deg,
                "max_lat": envelope.max_lat_deg,
                "max_lon": envelope.max_lon_deg,
            },
        },
        "filters": {
            "min_height_m": float(min_building_height_m),
            "storey_height_m": DEFAULT_STOREY_HEIGHT_M,
        },
        "buildings": derived_buildings,
    }


def select_tile_paths(citygml_dir: Path, requested_tiles: Sequence[str]) -> tuple[Path, ...]:
    all_tiles = {path.name: path for path in citygml_dir.glob("*.gml")}
    all_tiles.update({path.stem: path for path in citygml_dir.glob("*.gml")})
    if not requested_tiles:
        return tuple(sorted(citygml_dir.glob("*.gml")))
    selected: list[Path] = []
    for tile in requested_tiles:
        path = all_tiles.get(tile)
        if path is None:
            raise ValueError(f"Tile not found: {tile}")
        selected.append(path)
    return tuple(selected)


def write_tile_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def convert_tile_job(
    path_str: str,
    min_building_height_m: float,
    city_code: str | None,
    source_root_str: str,
) -> tuple[str, dict[str, object] | None]:
    path = Path(path_str)
    payload = convert_tile(
        path,
        min_building_height_m=min_building_height_m,
        city_code=city_code,
        source_root=Path(source_root_str),
    )
    return (path.name, payload)


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    citygml_dir = args.citygml_dir
    output_dir = args.output_dir
    city_code = detect_city_code(citygml_dir)
    tile_paths = select_tile_paths(citygml_dir, args.tile)
    workers = max(1, int(args.workers))

    written = 0
    skipped = 0
    if workers == 1:
        results = (
            convert_tile_job(
                str(path),
                min_building_height_m=float(args.min_building_height_m),
                city_code=city_code,
                source_root_str=str(citygml_dir),
            )
            for path in tile_paths
        )
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = executor.map(
                convert_tile_job,
                (str(path) for path in tile_paths),
                [float(args.min_building_height_m)] * len(tile_paths),
                [city_code] * len(tile_paths),
                [str(citygml_dir)] * len(tile_paths),
            )
            for tile_name, payload in results:
                if payload is None:
                    skipped += 1
                    continue
                destination = output_dir / f"{Path(tile_name).stem}.json"
                write_tile_json(destination, payload)
                written += 1
                print(f"[ok] {tile_name} -> {destination}")
        print(f"[ok] written={written} skipped={skipped}")
        return 0

    for tile_name, payload in results:
        if payload is None:
            skipped += 1
            continue
        destination = output_dir / f"{Path(tile_name).stem}.json"
        write_tile_json(destination, payload)
        written += 1
        print(f"[ok] {tile_name} -> {destination}")
    print(f"[ok] written={written} skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
