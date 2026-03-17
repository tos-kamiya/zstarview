#!/usr/bin/env python3
"""Build scan-tile candidates from Scraperbase skyscraper CSV."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

SCRAPERBASE_WORLD_CSV_URL = "https://www.scraperbase.com/export/world.csv"


@dataclass(frozen=True)
class BuildingRecord:
    name: str
    city: str
    country: str
    height_m: int
    lat_deg: float
    lon_deg: float
    status: str


@dataclass(frozen=True)
class ScanTile:
    zoom: int
    x: int
    y: int


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download or read Scraperbase's world CSV, keep skyscrapers at or above "
            "a given height, and group them into Web Mercator z/x/y scan tiles that "
            "can be converted into Overture bbox download regions."
        )
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=None,
        help="Optional local Scraperbase CSV path. If omitted, the script downloads the CSV.",
    )
    parser.add_argument(
        "--height-m",
        type=int,
        default=150,
        help="Minimum building height in metres to keep (default: 150).",
    )
    parser.add_argument(
        "--zoom",
        type=int,
        default=12,
        help="Web Mercator tile zoom used to cluster buildings into scan tiles (default: 12).",
    )
    parser.add_argument(
        "--status",
        nargs="*",
        default=("completed",),
        help="Statuses to keep, case-insensitive (default: completed).",
    )
    parser.add_argument(
        "--top-cities",
        type=int,
        default=30,
        help="How many cities to show in the summary (default: 30).",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional path to write the full tile summary JSON.",
    )
    return parser


def fetch_csv_text(input_csv: Path | None) -> str:
    if input_csv is not None:
        raw = input_csv.read_bytes()
    else:
        with urllib.request.urlopen(SCRAPERBASE_WORLD_CSV_URL) as response:
            raw = response.read()
    for encoding in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError("build_scraperbase_scan_tiles", b"", 0, 1, "unsupported CSV encoding")


def parse_int(value: str) -> int | None:
    value = value.strip()
    if not value or value == "-":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_float(value: str) -> float | None:
    value = value.strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def iter_buildings(csv_text: str) -> Iterable[BuildingRecord]:
    reader = csv.DictReader(csv_text.splitlines())
    for row in reader:
        height_m = parse_int(row.get("Height (whole metres)", ""))
        lat_deg = parse_float(row.get("Latitude", ""))
        lon_deg = parse_float(row.get("Longitude", ""))
        if height_m is None or lat_deg is None or lon_deg is None:
            continue
        yield BuildingRecord(
            name=(row.get("Building Name") or "").strip(),
            city=(row.get("City") or "").strip(),
            country=(row.get("Country") or "").strip(),
            height_m=height_m,
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            status=(row.get("Status") or "").strip(),
        )


def lonlat_to_tile(lon_deg: float, lat_deg: float, zoom: int) -> ScanTile:
    lat_deg = max(min(lat_deg, 85.05112878), -85.05112878)
    n = 2**zoom
    x = int((lon_deg + 180.0) / 360.0 * n)
    lat_rad = math.radians(lat_deg)
    y = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
    return ScanTile(zoom=zoom, x=max(0, min(n - 1, x)), y=max(0, min(n - 1, y)))


def tile_bounds(tile: ScanTile) -> tuple[float, float, float, float]:
    n = 2**tile.zoom
    west = tile.x / n * 360.0 - 180.0
    east = (tile.x + 1) / n * 360.0 - 180.0

    def tile_y_to_lat(y: int) -> float:
        merc_n = math.pi * (1.0 - 2.0 * y / n)
        return math.degrees(math.atan(math.sinh(merc_n)))

    north = tile_y_to_lat(tile.y)
    south = tile_y_to_lat(tile.y + 1)
    return (west, south, east, north)


def summarize_tiles(
    *,
    buildings: Sequence[BuildingRecord],
    zoom: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    city_summary: dict[tuple[str, str], dict[str, object]] = {}
    tile_summary: dict[ScanTile, dict[str, object]] = {}

    for building in buildings:
        city_key = (building.city, building.country)
        city_entry = city_summary.setdefault(
            city_key,
            {
                "city": building.city,
                "country": building.country,
                "building_count": 0,
                "max_height_m": 0,
                "sample_building": building.name,
            },
        )
        city_entry["building_count"] = int(city_entry["building_count"]) + 1
        if building.height_m > int(city_entry["max_height_m"]):
            city_entry["max_height_m"] = building.height_m
            city_entry["sample_building"] = building.name

        tile = lonlat_to_tile(building.lon_deg, building.lat_deg, zoom)
        tile_entry = tile_summary.setdefault(
            tile,
            {
                "z": tile.zoom,
                "x": tile.x,
                "y": tile.y,
                "building_count": 0,
                "max_height_m": 0,
                "cities": set(),
                "sample_buildings": [],
            },
        )
        tile_entry["building_count"] = int(tile_entry["building_count"]) + 1
        if building.height_m > int(tile_entry["max_height_m"]):
            tile_entry["max_height_m"] = building.height_m
        cast_cities = tile_entry["cities"]
        assert isinstance(cast_cities, set)
        cast_cities.add(f"{building.city}, {building.country}")
        cast_samples = tile_entry["sample_buildings"]
        assert isinstance(cast_samples, list)
        if len(cast_samples) < 5:
            cast_samples.append(
                {
                    "name": building.name,
                    "city": building.city,
                    "country": building.country,
                    "height_m": building.height_m,
                }
            )

    cities = sorted(
        city_summary.values(),
        key=lambda item: (-int(item["building_count"]), -int(item["max_height_m"]), str(item["city"])),
    )
    tiles: list[dict[str, object]] = []
    for tile, entry in sorted(
        tile_summary.items(),
        key=lambda item: (-int(item[1]["building_count"]), -int(item[1]["max_height_m"]), item[0].x, item[0].y),
    ):
        west, south, east, north = tile_bounds(tile)
        tiles.append(
            {
                "z": tile.zoom,
                "x": tile.x,
                "y": tile.y,
                "bbox": {
                    "west": west,
                    "south": south,
                    "east": east,
                    "north": north,
                },
                "building_count": int(entry["building_count"]),
                "max_height_m": int(entry["max_height_m"]),
                "city_count": len(entry["cities"]),
                "cities": sorted(entry["cities"]),
                "sample_buildings": entry["sample_buildings"],
            }
        )
    return cities, tiles


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.zoom < 0 or args.zoom > 22:
        parser.error("--zoom must be between 0 and 22.")
    if args.height_m <= 0:
        parser.error("--height-m must be positive.")

    csv_text = fetch_csv_text(args.input_csv)
    allowed_statuses = {status.strip().lower() for status in args.status if status.strip()}
    buildings = [
        building
        for building in iter_buildings(csv_text)
        if building.height_m >= args.height_m
        and (not allowed_statuses or building.status.lower() in allowed_statuses)
    ]
    cities, tiles = summarize_tiles(buildings=buildings, zoom=args.zoom)

    payload = {
        "source": args.input_csv.as_posix() if args.input_csv is not None else SCRAPERBASE_WORLD_CSV_URL,
        "height_m": args.height_m,
        "zoom": args.zoom,
        "status": sorted(allowed_statuses),
        "building_count": len(buildings),
        "city_count": len(cities),
        "tile_count": len(tiles),
        "top_cities": cities[: args.top_cities],
        "tiles": tiles,
    }

    print(f"source={payload['source']}")
    print(f"height_m>={args.height_m}")
    print(f"zoom={args.zoom}")
    print(f"statuses={','.join(sorted(allowed_statuses)) if allowed_statuses else 'all'}")
    print(f"building_count={len(buildings)}")
    print(f"city_count={len(cities)}")
    print(f"tile_count={len(tiles)}")
    print("")
    print(f"Top cities by building_count (top {min(args.top_cities, len(cities))})")
    print("rank city country building_count max_height_m sample_building")
    for index, city in enumerate(cities[: args.top_cities], start=1):
        print(
            f"{index:4d} {str(city['city'])[:28]:28s} {str(city['country'])[:20]:20s} "
            f"{int(city['building_count']):14d} {int(city['max_height_m']):12d} {city['sample_building']}"
        )

    print("")
    print("Top scan tiles by building_count")
    print("rank z/x/y building_count max_height_m city_count bbox_west_south_east_north")
    for index, tile in enumerate(tiles[: min(20, len(tiles))], start=1):
        bbox = tile["bbox"]
        assert isinstance(bbox, dict)
        print(
            f"{index:4d} {int(tile['z'])}/{int(tile['x'])}/{int(tile['y'])} "
            f"{int(tile['building_count']):14d} {int(tile['max_height_m']):12d} "
            f"{int(tile['city_count']):10d} "
            f"{bbox['west']:.6f},{bbox['south']:.6f},{bbox['east']:.6f},{bbox['north']:.6f}"
        )

    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        print("")
        print(f"output_json={args.output_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
