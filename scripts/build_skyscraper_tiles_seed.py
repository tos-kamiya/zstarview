#!/usr/bin/env python3
"""Build the curated skyscraper scan-tile seed JSON for zstarview."""

from __future__ import annotations

import argparse
import csv
import json
import math
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

SCRAPERBASE_WORLD_CSV_URL = "https://www.scraperbase.com/export/world.csv"
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CITY_SEED_PATH = REPO_ROOT / "dev-samples" / "skyscraper_candidate_cities.seed.json"
DEFAULT_OUTPUT_PATH = REPO_ROOT / "src" / "zstarview" / "data" / "skyscraper_tiles_z14.json"


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
            "Generate the curated skyscraper scan-tile seed JSON used by zstarview "
            "from the curated candidate-city list and Scraperbase's world CSV."
        )
    )
    parser.add_argument(
        "--city-seed-json",
        type=Path,
        default=DEFAULT_CITY_SEED_PATH,
        help=f"Curated city seed JSON (default: {DEFAULT_CITY_SEED_PATH}).",
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=None,
        help="Optional local Scraperbase CSV path. If omitted, the CSV is downloaded.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"Output path for the generated tile seed JSON (default: {DEFAULT_OUTPUT_PATH}).",
    )
    parser.add_argument(
        "--height-m",
        type=int,
        default=150,
        help="Minimum completed building height in metres to keep (default: 150).",
    )
    parser.add_argument(
        "--zoom",
        type=int,
        default=14,
        help="Web Mercator z/x/y zoom used to define scan tiles (default: 14).",
    )
    parser.add_argument(
        "--status",
        nargs="*",
        default=("completed",),
        help="Statuses to keep, case-insensitive (default: completed).",
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
    raise UnicodeDecodeError("build_skyscraper_tiles_seed", b"", 0, 1, "unsupported CSV encoding")


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


def normalize_city_key(city: str) -> str:
    return " ".join(city.strip().lower().split())


def load_city_seed(path: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    raw_cities = payload.get("cities")
    if not isinstance(raw_cities, list):
        raise ValueError(f"seed JSON does not contain a list at 'cities': {path}")
    cities: list[dict[str, object]] = []
    for item in raw_cities:
        if isinstance(item, str):
            cities.append({"city": item, "aliases": []})
            continue
        if isinstance(item, dict) and isinstance(item.get("city"), str):
            aliases = item.get("aliases", [])
            if not isinstance(aliases, list) or not all(isinstance(alias, str) for alias in aliases):
                raise ValueError(f"seed JSON aliases must be a string list: {path}")
            cities.append({"city": item["city"], "aliases": aliases})
            continue
        raise ValueError(f"unsupported city entry in seed JSON: {item!r}")
    return cities, payload


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


def filter_buildings_for_seed(
    *,
    buildings: Sequence[BuildingRecord],
    seed_cities: Sequence[dict[str, object]],
    min_height_m: int,
    allowed_statuses: set[str],
) -> list[BuildingRecord]:
    allowed_names: set[str] = set()
    for entry in seed_cities:
        city = str(entry["city"])
        allowed_names.add(normalize_city_key(city))
        for alias in entry.get("aliases", []):
            allowed_names.add(normalize_city_key(str(alias)))
    return [
        building
        for building in buildings
        if building.height_m >= min_height_m
        and (not allowed_statuses or building.status.lower() in allowed_statuses)
        and normalize_city_key(building.city) in allowed_names
    ]


def build_payload(
    *,
    city_seed_payload: dict[str, object],
    filtered_buildings: Sequence[BuildingRecord],
    zoom: int,
    min_height_m: int,
    allowed_statuses: set[str],
    source: str,
) -> dict[str, object]:
    tile_entries: dict[ScanTile, dict[str, object]] = {}
    city_summaries: dict[str, dict[str, object]] = {}

    for building in filtered_buildings:
        city_key = normalize_city_key(building.city)
        city_summary = city_summaries.setdefault(
            city_key,
            {
                "city": building.city,
                "country": building.country,
                "building_count": 0,
                "max_height_m": 0,
                "scan_tiles": set(),
                "sample_buildings": [],
            },
        )
        city_summary["building_count"] = int(city_summary["building_count"]) + 1
        city_summary["max_height_m"] = max(int(city_summary["max_height_m"]), building.height_m)
        cast_city_tiles = city_summary["scan_tiles"]
        cast_city_samples = city_summary["sample_buildings"]
        assert isinstance(cast_city_tiles, set)
        assert isinstance(cast_city_samples, list)
        tile = lonlat_to_tile(building.lon_deg, building.lat_deg, zoom)
        cast_city_tiles.add(tile)
        if len(cast_city_samples) < 5:
            cast_city_samples.append(
                {
                    "name": building.name,
                    "height_m": building.height_m,
                    "lat_deg": building.lat_deg,
                    "lon_deg": building.lon_deg,
                }
            )

        tile_entry = tile_entries.setdefault(
            tile,
            {
                "tile": {"z": tile.zoom, "x": tile.x, "y": tile.y},
                "bbox": {},
                "city_names": set(),
                "country_names": set(),
                "building_count": 0,
                "max_height_m": 0,
                "sample_buildings": [],
            },
        )
        tile_entry["building_count"] = int(tile_entry["building_count"]) + 1
        tile_entry["max_height_m"] = max(int(tile_entry["max_height_m"]), building.height_m)
        cast_tile_cities = tile_entry["city_names"]
        cast_tile_countries = tile_entry["country_names"]
        cast_tile_samples = tile_entry["sample_buildings"]
        assert isinstance(cast_tile_cities, set)
        assert isinstance(cast_tile_countries, set)
        assert isinstance(cast_tile_samples, list)
        cast_tile_cities.add(building.city)
        cast_tile_countries.add(building.country)
        if len(cast_tile_samples) < 5:
            cast_tile_samples.append(
                {
                    "name": building.name,
                    "city": building.city,
                    "country": building.country,
                    "height_m": building.height_m,
                    "lat_deg": building.lat_deg,
                    "lon_deg": building.lon_deg,
                }
            )

    tile_list: list[dict[str, object]] = []
    for tile, entry in sorted(
        tile_entries.items(),
        key=lambda item: (-int(item[1]["building_count"]), -int(item[1]["max_height_m"]), item[0].x, item[0].y),
    ):
        west, south, east, north = tile_bounds(tile)
        entry["bbox"] = {
            "west": west,
            "south": south,
            "east": east,
            "north": north,
        }
        tile_list.append(
            {
                "tile": entry["tile"],
                "bbox": entry["bbox"],
                "building_count": int(entry["building_count"]),
                "max_height_m": int(entry["max_height_m"]),
                "city_names": sorted(entry["city_names"]),
                "country_names": sorted(entry["country_names"]),
                "sample_buildings": entry["sample_buildings"],
            }
        )

    city_list: list[dict[str, object]] = []
    for city_summary in sorted(
        city_summaries.values(),
        key=lambda item: (-int(item["building_count"]), -int(item["max_height_m"]), str(item["city"])),
    ):
        scan_tiles = city_summary["scan_tiles"]
        assert isinstance(scan_tiles, set)
        city_list.append(
            {
                "city": city_summary["city"],
                "country": city_summary["country"],
                "building_count": int(city_summary["building_count"]),
                "max_height_m": int(city_summary["max_height_m"]),
                "scan_tile_count": len(scan_tiles),
                "sample_buildings": city_summary["sample_buildings"],
            }
        )

    return {
        "schema_version": 1,
        "description": (
            "Curated skyscraper scan tiles for the optional far-range urban-outline layer. "
            "Tiles are Web Mercator z/x/y helpers and must be converted back to bbox "
            "before calling overturemaps download."
        ),
        "selection": {
            "height_m_gte": min_height_m,
            "status": sorted(allowed_statuses),
            "zoom": zoom,
        },
        "sources": {
            "city_seed": city_seed_payload,
            "buildings": source,
        },
        "summary": {
            "city_count": len(city_list),
            "building_count": len(filtered_buildings),
            "tile_count": len(tile_list),
        },
        "cities": city_list,
        "tiles": tile_list,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.height_m <= 0:
        parser.error("--height-m must be positive.")
    if args.zoom < 0 or args.zoom > 22:
        parser.error("--zoom must be between 0 and 22.")

    city_seed, city_seed_payload = load_city_seed(args.city_seed_json)
    csv_text = fetch_csv_text(args.input_csv)
    buildings = list(iter_buildings(csv_text))
    allowed_statuses = {status.strip().lower() for status in args.status if status.strip()}
    filtered_buildings = filter_buildings_for_seed(
        buildings=buildings,
        seed_cities=city_seed,
        min_height_m=int(args.height_m),
        allowed_statuses=allowed_statuses,
    )
    payload = build_payload(
        city_seed_payload=city_seed_payload,
        filtered_buildings=filtered_buildings,
        zoom=int(args.zoom),
        min_height_m=int(args.height_m),
        allowed_statuses=allowed_statuses,
        source=args.input_csv.as_posix() if args.input_csv is not None else SCRAPERBASE_WORLD_CSV_URL,
    )

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"output_json={args.output_json}")
    print(f"city_count={payload['summary']['city_count']}")
    print(f"building_count={payload['summary']['building_count']}")
    print(f"tile_count={payload['summary']['tile_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
