#!/usr/bin/env python3
"""Match curated skyscraper-city seeds against Scraperbase data."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_SEED_PATH = REPO_ROOT / "dev-samples" / "skyscraper_candidate_cities.seed.json"
HELPER_PATH = REPO_ROOT / "dev-samples" / "build_scraperbase_scan_tiles.py"


def load_helper_module():
    spec = importlib.util.spec_from_file_location("build_scraperbase_scan_tiles", HELPER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load helper module from {HELPER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_helper = load_helper_module()
BuildingRecord = _helper.BuildingRecord
fetch_csv_text = _helper.fetch_csv_text
iter_buildings = _helper.iter_buildings
lonlat_to_tile = _helper.lonlat_to_tile


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Load a curated candidate-city list, match it against Scraperbase's "
            "world CSV, and report 150m+ building counts plus scan-tile counts "
            "for each candidate city."
        )
    )
    parser.add_argument(
        "--seed-json",
        type=Path,
        default=DEFAULT_SEED_PATH,
        help=f"Candidate city seed JSON (default: {DEFAULT_SEED_PATH}).",
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
        help="Web Mercator tile zoom used to count scan tiles (default: 12).",
    )
    parser.add_argument(
        "--status",
        nargs="*",
        default=("completed",),
        help="Statuses to keep, case-insensitive (default: completed).",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional path to write the matched result JSON.",
    )
    return parser


def load_seed_cities(path: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    cities = payload.get("cities")
    if not isinstance(cities, list):
        raise ValueError(f"seed JSON does not contain a list at 'cities': {path}")
    normalized: list[dict[str, object]] = []
    for item in cities:
        if isinstance(item, str):
            normalized.append({"city": item, "aliases": []})
            continue
        if isinstance(item, dict) and isinstance(item.get("city"), str):
            aliases = item.get("aliases", [])
            if not isinstance(aliases, list) or not all(isinstance(alias, str) for alias in aliases):
                raise ValueError(f"seed JSON aliases must be a string list: {path}")
            normalized.append({"city": item["city"], "aliases": aliases})
            continue
        raise ValueError(f"unsupported city entry in seed JSON: {item!r}")
    return normalized, payload


def normalize_city_key(city: str) -> str:
    return " ".join(city.strip().lower().split())


def filter_buildings(
    *,
    csv_text: str,
    min_height_m: int,
    allowed_statuses: set[str],
) -> list[BuildingRecord]:
    return [
        building
        for building in iter_buildings(csv_text)
        if building.height_m >= min_height_m
        and (not allowed_statuses or building.status.lower() in allowed_statuses)
    ]


def summarize_candidate_cities(
    *,
    seed_cities: Sequence[dict[str, object]],
    buildings: Sequence[BuildingRecord],
    zoom: int,
) -> list[dict[str, object]]:
    buckets: dict[str, dict[str, object]] = {}
    for entry in seed_cities:
        seed_city = str(entry["city"])
        aliases = [seed_city, *[str(alias) for alias in entry.get("aliases", [])]]
        bucket = {
            "seed_city": seed_city,
            "aliases": aliases,
            "matched_city_names": set(),
            "matched_country_names": set(),
            "building_count": 0,
            "max_height_m": 0,
            "sample_buildings": [],
            "scan_tiles": set(),
        }
        for alias in aliases:
            buckets[normalize_city_key(alias)] = bucket

    for building in buildings:
        key = normalize_city_key(building.city)
        bucket = buckets.get(key)
        if bucket is None:
            continue
        bucket["building_count"] = int(bucket["building_count"]) + 1
        if building.height_m > int(bucket["max_height_m"]):
            bucket["max_height_m"] = building.height_m
        matched_cities = bucket["matched_city_names"]
        matched_countries = bucket["matched_country_names"]
        scan_tiles = bucket["scan_tiles"]
        samples = bucket["sample_buildings"]
        assert isinstance(matched_cities, set)
        assert isinstance(matched_countries, set)
        assert isinstance(scan_tiles, set)
        assert isinstance(samples, list)
        matched_cities.add(building.city)
        matched_countries.add(building.country)
        scan_tiles.add(lonlat_to_tile(building.lon_deg, building.lat_deg, zoom))
        if len(samples) < 5:
            samples.append(
                {
                    "name": building.name,
                    "height_m": building.height_m,
                    "city": building.city,
                    "country": building.country,
                }
            )

    results: list[dict[str, object]] = []
    seen_seed_cities: set[str] = set()
    for entry in seed_cities:
        seed_city = str(entry["city"])
        if seed_city in seen_seed_cities:
            continue
        seen_seed_cities.add(seed_city)
        bucket = buckets[normalize_city_key(seed_city)]
        tiles = sorted(
            bucket["scan_tiles"],
            key=lambda tile: (tile.zoom, tile.x, tile.y),
        )
        results.append(
            {
                "seed_city": bucket["seed_city"],
                "aliases": [alias for alias in bucket["aliases"] if alias != bucket["seed_city"]],
                "building_count": int(bucket["building_count"]),
                "max_height_m": int(bucket["max_height_m"]),
                "matched_city_names": sorted(bucket["matched_city_names"]),
                "matched_country_names": sorted(bucket["matched_country_names"]),
                "scan_tile_count": len(tiles),
                "scan_tiles": [
                    {"z": tile.zoom, "x": tile.x, "y": tile.y}
                    for tile in tiles
                ],
                "sample_buildings": bucket["sample_buildings"],
            }
        )
    results.sort(key=lambda item: (-int(item["building_count"]), -int(item["max_height_m"]), str(item["seed_city"])))
    return results


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.height_m <= 0:
        parser.error("--height-m must be positive.")
    if args.zoom < 0 or args.zoom > 22:
        parser.error("--zoom must be between 0 and 22.")

    seed_cities, seed_payload = load_seed_cities(args.seed_json)
    allowed_statuses = {status.strip().lower() for status in args.status if status.strip()}
    csv_text = fetch_csv_text(args.input_csv)
    buildings = filter_buildings(
        csv_text=csv_text,
        min_height_m=int(args.height_m),
        allowed_statuses=allowed_statuses,
    )
    results = summarize_candidate_cities(
        seed_cities=seed_cities,
        buildings=buildings,
        zoom=int(args.zoom),
    )

    matched_count = sum(1 for item in results if int(item["building_count"]) > 0)
    print(f"seed_city_count={len(seed_cities)}")
    print(f"matched_city_count={matched_count}")
    print(f"height_m>={args.height_m}")
    print(f"zoom={args.zoom}")
    print(f"statuses={','.join(sorted(allowed_statuses)) if allowed_statuses else 'all'}")
    print("")
    print("rank seed_city building_count max_height_m scan_tile_count matched_country")
    for index, item in enumerate(results, start=1):
        countries = ",".join(item["matched_country_names"]) if item["matched_country_names"] else "-"
        print(
            f"{index:4d} {str(item['seed_city'])[:24]:24s} {int(item['building_count']):14d} "
            f"{int(item['max_height_m']):12d} {int(item['scan_tile_count']):15d} {countries}"
        )

    payload = {
        "seed": seed_payload,
        "height_m": args.height_m,
        "zoom": args.zoom,
        "status": sorted(allowed_statuses),
        "matched_city_count": matched_count,
        "results": results,
    }
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        print("")
        print(f"output_json={args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
