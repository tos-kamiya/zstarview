#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from zstarview.mountain_viewpoints import (
    _ascii_fallback_name as mountain_ascii_fallback_name,
    load_mountain_viewpoints,
)
from zstarview.paths import CITY_COORD_FILE
from zstarview.tower_viewpoints import (
    _ascii_fallback_name as tower_ascii_fallback_name,
    load_tower_viewpoints,
)
from zstarview.utils.resolve_city import _variants, CityRec


@dataclass(frozen=True)
class NameUse:
    kind: str
    canonical_name: str
    detail: str


def _iter_city_name_uses(cities_path: str) -> Iterable[tuple[str, NameUse]]:
    with open(cities_path, encoding="utf-8") as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 19:
                continue
            city = CityRec.from_cols(cols)
            detail = f"{city.cc}/{city.name} (geonameid={city.geonameid}, pop={city.pop})"
            raw_names = [city.name, city.asciiname]
            if cols[3]:
                raw_names.extend(name for name in cols[3].split(",") if name)
            seen_keys: set[str] = set()
            for raw_name in raw_names:
                for key in _variants(raw_name):
                    if key in seen_keys:
                        continue
                    seen_keys.add(key)
                    yield key, NameUse("city", city.name, detail)


def _iter_tower_name_uses() -> Iterable[tuple[str, NameUse]]:
    for tower in load_tower_viewpoints():
        detail = f"{tower.name} [{tower.persistent_key}]"
        raw_names = [
            tower.name,
            *tower.names,
            *tower.labels.values(),
        ]
        ascii_names = [
            tower.ascii_name,
            *(tower_ascii_fallback_name(name) for name in tower.names),
            *(tower_ascii_fallback_name(name) for name in tower.labels.values()),
        ]
        seen_keys: set[str] = set()
        for raw_name in [*raw_names, *ascii_names]:
            if not raw_name:
                continue
            for key in _variants(raw_name):
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                yield key, NameUse("tower", tower.name, detail)


def _iter_mountain_name_uses() -> Iterable[tuple[str, NameUse]]:
    for mountain in load_mountain_viewpoints():
        detail = f"{mountain.name} [{mountain.persistent_key}]"
        raw_names = [
            mountain.name,
            *mountain.names,
            *mountain.labels.values(),
        ]
        ascii_names = [
            mountain.ascii_name,
            *(mountain_ascii_fallback_name(name) for name in mountain.names),
            *(mountain_ascii_fallback_name(name) for name in mountain.labels.values()),
        ]
        seen_keys: set[str] = set()
        for raw_name in [*raw_names, *ascii_names]:
            if not raw_name:
                continue
            for key in _variants(raw_name):
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                yield key, NameUse("mountain", mountain.name, detail)


def build_collision_index(cities_path: str) -> dict[str, list[NameUse]]:
    index: dict[str, list[NameUse]] = defaultdict(list)
    for key, use in _iter_city_name_uses(cities_path):
        index[key].append(use)
    for key, use in _iter_tower_name_uses():
        index[key].append(use)
    for key, use in _iter_mountain_name_uses():
        index[key].append(use)
    return index


def _dedupe_uses(uses: list[NameUse]) -> list[NameUse]:
    seen: set[tuple[str, str, str]] = set()
    out: list[NameUse] = []
    for use in uses:
        key = (use.kind, use.canonical_name, use.detail)
        if key in seen:
            continue
        seen.add(key)
        out.append(use)
    return out


def find_cross_kind_collisions(cities_path: str) -> list[tuple[str, list[NameUse]]]:
    collisions: list[tuple[str, list[NameUse]]] = []
    for key, uses in build_collision_index(cities_path).items():
        deduped = _dedupe_uses(uses)
        kinds = {use.kind for use in deduped}
        if len(kinds) < 2:
            continue
        collisions.append((key, sorted(deduped, key=lambda use: (use.kind, use.canonical_name, use.detail))))
    collisions.sort(key=lambda item: (len({use.kind for use in item[1]}), item[0]), reverse=True)
    return collisions


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check name collisions between GeoNames cities and bundled tower/mountain viewpoints."
    )
    parser.add_argument(
        "--cities-file",
        default=CITY_COORD_FILE,
        help="Path to cities1000.txt (default: packaged cities1000.txt).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=50,
        help="Maximum number of collision groups to print (default: 50).",
    )
    args = parser.parse_args()

    cities_path = str(Path(args.cities_file))
    collisions = find_cross_kind_collisions(cities_path)
    print(f"cross-kind collision groups: {len(collisions)}")
    for index, (key, uses) in enumerate(collisions[: max(0, args.limit)], start=1):
        kinds = ", ".join(sorted({use.kind for use in uses}))
        print(f"\n[{index}] key={key!r} kinds={kinds}")
        for use in uses:
            print(f"  - {use.kind}: {use.detail}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
