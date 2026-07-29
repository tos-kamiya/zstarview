#!/usr/bin/env python3

"""
List duplicates in cities1000.txt where (country_code, name) is shared by multiple distinct places.
Now prints human-readable admin1 (state/province) name if available from admin1CodesASCII.txt.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass


@dataclass
class CityRec:
    geonameid: int
    name: str
    asciiname: str
    alts: list[str]
    lat: float
    lon: float
    cc: str
    admin1: str
    pop: int
    tz: str

    @classmethod
    def from_cols(cls, cols: list[str]) -> CityRec:
        return cls(
            geonameid=int(cols[0]),
            name=cols[1],
            asciiname=cols[2],
            alts=cols[3].split(",") if cols[3] else [],
            lat=float(cols[4]),
            lon=float(cols[5]),
            cc=cols[8],
            admin1=cols[10],
            pop=int(cols[14]) if cols[14] else 0,
            tz=cols[17],
        )


def load_admin1_names(path: str) -> dict[tuple[str, str], str]:
    """Load admin1CodesASCII.txt as (cc, admin1_code) -> state name."""
    mapping: dict[tuple[str, str], str] = {}
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2 or "." not in parts[0]:
                    continue
                cc, adm1 = parts[0].split(".", 1)
                name = parts[1]
                mapping[(cc, adm1)] = name
    except FileNotFoundError:
        print(f"[warn] admin1CodesASCII.txt not found: {path}", file=sys.stderr)
    return mapping


def find_name_duplicates(path: str, case_insensitive: bool = False) -> dict[tuple[str, str], list[CityRec]]:
    groups: dict[tuple[str, str], list[CityRec]] = {}
    with open(path, encoding="utf-8") as f:
        for ln, line in enumerate(f, 1):
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 19:
                continue
            rec = CityRec.from_cols(cols)
            cc, name = rec.cc, rec.name
            if case_insensitive:
                cc, name = cc.casefold(), name.casefold()
            key = (cc, name)
            groups.setdefault(key, []).append(rec)
    return {k: v for k, v in groups.items() if len(v) > 1}


def main():
    ap = argparse.ArgumentParser(description="List (CC, name) duplicates from cities1000.txt with state names")
    ap.add_argument("--file", "-f", default="cities1000.txt", help="Path to cities1000.txt")
    ap.add_argument("--admin1", default="admin1CodesASCII.txt", help="Path to admin1CodesASCII.txt")
    ap.add_argument("--case-insensitive", action="store_true", help="Case-insensitive comparison on city name")
    args = ap.parse_args()

    admin1_map = load_admin1_names(args.admin1)
    dups = find_name_duplicates(args.file, case_insensitive=args.case_insensitive)

    total_groups = 0
    for cc, name in sorted(dups.keys(), key=lambda x: (x[0], x[1])):
        records = sorted(dups[(cc, name)], key=lambda r: (-r.pop, r.geonameid))
        print(f"{cc}/{name}  (matches: {len(records)})")
        for r in records:
            state_name = admin1_map.get((r.cc, r.admin1), r.admin1 or "?")
            print(f"  - geonameid={r.geonameid}  state={state_name}  " f"lat={r.lat:.6f}  lon={r.lon:.6f}  pop={r.pop}  tz={r.tz}")
        total_groups += 1

    if total_groups == 0:
        print("No duplicates found.")
    else:
        print(f"\nTotal duplicate groups: {total_groups}")


if __name__ == "__main__":
    main()
