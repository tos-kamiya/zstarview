# -*- coding: utf-8 -*-
"""
City name resolution utility.

This module provides functions to find and parse city data from the GeoNames
cities1000.txt file. It supports resolving cities by name, country code + name,
or geonameid. It includes a CLI for testing and data exploration.
"""

import argparse
import re
import sys
import unicodedata
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

# ---------- Data Model ----------


@dataclass(frozen=True)
class CityRec:
    """A record representing a city from the GeoNames database."""

    geonameid: int
    name: str  # The name of the city in its local language.
    asciiname: str  # ASCII-only version of the city name.
    lat: float
    lon: float
    cc: str  # ISO 3166-1 alpha-2 country code (e.g., 'JP').
    admin1_code: str  # Code for the first-order administrative division (e.g., a state or province).
    admin1_name: Optional[str]  # Resolved ASCII name of the admin1 division.
    pop: int  # Population.
    tz: str  # IANA time zone identifier (e.g., 'Asia/Tokyo').

    @classmethod
    def from_cols(cls, cols: List[str], admin1_name: Optional[str] = None) -> "CityRec":
        """
        Creates a CityRec from a list of columns from a cities1000.txt line.

        Args:
            cols: A list of strings from a tab-separated line.
            admin1_name: The resolved name of the admin1 division, if available.

        Returns:
            A new CityRec instance.
        """
        return cls(
            geonameid=int(cols[0]),
            name=cols[1],
            asciiname=cols[2],
            lat=float(cols[4]),
            lon=float(cols[5]),
            cc=cols[8],
            admin1_code=cols[10],
            admin1_name=admin1_name,
            pop=int(cols[14]) if cols[14] else 0,
            tz=cols[17],
        )


# -------------------------
# Normalization Helpers
# -------------------------


def _nfkc_casefold(s: str) -> str:
    """Applies NFKC normalization and case-folding."""
    return unicodedata.normalize("NFKC", s).casefold().strip()


def _strip_diacritics(s: str) -> str:
    """Removes diacritical marks (accents) from a string."""
    return "".join(ch for ch in unicodedata.normalize("NFKD", s) if unicodedata.category(ch) != "Mn")


def _norm(s: str) -> str:
    """
    Fully normalizes a string for robust matching.

    The process includes case-folding, removing diacritics, stripping extra
    whitespace, and standardizing common symbols like hyphens and apostrophes.
    """
    s = _nfkc_casefold(s)
    s = _strip_diacritics(s)
    s = s.replace("’", "").replace("'", "")
    s = s.replace("–", "-").replace("—", "-")
    s = " ".join(s.split())
    return s


def _variants(s: str) -> List[str]:
    """
    Generates variations of a string to handle hyphens and spaces.

    For a given string, this creates versions with hyphens/spaces swapped or removed
    to match different spelling conventions (e.g., "Saint-Etienne", "Saint Etienne").
    """
    v = {_norm(s)}
    base = next(iter(v))
    if "-" in base:
        v.add(base.replace("-", " "))
        v.add(base.replace("-", ""))
    if " " in base:
        v.add(base.replace(" ", "-"))
        v.add(base.replace(" ", ""))
    return list(v)


def load_admin1_names(path: str) -> Dict[Tuple[str, str], str]:
    """
    Loads admin1CodesASCII.txt into a mapping.

    The mapping is from (country_code, admin1_code) to an ASCII admin1 name.
    This is used to enrich CityRec objects with human-readable state/province names.

    Args:
        path: The file path to admin1CodesASCII.txt.

    Returns:
        A dictionary mapping (cc, admin1_code) to the admin1 name.
    """
    mapping: Dict[Tuple[str, str], str] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or "." not in parts[0]:
                continue
            cc, adm1 = parts[0].split(".", 1)
            name = parts[2] or parts[1]  # Prefer ASCII Name over Name
            mapping[(cc, adm1)] = name
    return mapping


def _row_matches(cols: List[str], query_variants: set[str]) -> bool:
    """
    Checks if a city row matches any of the query variants.

    It checks the `name`, `asciiname`, and `alternatenames` fields.

    Args:
        cols: The columns of a line from the cities file.
        query_variants: A set of normalized string variants of the query name.

    Returns:
        True if a match is found.
    """
    # Check the 'name' field
    for k in _variants(cols[1]):
        if k in query_variants:
            return True
    # Check the 'asciiname' field
    if cols[2]:
        for k in _variants(cols[2]):
            if k in query_variants:
                return True
    # Check the 'alternatenames' field (comma-separated)
    if cols[3]:
        for alt in cols[3].split(","):
            if not alt:
                continue
            for k in _variants(alt):
                if k in query_variants:
                    return True
    return False


# -------------------------
# Core Resolution Functions
# -------------------------


def resolve_city(
    cc_slash_name: str,
    cities_path: str = "cities1000.txt",
    admin1_map: Optional[Dict[Tuple[str, str], str]] = None,
) -> List[CityRec]:
    """
    Resolves a 'CC/CityName' string against the GeoNames database.

    This function searches for cities within a specific country. It matches against
    the `name`, `asciiname`, and `alternatenames` fields.

    Args:
        cc_slash_name: The query string, e.g., 'ES/Zaragoza'.
        cities_path: The path to the cities1000.txt file.
        admin1_map: An optional pre-loaded map of admin1 names.

    Returns:
        A list of matching CityRec objects, sorted by population (descending)
        and then geonameid (ascending).
    """
    if "/" not in cc_slash_name:
        raise ValueError("Input must be in 'CC/CityName' format (e.g., 'ES/Zaragoza').")
    cc_in, city_in = cc_slash_name.split("/", 1)
    cc_key = _nfkc_casefold(cc_in)
    query_variants = set(_variants(city_in))

    matches: List[CityRec] = []
    with open(cities_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 19:
                continue
            # Filter by country code first for efficiency
            if _nfkc_casefold(cols[8]) != cc_key:
                continue
            if not _row_matches(cols, query_variants):
                continue

            adm1_name = admin1_map.get((cols[8], cols[10])) if admin1_map else None
            matches.append(CityRec.from_cols(cols, admin1_name=adm1_name))

    # Sort matches by population (most populous first) as the primary key.
    matches.sort(key=lambda r: (-r.pop, r.geonameid))
    return matches


def resolve_city_by_name(
    name: str,
    cities_path: str = "cities1000.txt",
    admin1_map: Optional[Dict[Tuple[str, str], str]] = None,
) -> List[CityRec]:
    """
    Resolves a city by its name, searching worldwide.

    This function does not filter by country and is useful when the country
    is unknown. It searches the `name`, `asciiname`, and `alternatenames` fields.

    Args:
        name: The city name or alias to search for.
        cities_path: The path to the cities1000.txt file.
        admin1_map: An optional pre-loaded map of admin1 names.

    Returns:
        A list of matching CityRec objects, sorted by population (descending)
        and then geonameid (ascending).
    """
    query_variants = set(_variants(name))
    matches: List[CityRec] = []

    with open(cities_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 19:
                continue
            if not _row_matches(cols, query_variants):
                continue

            adm1_name = admin1_map.get((cols[8], cols[10])) if admin1_map else None
            matches.append(CityRec.from_cols(cols, admin1_name=adm1_name))

    matches.sort(key=lambda r: (-r.pop, r.geonameid))
    return matches


def resolve_city_by_geonameid(
    prefer_geonameid: int,
    cities_path: str = "cities1000.txt",
    admin1_map: Optional[Dict[Tuple[str, str], str]] = None,
) -> Optional[CityRec]:
    """
    Finds a single city by its unique geonameid.

    Args:
        prefer_geonameid: The integer geonameid to find.
        cities_path: The path to the cities1000.txt file.
        admin1_map: An optional pre-loaded map of admin1 names.

    Returns:
        The matching CityRec object, or None if not found.
    """
    with open(cities_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 19:
                continue
            gid = int(cols[0])
            if gid == prefer_geonameid:
                adm1_name = admin1_map.get((cols[8], cols[10])) if admin1_map else None
                return CityRec.from_cols(cols, admin1_name=adm1_name)
    return None


# ---------- CLI ----------


def main():
    """Command-line interface for resolving city names."""
    ap = argparse.ArgumentParser(description="Resolve city names from cities1000.txt. Finds all matches and sorts by population.")
    ap.add_argument("city", help="Query string: 'CC/CityName' (e.g., 'ES/Zaragoza'), a city name, or a geonameid.")
    ap.add_argument("-f", "--file", default="cities1000.txt", help="Path to cities1000.txt")
    ap.add_argument("--admin1", default="admin1CodesASCII.txt", help="Path to admin1CodesASCII.txt for state/province names")
    args = ap.parse_args()

    try:
        admin1_map = load_admin1_names(args.admin1)
    except FileNotFoundError:
        print(f"Warning: Could not find admin1 file at '{args.admin1}'. State/province names will not be resolved.", file=sys.stderr)
        admin1_map = {}

    query = args.city
    recs: Optional[List[CityRec]] = None

    # Determine resolution strategy based on query format
    if re.match(r"^\d+$", query):
        # --- GeonameID lookup ---
        print(f"Resolving by geonameid: {query}...")
        rec = resolve_city_by_geonameid(int(query), args.file, admin1_map)
        recs = [rec] if rec else []
    elif "/" in query:
        # --- Country-specific name lookup ---
        print(f"Resolving by CC/Name: {query}...")
        recs = resolve_city(query, args.file, admin1_map)
    else:
        # --- Global name lookup ---
        print(f"Resolving by name (worldwide): {query}...")
        recs = resolve_city_by_name(query, args.file, admin1_map)

    # --- Print results ---
    if recs:
        print(f"Found {len(recs)} match(es) for '{query}':")
        for rec in recs:
            admin_info = f", {rec.admin1_name}" if rec.admin1_name else ""
            print(f"  - {rec.cc}/{rec.name}{admin_info} | " f"lat: {rec.lat:.6f}, lon: {rec.lon:.6f}, tz: {rec.tz}, pop: {rec.pop} " f"(geonameid={rec.geonameid})")
    else:
        print(f"No match found for '{query}'", file=sys.stderr)


if __name__ == "__main__":
    main()
