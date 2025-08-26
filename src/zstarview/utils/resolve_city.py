import argparse
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import unicodedata
import re

# ---------- data model ----------


@dataclass(frozen=True)
class CityRec:
    geonameid: int
    name: str  # The name of the city in the local language, from cities1000.txt.
    asciiname: str
    lat: float
    lon: float
    cc: str  # Country code (e.g., 'JP').
    admin1_code: str  # Admin1 code (e.g., for a state or province).
    admin1_name: Optional[str]  # Admin1 name (e.g., 'Saitama'); None if not resolved.
    pop: int
    tz: str

    @classmethod
    def from_cols(cls, cols: List[str], admin1_name: Optional[str] = None) -> "CityRec":
        """Creates a CityRec from a split line of cities1000.txt. The admin1_name can be added separately."""
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
# Helpers
# -------------------------


def _nfkc_casefold(s: str) -> str:
    return unicodedata.normalize("NFKC", s).casefold().strip()


def _strip_diacritics(s: str) -> str:
    """Removes diacritics from a string."""
    return "".join(ch for ch in unicodedata.normalize("NFKD", s) if unicodedata.category(ch) != "Mn")


def _norm(s: str) -> str:
    """Normalizes a string by case-folding, removing diacritics, stripping extra whitespace, and standardizing common symbols."""
    s = _nfkc_casefold(s)
    s = _strip_diacritics(s)
    s = s.replace("’", "").replace("'", "")
    s = s.replace("–", "-").replace("—", "-")
    s = " ".join(s.split())
    return s


def _variants(s: str) -> List[str]:
    """Generates variations of a string to account for hyphens and spaces."""
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
    """Loads admin1CodesASCII.txt into a mapping from (country_code, admin1_code) to an ASCII admin1 name."""
    mapping: Dict[Tuple[str, str], str] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or "." not in parts[0]:
                continue
            cc, adm1 = parts[0].split(".", 1)
            name = parts[2] or parts[1]  # Prefer ASCIIName
            mapping[(cc, adm1)] = name
    return mapping


def _row_matches(cols: List[str], query_variants: set[str]) -> bool:
    """Checks if the name, asciiname, or alternatenames in a row match any of the query variants."""
    # name
    for k in _variants(cols[1]):
        if k in query_variants:
            return True
    # asciiname
    if cols[2]:
        for k in _variants(cols[2]):
            if k in query_variants:
                return True
    # alternatenames (comma-separated)
    if cols[3]:
        for alt in cols[3].split(","):
            if not alt:
                continue
            for k in _variants(alt):
                if k in query_variants:
                    return True
    return False


# -------------------------
# Core utilities
# -------------------------


def resolve_city(
    cc_slash_name: str,
    cities_path: str = "cities1000.txt",
    admin1_map: Optional[Dict[Tuple[str, str], str]] = None,
) -> List[CityRec]:
    """Resolves a 'CC/CityNameOrAlias' string by searching the name, asciiname, and alternatenames fields in cities1000.txt.

    Returns a list of matching CityRec objects, sorted by population (descending)
    and then geonameid (ascending). Matching is exact (after normalization and
    variation handling). This function does not print warnings.
    """
    if "/" not in cc_slash_name:
        raise ValueError("Input must be 'CC/CityName' (e.g., 'ES/Zaragoza').")
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
            if _nfkc_casefold(cols[8]) != cc_key:
                continue
            if not _row_matches(cols, query_variants):
                continue
            adm1_name = admin1_map.get((cols[8], cols[10])) if admin1_map else None
            matches.append(CityRec.from_cols(cols, admin1_name=adm1_name))

    matches.sort(key=lambda r: (-r.pop, r.geonameid))
    return matches


def resolve_city_by_name(
    name: str,
    cities_path: str = "cities1000.txt",
    admin1_map: Optional[Dict[Tuple[str, str], str]] = None,
) -> List[CityRec]:
    """Resolves a city by its name or alias (searching name, asciiname, and alternatenames).

    Returns a list of matching CityRec objects, sorted by population (descending)
    and then geonameid (ascending). Does not filter by country code (searches
    worldwide).
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
    """Resolves a city by its geonameid."""
    with open(cities_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 19:
                continue
            gid = int(cols[0])
            if gid == prefer_geonameid:
                adm1_name = None
                if admin1_map is not None:
                    adm1_name = admin1_map.get((cols[8], cols[10]))
                rec = CityRec.from_cols(cols, admin1_name=adm1_name)
                return rec
    return None


# ---------- CLI ----------


def main():
    """CLI for resolving city names."""
    import sys

    ap = argparse.ArgumentParser(description="Resolve 'CC/CityName' in cities1000.txt, warn on multis, pick largest by population.")
    ap.add_argument("city", help="City string like 'ES/Zaragoza' (case-insensitive), OR just a geonameid (e.g., '3128760')")
    ap.add_argument("-f", "--file", default="cities1000.txt", help="Path to cities1000.txt")
    ap.add_argument("--admin1", default="admin1CodesASCII.txt", help="Path to admin1CodesASCII.txt (for state names)")
    args = ap.parse_args()

    admin1_map = load_admin1_names(args.admin1)
    c = args.city
    if re.match(r"^\d+$", c):
        # If input is just a geonameid, resolve it directly
        rec = resolve_city_by_geonameid(int(c), args.file)
        if rec:
            print(f"Resolved geonameid {c} to {rec.cc}/{rec.name}, lat: {rec.lat:.6f}, lon: {rec.lon:.6f}, tz: {rec.tz}")
        else:
            print(f"No city found for geonameid {c}", file=sys.stderr)
    elif not "/" in c:
        # If input is just a city name, resolve it by name
        recs = resolve_city_by_name(c, args.file, admin1_map)
        if not recs:
            print(f"No match for '{c}'", file=sys.stderr)
        else:
            print(f"Found {len(recs)} match(es) for '{c}':")
            for rec in recs:
                print(f"{rec.cc}/{rec.name}, lat: {rec.lat:.6f}, lon: {rec.lon:.6f}, tz: {rec.tz}  (geonameid={rec.geonameid})")
    else:
        recs = do_resolve_city(args.city, args.file, admin1_map)
        if recs:
            print(f"Found {len(recs)} match(es) for '{args.city}':")
            for rec in recs:
                print(f"{rec.cc}/{rec.name}, lat: {rec.lat:.6f}, lon: {rec.lon:.6f}, tz: {rec.tz}  (geonameid={rec.geonameid})")
        else:
            print(f"No match for '{args.city}'", file=sys.stderr)


if __name__ == "__main__":
    main()
