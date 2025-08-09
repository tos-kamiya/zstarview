import csv
from typing import Dict, List, Optional, Tuple

from astropy.coordinates import SkyCoord
import astropy.units as u

from .types import StarRecord


def load_city_coords(filename: str) -> Dict[str, Tuple[float, float, str]]:
    """Loads city coordinates and timezone from the data file.

    Returns a dict keyed by "{country}/{name}" (lowercase) -> (lat, lon, tz).
    """
    city_table: Dict[str, Tuple[float, float, str]] = {}
    with open(filename, encoding="utf-8") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 18:
                continue
            name = cols[1]
            lat = float(cols[4])
            lon = float(cols[5])
            country = cols[8]
            timezone_name = cols[17]
            key = f"{country.lower()}/{name.lower()}"
            city_table[key] = (lat, lon, timezone_name)
    return city_table


def load_star_catalog(filename: str, vmag_threshold: Optional[float] = 7.0) -> List[StarRecord]:
    """Loads the star catalog from a CSV file.

    Each row contains: name, SkyCoord, Vmag, B-V.
    If vmag_threshold is not None, keeps only rows with Vmag <= threshold.
    """
    result: List[StarRecord] = []
    with open(filename, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                v_raw = row.get("Vmag")
                if v_raw is None or v_raw == "":
                    continue
                v = float(v_raw)
                if (vmag_threshold is not None) and (v > vmag_threshold):
                    continue

                ra_h = float(row["RAh"])          # RAh（時間）想定
                dec = float(row["Dec"])
                bv = float(row["B-V"]) if row.get("B-V") else float("nan")
                name = row.get("Name", "")
                coord = SkyCoord(ra=(ra_h * 15.0) * u.deg, dec=dec * u.deg, frame="icrs")
                result.append(StarRecord(name=name, coord=coord, vmag=v, bv=bv))
            except Exception:
                continue
    return result

