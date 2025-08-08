import csv
from typing import Any, Dict, List, Tuple

from astropy.coordinates import SkyCoord
import astropy.units as u


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


def load_star_catalog(filename: str) -> List[Dict[str, Any]]:
    """Loads the star catalog from a CSV file.

    Each row contains: name, SkyCoord, Vmag, B-V.
    """
    star_catalog: List[Dict[str, Any]] = []
    with open(filename, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                ra = float(row["RA"]) * 15
                name = row["Name"]
                dec = float(row["Dec"])
                vmag = float(row["Vmag"])
                bv = float(row["B-V"])
                coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
                star_catalog.append({"name": name, "coord": coord, "vmag": vmag, "bv": bv})
            except Exception:
                continue
    return star_catalog
