from typing import Dict, List, Optional, Tuple

import polars as pl


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


def load_star_catalog(filename: str, vmag_threshold: Optional[float] = 7.0) -> pl.DataFrame:
    """Loads the star catalog from a CSV file using Polars.

    If vmag_threshold is not None, keeps only rows with Vmag <= threshold.
    Returns a Polars DataFrame.
    """
    # Use fill_null to handle empty strings for name, etc.
    df = pl.read_csv(filename, try_parse_dates=False, null_values="").fill_null("")
    if vmag_threshold is not None:
        # Vmag can be empty string, cast to float handles this (becomes null)
        # then filter out nulls and values > threshold
        vmag_col = pl.col("Vmag").cast(pl.Float64, strict=False)
        df = df.filter((vmag_col.is_not_null()) & (vmag_col <= vmag_threshold))
    return df
