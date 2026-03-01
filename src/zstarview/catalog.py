import math
from pathlib import Path
from typing import List, Optional

import polars as pl


def _split_files_for_threshold(filename: str, vmag_threshold: Optional[float]) -> Optional[List[Path]]:
    """Return split catalog files for threshold-based loading, if available."""
    if vmag_threshold is None or not math.isfinite(float(vmag_threshold)):
        return None

    base_csv = Path(filename)
    # Only apply split-loading policy to canonical stars.csv.
    if base_csv.name != "stars.csv":
        return None

    split_dir = base_csv.parent / "stars"
    files: List[Path] = [split_dir / "stars_base.csv"]

    t = float(vmag_threshold)
    if t > 6.0:
        files.append(split_dir / "stars_extra7.csv")
    if t > 7.0:
        files.append(split_dir / "stars_extra8.csv")
    if t > 8.0:
        files.append(split_dir / "stars_extra9.csv")
    if t > 9.0:
        files.append(split_dir / "stars_extra10.csv")

    # For >10.0 we currently rely on the unified stars.csv.
    if t > 10.0:
        return None

    return files if all(p.exists() for p in files) else None


def load_star_catalog(filename: str, vmag_threshold: Optional[float] = 7.0) -> pl.DataFrame:
    """Loads the star catalog from a CSV file using Polars.

    If vmag_threshold is not None, keeps only rows with Vmag <= threshold.
    Returns a Polars DataFrame.
    """
    split_files = _split_files_for_threshold(filename, vmag_threshold)
    if split_files is not None:
        parts = [pl.read_csv(str(p), try_parse_dates=False, null_values="").fill_null("") for p in split_files]
        df = pl.concat(parts, how="vertical_relaxed")
    else:
        # Use fill_null to handle empty strings for name, etc.
        df = pl.read_csv(filename, try_parse_dates=False, null_values="").fill_null("")
    if vmag_threshold is not None:
        # Vmag can be empty string, cast to float handles this (becomes null)
        # then filter out nulls and values > threshold
        vmag_col = pl.col("Vmag").cast(pl.Float64, strict=False)
        df = df.filter((vmag_col.is_not_null()) & (vmag_col <= vmag_threshold))
    return df
