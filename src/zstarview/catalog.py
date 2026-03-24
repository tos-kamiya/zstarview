import math
from pathlib import Path
from typing import List, Optional

import polars as pl

from .asterisms import ASTERISM_REQUIRED_SOURCE_IDS


def _resolve_split_dir(filename: str) -> Optional[Path]:
    path = Path(filename)
    if path.name == "stars.csv":
        return path.parent / "stars"
    if path.name == "stars_base.csv":
        return path.parent
    if path.is_dir():
        return path
    return None


def _split_files_for_threshold(filename: str, vmag_threshold: Optional[float]) -> Optional[List[Path]]:
    """Return split catalog files for threshold-based loading, if available."""
    split_dir = _resolve_split_dir(filename)
    if split_dir is None:
        return None

    base = split_dir / "stars_base.csv"
    extra7 = split_dir / "stars_extra7.csv"
    extra8 = split_dir / "stars_extra8.csv"
    extra9 = split_dir / "stars_extra9.csv"
    extra10 = split_dir / "stars_extra10.csv"
    extra_faint = split_dir / "stars_extra_faint.csv"

    if vmag_threshold is None or not math.isfinite(float(vmag_threshold)):
        all_files = [p for p in (base, extra7, extra8, extra9, extra10, extra_faint) if p.exists()]
        return all_files if all_files else None

    t = float(vmag_threshold)
    selected_files: List[Path] = [base]
    required_buckets: List[Path] = []
    if t > 6.0:
        required_buckets.append(extra7)
    if t > 7.0:
        required_buckets.append(extra8)
    if t > 8.0:
        required_buckets.append(extra9)
    if t > 9.0:
        required_buckets.append(extra10)

    if not base.exists() or any(not p.exists() for p in required_buckets):
        return None

    selected_files.extend(required_buckets)
    if t > 10.0 and extra_faint.exists():
        selected_files.append(extra_faint)

    return selected_files


def _all_split_files(filename: str) -> List[Path]:
    split_dir = _resolve_split_dir(filename)
    if split_dir is None:
        return []
    files = [
        split_dir / "stars_base.csv",
        split_dir / "stars_extra7.csv",
        split_dir / "stars_extra8.csv",
        split_dir / "stars_extra9.csv",
        split_dir / "stars_extra10.csv",
        split_dir / "stars_extra_faint.csv",
    ]
    return [p for p in files if p.exists()]


def _read_star_csv(path: Path) -> pl.DataFrame:
    return pl.read_csv(str(path), try_parse_dates=False, null_values="").fill_null("")


def _append_required_asterism_rows(df: pl.DataFrame, filename: str, selected_files: List[Path]) -> pl.DataFrame:
    if not ASTERISM_REQUIRED_SOURCE_IDS:
        return df
    present_ids = set(df.get_column("SourceId").to_list()) if "SourceId" in df.columns else set()
    needed_ids = set(ASTERISM_REQUIRED_SOURCE_IDS) - present_ids
    if not needed_ids:
        return df

    extra_parts: List[pl.DataFrame] = []
    for path in _all_split_files(filename):
        if path in selected_files or not needed_ids:
            continue
        part = _read_star_csv(path).filter(pl.col("SourceId").is_in(sorted(needed_ids)))
        if part.height == 0:
            continue
        extra_parts.append(part)
        needed_ids -= set(part.get_column("SourceId").to_list())

    if not extra_parts:
        return df
    return pl.concat([df, *extra_parts], how="vertical_relaxed")


def load_star_catalog(filename: str, vmag_threshold: Optional[float] = 7.0) -> pl.DataFrame:
    """Loads the star catalog from a CSV file using Polars.

    If vmag_threshold is not None, keeps only rows with Vmag <= threshold.
    Returns a Polars DataFrame.
    """
    split_files = _split_files_for_threshold(filename, vmag_threshold)
    if split_files is not None:
        parts = [_read_star_csv(p) for p in split_files]
        df = pl.concat(parts, how="vertical_relaxed")
        df = _append_required_asterism_rows(df, filename, split_files)
    else:
        # Use fill_null to handle empty strings for name, etc.
        df = pl.read_csv(filename, try_parse_dates=False, null_values="").fill_null("")
    if vmag_threshold is not None:
        # Vmag can be empty string, cast to float handles this (becomes null)
        # then filter out nulls and values > threshold
        vmag_col = pl.col("Vmag").cast(pl.Float64, strict=False)
        required_mask = pl.col("SourceId").is_in(sorted(ASTERISM_REQUIRED_SOURCE_IDS)) if ASTERISM_REQUIRED_SOURCE_IDS else pl.lit(False)
        df = df.filter((vmag_col.is_not_null()) & ((vmag_col <= vmag_threshold) | required_mask))
    return df


def load_dso_catalog(filename: str) -> pl.DataFrame:
    """Load DSO catalog CSV and keep rows that have valid RA/Dec."""
    df = pl.read_csv(filename, try_parse_dates=False, null_values="").fill_null("")
    ra = pl.col("RAh").cast(pl.Float64, strict=False)
    dec = pl.col("Dec").cast(pl.Float64, strict=False)
    return df.filter(ra.is_not_null() & dec.is_not_null())
