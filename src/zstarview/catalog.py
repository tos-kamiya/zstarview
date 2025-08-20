from typing import Optional

import polars as pl


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
