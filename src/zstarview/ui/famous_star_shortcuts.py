# -*- coding: utf-8 -*-
"""Helpers for famous-star jump shortcuts."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import polars as pl


DEC_BAND_NORTH = "north"
DEC_BAND_EQUATOR = "equatorial"
DEC_BAND_SOUTH = "south"
DEC_BANDS = (DEC_BAND_NORTH, DEC_BAND_EQUATOR, DEC_BAND_SOUTH)


@dataclass(frozen=True)
class FamousStarShortcut:
    name: str
    ra_hours: float
    dec_deg: float
    vmag: float
    band: str


def classify_declination_band(dec_deg: float) -> str:
    if dec_deg >= 20.0:
        return DEC_BAND_NORTH
    if dec_deg <= -20.0:
        return DEC_BAND_SOUTH
    return DEC_BAND_EQUATOR


def build_famous_star_shortcuts(star_catalog: pl.DataFrame, max_vmag: float = 2.0) -> Dict[str, List[FamousStarShortcut]]:
    """Build famous-star candidates grouped by declination band.

    Rules:
    - Name must be non-empty.
    - Vmag must be finite and <= max_vmag.
    - For duplicate names, keep the brightest entry (lowest Vmag).
    - Sort each band by Vmag asc, then name asc.
    """
    bands: Dict[str, List[FamousStarShortcut]] = {key: [] for key in DEC_BANDS}
    best_by_name: dict[str, FamousStarShortcut] = {}

    rows = (
        star_catalog.select(["Name", "RAh", "Dec", "Vmag"])
        .iter_rows(named=True)
    )
    for row in rows:
        name_raw = row.get("Name")
        name = str(name_raw).strip() if name_raw is not None else ""
        if not name:
            continue

        try:
            ra_h = float(row["RAh"])
            dec = float(row["Dec"])
            vmag = float(row["Vmag"])
        except (TypeError, ValueError):
            continue
        if vmag > float(max_vmag):
            continue

        band = classify_declination_band(dec)
        candidate = FamousStarShortcut(
            name=name,
            ra_hours=ra_h,
            dec_deg=dec,
            vmag=vmag,
            band=band,
        )
        current = best_by_name.get(name)
        if current is None or candidate.vmag < current.vmag:
            best_by_name[name] = candidate

    for star in best_by_name.values():
        bands[star.band].append(star)

    for key in DEC_BANDS:
        bands[key].sort(key=lambda s: (s.vmag, s.name.casefold()))
    return bands
