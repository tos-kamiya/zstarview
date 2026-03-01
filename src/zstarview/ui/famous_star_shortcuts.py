# -*- coding: utf-8 -*-
"""Helpers for named-star jump shortcuts."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import polars as pl


DEC_BAND_NORTH = "north"
DEC_BAND_EQUATOR = "equatorial"
DEC_BAND_SOUTH = "south"
DEC_BANDS = (DEC_BAND_NORTH, DEC_BAND_EQUATOR, DEC_BAND_SOUTH)


@dataclass(frozen=True)
class NamedStarShortcut:
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


def build_named_star_shortcuts(star_catalog: pl.DataFrame, max_vmag: Optional[float] = 2.0) -> Dict[str, List[NamedStarShortcut]]:
    """Build named-star candidates grouped by declination band.

    Rules:
    - Name must be non-empty.
    - Vmag must be finite and <= max_vmag.
    - For duplicate names, keep the brightest entry (lowest Vmag).
    - Sort each band by Vmag asc, then name asc.
    """
    bands: Dict[str, List[NamedStarShortcut]] = {key: [] for key in DEC_BANDS}
    best_by_name: dict[str, NamedStarShortcut] = {}

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
        if max_vmag is not None and vmag > float(max_vmag):
            continue

        band = classify_declination_band(dec)
        candidate = NamedStarShortcut(
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


def flatten_named_star_shortcuts(stars_by_band: Dict[str, List[NamedStarShortcut]]) -> List[NamedStarShortcut]:
    """Flatten grouped shortcuts to a single list sorted by brightness then name."""
    out: List[NamedStarShortcut] = []
    for key in DEC_BANDS:
        out.extend(stars_by_band.get(key, []))
    out.sort(key=lambda s: (s.vmag, s.name.casefold()))
    return out
