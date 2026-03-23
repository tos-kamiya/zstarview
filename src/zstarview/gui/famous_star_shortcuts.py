# -*- coding: utf-8 -*-
"""Helpers for named-star jump shortcuts."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, List, Optional

import polars as pl

from ..asterisms import ASTERISMS


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
    kind: str = "star"
    subtitle: str = ""


@dataclass(frozen=True)
class SearchJumpTarget:
    label: str
    ra_hours: float
    dec_deg: float
    kind: str
    sort_key: tuple[float, str]
    subtitle: str = ""
    object_key: str = ""


SATELLITE_JUMP_SHORTCUTS = (
    NamedStarShortcut(
        name="ISS",
        ra_hours=0.0,
        dec_deg=0.0,
        vmag=99.0,
        band=DEC_BAND_EQUATOR,
        kind="satellite",
        subtitle="Satellite",
    ),
)


def classify_declination_band(dec_deg: float) -> str:
    if dec_deg >= 20.0:
        return DEC_BAND_NORTH
    if dec_deg <= -20.0:
        return DEC_BAND_SOUTH
    return DEC_BAND_EQUATOR


def build_named_star_shortcuts(
    star_catalog: pl.DataFrame,
    max_vmag: Optional[float] = 2.0,
    *,
    include_satellites: bool = False,
) -> Dict[str, List[NamedStarShortcut]]:
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
    if include_satellites:
        for shortcut in SATELLITE_JUMP_SHORTCUTS:
            bands[shortcut.band].append(shortcut)

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


def _circular_mean_hours(values: List[float]) -> float:
    if not values:
        return 0.0
    angles = [math.tau * (value / 24.0) for value in values]
    x = sum(math.cos(angle) for angle in angles)
    y = sum(math.sin(angle) for angle in angles)
    if x == 0.0 and y == 0.0:
        return values[0] % 24.0
    angle = math.atan2(y, x)
    if angle < 0.0:
        angle += math.tau
    return (angle / math.tau) * 24.0


def build_search_jump_targets(
    star_catalog: pl.DataFrame,
    *,
    include_satellites: bool = False,
) -> List[SearchJumpTarget]:
    """Build search targets for both named stars and asterisms."""
    targets: List[SearchJumpTarget] = []

    for star in flatten_named_star_shortcuts(
        build_named_star_shortcuts(star_catalog, max_vmag=None, include_satellites=include_satellites)
    ):
        targets.append(
            SearchJumpTarget(
                label=star.name,
                ra_hours=star.ra_hours,
                dec_deg=star.dec_deg,
                kind=star.kind,
                subtitle=star.subtitle or f"Vmag {star.vmag:.2f}",
                sort_key=(star.vmag, star.name.casefold()),
                object_key=star.name if star.kind == "satellite" else "",
            )
        )

    rows = star_catalog.select(["SourceId", "RAh", "Dec"]).iter_rows(named=True)
    by_source_id: dict[str, tuple[float, float]] = {}
    for row in rows:
        source_id_raw = row.get("SourceId")
        source_id = str(source_id_raw).strip() if source_id_raw is not None else ""
        if not source_id:
            continue
        try:
            by_source_id[source_id] = (float(row["RAh"]), float(row["Dec"]))
        except (TypeError, ValueError):
            continue

    for asterism in ASTERISMS:
        coords = [by_source_id[source_id] for source_id in asterism.source_ids if source_id in by_source_id]
        if not coords:
            continue
        ra_values = [ra for ra, _ in coords]
        dec_values = [dec for _, dec in coords]
        targets.append(
            SearchJumpTarget(
                label=asterism.name,
                ra_hours=_circular_mean_hours(ra_values),
                dec_deg=sum(dec_values) / float(len(dec_values)),
                kind="asterism",
                subtitle="Asterism",
                sort_key=(99.0, asterism.name.casefold()),
            )
        )

    targets.sort(key=lambda t: t.sort_key)
    return targets
