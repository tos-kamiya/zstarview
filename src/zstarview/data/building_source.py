"""Select prepared building data sources for runtime consumers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from zstarview.data.plateau_building_cache import find_plateau_building_dirs
from zstarview.paths import PLATEAU_DERIVED_ROOT_DIR


@dataclass(frozen=True)
class BuildingSourceSelection:
    source: str
    derived_dirs: tuple[Path, ...]


def select_prepared_building_source(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float,
    plateau_root_dir: str | Path = PLATEAU_DERIVED_ROOT_DIR,
) -> BuildingSourceSelection:
    plateau_dirs = find_plateau_building_dirs(
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        radius_km=radius_km,
        root_dir=plateau_root_dir,
    )
    if plateau_dirs:
        return BuildingSourceSelection(source="plateau", derived_dirs=plateau_dirs)
    return BuildingSourceSelection(source="overture", derived_dirs=())
