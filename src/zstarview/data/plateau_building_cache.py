"""Discover and validate prepared PLATEAU building caches."""

from __future__ import annotations

import json
from pathlib import Path

from zstarview.data.plateau_buildings import (
    CACHE_METADATA_SCHEMA_VERSION,
    DERIVED_TILE_SCHEMA_VERSION,
)
from zstarview.paths import PLATEAU_DERIVED_ROOT_DIR

CACHE_METADATA_FILENAME = "cache_meta.json"


def read_plateau_cache_metadata(dataset_dir: Path) -> dict[str, object] | None:
    path = dataset_dir / CACHE_METADATA_FILENAME
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def is_valid_plateau_cache(dataset_dir: Path) -> bool:
    from zstarview.data.derived_tile_cache import load_derived_tile_index

    metadata = read_plateau_cache_metadata(dataset_dir)
    if metadata is None:
        return False
    if metadata.get("status") != "complete":
        return False
    if metadata.get("source") != "PLATEAU-CityGML":
        return False
    if metadata.get("metadata_schema_version") != CACHE_METADATA_SCHEMA_VERSION:
        return False
    if metadata.get("derived_tile_schema_version") != DERIVED_TILE_SCHEMA_VERSION:
        return False
    derived_dir = dataset_dir / "bldg"
    indexed = load_derived_tile_index(derived_dir)
    if indexed is None:
        return False
    _city_envelope, tile_envelopes = indexed
    return bool(tile_envelopes) and all(
        envelope.path.is_file() for envelope in tile_envelopes
    )


def find_plateau_building_dirs(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float,
    root_dir: str | Path = PLATEAU_DERIVED_ROOT_DIR,
) -> tuple[Path, ...]:
    from zstarview.data.derived_tile_cache import select_derived_tile_envelopes

    root = Path(root_dir)
    if not root.is_dir():
        return ()
    matches: list[Path] = []
    for dataset_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        if not is_valid_plateau_cache(dataset_dir):
            continue
        derived_dir = dataset_dir / "bldg"
        try:
            select_derived_tile_envelopes(
                derived_dir,
                observer_lat_deg=float(observer_lat_deg),
                observer_lon_deg=float(observer_lon_deg),
                radius_km=float(radius_km),
            )
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        matches.append(derived_dir)
    return tuple(matches)
