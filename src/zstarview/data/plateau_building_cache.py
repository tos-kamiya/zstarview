"""Discover and validate prepared PLATEAU building caches.

Cache layout under the PLATEAU root directory:

- ``<code>/``: current layout, one generation per five-digit municipality code.
- ``<code>_<yyyy>/``: legacy layout.  Read only when ``<code>/`` is not a valid
  cache, and then only the valid one with the largest year.
- Anything else that belongs to a municipality (``.<code>-*`` staging or
  retired directories, ``<code>_<year>.outdated-*``, ``.replaced-*``,
  ``.incomplete-*``, ``<code>_latest``, unselected legacy years, invalid
  caches) is never read and can be removed.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

from zstarview.data.plateau_buildings import (
    CACHE_METADATA_SCHEMA_VERSION,
    DERIVED_TILE_SCHEMA_VERSION,
)
from zstarview.paths import PLATEAU_DERIVED_ROOT_DIR

CACHE_METADATA_FILENAME = "cache_meta.json"

_CURRENT_DIR_RE = re.compile(r"^(\d{5})$")
_LEGACY_DIR_RE = re.compile(r"^(\d{5})_(\d{4})$")
_LEGACY_OTHER_RE = re.compile(r"^(\d{5})_")
_DOT_DIR_RE = re.compile(r"^\.(\d{5})-")

logger = logging.getLogger(__name__)
_logged_messages: set[str] = set()


@dataclass(frozen=True)
class PlateauCacheEntry:
    """One municipality's selected cache and its removable leftovers."""

    city_code: str
    selected: Path | None
    legacy: bool
    removable: tuple[Path, ...] = field(default_factory=tuple)


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
    try:
        indexed = load_derived_tile_index(derived_dir)
    except (OSError, ValueError):
        return False
    if indexed is None:
        return False
    _city_envelope, tile_envelopes = indexed
    return bool(tile_envelopes) and all(
        envelope.path.is_file() for envelope in tile_envelopes
    )


def plateau_dir_city_code(name: str) -> str | None:
    """Return the municipality code a cache-root entry name belongs to."""

    for pattern in (_CURRENT_DIR_RE, _LEGACY_OTHER_RE, _DOT_DIR_RE):
        match = pattern.match(name)
        if match is not None:
            return match.group(1)
    return None


def plateau_city_dirs(root_dir: str | Path, city_code: str) -> tuple[Path, ...]:
    """List every directory under the root that belongs to ``city_code``."""

    root = Path(root_dir)
    if not root.is_dir():
        return ()
    return tuple(
        sorted(
            path
            for path in root.iterdir()
            if path.is_dir() and plateau_dir_city_code(path.name) == city_code
        )
    )


def classify_plateau_caches(
    root_dir: str | Path = PLATEAU_DERIVED_ROOT_DIR,
) -> tuple[PlateauCacheEntry, ...]:
    """Select at most one readable cache per municipality code."""

    root = Path(root_dir)
    if not root.is_dir():
        return ()
    grouped: dict[str, list[Path]] = {}
    for path in root.iterdir():
        if not path.is_dir():
            continue
        city_code = plateau_dir_city_code(path.name)
        if city_code is not None:
            grouped.setdefault(city_code, []).append(path)
    entries: list[PlateauCacheEntry] = []
    for city_code in sorted(grouped):
        paths = sorted(grouped[city_code])
        selected: Path | None = None
        legacy = False
        current = root / city_code
        if current in paths and is_valid_plateau_cache(current):
            selected = current
        else:
            legacy_dirs = sorted(
                (
                    path
                    for path in paths
                    if _LEGACY_DIR_RE.match(path.name) is not None
                ),
                key=lambda path: path.name,
                reverse=True,
            )
            for path in legacy_dirs:
                if is_valid_plateau_cache(path):
                    selected = path
                    legacy = True
                    break
        removable = tuple(path for path in paths if path != selected)
        entries.append(
            PlateauCacheEntry(
                city_code=city_code,
                selected=selected,
                legacy=legacy,
                removable=removable,
            )
        )
    return tuple(entries)


def _log_once(message: str) -> None:
    if message in _logged_messages:
        return
    _logged_messages.add(message)
    logger.info(message)


def find_plateau_building_dirs(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float,
    root_dir: str | Path = PLATEAU_DERIVED_ROOT_DIR,
) -> tuple[Path, ...]:
    from zstarview.data.derived_tile_cache import select_derived_tile_envelopes

    matches: list[Path] = []
    for entry in classify_plateau_caches(root_dir):
        if entry.removable:
            _log_once(
                f"PLATEAU {entry.city_code}: {len(entry.removable)} unused cache "
                "director"
                + ("y" if len(entry.removable) == 1 else "ies")
                + " can be removed; see 'zstarview-download-plateau-buildings list'"
            )
        if entry.selected is None:
            continue
        derived_dir = entry.selected / "bldg"
        try:
            select_derived_tile_envelopes(
                derived_dir,
                observer_lat_deg=float(observer_lat_deg),
                observer_lon_deg=float(observer_lon_deg),
                radius_km=float(radius_km),
            )
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        if entry.legacy:
            _log_once(
                f"PLATEAU {entry.city_code}: using legacy cache {entry.selected}; "
                "run 'zstarview-download-plateau-buildings prepare "
                f"{entry.city_code}' to migrate it"
            )
        matches.append(derived_dir)
    return tuple(matches)
