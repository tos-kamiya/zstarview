from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .derived_tile_cache import TileEnvelope, envelope_max_distance_km, envelope_min_distance_km
from ..paths import OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR, SKYSCRAPER_TILES_FILE


@dataclass(frozen=True)
class SkyscraperSeedTile:
    zoom: int
    x: int
    y: int
    envelope: TileEnvelope
    building_count: int
    max_height_m: float

    @property
    def cache_key(self) -> str:
        return f"z{self.zoom}_{self.x}_{self.y}"

    def derived_dir(self, root_dir: Path) -> Path:
        return root_dir / self.cache_key / "bldg"


@lru_cache(maxsize=1)
def load_skyscraper_seed_tiles(
    seed_file: str | Path = SKYSCRAPER_TILES_FILE,
) -> tuple[SkyscraperSeedTile, ...]:
    payload = json.loads(Path(seed_file).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return ()
    raw_tiles = payload.get("tiles")
    if not isinstance(raw_tiles, list):
        return ()

    tiles: list[SkyscraperSeedTile] = []
    for row in raw_tiles:
        if not isinstance(row, dict):
            continue
        raw_tile = row.get("tile")
        raw_bbox = row.get("bbox")
        if not isinstance(raw_tile, dict) or not isinstance(raw_bbox, dict):
            continue
        try:
            zoom = int(raw_tile["z"])
            x = int(raw_tile["x"])
            y = int(raw_tile["y"])
            west = float(raw_bbox["west"])
            south = float(raw_bbox["south"])
            east = float(raw_bbox["east"])
            north = float(raw_bbox["north"])
        except (KeyError, TypeError, ValueError):
            continue
        tiles.append(
            SkyscraperSeedTile(
                zoom=zoom,
                x=x,
                y=y,
                envelope=TileEnvelope(
                    path=Path(f"{zoom}/{x}/{y}"),
                    min_lat_deg=south,
                    min_lon_deg=west,
                    max_lat_deg=north,
                    max_lon_deg=east,
                ),
                building_count=int(row.get("building_count", 0) or 0),
                max_height_m=float(row.get("max_height_m", 0.0) or 0.0),
            )
        )
    return tuple(tiles)


def select_skyscraper_seed_tiles_for_viewer(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    inner_radius_km: float,
    outer_radius_km: float,
    seed_file: str | Path = SKYSCRAPER_TILES_FILE,
) -> tuple[SkyscraperSeedTile, ...]:
    if inner_radius_km < 0.0:
        raise ValueError("inner_radius_km must be non-negative")
    if outer_radius_km <= inner_radius_km:
        raise ValueError("outer_radius_km must be greater than inner_radius_km")

    selected: list[SkyscraperSeedTile] = []
    for tile in load_skyscraper_seed_tiles(seed_file):
        min_distance_km = envelope_min_distance_km(
            tile.envelope,
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
        )
        if min_distance_km > outer_radius_km:
            continue
        max_distance_km = envelope_max_distance_km(
            tile.envelope,
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
        )
        if max_distance_km <= inner_radius_km:
            continue
        selected.append(tile)
    return tuple(selected)


def skyscraper_tile_derived_dir(
    tile: SkyscraperSeedTile,
    *,
    derived_root_dir: str | Path = OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR,
) -> Path:
    return tile.derived_dir(Path(derived_root_dir))
