from __future__ import annotations

import json
from pathlib import Path

from zstarview.data.skyscraper_tiles import (
    load_skyscraper_seed_tiles,
    select_skyscraper_seed_tiles_for_viewer,
    skyscraper_tile_derived_dir,
)


def test_select_skyscraper_seed_tiles_for_viewer_respects_inner_and_outer_rings(tmp_path: Path) -> None:
    seed_path = tmp_path / "skyscraper_tiles.json"
    seed_path.write_text(
        json.dumps(
            {
                "tiles": [
                    {
                        "tile": {"z": 14, "x": 1, "y": 1},
                        "bbox": {"west": 139.70, "south": 35.60, "east": 139.72, "north": 35.62},
                        "building_count": 3,
                        "max_height_m": 220,
                    },
                    {
                        "tile": {"z": 14, "x": 2, "y": 2},
                        "bbox": {"west": 139.78, "south": 35.67, "east": 139.80, "north": 35.69},
                        "building_count": 4,
                        "max_height_m": 260,
                    },
                    {
                        "tile": {"z": 14, "x": 3, "y": 3},
                        "bbox": {"west": 140.20, "south": 35.95, "east": 140.22, "north": 35.97},
                        "building_count": 2,
                        "max_height_m": 210,
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    load_skyscraper_seed_tiles.cache_clear()
    got = select_skyscraper_seed_tiles_for_viewer(
        observer_lat_deg=35.68,
        observer_lon_deg=139.79,
        inner_radius_km=2.5,
        outer_radius_km=10.0,
        seed_file=seed_path,
    )

    assert [(tile.x, tile.y) for tile in got] == [(1, 1)]


def test_skyscraper_tile_derived_dir_builds_expected_cache_path(tmp_path: Path) -> None:
    seed_path = tmp_path / "skyscraper_tiles.json"
    seed_path.write_text(
        json.dumps(
            {
                "tiles": [
                    {
                        "tile": {"z": 14, "x": 14525, "y": 6453},
                        "bbox": {"west": 139.0, "south": 35.0, "east": 139.1, "north": 35.1},
                        "building_count": 1,
                        "max_height_m": 200,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    load_skyscraper_seed_tiles.cache_clear()
    tile = load_skyscraper_seed_tiles(seed_path)[0]

    got = skyscraper_tile_derived_dir(tile, derived_root_dir=tmp_path / "cache")

    assert got == tmp_path / "cache" / "z14_14525_6453" / "bldg"
