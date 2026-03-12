from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "plateau_derived_tiles.py"
    spec = importlib.util.spec_from_file_location("plateau_derived_tiles", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_derived_tile_envelope_reads_bbox(tmp_path: Path) -> None:
    mod = _load_module()
    path = tmp_path / "tile.json"
    path.write_text(
        json.dumps(
            {
                "tile": {
                    "id": "tile",
                    "bbox": {
                        "min_lat": 35.0,
                        "min_lon": 139.0,
                        "max_lat": 35.1,
                        "max_lon": 139.1,
                    },
                },
                "buildings": [],
            }
        ),
        encoding="utf-8",
    )

    envelope = mod.load_derived_tile_envelope(path)

    assert envelope is not None
    assert envelope.min_lat_deg == 35.0
    assert envelope.max_lon_deg == 139.1


def test_parse_derived_tile_buildings_filters_low_buildings(tmp_path: Path) -> None:
    mod = _load_module()
    path = tmp_path / "tile.json"
    path.write_text(
        json.dumps(
            {
                "buildings": [
                    {
                        "id": "high",
                        "height_m": 80.0,
                        "rings": [
                            [
                                [139.0, 35.0],
                                [139.001, 35.0],
                                [139.001, 35.001],
                                [139.0, 35.001],
                                [139.0, 35.0],
                            ]
                        ],
                    },
                    {
                        "id": "low",
                        "height_m": 20.0,
                        "rings": [
                            [
                                [139.0, 35.0],
                                [139.001, 35.0],
                                [139.001, 35.001],
                                [139.0, 35.001],
                                [139.0, 35.0],
                            ]
                        ],
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    buildings = mod.parse_derived_tile_buildings(path, min_building_height_m=40.0)

    assert len(buildings) == 1
    assert buildings[0].building_id == "high"
    assert buildings[0].height_m == 80.0
