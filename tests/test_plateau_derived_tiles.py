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


def test_parse_derived_tile_buildings_defaults_to_no_height_filter(tmp_path: Path) -> None:
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

    buildings = mod.parse_derived_tile_buildings(path)

    assert len(buildings) == 2
    assert [building.building_id for building in buildings] == ["high", "low"]


def test_parse_derived_tile_buildings_can_filter_low_buildings(tmp_path: Path) -> None:
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


def test_select_derived_tile_envelopes_prefers_tile_index(tmp_path: Path) -> None:
    mod = _load_module()
    derived_dir = tmp_path / "13100_tokyo23" / "bldg"
    derived_dir.mkdir(parents=True)
    (derived_dir / "tile_index.json").write_text(
        json.dumps(
            {
                "bbox": {
                    "min_lat": 35.0,
                    "min_lon": 139.0,
                    "max_lat": 35.2,
                    "max_lon": 139.2,
                },
                "tiles": [
                    {
                        "path": "near.json",
                        "bbox": {
                            "min_lat": 35.0,
                            "min_lon": 139.0,
                            "max_lat": 35.01,
                            "max_lon": 139.01,
                        },
                    },
                    {
                        "path": "far.json",
                        "bbox": {
                            "min_lat": 35.2,
                            "min_lon": 139.2,
                            "max_lat": 35.21,
                            "max_lon": 139.21,
                        },
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    envelopes = mod.select_derived_tile_envelopes(
        derived_dir,
        observer_lat_deg=35.005,
        observer_lon_deg=139.005,
        radius_km=3.0,
    )

    assert [envelope.path.name for envelope in envelopes] == ["near.json"]


def test_select_derived_tile_envelopes_uses_city_bbox_to_reject(tmp_path: Path) -> None:
    mod = _load_module()
    derived_dir = tmp_path / "26100_kyoto" / "bldg"
    derived_dir.mkdir(parents=True)
    (derived_dir / "tile_index.json").write_text(
        json.dumps(
            {
                "bbox": {
                    "min_lat": 35.0,
                    "min_lon": 135.0,
                    "max_lat": 35.1,
                    "max_lon": 135.1,
                },
                "tiles": [],
            }
        ),
        encoding="utf-8",
    )

    try:
        mod.select_derived_tile_envelopes(
            derived_dir,
            observer_lat_deg=34.69,
            observer_lon_deg=135.50,
            radius_km=3.0,
        )
    except ValueError as exc:
        assert "No derived building tiles found" in str(exc)
    else:
        raise AssertionError("Expected ValueError when observer is outside indexed city bbox")
