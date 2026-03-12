from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "build_plateau_tile_index.py"
    spec = importlib.util.spec_from_file_location("build_plateau_tile_index", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_tile_index_payload_collects_bbox_and_tiles(tmp_path: Path) -> None:
    mod = _load_module()
    derived_dir = tmp_path / "13100_tokyo23" / "bldg"
    derived_dir.mkdir(parents=True)
    for name, bbox in (
        ("a.json", {"min_lat": 35.0, "min_lon": 139.0, "max_lat": 35.1, "max_lon": 139.1}),
        ("b.json", {"min_lat": 35.2, "min_lon": 139.2, "max_lat": 35.3, "max_lon": 139.3}),
    ):
        (derived_dir / name).write_text(
            json.dumps(
                {
                    "tile": {"id": name[:-5], "bbox": bbox},
                    "filters": {"min_height_m": 40.0},
                    "buildings": [{"id": "x"}],
                }
            ),
            encoding="utf-8",
        )

    payload = mod.build_tile_index_payload(derived_dir)

    assert payload["schema_version"] == 1
    assert payload["city_code"] == "13100"
    assert payload["city_name"] == "tokyo23"
    assert payload["tile_count"] == 2
    assert payload["bbox"] == {
        "min_lat": 35.0,
        "min_lon": 139.0,
        "max_lat": 35.3,
        "max_lon": 139.3,
    }
    assert payload["tiles"][0]["path"] == "a.json"


def test_main_writes_tile_index_json(tmp_path: Path) -> None:
    mod = _load_module()
    derived_dir = tmp_path / "26100_kyoto" / "bldg"
    derived_dir.mkdir(parents=True)
    (derived_dir / "tile.json").write_text(
        json.dumps(
            {
                "tile": {
                    "id": "tile",
                    "bbox": {
                        "min_lat": 35.0,
                        "min_lon": 135.0,
                        "max_lat": 35.1,
                        "max_lon": 135.1,
                    },
                },
                "filters": {"min_height_m": 40.0},
                "buildings": [{"id": "b1"}, {"id": "b2"}],
            }
        ),
        encoding="utf-8",
    )

    rc = mod.main(["--derived-dir", str(derived_dir)])

    assert rc == 0
    payload = json.loads((derived_dir / "tile_index.json").read_text(encoding="utf-8"))
    assert payload["city_code"] == "26100"
    assert payload["tile_count"] == 1
    assert payload["tiles"][0]["building_count"] == 2
