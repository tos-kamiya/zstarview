from __future__ import annotations

import importlib.util
import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[2]
    mod_path = root / "src" / "zstarview" / "data" / "import_overture_buildings.py"
    spec = importlib.util.spec_from_file_location("import_overture_buildings", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_bbox_from_point_covers_requested_radius() -> None:
    mod = _load_module()

    west, south, east, north = mod.bbox_from_point(35.6896, 139.7006, 3.0)

    assert west < 139.7006 < east
    assert south < 35.6896 < north
    assert round(east - west, 4) > 0.05


def test_convert_feature_to_building_uses_height_and_geometry() -> None:
    mod = _load_module()

    feature = {
        "id": "building-1",
        "geometry": {
            "type": "Polygon",
            "coordinates": [
                [
                    [139.0, 35.0],
                    [139.001, 35.0],
                    [139.001, 35.001],
                    [139.0, 35.001],
                    [139.0, 35.0],
                ]
            ],
        },
        "properties": {
            "height": 55.0,
        },
    }

    got = mod.convert_feature_to_building(feature, min_building_height_m=40.0)

    assert got is not None
    assert got["id"] == "building-1"
    assert got["height_m"] == 55.0
    assert got["height_source"] == "height"
    assert got["parent_building_id"] is None


def test_convert_feature_to_building_can_fall_back_to_num_floors() -> None:
    mod = _load_module()

    feature = {
        "id": "building-2",
        "geometry": {
            "type": "Polygon",
            "coordinates": [
                [
                    [139.0, 35.0],
                    [139.001, 35.0],
                    [139.001, 35.001],
                    [139.0, 35.001],
                    [139.0, 35.0],
                ]
            ],
        },
        "properties": {
            "num_floors": 10,
        },
    }

    got = mod.convert_feature_to_building(feature, min_building_height_m=0.0)

    assert got is not None
    assert got["height_m"] == 35.0
    assert got["height_source"] == "num_floors*3.5"


def test_convert_feature_to_building_keeps_parent_building_id_for_parts() -> None:
    mod = _load_module()

    feature = {
        "id": "part-1",
        "geometry": {
            "type": "Polygon",
            "coordinates": [
                [
                    [139.0, 35.0],
                    [139.0005, 35.0],
                    [139.0005, 35.0005],
                    [139.0, 35.0005],
                    [139.0, 35.0],
                ]
            ],
        },
        "properties": {
            "height": 30.0,
            "building_id": "building-1",
        },
    }

    got = mod.convert_feature_to_building(feature, min_building_height_m=0.0)

    assert got is not None
    assert got["parent_building_id"] == "building-1"


def test_main_imports_geojsonseq_download_into_derived_dir(tmp_path: Path, monkeypatch) -> None:
    mod = _load_module()
    derived_root = tmp_path / "derived-root"
    raw_download = tmp_path / "raw-download.geojsonseq"

    def fake_run(_command, check=False):
        assert check is False
        raw_download.write_text(
            json.dumps(
                {
                    "id": "building-1",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [139.0, 35.0],
                                [139.001, 35.0],
                                [139.001, 35.001],
                                [139.0, 35.001],
                                [139.0, 35.0],
                            ]
                        ],
                    },
                    "properties": {
                        "height": 55.0,
                    },
                }
            )
            + "\n",
            encoding="utf-8",
        )
        download_path = Path(_command[-1])
        download_path.write_text(raw_download.read_text(encoding="utf-8"), encoding="utf-8")
        return type("CompletedProcess", (), {"returncode": 0})()

    monkeypatch.setattr(mod.shutil, "which", lambda _name: "/usr/bin/overturemaps")
    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    rc = mod.main(
        [
            "--lat",
            "35.6896",
            "--lon",
            "139.7006",
            "--radius-km",
            "3",
            "--derived-root-dir",
            str(derived_root),
            "--dataset-name",
            "shinjuku_test",
        ]
    )

    assert rc == 0
    tile_path = derived_root / "shinjuku_test" / "bldg" / "shinjuku_test.json"
    index_path = derived_root / "shinjuku_test" / "bldg" / "tile_index.json"
    payload = json.loads(tile_path.read_text(encoding="utf-8"))
    index_payload = json.loads(index_path.read_text(encoding="utf-8"))
    metadata_payload = json.loads((derived_root / "shinjuku_test" / "cache_meta.json").read_text(encoding="utf-8"))
    assert payload["source"]["format"] == "Overture-Buildings"
    assert payload["buildings"][0]["height_m"] == 55.0
    assert index_payload["tile_count"] == 1
    assert metadata_payload["dataset_name"] == "shinjuku_test"
    assert "fetched_at_utc" in metadata_payload


def test_read_derived_dataset_fetched_at_utc_migrates_legacy_cache(tmp_path: Path) -> None:
    mod = _load_module()
    derived_dir = tmp_path / "dataset" / "bldg"
    derived_dir.mkdir(parents=True)
    now = datetime(2026, 3, 27, 1, 30, tzinfo=timezone.utc)

    fetched_at_utc = mod.read_derived_dataset_fetched_at_utc(derived_dir, now_utc=now)

    assert fetched_at_utc == now
    metadata_payload = json.loads((tmp_path / "dataset" / "cache_meta.json").read_text(encoding="utf-8"))
    assert metadata_payload["dataset_name"] == "dataset"
    assert metadata_payload["fetched_at_utc"] == now.isoformat()


def test_is_derived_dataset_stale_uses_fetched_at_utc(tmp_path: Path) -> None:
    mod = _load_module()
    derived_dir = tmp_path / "dataset" / "bldg"
    derived_dir.mkdir(parents=True)
    now = datetime(2026, 3, 27, 1, 30, tzinfo=timezone.utc)
    (tmp_path / "dataset" / "cache_meta.json").write_text(
        json.dumps(
            {
                "dataset_name": "dataset",
                "fetched_at_utc": (now - timedelta(days=31)).isoformat(),
            }
        ),
        encoding="utf-8",
    )

    assert mod.is_derived_dataset_stale(derived_dir, now_utc=now) is True
    assert mod.is_derived_dataset_stale(derived_dir, ttl_days=40, now_utc=now) is False
