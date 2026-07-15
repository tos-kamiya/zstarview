from __future__ import annotations

import json
from pathlib import Path

from zstarview.data.building_source import select_prepared_building_source
from zstarview.data.plateau_building_cache import (
    find_plateau_building_dirs,
    is_valid_plateau_cache,
)


def _write_cache(root: Path, *, status: str = "complete") -> Path:
    dataset_dir = root / "32201_2024"
    derived_dir = dataset_dir / "bldg"
    derived_dir.mkdir(parents=True)
    tile_name = "53394500_bldg_6697.json"
    tile_payload = {
        "schema_version": 1,
        "tile": {
            "id": "53394500_bldg_6697",
            "bbox": {
                "min_lat": 35.45,
                "min_lon": 133.05,
                "max_lat": 35.46,
                "max_lon": 133.06,
            },
        },
        "buildings": [],
    }
    (derived_dir / tile_name).write_text(json.dumps(tile_payload), encoding="utf-8")
    (derived_dir / "tile_index.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "bbox": {
                    "min_lat": 35.45,
                    "min_lon": 133.05,
                    "max_lat": 35.46,
                    "max_lon": 133.06,
                },
                "tiles": [
                    {
                        "path": tile_name,
                        "bbox": {
                            "min_lat": 35.45,
                            "min_lon": 133.05,
                            "max_lat": 35.46,
                            "max_lon": 133.06,
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    (dataset_dir / "cache_meta.json").write_text(
        json.dumps(
            {
                "metadata_schema_version": 1,
                "derived_tile_schema_version": 3,
                "source": "PLATEAU-CityGML",
                "status": status,
            }
        ),
        encoding="utf-8",
    )
    return dataset_dir


def test_find_plateau_building_dirs_selects_cache_covering_observer(
    tmp_path: Path,
) -> None:
    dataset_dir = _write_cache(tmp_path)

    assert is_valid_plateau_cache(dataset_dir)
    assert find_plateau_building_dirs(
        observer_lat_deg=35.455,
        observer_lon_deg=133.055,
        radius_km=2.0,
        root_dir=tmp_path,
    ) == (dataset_dir / "bldg",)


def test_find_plateau_building_dirs_ignores_incomplete_and_outside_cache(
    tmp_path: Path,
) -> None:
    incomplete = _write_cache(tmp_path / "incomplete", status="running")
    _write_cache(tmp_path / "complete")

    assert not is_valid_plateau_cache(incomplete)
    assert (
        find_plateau_building_dirs(
            observer_lat_deg=36.0,
            observer_lon_deg=134.0,
            radius_km=1.0,
            root_dir=tmp_path,
        )
        == ()
    )


def test_select_prepared_building_source_falls_back_to_overture(tmp_path: Path) -> None:
    selection = select_prepared_building_source(
        observer_lat_deg=35.455,
        observer_lon_deg=133.055,
        radius_km=2.0,
        plateau_root_dir=tmp_path / "missing",
    )

    assert selection.source == "overture"
    assert selection.derived_dirs == ()


def test_select_prepared_building_source_prefers_plateau(tmp_path: Path) -> None:
    dataset_dir = _write_cache(tmp_path)

    selection = select_prepared_building_source(
        observer_lat_deg=35.455,
        observer_lon_deg=133.055,
        radius_km=2.0,
        plateau_root_dir=tmp_path,
    )

    assert selection.source == "plateau"
    assert selection.derived_dirs == (dataset_dir / "bldg",)
