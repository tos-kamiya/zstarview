from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "urban_outline_from_buildings.py"
    spec = importlib.util.spec_from_file_location("urban_outline_from_buildings", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_compute_debug_outlines_emits_polylines_for_nearby_building() -> None:
    mod = _load_module()
    tower = mod.resolve_tower_viewpoint("Tokyo Skytree")
    assert tower is not None
    building = mod.BuildingFootprint(
        building_id="test",
        height_m=60.0,
        rings_lonlat=(
            (
                (139.8112, 35.7102),
                (139.8114, 35.7102),
                (139.8114, 35.7104),
                (139.8112, 35.7104),
                (139.8112, 35.7102),
            ),
        ),
    )

    result = mod.compute_debug_outlines(
        tower,
        (building,),
        radius_km=5.0,
        edge_sample_step_m=10.0,
    )

    assert result.buildings_considered == 1
    assert result.outlines_emitted >= 1
    assert len(result.outlines[0]) >= 2


def test_compute_debug_outlines_excludes_building_containing_observer() -> None:
    mod = _load_module()
    tower = mod.resolve_tower_viewpoint("Tokyo Tower")
    assert tower is not None
    building = mod.BuildingFootprint(
        building_id="observer-building",
        height_m=350.0,
        rings_lonlat=(
            (
                (tower.longitude_deg - 0.0002, tower.latitude_deg - 0.0002),
                (tower.longitude_deg + 0.0002, tower.latitude_deg - 0.0002),
                (tower.longitude_deg + 0.0002, tower.latitude_deg + 0.0002),
                (tower.longitude_deg - 0.0002, tower.latitude_deg + 0.0002),
                (tower.longitude_deg - 0.0002, tower.latitude_deg - 0.0002),
            ),
        ),
    )

    result = mod.compute_debug_outlines(
        tower,
        (building,),
        radius_km=3.0,
        edge_sample_step_m=10.0,
    )

    assert result.buildings_considered == 0
    assert result.outlines_emitted == 0


def test_compute_debug_outlines_uses_legacy_observer_height_attr() -> None:
    mod = _load_module()
    tower = SimpleNamespace(
        latitude_deg=35.710055555,
        longitude_deg=139.810722222,
        viewpoint_height_m=None,
        observer_height_m=500.0,
    )
    building = mod.BuildingFootprint(
        building_id="legacy-height",
        height_m=60.0,
        rings_lonlat=(
            (
                (139.8112, 35.7102),
                (139.8114, 35.7102),
                (139.8114, 35.7104),
                (139.8112, 35.7104),
                (139.8112, 35.7102),
            ),
        ),
    )

    result = mod.compute_debug_outlines(
        tower,
        (building,),
        radius_km=5.0,
        edge_sample_step_m=10.0,
    )

    assert result.buildings_considered == 1
    assert result.outlines_emitted >= 1
