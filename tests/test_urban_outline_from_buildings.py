from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


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


def test_compute_urban_outlines_emits_polylines_for_nearby_building() -> None:
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

    result = mod.compute_urban_outlines(
        tower,
        (building,),
        radius_km=5.0,
        edge_sample_step_m=10.0,
    )

    assert result.buildings_considered == 1
    assert result.outlines_emitted >= 1
    assert result.outlines[0].height_m == 60.0
    assert result.outlines[0].distance_km > 0.0
    assert len(result.outlines[0].points) >= 2


def test_compute_urban_outlines_uses_lod2_roof_surface_elevations() -> None:
    mod = _load_module()
    tower = mod.resolve_tower_viewpoint("Tokyo Skytree")
    assert tower is not None
    footprint = (
        (139.8112, 35.7102),
        (139.8114, 35.7102),
        (139.8114, 35.7104),
        (139.8112, 35.7104),
        (139.8112, 35.7102),
    )
    roof_surface = tuple(
        (lon, lat, elevation)
        for (lon, lat), elevation in zip(
            footprint,
            (100.0, 100.0, 140.0, 140.0, 100.0),
            strict=True,
        )
    )
    building = mod.BuildingFootprint(
        building_id="lod2-roof",
        height_m=60.0,
        rings_lonlat=(footprint,),
        geometry_lod=2,
        roof_surfaces_lonlat=(roof_surface,),
    )

    result = mod.compute_urban_outlines(
        tower,
        (building,),
        radius_km=5.0,
        edge_sample_step_m=10.0,
    )

    altitudes = [
        point.altitude_deg for outline in result.outlines for point in outline.points
    ]
    assert result.outlines_emitted >= 1
    assert max(altitudes) > min(altitudes)


def test_compute_urban_outlines_keeps_surfaces_when_tolerance_is_negative() -> None:
    mod = _load_module()
    tower = mod.resolve_tower_viewpoint("Tokyo Skytree")
    assert tower is not None
    outer = (
        (139.8112, 35.7102, 100.0),
        (139.8114, 35.7102, 100.0),
        (139.8114, 35.7104, 100.0),
        (139.8112, 35.7104, 100.0),
        (139.8112, 35.7102, 100.0),
    )
    inner = (
        (139.81125, 35.71025, 100.0),
        (139.81135, 35.71025, 100.0),
        (139.81135, 35.71035, 100.0),
        (139.81125, 35.71035, 100.0),
        (139.81125, 35.71025, 100.0),
    )
    footprint = tuple((lon, lat) for lon, lat, _elevation in outer)
    building = mod.BuildingFootprint(
        building_id="lod2-contained-roof",
        height_m=60.0,
        rings_lonlat=(footprint,),
        geometry_lod=2,
        roof_surfaces_lonlat=(outer, inner),
    )

    result = mod.compute_urban_outlines(
        tower,
        (building,),
        radius_km=5.0,
        edge_sample_step_m=10.0,
    )

    assert result.outlines_emitted == 2


def test_merge_projected_roof_surface_group_edges() -> None:
    mod = _load_module()
    left = np.array(
        [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [0.0, 0.0]]
    )
    right = np.array(
        [[10.0, 0.0], [20.0, 0.0], [20.0, 10.0], [10.0, 10.0], [10.0, 0.0]]
    )
    elevations = np.full(5, 100.0)

    merged = mod._merge_projected_roof_surface_groups(
        ((left, elevations), (right, elevations)),
        elevation_tolerance_m=0.0,
    )

    assert len(merged) == 1
    assert len(merged[0][0]) == 7
    assert abs(mod._ring_area_xy(merged[0][0])) == 200.0


def test_roof_surface_elevation_tolerance_increases_with_distance() -> None:
    mod = _load_module()

    assert mod._roof_surface_elevation_tolerance(0.0) == -1.0
    assert mod._roof_surface_elevation_tolerance(1000.0) == 0.0
    assert mod._roof_surface_elevation_tolerance(2000.0) == 1.0
    assert mod._roof_surface_elevation_tolerance(3000.0) == 2.0
    assert mod._roof_surface_elevation_tolerance(4000.0) == 3.0
    assert mod._roof_surface_elevation_tolerance(5000.0) == 3.0


def test_negative_roof_surface_tolerance_skips_merging() -> None:
    mod = _load_module()
    ring = np.array(
        [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [0.0, 0.0]]
    )
    elevations = np.full(5, 100.0)

    merged = mod._merge_projected_roof_surface_groups(
        ((ring, elevations), (ring, elevations)),
        elevation_tolerance_m=-0.1,
    )

    assert len(merged) == 2


def test_compute_urban_outlines_excludes_building_containing_observer() -> None:
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

    result = mod.compute_urban_outlines(
        tower,
        (building,),
        radius_km=3.0,
        edge_sample_step_m=10.0,
    )

    assert result.buildings_considered == 0
    assert result.outlines_emitted == 0


def test_compute_urban_outlines_uses_legacy_observer_height_attr() -> None:
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

    result = mod.compute_urban_outlines(
        tower,
        (building,),
        radius_km=5.0,
        edge_sample_step_m=10.0,
    )

    assert result.buildings_considered == 1
    assert result.outlines_emitted >= 1


def test_compute_urban_outlines_uses_ground_elevation_difference() -> None:
    mod = _load_module()
    tower = SimpleNamespace(
        latitude_deg=35.710055555,
        longitude_deg=139.810722222,
        viewpoint_height_m=20.0,
    )
    building = mod.BuildingFootprint(
        building_id="elevated-site",
        height_m=60.0,
        ground_elevation_m=100.0,
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

    result = mod.compute_urban_outlines(
        tower,
        (building,),
        radius_km=5.0,
        observer_ground_elevation_m=10.0,
        edge_sample_step_m=10.0,
    )

    assert result.outlines_emitted >= 1
    assert max(point.altitude_deg for point in result.outlines[0].points) > 0.0


def test_compute_urban_outlines_min_height_does_not_change_top_height() -> None:
    mod = _load_module()
    tower = SimpleNamespace(
        latitude_deg=35.710055555,
        longitude_deg=139.810722222,
        viewpoint_height_m=20.0,
    )
    low_part = mod.BuildingFootprint(
        building_id="part-low",
        height_m=60.0,
        min_height_m=0.0,
        ground_elevation_m=100.0,
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
    floating_part = mod.BuildingFootprint(
        building_id="part-floating",
        height_m=60.0,
        min_height_m=12.0,
        ground_elevation_m=100.0,
        rings_lonlat=low_part.rings_lonlat,
    )

    low_result = mod.compute_urban_outlines(
        tower,
        (low_part,),
        radius_km=5.0,
        observer_ground_elevation_m=10.0,
        edge_sample_step_m=10.0,
    )
    floating_result = mod.compute_urban_outlines(
        tower,
        (floating_part,),
        radius_km=5.0,
        observer_ground_elevation_m=10.0,
        edge_sample_step_m=10.0,
    )

    assert low_result.outlines_emitted >= 1
    assert floating_result.outlines_emitted >= 1
    assert low_result.outlines[0].points == floating_result.outlines[0].points


def test_compute_urban_outlines_skips_small_hole_rings_for_far_building() -> None:
    mod = _load_module()
    tower = SimpleNamespace(
        latitude_deg=35.0,
        longitude_deg=139.0,
        viewpoint_height_m=0.0,
    )
    building = mod.BuildingFootprint(
        building_id="far-hole",
        height_m=80.0,
        rings_lonlat=(
            (
                (139.0280, 35.0000),
                (139.0320, 35.0000),
                (139.0320, 35.0040),
                (139.0280, 35.0040),
                (139.0280, 35.0000),
            ),
            (
                (139.0299, 35.0019),
                (139.0301, 35.0019),
                (139.0301, 35.0021),
                (139.0299, 35.0021),
                (139.0299, 35.0019),
            ),
        ),
    )

    result = mod.compute_urban_outlines(
        tower,
        (building,),
        radius_km=5.0,
        edge_sample_step_m=10.0,
    )

    assert result.buildings_considered == 1
    assert result.outlines_emitted == 1


def test_compute_urban_outlines_keeps_hole_rings_for_near_building() -> None:
    mod = _load_module()
    tower = SimpleNamespace(
        latitude_deg=35.0,
        longitude_deg=139.0,
        viewpoint_height_m=0.0,
    )
    building = mod.BuildingFootprint(
        building_id="near-hole",
        height_m=80.0,
        rings_lonlat=(
            (
                (139.0030, 35.0000),
                (139.0070, 35.0000),
                (139.0070, 35.0040),
                (139.0030, 35.0040),
                (139.0030, 35.0000),
            ),
            (
                (139.0049, 35.0019),
                (139.0051, 35.0019),
                (139.0051, 35.0021),
                (139.0049, 35.0021),
                (139.0049, 35.0019),
            ),
        ),
    )

    result = mod.compute_urban_outlines(
        tower,
        (building,),
        radius_km=5.0,
        edge_sample_step_m=10.0,
    )

    assert result.buildings_considered == 1
    assert result.outlines_emitted >= 2


def test_maybe_linearize_run_points_collapses_vertical_thin_run() -> None:
    mod = _load_module()
    run_points = [
        mod.UrbanPolylinePoint(azimuth_deg=10.0, altitude_deg=1.0),
        mod.UrbanPolylinePoint(azimuth_deg=20.0, altitude_deg=1.0),
        mod.UrbanPolylinePoint(azimuth_deg=30.0, altitude_deg=1.0),
    ]

    result = mod._maybe_linearize_run_points(
        run_points,
        view_center=(1.0, 20.0),
        edge_fov_deg=90.0,
    )

    assert len(result) == 2
    assert result[0].altitude_deg == result[1].altitude_deg
    assert len({point.azimuth_deg for point in result}) == 2
    assert {point.azimuth_deg for point in result} <= {10.0, 20.0, 30.0}


def test_maybe_linearize_run_points_keeps_taller_run() -> None:
    mod = _load_module()
    run_points = [
        mod.UrbanPolylinePoint(azimuth_deg=10.0, altitude_deg=0.0),
        mod.UrbanPolylinePoint(azimuth_deg=10.0, altitude_deg=2.5),
        mod.UrbanPolylinePoint(azimuth_deg=10.0, altitude_deg=5.0),
    ]

    result = mod._maybe_linearize_run_points(
        run_points,
        view_center=(0.0, 10.0),
        edge_fov_deg=90.0,
    )

    assert result == run_points


def test_maybe_linearize_run_points_keeps_high_altitude_run() -> None:
    mod = _load_module()
    run_points = [
        mod.UrbanPolylinePoint(azimuth_deg=10.0, altitude_deg=5.1),
        mod.UrbanPolylinePoint(azimuth_deg=20.0, altitude_deg=5.2),
        mod.UrbanPolylinePoint(azimuth_deg=30.0, altitude_deg=5.3),
    ]

    result = mod._maybe_linearize_run_points(
        run_points,
        view_center=(1.0, 20.0),
        edge_fov_deg=90.0,
    )

    assert result == run_points


def test_maybe_linearize_run_points_keeps_distinct_endpoints_for_tight_vertical_run() -> None:
    mod = _load_module()
    run_points = [
        mod.UrbanPolylinePoint(azimuth_deg=22.0, altitude_deg=30.0),
        mod.UrbanPolylinePoint(azimuth_deg=22.0, altitude_deg=30.05),
        mod.UrbanPolylinePoint(azimuth_deg=22.0, altitude_deg=30.1),
    ]

    result = mod._maybe_linearize_run_points(
        run_points,
        view_center=(5.0, 0.0),
        edge_fov_deg=95.0,
    )

    assert result == run_points


def test_urban_outline_ring_width_angle_uses_perpendicular_projection() -> None:
    mod = _load_module()
    ring_xy = np.array(
        [
            [100.0, -25.0],
            [100.0, 25.0],
            [120.0, 25.0],
            [120.0, -25.0],
            [100.0, -25.0],
        ],
        dtype=np.float64,
    )

    width_deg = mod._urban_outline_ring_width_angle_deg(
        ring_xy,
        building_distance_m=1000.0,
    )

    assert width_deg > 0.0


def test_compute_urban_outlines_limits_candidate_sampling_before_expensive_projection(monkeypatch) -> None:
    mod = _load_module()
    sampled_rings: list[np.ndarray] = []

    def fake_sample_ring_points_xy(ring_xy: np.ndarray, *, sample_step_m: float) -> np.ndarray:
        sampled_rings.append(np.array(ring_xy, copy=True))
        return np.array(ring_xy, copy=True)

    monkeypatch.setattr(mod, "sample_ring_points_xy", fake_sample_ring_points_xy)
    monkeypatch.setattr(mod, "MAX_URBAN_OUTLINE_CANDIDATES", 1)

    tower = SimpleNamespace(
        latitude_deg=35.710055555,
        longitude_deg=139.810722222,
        viewpoint_height_m=10.0,
    )
    low = mod.BuildingFootprint(
        building_id="low",
        height_m=20.0,
        rings_lonlat=(
            (
                (139.8120, 35.7100),
                (139.8121, 35.7100),
                (139.8121, 35.7101),
                (139.8120, 35.7101),
                (139.8120, 35.7100),
            ),
        ),
    )
    high = mod.BuildingFootprint(
        building_id="high",
        height_m=5000.0,
        rings_lonlat=(
            (
                (139.8130, 35.7100),
                (139.8131, 35.7100),
                (139.8131, 35.7101),
                (139.8130, 35.7101),
                (139.8130, 35.7100),
            ),
        ),
    )

    result = mod.compute_urban_outlines(
        tower,
        (low, high),
        radius_km=5.0,
        edge_sample_step_m=10.0,
        max_candidates=1,
    )

    assert len(sampled_rings) == 1
    assert result.outlines_emitted == 1
