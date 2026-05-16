from __future__ import annotations

import numpy as np
from rasterio.transform import Affine

import zstarview.water_mask_interface as mod
from zstarview.water_overlay import DEFAULT_WATER_SAMPLE_GROWTH_FACTOR
from zstarview.water_mask_interface import WaterSurfaceBandStats


def test_extract_lonlat_points_from_mask_projects_all_water_pixels() -> None:
    mask = np.array(
        [
            [0, 0, 0, 0, 0],
            [0, 1, 1, 1, 0],
            [0, 1, 1, 1, 0],
            [0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0],
        ],
        dtype=np.uint8,
    )

    transform = Affine.translation(10.0, 20.0) * Affine.scale(0.5, -0.5)

    points = mod._extract_lonlat_points_from_mask(mask, transform=transform, stride=1)  # noqa: SLF001

    assert len(points) == 9
    assert (10.75, 19.25) in points


def test_collapse_tile_points_for_500m_root_picks_one_point() -> None:
    boundary_points = ((0.0, 0.0), (2.0, 1.0), (9.0, 9.0))
    got = mod._collapse_tile_points_for_root(  # noqa: SLF001
        boundary_points,
        tile_root=mod.DEFAULT_WATER_TILES_ROOT_500M,
        tile_bounds=(0.0, 0.0, 10.0, 10.0),
    )

    assert got == ((2.0, 1.0),)


def test_collapse_tile_points_for_non_sparse_root_keeps_points() -> None:
    boundary_points = ((0.0, 0.0), (2.0, 1.0))
    got = mod._collapse_tile_points_for_root(  # noqa: SLF001
        boundary_points,
        tile_root=mod.DEFAULT_WATER_TILES_ROOT_125M,
        tile_bounds=(0.0, 0.0, 10.0, 10.0),
    )

    assert got == boundary_points


def test_collapse_tile_points_for_250m_root_picks_one_point() -> None:
    boundary_points = ((0.0, 0.0), (2.0, 1.0), (9.0, 9.0))
    got = mod._collapse_tile_points_for_root(  # noqa: SLF001
        boundary_points,
        tile_root=mod.DEFAULT_WATER_TILES_ROOT_250M,
        tile_bounds=(0.0, 0.0, 10.0, 10.0),
    )

    assert got == ((2.0, 1.0),)


def test_sample_water_surface_interface_points_labels_tile_bands(monkeypatch) -> None:
    def _fake_load(*, tile_root, **_kwargs):
        if tile_root == mod.DEFAULT_WATER_TILES_ROOT_125M:
            return (
                (mod.WaterOverlayPoint("water-mask", 1.0, 2.0, 3.0, water_category="sea-125"),),
                WaterSurfaceBandStats("125m", 0, 0, 0, 1),
            )
        if tile_root == mod.DEFAULT_WATER_TILES_ROOT_250M:
            return (
                (mod.WaterOverlayPoint("water-mask", 1.0, 2.0, 3.0, water_category="sea-250"),),
                WaterSurfaceBandStats("250m", 0, 0, 0, 1),
            )
        if tile_root == mod.DEFAULT_WATER_TILES_ROOT_500M:
            return (
                (mod.WaterOverlayPoint("water-mask", 1.0, 2.0, 3.0, water_category="sea-500"),),
                WaterSurfaceBandStats("500m", 0, 0, 0, 1),
            )
        return (), WaterSurfaceBandStats("other", 0, 0, 0, 0)

    monkeypatch.setattr(
        mod,
        "_sample_water_surface_interface_ray_points_for_root_with_stats",
        _fake_load,
    )

    points, band_stats = mod.sample_water_surface_interface_points_with_stats(  # noqa: SLF001
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        observer_height_m=1.7,
        max_distance_km=10.0,
    )

    assert [point.water_category for point in points] == ["sea-125", "sea-250", "sea-500"]
    assert band_stats == (
        WaterSurfaceBandStats("125m", 0, 0, 0, 1),
        WaterSurfaceBandStats("250m", 0, 0, 0, 1),
        WaterSurfaceBandStats("500m", 0, 0, 0, 1),
    )


def test_sample_water_surface_interface_points_can_use_ground_sampler(monkeypatch) -> None:
    captured: list[float] = []

    def _fake_load(*, tile_root, target_ground_elevation_m_sampler=None, **_kwargs):
        if tile_root == mod.DEFAULT_WATER_TILES_ROOT_125M:
            if target_ground_elevation_m_sampler is not None:
                captured.append(float(target_ground_elevation_m_sampler(35.0, 139.0)))
            return (
                (mod.WaterOverlayPoint("water-mask", 1.0, 2.0, 3.0, water_category="sea-125"),),
                WaterSurfaceBandStats("125m", 0, 0, 0, 1),
            )
        return (), WaterSurfaceBandStats("other", 0, 0, 0, 0)

    monkeypatch.setattr(
        mod,
        "_sample_water_surface_interface_ray_points_for_root_with_stats",
        _fake_load,
    )

    points, _band_stats = mod.sample_water_surface_interface_points_with_stats(  # noqa: SLF001
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        observer_height_m=1.7,
        max_distance_km=10.0,
        target_ground_elevation_m_sampler=lambda _lat, _lon: 42.0,
    )

    assert points[0].water_category == "sea-125"
    assert captured == [42.0]


def test_sample_water_surface_interface_ray_points_start_at_125m(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def _fake_build_geometric_distance_samples(
        max_distance_km,
        sample_start_m,
        *,
        growth_factor=DEFAULT_WATER_SAMPLE_GROWTH_FACTOR,
    ):
        captured["args"] = (float(max_distance_km), float(sample_start_m), float(growth_factor))
        return np.asarray([37.5, 114.0, 349.0], dtype=np.float64)

    def _fake_scan(*, distance_samples_m, **_kwargs):
        captured["distances"] = tuple(float(value) for value in distance_samples_m.tolist())
        return type(
            "Scan",
            (),
            {
                "azimuths_deg": [0.0],
                "distance_grid_m": np.asarray([[349.0]], dtype=np.float64),
                "ray_lon_deg": np.asarray([[139.001]], dtype=np.float64),
                "ray_lat_deg": np.asarray([[35.001]], dtype=np.float64),
            },
        )()

    monkeypatch.setattr(mod, "build_geometric_distance_samples", _fake_build_geometric_distance_samples)
    monkeypatch.setattr(mod, "build_ray_scan_grid", _fake_scan)
    monkeypatch.setattr(mod, "_sample_water_mask_for_lonlat_points", lambda lonlat_points, **_kwargs: [True] * len(lonlat_points))
    monkeypatch.setattr(
        mod,
        "project_place_targets_to_altaz",
        lambda **_kwargs: (
            type("Projection", (), {"alt_deg": 1.0, "az_deg": 2.0, "distance_km": 3.0})(),
        ),
    )

    points, _stats = mod._sample_water_surface_interface_ray_points_for_root_with_stats(  # noqa: SLF001
        center_lat_deg=35.0,
        center_lon_deg=139.0,
        observer_height_m=1.7,
        radius_km=10.0,
        tile_root=mod.DEFAULT_WATER_TILES_ROOT_125M,
    )

    assert captured["args"] == (
        10.0,
        float(mod.DEFAULT_WATER_SAMPLE_STEP_M),
        float(DEFAULT_WATER_SAMPLE_GROWTH_FACTOR),
    )
    assert captured["distances"] == (349.0,)
    assert [point.scan_distance_index for point in points] == [0]


def test_sample_water_surface_horizon_points_uses_horizon_altaz(monkeypatch) -> None:
    monkeypatch.setattr(
        mod,
        "_sample_water_mask_for_lonlat_points",
        lambda lonlat_points, **_kwargs: [True, False],
    )

    points = mod.sample_water_surface_horizon_points(  # noqa: SLF001
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        horizon_profile_altaz=[(1.5, 10.0), (2.0, 20.0)],
        horizon_profile_distances_m=[500.0, 1000.0],
    )

    assert len(points) == 1
    assert points[0].alt_deg == 1.5
    assert points[0].az_deg == 10.0
    assert points[0].distance_km == 0.5
    assert points[0].water_category == "sea"


def test_sample_water_surface_horizon_points_uses_near_mid_and_far_tile_roots(monkeypatch) -> None:
    seen: list[object] = []

    def _fake_sample(lonlat_points, *, tile_root=None):
        seen.append(tile_root)
        return [True] * len(lonlat_points)

    monkeypatch.setattr(mod, "_sample_water_mask_for_lonlat_points", _fake_sample)

    points = mod.sample_water_surface_horizon_points(  # noqa: SLF001
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        horizon_profile_altaz=[(1.5, 10.0), (2.0, 20.0), (3.0, 30.0)],
        horizon_profile_distances_m=[500.0, 4000.0, 9000.0],
    )

    assert len(points) == 3
    assert seen == [
        mod.DEFAULT_WATER_TILES_ROOT_125M,
        mod.DEFAULT_WATER_TILES_ROOT_250M,
        mod.DEFAULT_WATER_TILES_ROOT_500M,
    ]


def test_sample_water_surface_horizon_layers_points_combines_secondary_layers(monkeypatch) -> None:
    def _fake_sample(lonlat_points, *, tile_root=None):
        return [True] * len(lonlat_points)

    monkeypatch.setattr(mod, "_sample_water_mask_for_lonlat_points", _fake_sample)

    points = mod.sample_water_surface_horizon_layers_points(  # noqa: SLF001
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        horizon_profile_altaz=[(1.5, 10.0)],
        horizon_profile_distances_m=[500.0],
        secondary_horizon_profile_altaz_layers=[[(1.0, 20.0)]],
        secondary_horizon_profile_distances_m_layers=[[1000.0]],
    )

    assert len(points) == 2
    assert {point.az_deg for point in points} == {10.0, 20.0}


def test_sample_water_mask_for_lonlat_points_uses_marker_tiles(tmp_path) -> None:
    water_tile = tmp_path / "tile_y1_x4.1"
    land_tile = tmp_path / "tile_y1_x5.0"
    water_tile.write_bytes(b"")
    land_tile.write_bytes(b"")

    flags = mod._sample_water_mask_for_lonlat_points(  # noqa: SLF001
        [(10.0, 10.0), (60.0, 10.0)],
        tile_root=tmp_path,
    )

    assert flags == [True, False]


def test_load_water_surface_interface_lonlat_points_skips_marker_tiles(
    tmp_path, monkeypatch
) -> None:
    water_tile = tmp_path / "tile_y1_x4.1"
    water_tile.write_bytes(b"")

    import rasterio

    monkeypatch.setattr(rasterio, "open", lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("should not open marker tiles")))

    points = mod.load_water_surface_interface_lonlat_points(  # noqa: SLF001
        center_lat_deg=10.0,
        center_lon_deg=10.0,
        radius_km=20.0,
        tile_root=tmp_path,
    )

    assert points == ()


def test_load_water_surface_interface_lonlat_points_combines_near_and_far_roots(
    monkeypatch,
) -> None:
    seen: list[tuple[float, object]] = []

    def _fake_load(*, center_lat_deg, center_lon_deg, radius_km, tile_root, bbox_scale, stride, **_kwargs):
        seen.append((float(radius_km), tile_root))
        if tile_root == mod.DEFAULT_WATER_TILES_ROOT_125M:
            return ((1.0, 1.0),)
        if tile_root == mod.DEFAULT_WATER_TILES_ROOT_250M:
            return ((1.5, 1.5),)
        if tile_root == mod.DEFAULT_WATER_TILES_ROOT_500M:
            return ((2.0, 2.0),)
        return ()

    monkeypatch.setattr(
        mod,
        "_load_water_surface_interface_lonlat_points_for_root",
        _fake_load,
    )

    points = mod.load_water_surface_interface_lonlat_points(  # noqa: SLF001
        center_lat_deg=10.0,
        center_lon_deg=10.0,
        radius_km=20.0,
    )

    assert points == ((1.0, 1.0), (1.5, 1.5), (2.0, 2.0))
    assert seen[0][0] == 2.0
    assert seen[0][1] == mod.DEFAULT_WATER_TILES_ROOT_125M
    assert seen[1][0] == 6.0
    assert seen[1][1] == mod.DEFAULT_WATER_TILES_ROOT_250M
    assert seen[2][0] == 12.0
    assert seen[2][1] == mod.DEFAULT_WATER_TILES_ROOT_500M
    assert seen[3][0] == 18.0
    assert seen[3][1] == mod.DEFAULT_WATER_TILES_ROOT_500M
    assert seen[4][0] == 20.0
    assert seen[4][1] == mod.DEFAULT_WATER_TILES_ROOT_500M
