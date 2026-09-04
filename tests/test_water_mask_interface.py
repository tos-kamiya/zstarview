from __future__ import annotations

import threading

import numpy as np
import pytest
from rasterio.transform import Affine

import zstarview.water_mask_interface as mod
from zstarview.clouddisc.types import DownloadCancelledError
from zstarview.water_mask_interface import WaterSurfaceBandStats
from zstarview.water_overlay import DEFAULT_WATER_SAMPLE_GROWTH_FACTOR


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

    points = mod._extract_lonlat_points_from_mask(mask, transform=transform, stride=1)

    assert len(points) == 9
    assert (10.75, 19.25) in points


def test_collapse_tile_points_for_500m_root_picks_one_point() -> None:
    boundary_points = ((0.0, 0.0), (2.0, 1.0), (9.0, 9.0))
    got = mod._collapse_tile_points_for_root(
        boundary_points,
        tile_root=mod.DEFAULT_WATER_TILES_ROOT_500M,
        tile_bounds=(0.0, 0.0, 10.0, 10.0),
    )

    assert got == ((2.0, 1.0),)


def test_collapse_tile_points_for_non_sparse_root_keeps_points() -> None:
    boundary_points = ((0.0, 0.0), (2.0, 1.0))
    got = mod._collapse_tile_points_for_root(
        boundary_points,
        tile_root=mod.DEFAULT_WATER_TILES_ROOT_125M,
        tile_bounds=(0.0, 0.0, 10.0, 10.0),
    )

    assert got == boundary_points


def test_collapse_tile_points_for_250m_root_picks_one_point() -> None:
    boundary_points = ((0.0, 0.0), (2.0, 1.0), (9.0, 9.0))
    got = mod._collapse_tile_points_for_root(
        boundary_points,
        tile_root=mod.DEFAULT_WATER_TILES_ROOT_250M,
        tile_bounds=(0.0, 0.0, 10.0, 10.0),
    )

    assert got == ((2.0, 1.0),)


def test_tile_key_for_lonlat_uses_resolution_specific_grids() -> None:
    assert mod._tile_key_for_lonlat(135.0, 0.0, tile_root=mod.DEFAULT_WATER_TILES_ROOT_125M) == (8, 28)
    assert mod._tile_key_for_lonlat(135.0, 0.0, tile_root=mod.DEFAULT_WATER_TILES_ROOT_250M) == (4, 14)
    assert mod._tile_key_for_lonlat(135.0, 0.0, tile_root=mod.DEFAULT_WATER_TILES_ROOT_500M) == (2, 7)


def test_water_bands_use_25m_near_cache_when_available(monkeypatch) -> None:
    fake_25m_root = object()
    monkeypatch.setattr(mod._ZipWaterMaskRoot, "from_cache", staticmethod(lambda: fake_25m_root))

    specs = mod._water_band_specs(tile_root=None, max_distance_km=10.0)

    assert specs[0][0] is fake_25m_root
    assert specs[0][1:3] == (0.0, 0.25)
    assert specs[1][0] == mod.DEFAULT_WATER_TILES_ROOT_125M
    assert specs[1][1:3] == (0.25, 2.0)
    assert specs[2][0] == mod.DEFAULT_WATER_TILES_ROOT_250M


def test_zip_water_mask_accepts_release_archive_directory_prefix(tmp_path, monkeypatch) -> None:
    import hashlib
    import json
    import zipfile

    cache_root = tmp_path / "resolution-25m"
    cache_root.mkdir()
    archive_path = cache_root / mod.WATER_MASK_ASSET_NAME
    with zipfile.ZipFile(archive_path, "w") as archive:
        for row in range(16):
            for column in range(32):
                archive.writestr(
                    f"raw-data/water-tiles-25m-global-20260725/tile_y{row:02d}_x{column:02d}.0",
                    b"",
                )
    manifest = {
        "schema": 1,
        "data_source_date": "2026-07-25",
        "coverage": {"columns": 32, "latitude_rows": 16},
        "raster": {"resolution_m": 25},
        "assets": [
            {
                "name": mod.WATER_MASK_ASSET_NAME,
                "bytes": archive_path.stat().st_size,
                "sha256": hashlib.sha256(archive_path.read_bytes()).hexdigest(),
            }
        ],
    }
    (cache_root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (cache_root / "READY").write_text("ready\n", encoding="ascii")
    monkeypatch.setattr(mod, "water_mask_dataset_root", lambda: cache_root)

    root = mod._ZipWaterMaskRoot.from_cache()
    assert root is not None
    assert len(root.members) == 512
    root.close()


def test_sample_water_surface_interface_points_labels_tile_bands(monkeypatch) -> None:
    monkeypatch.setattr(mod._ZipWaterMaskRoot, "from_cache", staticmethod(lambda: None))

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

    points, band_stats = mod.sample_water_surface_interface_points_with_stats(
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


def test_sample_water_surface_interface_points_keeps_sea_mask_at_zero_m(monkeypatch) -> None:
    projected_heights: list[list[float]] = []

    def _fake_build_ray_scan_grid(**_kwargs):
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

    def _fake_sample_water_mask_for_lonlat_points_with_stats(lonlat_points, **_kwargs):
        return [True] * len(lonlat_points), 1

    def _fake_project_place_targets_to_altaz(**kwargs):
        projected_heights.append([float(value) for value in kwargs["target_height_m"]])
        return (
            type("Projection", (), {"alt_deg": 1.0, "az_deg": 2.0, "distance_km": 3.0})(),
        )

    monkeypatch.setattr(
        mod,
        "build_ray_scan_grid",
        _fake_build_ray_scan_grid,
    )
    monkeypatch.setattr(
        mod,
        "_sample_water_mask_for_lonlat_points_with_stats",
        _fake_sample_water_mask_for_lonlat_points_with_stats,
    )
    monkeypatch.setattr(mod, "project_place_targets_to_altaz", _fake_project_place_targets_to_altaz)

    points, _band_stats = mod._sample_water_surface_interface_ray_points_for_root_with_stats(
        center_lat_deg=35.0,
        center_lon_deg=139.0,
        observer_height_m=1.7,
        radius_km=10.0,
        tile_root=mod.DEFAULT_WATER_TILES_ROOT_125M,
    )

    assert points[0].water_category == "sea-125"
    assert projected_heights == [[0.0]]


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
    monkeypatch.setattr(
        mod,
        "_sample_water_mask_for_lonlat_points_with_stats",
        lambda lonlat_points, **_kwargs: ([True] * len(lonlat_points), 1),
    )
    monkeypatch.setattr(
        mod,
        "project_place_targets_to_altaz",
        lambda **_kwargs: (
            type("Projection", (), {"alt_deg": 1.0, "az_deg": 2.0, "distance_km": 3.0})(),
        ),
    )

    points, _stats = mod._sample_water_surface_interface_ray_points_for_root_with_stats(
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


def test_sample_water_surface_interface_ray_points_can_be_cancelled(monkeypatch) -> None:
    abort_event = threading.Event()
    abort_event.set()

    monkeypatch.setattr(
        mod,
        "build_ray_scan_grid",
        lambda **_kwargs: type(
            "Scan",
            (),
            {
                "azimuths_deg": [0.0],
                "distance_grid_m": np.asarray([[349.0]], dtype=np.float64),
                "ray_lon_deg": np.asarray([[139.001]], dtype=np.float64),
                "ray_lat_deg": np.asarray([[35.001]], dtype=np.float64),
            },
        )(),
    )

    with pytest.raises(DownloadCancelledError):
        mod._sample_water_surface_interface_ray_points_for_root_with_stats(
            center_lat_deg=35.0,
            center_lon_deg=139.0,
            observer_height_m=1.7,
            radius_km=10.0,
            tile_root=mod.DEFAULT_WATER_TILES_ROOT_125M,
            abort_event=abort_event,
        )


def test_sample_water_surface_horizon_points_uses_horizon_altaz(monkeypatch) -> None:
    monkeypatch.setattr(
        mod,
        "_sample_water_mask_for_lonlat_points",
        lambda lonlat_points, **_kwargs: [True, False],
    )

    points = mod.sample_water_surface_horizon_points(
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

    points = mod.sample_water_surface_horizon_points(
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

    points = mod.sample_water_surface_horizon_layers_points(
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

    flags = mod._sample_water_mask_for_lonlat_points(
        [(10.0, 10.0), (60.0, 10.0)],
        tile_root=tmp_path,
    )

    assert flags == [True, False]


def test_sample_water_mask_for_lonlat_points_uses_resolution_specific_keys(
    tmp_path, monkeypatch
) -> None:
    root_125 = tmp_path / "125"
    root_250 = tmp_path / "250"
    root_500 = tmp_path / "500"
    root_125.mkdir()
    root_250.mkdir()
    root_500.mkdir()

    (root_125 / "tile_y8_x28.1").write_bytes(b"")
    (root_250 / "tile_y4_x14.1").write_bytes(b"")
    (root_500 / "tile_y2_x7.1").write_bytes(b"")

    monkeypatch.setattr(mod, "DEFAULT_WATER_TILES_ROOT_125M", root_125)
    monkeypatch.setattr(mod, "DEFAULT_WATER_TILES_ROOT_250M", root_250)
    monkeypatch.setattr(mod, "DEFAULT_WATER_TILES_ROOT_500M", root_500)

    assert mod._sample_water_mask_for_lonlat_points(
        [(135.0, 0.0)],
        tile_root=root_125,
    ) == [True]
    assert mod._sample_water_mask_for_lonlat_points(
        [(135.0, 0.0)],
        tile_root=root_250,
    ) == [True]
    assert mod._sample_water_mask_for_lonlat_points(
        [(135.0, 0.0)],
        tile_root=root_500,
    ) == [True]


def test_sample_water_mask_for_lonlat_points_with_stats_counts_open_tiles(
    tmp_path, monkeypatch
) -> None:
    root = tmp_path / "water"
    root.mkdir()
    (root / "tile_y2_x7.tif").write_bytes(b"")
    (root / "tile_y2_x6.0").write_bytes(b"")

    class _FakeDataset:
        bounds = type("Bounds", (), {"left": 135.0, "right": 180.0, "bottom": -45.0, "top": 0.0})()

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def sample(self, coords):
            return [np.asarray([1], dtype=np.uint8) for _coord in coords]

    import rasterio

    monkeypatch.setattr(rasterio, "open", lambda *_args, **_kwargs: _FakeDataset())

    flags, opened_tile_count = mod._sample_water_mask_for_lonlat_points_with_stats(
        [(135.0, 0.0), (134.0, 0.0)],
        tile_root=root,
    )

    assert flags == [True, False]
    assert opened_tile_count == 1


def test_load_water_surface_interface_lonlat_points_skips_marker_tiles(
    tmp_path, monkeypatch
) -> None:
    water_tile = tmp_path / "tile_y1_x4.1"
    water_tile.write_bytes(b"")

    import rasterio

    monkeypatch.setattr(rasterio, "open", lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("should not open marker tiles")))

    points = mod.load_water_surface_interface_lonlat_points(
        center_lat_deg=10.0,
        center_lon_deg=10.0,
        radius_km=20.0,
        tile_root=tmp_path,
    )

    assert points == ()


def test_load_water_surface_interface_lonlat_points_combines_near_and_far_roots(
    monkeypatch,
) -> None:
    monkeypatch.setattr(mod._ZipWaterMaskRoot, "from_cache", staticmethod(lambda: None))

    seen: list[tuple[float, object]] = []

    def _fake_load(*, center_lat_deg, center_lon_deg, radius_km, tile_root, stride, **_kwargs):
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

    points = mod.load_water_surface_interface_lonlat_points(
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
