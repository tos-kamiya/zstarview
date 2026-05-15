from __future__ import annotations

import numpy as np
from rasterio.transform import Affine

import zstarview.water_mask_interface as mod


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


def test_sample_water_surface_interface_points_wraps_as_sea(monkeypatch) -> None:
    monkeypatch.setattr(
        mod,
        "load_water_surface_interface_lonlat_points",
        lambda **_kwargs: ((139.0, 35.0),),
    )

    points = mod.sample_water_surface_interface_points(  # noqa: SLF001
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        observer_height_m=1.7,
        max_distance_km=10.0,
    )

    assert len(points) == 1
    assert points[0].water_category == "sea"


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


def test_sample_water_surface_horizon_points_uses_near_and_far_tile_roots(monkeypatch) -> None:
    seen: list[tuple[tuple[tuple[float, float], ...], object]] = []

    def _fake_sample(lonlat_points, *, tile_root=None):
        seen.append((tuple(lonlat_points), tile_root))
        return [True] * len(lonlat_points)

    monkeypatch.setattr(mod, "_sample_water_mask_for_lonlat_points", _fake_sample)

    points = mod.sample_water_surface_horizon_points(  # noqa: SLF001
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        horizon_profile_altaz=[(1.5, 10.0), (2.0, 20.0)],
        horizon_profile_distances_m=[500.0, 3000.0],
    )

    assert len(points) == 2
    assert len(seen) == 2
    assert seen[0][1] == mod.DEFAULT_WATER_TILES_ROOT_100M
    assert seen[1][1] == mod.DEFAULT_WATER_TILES_ROOT_500M


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

    def _fake_load(*, center_lat_deg, center_lon_deg, radius_km, tile_root, bbox_scale, stride):
        seen.append((float(radius_km), tile_root))
        if tile_root == mod.DEFAULT_WATER_TILES_ROOT_100M:
            return ((1.0, 1.0),)
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

    assert points == ((1.0, 1.0), (2.0, 2.0))
    assert seen[0][0] == 2.0
    assert seen[0][1] == mod.DEFAULT_WATER_TILES_ROOT_100M
    assert seen[1][0] == 20.0
    assert seen[1][1] == mod.DEFAULT_WATER_TILES_ROOT_500M
