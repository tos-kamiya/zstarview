from __future__ import annotations

import numpy as np
from rasterio.transform import Affine

from zstarview.terrain.dem import WGS84_GEOD
from zstarview.terrain.horizon import (
    DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM,
    ObserverLocation,
    compute_horizon_layers,
    _prune_secondary_peak_indices_by_visibility,
    _prune_secondary_peak_indices_by_main_profile,
    _select_distance_band_peak_index,
    _should_break_secondary_ridge,
    build_distance_samples,
    compute_apparent_altitudes,
)
from zstarview.render.terrain import _distance_band_alpha
from zstarview.terrain.dem import DemGrid, sample_ground_elevation


class _IdentityTransformer:
    def transform(self, lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return lon_deg, lat_deg


class _FlatDemGrid:
    def sample_lonlat(
        self,
        lon_deg: np.ndarray,
        lat_deg: np.ndarray,
        *,
        method: str = "bilinear",
    ) -> np.ndarray:
        return np.zeros_like(lon_deg, dtype=np.float64)


def test_build_distance_samples_covers_requested_distance() -> None:
    samples = build_distance_samples(1.0, 300.0)

    assert len(samples) == 4
    assert samples[0] == 300.0
    assert samples[-1] == 1000.0


def test_compute_apparent_altitudes_is_zero_for_equal_height_flat_short_range() -> None:
    got = compute_apparent_altitudes(
        observer_elevation_m=100.0,
        target_elevation_m=np.array([100.0]),
        surface_distance_m=np.array([1.0]),
        earth_radius_m=6_371_008.8,
        refraction_coefficient=0.0,
    )

    assert abs(float(got[0])) < 1e-3


def test_should_break_secondary_ridge_on_large_distance_jump() -> None:
    assert _should_break_secondary_ridge(
        3_000.0,
        9_500.0,
        max_distance_jump_ratio=0.25,
    )


def test_should_keep_secondary_ridge_on_small_distance_change() -> None:
    assert not _should_break_secondary_ridge(
        3_000.0,
        3_700.0,
        max_distance_jump_ratio=0.25,
    )


def test_prune_secondary_peak_indices_by_visibility_drops_farther_lower_peak() -> None:
    altitudes = np.array([0.2, 0.7, 0.6, 0.9], dtype=np.float64)

    got = _prune_secondary_peak_indices_by_visibility(altitudes, [0, 1, 2, 3])

    assert got == [0, 1, 3]


def test_prune_secondary_peak_indices_by_visibility_keeps_farther_higher_peak() -> None:
    altitudes = np.array([0.2, 0.1, 0.4, 0.3], dtype=np.float64)

    got = _prune_secondary_peak_indices_by_visibility(altitudes, [0, 1, 2, 3])

    assert got == [0, 2]


def test_prune_secondary_peak_indices_by_main_profile_drops_lower_than_main() -> None:
    altitudes = np.array([0.2, 0.7, 0.6, 0.9], dtype=np.float64)
    distances = np.array([100.0, 200.0, 300.0, 400.0], dtype=np.float64)

    got = _prune_secondary_peak_indices_by_main_profile(
        altitudes,
        distances,
        [0, 1, 2, 3],
        main_peak_distance_m=250.0,
    )

    assert got == [0, 1]


def test_prune_secondary_peak_indices_by_main_profile_keeps_higher_than_main() -> None:
    altitudes = np.array([0.2, 0.7, 0.6, 0.9], dtype=np.float64)
    distances = np.array([100.0, 200.0, 300.0, 400.0], dtype=np.float64)

    got = _prune_secondary_peak_indices_by_main_profile(
        altitudes,
        distances,
        [0, 1, 2, 3],
        main_peak_distance_m=450.0,
    )

    assert got == [0, 1, 2, 3]


def test_select_distance_band_peak_index_prefers_highest_in_band() -> None:
    altitudes = np.array([0.1, 0.6, 0.3, 0.9, 0.4], dtype=np.float64)
    distances = np.array([5_000.0, 20_000.0, 35_000.0, 42_000.0, 70_000.0], dtype=np.float64)

    got = _select_distance_band_peak_index(
        altitudes,
        distances,
        band_min_m=10_000.0,
        band_max_m=50_000.0,
        include_upper_bound=False,
    )

    assert got == 3


def test_select_distance_band_peak_index_includes_last_upper_bound() -> None:
    altitudes = np.array([0.1, 0.6, 0.3, 0.9, 0.4], dtype=np.float64)
    distances = np.array([5_000.0, 20_000.0, 35_000.0, 42_000.0, 120_000.0], dtype=np.float64)

    got = _select_distance_band_peak_index(
        altitudes,
        distances,
        band_min_m=80_000.0,
        band_max_m=120_000.0,
        include_upper_bound=True,
    )

    assert got == 4


def test_distance_band_alpha_drops_with_distance() -> None:
    near = _distance_band_alpha(0, 7, 0.38)
    far = _distance_band_alpha(6, 7, 0.38)

    assert near > far
    assert far < 0.1


def test_compute_horizon_layers_adds_nearest_secondary_band() -> None:
    layers = compute_horizon_layers(
        dem_grid=_FlatDemGrid(),
        geod=WGS84_GEOD,
        observer=ObserverLocation(
            latitude_deg=35.0,
            longitude_deg=139.0,
            observer_ground_m=0.0,
            observer_eye_m=0.0,
        ),
        azimuth_step_deg=180.0,
        distance_samples_m=build_distance_samples(128.0, 500.0),
        dem_resampling="bilinear",
        earth_radius_m=6_371_008.8,
        refraction_coefficient=0.0,
    )

    assert DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM[0] == 2.0
    assert len(layers.secondary_layers) == len(DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM)
    assert [point.distance_m for point in layers.secondary_layers[0]] == [500.0, 500.0]


def test_sample_ground_elevation_uses_default_elevation_for_nodata() -> None:
    grid = DemGrid(
        array=np.array([[np.nan]], dtype=np.float64),
        transform=Affine.translation(0.0, 1.0) * Affine.scale(1.0, -1.0),
        crs="EPSG:4326",
        default_elevation_m=0.0,
        _to_grid=_IdentityTransformer(),
    )

    got = sample_ground_elevation(
        grid,
        latitude_deg=0.5,
        longitude_deg=0.5,
    )

    assert got == 0.0
