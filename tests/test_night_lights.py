from __future__ import annotations

import numpy as np

from zstarview import night_lights


def test_terrain_visibility_threshold_curve_uses_distance_order() -> None:
    distances_m = np.asarray([500.0, 1000.0, 1500.0, 2500.0], dtype=np.float64)
    threshold = night_lights._terrain_visibility_threshold_curve(
        azimuth_deg=180.0,
        distances_m=distances_m,
        terrain_profile_altaz=[(1.0, 180.0), (2.0, 180.0)],
        terrain_profile_distances_m=[1000.0, 2000.0],
        terrain_secondary_ridges_altaz_layers=[
            [(3.0, 180.0)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1500.0],
        ],
    )

    assert threshold is not None
    assert np.allclose(
        threshold,
        np.asarray([-np.inf, 1.0, 3.0, 3.0], dtype=np.float64),
        equal_nan=False,
    )


def test_surface_point_apparent_altitudes_decrease_with_distance() -> None:
    distances_m = np.asarray([1000.0, 5000.0, 10000.0], dtype=np.float64)

    altitudes = night_lights._surface_point_apparent_altitudes(
        distances_m,
        observer_height_m=120.0,
        refraction_coefficient=0.13,
    )

    assert altitudes.shape == distances_m.shape
    assert np.all(np.isfinite(altitudes))
    assert altitudes[0] < altitudes[1] < altitudes[2] < 0.0


def test_apply_night_light_sample_floor_keeps_masked_samples_dark() -> None:
    samples = np.asarray([0.0, 1.0, 2.0, 0.0], dtype=np.float64)
    visibility_mask = np.asarray([True, True, False, True], dtype=bool)

    got = night_lights._apply_night_light_sample_floor(
        samples,
        visibility_mask,
        floor_value=0.2,
    )

    assert np.allclose(got, np.asarray([0.2, 1.2, 0.0, 0.2], dtype=np.float64))


def test_terrain_sample_edge_strength_rows_use_dem_height() -> None:
    terrain_sample_distances_m = np.asarray([1_000.0, 3_000.0], dtype=np.float64)
    terrain_sample_terrain_elevation_m = np.asarray(
        [
            [100.0, 250.0],
            [50.0, 75.0],
        ],
        dtype=np.float64,
    )
    source_distances_m = np.asarray([1_000.0, 2_000.0, 3_000.0], dtype=np.float64)

    rows = night_lights._terrain_sample_edge_strength_rows(
        terrain_sample_distances_m=terrain_sample_distances_m,
        terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
        source_distances_m=source_distances_m,
    )

    assert rows is not None
    assert rows.shape == (2, 3)
    assert np.all(rows >= 0.0)
    assert rows[0, 2] > rows[0, 0]
    assert rows[1, 2] > rows[1, 0]


def test_terrain_band_target_mask_uses_altitude_fade(monkeypatch) -> None:
    def fake_visibility_threshold_curve(**_kwargs) -> np.ndarray:
        return np.asarray([2.0], dtype=np.float64)

    monkeypatch.setattr(night_lights, "_terrain_visibility_threshold_curve", fake_visibility_threshold_curve)

    mask = night_lights._terrain_band_target_mask(
        az_grid=np.asarray([180.0, 180.0, 180.0], dtype=np.float64),
        target_altitudes=np.asarray([0.0, 1.0, 2.0], dtype=np.float64),
        band_distance_m=1_000.0,
        terrain_profile_altaz=[(2.0, 180.0)],
        terrain_profile_distances_m=[1_000.0],
        terrain_secondary_ridges_altaz_layers=[[(2.0, 180.0)]],
        terrain_secondary_ridges_distances_m_layers=[[1_000.0]],
    )

    assert mask.dtype == np.float64
    assert np.allclose(mask, np.asarray([0.0, 0.5, 1.0], dtype=np.float64))


def test_terrain_band_target_altaz_mask_uses_altitude_fade(monkeypatch) -> None:
    def fake_visibility_threshold_curve(**_kwargs) -> np.ndarray:
        return np.asarray([2.0], dtype=np.float64)

    monkeypatch.setattr(night_lights, "_terrain_visibility_threshold_curve", fake_visibility_threshold_curve)

    mask = night_lights._terrain_band_target_altaz_mask(
        az_grid=np.asarray([180.0, 190.0], dtype=np.float64),
        target_altitudes=np.asarray([0.0, 1.0, 2.0], dtype=np.float64),
        band_distance_m=1_000.0,
        terrain_profile_altaz=[(2.0, 180.0)],
        terrain_profile_distances_m=[1_000.0],
        terrain_secondary_ridges_altaz_layers=[[(2.0, 180.0)]],
        terrain_secondary_ridges_distances_m_layers=[[1_000.0]],
    )

    assert mask.dtype == np.float64
    assert mask.shape == (3, 2)
    assert np.allclose(mask[:, 0], np.asarray([0.0, 0.5, 1.0], dtype=np.float64))
    assert np.allclose(mask[:, 1], np.asarray([0.0, 0.5, 1.0], dtype=np.float64))


def test_night_light_distance_boost_grows_linearly() -> None:
    distances_m = np.asarray([0.0, 64_000.0, 128_000.0], dtype=np.float64)

    boost = night_lights._night_light_distance_boost(distances_m)

    assert np.allclose(boost, np.asarray([1.0, 1.5, 2.0], dtype=np.float64))


def test_night_light_strength_factor_uses_minus_nine_to_minus_four_blend() -> None:
    assert np.isclose(night_lights.night_light_strength_factor(-9.0), 1.0)
    assert np.isclose(night_lights.night_light_strength_factor(-6.5), 0.5)
    assert np.isclose(night_lights.night_light_strength_factor(-4.0), 0.0)
    assert night_lights.is_night_light_enabled(-4.1)
    assert not night_lights.is_night_light_enabled(-4.0)


def test_gaussian_weight_lut_uses_half_degree_bins() -> None:
    sigma = night_lights.NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG
    weights = night_lights._lookup_gaussian_weights(
        np.asarray([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float64),
        sigma_deg=sigma,
        step_deg=0.5,
    )

    assert weights.shape == (5,)
    assert np.isclose(weights[0], 1.0)
    assert np.isclose(weights[1], weights[0])
    assert np.isclose(weights[2], night_lights._gaussian_weight_lut(sigma, 0.5)[1])
    assert weights[2] > weights[3] >= weights[4]


def test_gaussian_weights_clip_beyond_azimuth_limit() -> None:
    sigma = night_lights.NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG
    weights = night_lights._lookup_gaussian_weights(
        np.asarray([24.5, 25.0, 25.5, 30.0], dtype=np.float64),
        sigma_deg=sigma,
        step_deg=0.5,
        max_delta_deg=night_lights.NIGHT_LIGHTS_NEIGHBORHOOD_MAX_AZ_DELTA_DEG,
    )

    assert weights[0] > 0.0
    assert weights[1] > 0.0
    assert weights[2] == 0.0
    assert weights[3] == 0.0


def test_target_altitude_bins_use_two_degree_step() -> None:
    bins = night_lights._target_altitude_bins()

    assert bins.size > 1
    assert np.isclose(bins[1] - bins[0], 2.0)
    assert np.isclose(bins[0], -90.0)


def test_sample_ray_brightness_curve_uses_linear_distance_boost(
    monkeypatch,
    tmp_path,
) -> None:
    tile_path = tmp_path / "BlackMarble_2016_C1_geo_gray.tif"
    tile_path.write_bytes(b"stub")
    tile_paths = {"C1": tile_path}

    monkeypatch.setattr(
        night_lights,
        "_sample_dataset_points",
        lambda _dataset, coords: np.full(len(coords), 2.0, dtype=np.float64),
    )
    monkeypatch.setattr(
        night_lights,
        "_open_dataset_cached",
        lambda *_args, **_kwargs: object(),
    )

    distances_m = np.asarray([1_000.0, 64_000.0, 128_000.0], dtype=np.float64)
    curve = night_lights._sample_ray_brightness_curve(
        tile_paths=tile_paths,
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        azimuth_deg=90.0,
        distances_m=distances_m,
    )

    assert np.allclose(
        curve,
        np.asarray([2.0, 4.0, 6.0], dtype=np.float64),
    )
