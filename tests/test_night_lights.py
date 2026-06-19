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
