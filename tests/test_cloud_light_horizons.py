from __future__ import annotations

import numpy as np
import pytest

from zstarview.cloud_night_lighting import (
    DIRECT_LIGHT_HORIZON_FADE_DEG,
    cloud_horizon_depression_deg,
    cloud_light_visibility,
    moon_cloud_strength,
)
from zstarview.render.cloud_voxels import (
    _altaz_to_directions,
    _layer_lighting_metadata,
    _night_lighting_inputs,
    _render,
)


def test_cloud_horizons_extend_farther_below_observer_horizon_at_height() -> None:
    assert cloud_horizon_depression_deg(0) == 0
    assert cloud_horizon_depression_deg(1) == pytest.approx(1.0, abs=0.06)
    assert cloud_horizon_depression_deg(3) == pytest.approx(1.8, abs=0.06)
    assert cloud_horizon_depression_deg(11) == pytest.approx(3.4, abs=0.06)
    assert cloud_light_visibility(-2, 3) == 0
    assert cloud_light_visibility(-2, 11) == 1
    assert cloud_light_visibility(-5, 11) == 0


def test_solar_horizon_fade_is_continuous_and_only_above_geometric_limit() -> None:
    cutoff = -cloud_horizon_depression_deg(3)
    assert cloud_light_visibility(cutoff - 0.001, 3) == 0
    assert cloud_light_visibility(cutoff, 3) == 0
    assert 0 < cloud_light_visibility(cutoff + 0.05, 3) < 1
    assert cloud_light_visibility(cutoff + DIRECT_LIGHT_HORIZON_FADE_DEG + 0.001, 3) == 1
    assert cloud_light_visibility(cutoff + 1e-6, 3) < 1e-8


def test_moon_phase_survives_below_observer_horizon_until_layer_evaluation() -> None:
    assert moon_cloud_strength(-2, 180) == 0
    assert moon_cloud_strength(-2, 180, 3) == 0
    assert moon_cloud_strength(-2, 180, 11) > moon_cloud_strength(-2, 90, 11) > 0
    assert moon_cloud_strength(-2, 0, 11) == 0
    assert moon_cloud_strength(-5, 180, 11) == 0
    _, peak, _ = _night_lighting_inputs(
        (-2, 90, 180), None, 0, 35, 135, (1, 1, 1), np.zeros(2),
        np.zeros((1, 2)), np.eye(2)[None],
    )
    assert peak == moon_cloud_strength(90, 180)


def test_source_metadata_exposes_actual_native_and_geo_layer_limits() -> None:
    native = _layer_lighting_metadata(-2, (-2, 90, 180), 2.5, 9)
    geo = _layer_lighting_metadata(-2, (-2, 90, 180), 0.5, 9)
    assert native["direct_light_horizon_alt_deg_by_layer"][0] == pytest.approx(-1.8, abs=0.06)
    assert geo["direct_light_horizon_alt_deg_by_layer"][0] == pytest.approx(-1.0, abs=0.06)
    for metadata in (native, geo):
        assert metadata["moon_light_strength"] == 0  # Observer reference only.
        assert metadata["moon_peak_light_strength"] > 0
        assert metadata["moon_light_strength_by_layer"][0] == 0
        assert metadata["moon_light_strength_by_layer"][-1] > 0
        assert metadata["sun_visibility_by_layer"][0] == 0
        assert metadata["sun_visibility_by_layer"][-1] == 1


def _single_layer(height_km, altitude, *, solar, blocker=0):
    density = np.array([[[0.4]], [[blocker]]])
    light = _altaz_to_directions([altitude], [90])[0]
    return _render(
        density, np.array([-1.0, 0.5, 0.5]), np.array([[1.0, 0.0, 0.0]]),
        light, 1.0 if solar else 0.0, np.array([0.36, 0.385, 0.44]),
        np.zeros((1, 3)), 1.0, np.zeros(2), np.eye(2)[None], np.zeros((1, 2)),
        np.eye(2)[None], light, 0.25 if not solar else 0.0, None,
        layer_base_km=height_km - 0.5,
    )


@pytest.mark.parametrize("solar", [True, False])
def test_below_horizon_light_reaches_high_layer_without_changing_view_transmission(solar) -> None:
    low_rgb, low_transmission = _single_layer(3, -2, solar=solar)
    high_rgb, high_transmission = _single_layer(11, -2, solar=solar)
    assert np.all(high_rgb > low_rgb)
    assert np.all(low_rgb > 0)  # Diffuse solar or weak night ambient remains.
    np.testing.assert_array_equal(high_transmission, low_transmission)
    deeply_set_low, _ = _single_layer(3, -5, solar=solar)
    deeply_set_high, _ = _single_layer(11, -5, solar=solar)
    np.testing.assert_array_equal(deeply_set_high, deeply_set_low)


@pytest.mark.parametrize("solar", [True, False])
def test_geo_and_native_low_layers_use_their_actual_height(solar) -> None:
    # At -1.3 deg, a 1-km Geo layer is shadowed but a 3-km native layer is lit.
    geo_rgb, geo_transmission = _single_layer(1, -1.3, solar=solar)
    native_rgb, native_transmission = _single_layer(3, -1.3, solar=solar)
    assert np.all(native_rgb > geo_rgb)
    np.testing.assert_array_equal(native_transmission, geo_transmission)


def test_lunar_visibility_boundary_is_continuous_and_zero_at_geometric_cutoff() -> None:
    cutoff = -cloud_horizon_depression_deg(11)
    assert moon_cloud_strength(cutoff, 180, 11) == 0
    assert moon_cloud_strength(cutoff - 0.01, 180, 11) == 0
    assert moon_cloud_strength(cutoff + 1e-6, 180, 11) < 1e-12
    assert moon_cloud_strength(cutoff + 0.1, 180, 11) > 0
