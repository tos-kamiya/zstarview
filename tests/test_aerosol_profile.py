from __future__ import annotations

import numpy as np

from zstarview.render.aerosol_profile import (
    DEFAULT_AOD550,
    bundled_aod550_or_default,
    load_bundled_climatology,
)
from zstarview.render.atmosphere import atmospheric_sky_samples
from zstarview.render.sky_disc import sky_color_near_solar_horizon, sky_color_samples


def test_bundled_climatology_has_global_monthly_grid() -> None:
    climatology = load_bundled_climatology()

    assert climatology.aod550.shape == (12, 241, 480)
    assert climatology.latitude[0] > climatology.latitude[-1]
    assert climatology.longitude[0] == 0.0
    assert np.all(np.isfinite(climatology.aod550))


def test_aod_lookup_wraps_longitude_and_selects_month() -> None:
    climatology = load_bundled_climatology()

    east = climatology.value_at(35.0, 139.0, 8)
    wrapped = climatology.value_at(35.0, -221.0, 8)

    assert east == wrapped
    assert east != climatology.value_at(35.0, 139.0, 1)


def test_invalid_asset_uses_fixed_default(monkeypatch) -> None:
    def fail_to_load():
        raise ValueError("test asset failure")

    monkeypatch.setattr(
        "zstarview.render.aerosol_profile.load_bundled_climatology",
        fail_to_load,
    )

    assert bundled_aod550_or_default(35.0, 139.0, 8) == DEFAULT_AOD550


def test_aod_changes_aerosol_scattering() -> None:
    altitudes = np.asarray([2.0, 20.0], dtype=np.float32)
    azimuths = np.asarray([0.0, 90.0], dtype=np.float32)
    clean = atmospheric_sky_samples(
        altitudes,
        azimuths,
        (2.0, 0.0),
        aerosol_optical_depth=0.02,
    )
    hazy = atmospheric_sky_samples(
        altitudes,
        azimuths,
        (2.0, 0.0),
        aerosol_optical_depth=1.0,
    )

    assert not np.array_equal(clean, hazy)


def test_aerosol_extinction_increases_sunset_redness() -> None:
    colors = atmospheric_sky_samples(
        np.asarray([0.0], dtype=np.float32),
        np.asarray([0.0], dtype=np.float32),
        (0.0, 0.0),
        aerosol_optical_depth=0.3,
    )[0]

    assert float(colors[0] / colors[2]) > 1.5


def test_daytime_rayleigh_color_keeps_green_below_blue() -> None:
    colors = atmospheric_sky_samples(
        np.asarray([20.0], dtype=np.float32),
        np.asarray([90.0], dtype=np.float32),
        (45.0, 0.0),
    )[0]

    assert float(colors[1]) < float(colors[2])


def test_solar_horizon_color_averages_zero_to_ten_degrees() -> None:
    sun_altaz = (2.0, 135.0)
    altitudes = np.linspace(0.0, 10.0, 6, dtype=np.float32)
    expected = sky_color_samples(
        altitudes,
        np.full_like(altitudes, 135.0),
        sun_altaz,
        aerosol_optical_depth=0.3,
    ).mean(axis=0)
    actual = sky_color_near_solar_horizon(
        sun_altaz,
        aerosol_optical_depth=0.3,
    )

    assert np.allclose(np.asarray(actual[:3]) / 255.0, expected, atol=1 / 255)
