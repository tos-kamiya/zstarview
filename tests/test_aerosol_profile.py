from __future__ import annotations

import numpy as np

from zstarview.render.aerosol_profile import (
    DEFAULT_AOD550,
    bundled_aod550_or_default,
    load_bundled_climatology,
)
from zstarview.render.atmosphere import atmospheric_sky_samples


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
