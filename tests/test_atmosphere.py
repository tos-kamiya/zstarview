from __future__ import annotations

import numpy as np

from zstarview.render import atmosphere
from zstarview.render.atmosphere import (
    atmospheric_sky_samples,
    direct_atmospheric_transmission_rgb,
)


def _samples(
    altitudes: list[float],
    azimuths: list[float],
    sun_altaz: tuple[float, float],
    *,
    view_steps: int = 32,
) -> np.ndarray:
    return atmospheric_sky_samples(
        np.asarray(altitudes, dtype=np.float32),
        np.asarray(azimuths, dtype=np.float32),
        sun_altaz,
        view_steps=view_steps,
        sun_steps=12,
    )


def test_atmospheric_samples_have_rgb_shape_and_display_range() -> None:
    colors = _samples([0.0, 10.0, 45.0], [0.0, 90.0, 180.0], (15.0, 0.0))

    assert colors.shape == (3, 3)
    assert colors.dtype == np.float32
    assert np.all(colors >= 0.0)
    assert np.all(colors <= 1.0)


def test_direct_solar_transmission_is_rgb_and_reduces_short_wavelengths() -> None:
    transmission = direct_atmospheric_transmission_rgb(45.0)

    assert transmission.shape == (3,)
    assert transmission.dtype == np.float32
    assert np.all((transmission > 0.0) & (transmission <= 1.0))
    assert float(transmission[0]) > float(transmission[2])


def test_direct_solar_transmission_becomes_redder_at_low_altitude() -> None:
    high = direct_atmospheric_transmission_rgb(60.0)
    low = direct_atmospheric_transmission_rgb(2.0, aerosol_optical_depth=0.3)

    assert float(low[0] / low[2]) > float(high[0] / high[2])


def test_direct_solar_transmission_rejects_sun_below_horizon() -> None:
    np.testing.assert_array_equal(
        direct_atmospheric_transmission_rgb(-1.0),
        np.zeros(3, dtype=np.float32),
    )


def test_view_rays_into_earth_produce_no_sky_radiance() -> None:
    colors = _samples([-1.0, -10.0], [0.0, 180.0], (30.0, 0.0))

    assert np.array_equal(colors, np.zeros((2, 3), dtype=np.float32))


def test_mie_scattering_is_stronger_near_the_solar_direction() -> None:
    colors = _samples([2.0, 2.0], [0.0, 90.0], (2.0, 0.0))

    assert float(colors[0].mean()) > float(colors[1].mean())


def test_high_atmosphere_remains_lit_after_sun_sets() -> None:
    colors = _samples([10.0, 20.0], [0.0, 90.0], (-3.0, 0.0))

    assert float(colors.max()) > 0.0


def test_ozone_shell_vertical_path_is_its_thickness() -> None:
    origin = np.asarray([[0.0, 0.0, atmosphere.EARTH_RADIUS_KM]], dtype=np.float32)
    direction = np.asarray([[0.0, 0.0, 1.0]], dtype=np.float32)

    path = atmosphere._shell_path_length(
        origin,
        direction,
        np.asarray([100.0], dtype=np.float32),
        atmosphere.EARTH_RADIUS_KM + atmosphere.OZONE_SHELL_BOTTOM_KM,
        atmosphere.EARTH_RADIUS_KM + atmosphere.OZONE_SHELL_TOP_KM,
    )

    np.testing.assert_allclose(
        path,
        [atmosphere.OZONE_SHELL_TOP_KM - atmosphere.OZONE_SHELL_BOTTOM_KM],
    )


def test_ozone_absorption_makes_twilight_relatively_bluer(monkeypatch) -> None:
    altitudes = np.asarray([20.0], dtype=np.float32)
    azimuths = np.asarray([90.0], dtype=np.float32)
    monkeypatch.setattr(
        atmosphere,
        "OZONE_EXTINCTION_RGB",
        np.asarray([0.01, 0.02, 0.00075], dtype=np.float32),
    )
    with_ozone = atmospheric_sky_samples(altitudes, azimuths, (-3.0, 0.0))[0]
    monkeypatch.setattr(
        atmosphere,
        "OZONE_EXTINCTION_RGB",
        np.zeros(3, dtype=np.float32),
    )
    without_ozone = atmospheric_sky_samples(altitudes, azimuths, (-3.0, 0.0))[0]

    assert float(with_ozone[0] / with_ozone[2]) < float(
        without_ozone[0] / without_ozone[2]
    )
    assert float(with_ozone[1] / with_ozone[2]) < float(
        without_ozone[1] / without_ozone[2]
    )


def test_twilight_multiple_scattering_lifts_the_upper_sky(monkeypatch) -> None:
    upper = atmospheric_sky_samples(
        np.asarray([75.0], dtype=np.float32),
        np.asarray([90.0], dtype=np.float32),
        (0.0, 0.0),
    )[0]
    lower = atmospheric_sky_samples(
        np.asarray([5.0], dtype=np.float32),
        np.asarray([90.0], dtype=np.float32),
        (0.0, 0.0),
    )[0]
    monkeypatch.setattr(
        atmosphere,
        "TWILIGHT_MULTIPLE_SCATTERING_RGB",
        np.zeros(3, dtype=np.float32),
    )
    upper_without = atmospheric_sky_samples(
        np.asarray([75.0], dtype=np.float32),
        np.asarray([90.0], dtype=np.float32),
        (0.0, 0.0),
    )[0]
    lower_without = atmospheric_sky_samples(
        np.asarray([5.0], dtype=np.float32),
        np.asarray([90.0], dtype=np.float32),
        (0.0, 0.0),
    )[0]

    assert float(upper[2]) > float(upper_without[2])
    assert float(upper[2] - upper_without[2]) > float(
        lower[2] - lower_without[2]
    )


def test_twilight_multiple_scattering_is_off_at_deep_night(monkeypatch) -> None:
    colors = atmospheric_sky_samples(
        np.asarray([75.0], dtype=np.float32),
        np.asarray([90.0], dtype=np.float32),
        (-20.0, 0.0),
    )
    monkeypatch.setattr(
        atmosphere,
        "TWILIGHT_MULTIPLE_SCATTERING_RGB",
        np.zeros(3, dtype=np.float32),
    )
    colors_without = atmospheric_sky_samples(
        np.asarray([75.0], dtype=np.float32),
        np.asarray([90.0], dtype=np.float32),
        (-20.0, 0.0),
    )

    np.testing.assert_allclose(colors, colors_without)


def test_sunset_scattering_balance_reduces_red_relative_to_green(monkeypatch) -> None:
    altitudes = np.asarray([10.0], dtype=np.float32)
    azimuths = np.asarray([90.0], dtype=np.float32)
    balanced = atmospheric_sky_samples(altitudes, azimuths, (-1.0, 270.0))[0]
    monkeypatch.setattr(
        atmosphere,
        "SCATTERING_RGB_BALANCE",
        np.ones(3, dtype=np.float32),
    )
    unbalanced = atmospheric_sky_samples(altitudes, azimuths, (-1.0, 270.0))[0]

    assert float(balanced[1] / balanced[0]) > float(
        unbalanced[1] / unbalanced[0]
    )


def test_sky_fades_as_sun_moves_below_horizon() -> None:
    colors = np.stack(
        [
            _samples([20.0], [90.0], (sun_alt, 0.0))[0]
            for sun_alt in (0.0, -3.0, -6.0)
        ]
    )

    assert float(colors[0].mean()) > float(colors[1].mean()) > float(colors[2].mean())


def test_thirty_two_view_steps_are_close_to_sixty_four() -> None:
    altitudes = [0.0, 5.0, 15.0, 45.0]
    azimuths = [0.0, 45.0, 90.0, 180.0]
    coarse = _samples(altitudes, azimuths, (0.0, 0.0), view_steps=32)
    fine = _samples(altitudes, azimuths, (0.0, 0.0), view_steps=64)

    assert float(np.abs(coarse - fine).mean()) < 0.02
    assert float(np.abs(coarse - fine).max()) < 0.05
