from __future__ import annotations

import numpy as np

from zstarview.render.atmosphere import atmospheric_sky_samples


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


def test_view_rays_into_earth_produce_no_sky_radiance() -> None:
    colors = _samples([-1.0, -10.0], [0.0, 180.0], (30.0, 0.0))

    assert np.array_equal(colors, np.zeros((2, 3), dtype=np.float32))


def test_mie_scattering_is_stronger_near_the_solar_direction() -> None:
    colors = _samples([2.0, 2.0], [0.0, 90.0], (2.0, 0.0))

    assert float(colors[0].mean()) > float(colors[1].mean())


def test_high_atmosphere_remains_lit_after_sun_sets() -> None:
    colors = _samples([10.0, 20.0], [0.0, 90.0], (-3.0, 0.0))

    assert float(colors.max()) > 0.0


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
