from __future__ import annotations

import numpy as np
import pytest

from zstarview.render.cloud_voxels import _render


def _render_with_neighbor(
    blocker_amount: float,
    sunlight_mix: float,
    solar_vertical: float,
    target_amount: float = 0.4,
    layer_count: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    # The observer sees only the first cell. Its light ray crosses a second
    # cell outside the viewing ray, so the second cell casts only a shadow.
    density = np.array([[[target_amount]], [[blocker_amount]]], dtype=np.float64)
    density = np.repeat(density / layer_count, layer_count, axis=2)
    sun = np.array([1.0, 0.0, solar_vertical])
    sun /= np.linalg.norm(sun)
    basis = np.repeat(np.eye(2)[None, :, :], layer_count, axis=0)
    return _render(
        density,
        np.array([0.5, 0.5, -0.5]),
        np.array([[0.0, 0.0, 1.0]]),
        sun,
        sunlight_mix,
        np.array([0.36, 0.385, 0.44]),
        np.zeros((1, 3)),
        1.0,
        np.zeros(2),
        basis,
        np.zeros((layer_count, 2)),
        basis,
    )


@pytest.mark.parametrize("sunlight_mix", [0.5, 1.0])
@pytest.mark.parametrize("solar_vertical", [-0.2, 0.2])
def test_neighbor_casts_shadow_without_changing_view_transmission(
    sunlight_mix: float, solar_vertical: float
) -> None:
    clear_rgb, clear_transmission = _render_with_neighbor(
        0.0, sunlight_mix, solar_vertical
    )
    shaded_rgb, shaded_transmission = _render_with_neighbor(
        4.0, sunlight_mix, solar_vertical
    )

    assert np.all(shaded_rgb < clear_rgb)
    assert np.all(shaded_rgb > 0.0)
    np.testing.assert_array_equal(shaded_transmission, clear_transmission)


@pytest.mark.parametrize("solar_vertical", [-0.2, 0.2])
def test_ambient_night_color_has_no_directional_shadow(solar_vertical: float) -> None:
    clear_rgb, clear_transmission = _render_with_neighbor(0.0, 0.0, -0.2)
    blocked_rgb, blocked_transmission = _render_with_neighbor(
        4.0, 0.0, solar_vertical
    )

    np.testing.assert_array_equal(blocked_rgb, clear_rgb)
    np.testing.assert_array_equal(blocked_transmission, clear_transmission)


def test_ambient_night_color_increases_with_cloud_amount() -> None:
    thin_rgb, thin_transmission = _render_with_neighbor(0.0, 0.0, -0.2, 0.1)
    thick_rgb, thick_transmission = _render_with_neighbor(0.0, 0.0, -0.2, 0.8)
    thin_alpha = 1.0 - thin_transmission[0]
    thick_alpha = 1.0 - thick_transmission[0]

    assert np.all(thick_rgb > thin_rgb)
    assert np.all(thick_rgb / thick_alpha > thin_rgb / thin_alpha)
    assert thick_transmission[0] < thin_transmission[0]


def test_night_color_preserves_cloud_amount_distributed_across_layers() -> None:
    single_rgb, single_transmission = _render_with_neighbor(
        0.0, 0.0, -0.2, 0.8, layer_count=1
    )
    layered_rgb, layered_transmission = _render_with_neighbor(
        0.0, 0.0, -0.2, 0.8, layer_count=9
    )

    np.testing.assert_allclose(layered_rgb, single_rgb)
    np.testing.assert_allclose(layered_transmission, single_transmission)


@pytest.mark.parametrize("blocker_amount", [0.0, 4.0])
def test_twilight_blends_directional_day_and_ambient_night(blocker_amount: float) -> None:
    day_rgb, _ = _render_with_neighbor(blocker_amount, 1.0, -0.2)
    night_rgb, _ = _render_with_neighbor(blocker_amount, 0.0, -0.2)
    mixed_rgb, _ = _render_with_neighbor(blocker_amount, 0.5, -0.2)

    np.testing.assert_allclose(mixed_rgb, 0.5 * (day_rgb + night_rgb))
