from __future__ import annotations

import numpy as np
import pytest

from zstarview.render.cloud_voxels import _refine_geo_cloud_amount, _render, _segments


@pytest.mark.parametrize("factor", [3, 4])
@pytest.mark.parametrize("sign", [-1, 1])
def test_mixed_cell_sizes_preserve_distance_and_extinction(
    factor: int, sign: int
) -> None:
    # A horizontal ray crosses coarse cells on both sides of a refined patch.
    coarse_shape = np.array([6, 4, 1])
    origin = np.array([-1.0 if sign > 0 else 7.0, 1.25, 0.5])
    ray = np.array([float(sign), 0.0, 0.0])
    basis = np.eye(2)[None]
    centers = np.zeros((1, 2))
    window = (2 * factor, 4 * factor, factor, 3 * factor)
    segments = _segments(
        origin,
        ray,
        coarse_shape * np.array([factor, factor, 1]),
        np.zeros(2),
        basis * factor,
        centers,
        window,
        factor,
    )
    assert sum(segment[3] for segment in segments) == pytest.approx(6.0)
    assert len(segments) == 4 + 2 * factor
    density = np.full((6 * factor, 4 * factor, 1), 0.2)
    rgb, transmission = _render(
        density,
        origin,
        ray[None],
        np.array([0.0, 0.0, 1.0]),
        0.0,
        np.ones(3),
        np.zeros((1, 3)),
        1.0,
        False,
        np.zeros(2),
        basis * factor,
        centers,
        basis / factor,
        window,
        factor,
    )
    expected = np.exp(-2.4 * 0.2 * 6.0)
    assert transmission[0] == pytest.approx(expected)
    np.testing.assert_allclose(rgb[0], 0.2 * (1.0 - expected))


def test_refined_view_receives_shadow_from_distant_coarse_cell() -> None:
    factor = 4
    basis = np.eye(2)[None]
    window = (0, factor, 0, factor)

    def render(blocker: float):
        density = np.full((3 * factor, factor, 1), 0.4)
        density[2 * factor :] = blocker
        sun = np.array([1.0, 0.0, 0.1])
        sun /= np.linalg.norm(sun)
        return _render(
            density,
            np.array([0.5, 0.5, -0.5]),
            np.array([[0.0, 0.0, 1.0]]),
            sun,
            1.0,
            np.ones(3),
            np.zeros((1, 3)),
            1.0,
            False,
            np.zeros(2),
            basis * factor,
            np.zeros((1, 2)),
            basis / factor,
            window,
            factor,
        )

    clear_rgb, clear_transmission = render(0.0)
    shaded_rgb, shaded_transmission = render(4.0)
    assert np.all(shaded_rgb < clear_rgb)
    np.testing.assert_allclose(shaded_transmission, clear_transmission)


def test_interpolation_refines_every_column_and_ignores_missing() -> None:
    amount = np.full((4, 4), 0.4, dtype=np.float32)
    valid = np.ones_like(amount, dtype=bool)
    amount[1, 1] = 0.0
    valid[1, 1] = False
    fine = _refine_geo_cloud_amount(amount, 3, valid)
    assert fine.shape == (12, 12)
    np.testing.assert_array_equal(fine[3:6, 3:6], 0.0)
    np.testing.assert_allclose(fine[6:9, 3:6], 0.4)
    assert fine[0, 0] == pytest.approx(0.4)
    assert fine[-1, -1] == pytest.approx(0.4)


def test_interpolation_reduces_jump_at_nearby_original_cell_boundary() -> None:
    amount = np.tile(np.linspace(0.0, 1.0, 4), (4, 1))
    fine = _refine_geo_cloud_amount(amount, 3, np.ones_like(amount, dtype=bool))
    coarse_jump = amount[1, 2] - amount[1, 1]
    fine_jump = fine[4, 6] - fine[4, 5]
    assert fine_jump == pytest.approx(coarse_jump / 3)


def test_full_refinement_preserves_distance_based_extinction() -> None:
    factor = 3
    density = np.full((6 * factor, factor, 1), 0.2)
    window = (0, density.shape[0], 0, density.shape[1])
    basis = np.eye(2)[None]
    origin = np.array([-1.0, 0.5, 0.5])
    ray = np.array([[1.0, 0.0, 0.0]])
    segments = _segments(
        origin,
        ray[0],
        np.array(density.shape),
        np.zeros(2),
        basis * factor,
        np.zeros((1, 2)),
        window,
        factor,
    )
    assert len(segments) == 6 * factor
    assert sum(segment[3] for segment in segments) == pytest.approx(6.0)
    _, transmission = _render(
        density,
        origin,
        ray,
        np.array([0.0, 0.0, 1.0]),
        0.0,
        np.ones(3),
        np.zeros((1, 3)),
        1.0,
        False,
        np.zeros(2),
        basis * factor,
        np.zeros((1, 2)),
        basis / factor,
        window,
        factor,
    )
    assert transmission[0] == pytest.approx(np.exp(-2.4 * 0.2 * 6.0))
