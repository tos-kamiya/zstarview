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


def test_interpolation_refines_only_nearest_thirty_six_columns_and_ignores_missing() -> (
    None
):
    amount = np.full((10, 10), 0.4, dtype=np.float32)
    valid = np.ones_like(amount, dtype=bool)
    amount[4, 4] = 0.0
    valid[4, 4] = False
    fine, window = _refine_geo_cloud_amount(amount, np.array([4.9, 4.9]), 4, valid)
    assert window == (8, 32, 8, 32)
    np.testing.assert_array_equal(fine[16:20, 16:20], 0.0)
    np.testing.assert_allclose(fine[20:24, 16:24], 0.4)
    np.testing.assert_array_equal(
        fine[:8], np.repeat(amount[:2], 4, axis=0).repeat(4, axis=1)
    )


def test_interpolation_reduces_jump_at_nearby_original_cell_boundary() -> None:
    amount = np.tile(np.linspace(0.0, 1.0, 6), (6, 1))
    fine, _ = _refine_geo_cloud_amount(amount, np.array([3.0, 3.0]), 4)
    coarse_jump = amount[2, 3] - amount[2, 2]
    fine_jump = fine[11, 12] - fine[11, 11]
    assert fine_jump == pytest.approx(coarse_jump / 4)
