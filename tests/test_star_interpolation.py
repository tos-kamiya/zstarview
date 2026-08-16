from __future__ import annotations

import numpy as np

from zstarview.render.star_interpolation import (
    STAR_INTERPOLATION_COVERAGE,
    STAR_INTERPOLATION_MAX_UPDATE_INTERVAL_SECONDS,
    build_star_interpolation_mesh,
    should_interpolate_stars,
)


def test_star_interpolation_coverage_is_configurable() -> None:
    assert 0.0 < STAR_INTERPOLATION_COVERAGE <= 1.0


def test_star_interpolation_is_disabled_above_ninety_second_update_interval() -> None:
    assert STAR_INTERPOLATION_MAX_UPDATE_INTERVAL_SECONDS == 90
    assert should_interpolate_stars(90)
    assert not should_interpolate_stars(90.1)


def test_star_interpolation_mesh_uses_square_100px_cells() -> None:
    source, target = build_star_interpolation_mesh(
        width_px=1600,
        height_px=900,
        geometry_center=(800.0, 450.0),
        geometry_radius=450.0,
        view_center_altaz_deg=(45.0, 0.0),
        observer_lat_deg=35.0,
        edge_fov_deg=90.0,
        elapsed_seconds=30.0,
    )

    assert source.shape == (17 * 10, 2)
    assert target.shape == source.shape
    assert np.all(np.isfinite(target))
    assert np.allclose(source[0], (0.0, 0.0))
    assert np.allclose(source[-1], (1600.0, 900.0))
