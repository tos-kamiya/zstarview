from __future__ import annotations

import numpy as np

from zstarview.clouddisc.projectors.az import (
    az_project_lonlat_grid,
    build_projection_context,
    project_lonlat_grid_from_context,
)


def test_projection_context_matches_public_projection_api() -> None:
    context = build_projection_context(
        lat0_deg=35.0,
        lon0_deg=139.0,
        alt0_deg=45.0,
        az0_deg=180.0,
        radius_px=4,
        alt_min_deg=0.0,
        mask_fov_deg=90.0,
        edge_fov_deg=90.0,
    )

    expected = az_project_lonlat_grid(
        lat0_deg=35.0,
        lon0_deg=139.0,
        alt0_deg=45.0,
        az0_deg=180.0,
        radius_px=4,
        cloud_shell_km=6374.0,
        alt_min_deg=0.0,
        mask_fov_deg=90.0,
        edge_fov_deg=90.0,
    )
    actual = project_lonlat_grid_from_context(context, 6374.0)

    for got, want in zip(actual, expected, strict=True):
        np.testing.assert_allclose(got, want, equal_nan=True)


def test_projection_context_reuses_shell_independent_mask() -> None:
    context = build_projection_context(
        lat0_deg=35.0,
        lon0_deg=139.0,
        alt0_deg=45.0,
        az0_deg=180.0,
        radius_px=4,
        alt_min_deg=0.0,
        mask_fov_deg=90.0,
        edge_fov_deg=90.0,
    )

    _, _, mask_low = project_lonlat_grid_from_context(context, 6374.0)
    _, _, mask_high = project_lonlat_grid_from_context(context, 6378.0)

    assert np.array_equal(mask_low, mask_high)
