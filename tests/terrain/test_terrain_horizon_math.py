from __future__ import annotations

import numpy as np
from rasterio.transform import Affine

from zstarview.terrain.horizon import build_distance_samples, compute_apparent_altitudes
from zstarview.terrain.dem import DemGrid, sample_ground_elevation


class _IdentityTransformer:
    def transform(self, lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return lon_deg, lat_deg


def test_build_distance_samples_covers_requested_distance() -> None:
    samples = build_distance_samples(1.0, 300.0)

    assert len(samples) == 4
    assert samples[0] == 300.0
    assert samples[-1] == 1000.0


def test_compute_apparent_altitudes_is_zero_for_equal_height_flat_short_range() -> None:
    got = compute_apparent_altitudes(
        observer_elevation_m=100.0,
        target_elevation_m=np.array([100.0]),
        surface_distance_m=np.array([1.0]),
        earth_radius_m=6_371_008.8,
        refraction_coefficient=0.0,
    )

    assert abs(float(got[0])) < 1e-3


def test_sample_ground_elevation_uses_default_elevation_for_nodata() -> None:
    grid = DemGrid(
        array=np.array([[np.nan]], dtype=np.float64),
        transform=Affine.translation(0.0, 1.0) * Affine.scale(1.0, -1.0),
        crs="EPSG:4326",
        default_elevation_m=0.0,
        _to_grid=_IdentityTransformer(),
    )

    got = sample_ground_elevation(
        grid,
        latitude_deg=0.5,
        longitude_deg=0.5,
    )

    assert got == 0.0
