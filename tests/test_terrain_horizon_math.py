from __future__ import annotations

import numpy as np

from zstarview.terrain.horizon import build_distance_samples, compute_apparent_altitudes


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
