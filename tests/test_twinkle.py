import math
import random
from typing import cast

import numpy as np

from zstarview.render.twinkle import (
    TWINKLE_MAX_DISTANCE_DEG,
    TWINKLE_TARGET_COUNT,
    nearest_twinkle_star_index,
    sample_twinkle_direction,
    twinkle_alpha,
    twinkle_strength,
    spherical_distance_deg,
)
from zstarview.types import StarsTable


def _stars() -> StarsTable:
    return cast(
        StarsTable,
        {
            "star_index": np.array([10, 11, 12], dtype=np.int32),
            "alt": np.array([45.0, 45.4, -2.0]),
            "az": np.array([180.0, 180.2, 180.0]),
            "vmag": np.array([4.1, 2.0, 5.0]),
            "bv": np.array([0.0, 0.0, 0.0]),
            "size_factor": np.array([0.5, 0.5, 0.5]),
            "color_factor_base": np.array([0.5, 0.5, 0.5]),
        },
    )


def test_twinkle_strength_changes_linearly_with_altitude() -> None:
    assert math.isclose(twinkle_strength(0.1), 1.0)
    assert math.isclose(twinkle_strength(10.0), 1.0)
    assert math.isclose(twinkle_strength(50.0), 0.5)
    assert math.isclose(twinkle_strength(90.0), 0.0)
    assert math.isclose(twinkle_strength(0.0), 1.0)


def test_twinkle_alpha_interpolates_from_point_five_to_zero() -> None:
    assert math.isclose(twinkle_alpha(5.0), 0.5)
    assert math.isclose(twinkle_alpha(0.0), 0.5)
    assert math.isclose(twinkle_alpha(10.0), 0.5)
    assert math.isclose(twinkle_alpha(50.0), 0.25)
    assert math.isclose(twinkle_alpha(90.0), 0.0)
    assert twinkle_alpha(-1.0) == 0.0


def test_spherical_distance_handles_azimuth_wrap() -> None:
    assert math.isclose(spherical_distance_deg(45.0, 359.5, 45.0, 0.5), 0.7071, rel_tol=1e-3)


def test_nearest_twinkle_star_excludes_bright_and_below_horizon() -> None:
    assert nearest_twinkle_star_index(
        _stars(),
        target_alt_deg=45.0,
        target_az_deg=180.0,
        view_center=(45.0, 180.0),
        content_fov_deg=90.0,
        vmag_limit=7.0,
    ) == 10


def test_nearest_twinkle_star_rejects_distance_over_three_degrees() -> None:
    assert nearest_twinkle_star_index(
        _stars(),
        target_alt_deg=60.0,
        target_az_deg=180.0,
        view_center=(45.0, 180.0),
        content_fov_deg=90.0,
        vmag_limit=7.0,
    ) is None


def test_twinkle_search_radius_is_three_degrees() -> None:
    assert TWINKLE_MAX_DISTANCE_DEG == 3.0


def test_twinkle_target_count_is_thirty() -> None:
    assert TWINKLE_TARGET_COUNT == 30


def test_sample_twinkle_direction_uses_altitude_distribution() -> None:
    alt, az = sample_twinkle_direction(rng=random.Random(3))
    assert 10.0 <= alt <= 90.0
    assert 0.0 <= az < 360.0
