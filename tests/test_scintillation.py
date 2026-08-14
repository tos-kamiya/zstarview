import math
import random
from typing import cast

import numpy as np

from zstarview.render.scintillation import (
    SCINTILLATION_MAX_DISTANCE_DEG,
    SCINTILLATION_TARGET_COUNT,
    nearest_scintillation_star_index,
    sample_scintillation_direction,
    scintillation_alpha,
    scintillation_strength,
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


def test_scintillation_strength_clamps_below_ten_degrees() -> None:
    assert math.isclose(scintillation_strength(0.1), 1.0)
    assert math.isclose(scintillation_strength(10.0), 1.0)
    assert scintillation_strength(45.0) < 1.0
    assert scintillation_strength(0.0) == 0.0


def test_scintillation_alpha_uses_point_six_at_ten_degrees() -> None:
    assert math.isclose(scintillation_alpha(10.0), 0.6)
    assert scintillation_alpha(-1.0) == 0.0


def test_spherical_distance_handles_azimuth_wrap() -> None:
    assert math.isclose(spherical_distance_deg(45.0, 359.5, 45.0, 0.5), 0.7071, rel_tol=1e-3)


def test_nearest_scintillation_star_excludes_bright_and_below_horizon() -> None:
    assert nearest_scintillation_star_index(
        _stars(),
        target_alt_deg=45.0,
        target_az_deg=180.0,
        view_center=(45.0, 180.0),
        content_fov_deg=90.0,
        vmag_limit=7.0,
    ) == 10


def test_nearest_scintillation_star_rejects_distance_over_three_degrees() -> None:
    assert nearest_scintillation_star_index(
        _stars(),
        target_alt_deg=60.0,
        target_az_deg=180.0,
        view_center=(45.0, 180.0),
        content_fov_deg=90.0,
        vmag_limit=7.0,
    ) is None


def test_scintillation_search_radius_is_three_degrees() -> None:
    assert SCINTILLATION_MAX_DISTANCE_DEG == 3.0


def test_scintillation_target_count_is_ten() -> None:
    assert SCINTILLATION_TARGET_COUNT == 10


def test_sample_scintillation_direction_is_in_reasonable_view() -> None:
    alt, az = sample_scintillation_direction(
        (45.0, 180.0), 90.0, rng=random.Random(3)
    )
    assert -45.0 <= alt <= 90.0
    assert 0.0 <= az < 360.0
