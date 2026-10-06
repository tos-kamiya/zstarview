from __future__ import annotations

import numpy as np
import pytest

from zstarview.render.cloud_voxels import (
    _altaz_to_directions,
    _validate_cloud_appearance,
)


def test_validate_cloud_appearance_clips_values_and_preserves_none() -> None:
    sunlight_mix, night_color, subtraction = _validate_cloud_appearance(
        1.5, [-0.2, 0.5, 1.2], None
    )

    assert sunlight_mix == 1.0
    np.testing.assert_array_equal(night_color, np.array([0.0, 0.5, 1.0]))
    assert night_color.dtype == np.float32
    assert subtraction is None


@pytest.mark.parametrize(
    ("sunlight_mix", "night_color_rgb", "subtraction", "message"),
    [
        (np.nan, [0.0, 0.0, 0.0], None, "sunlight_mix must be finite"),
        (0.5, [0.0, np.inf, 0.0], None, "night_color_rgb must contain"),
        (0.5, [0.0, 0.0, 0.0], -0.1, "cloud_amount_subtract must be finite"),
        (0.5, [0.0, 0.0, 0.0], 1.1, "cloud_amount_subtract must be finite"),
    ],
)
def test_validate_cloud_appearance_rejects_invalid_values(
    sunlight_mix: float,
    night_color_rgb: list[float],
    subtraction: float | None,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        _validate_cloud_appearance(sunlight_mix, night_color_rgb, subtraction)


def test_altaz_to_directions_handles_broadcast_angles_and_cardinal_axes() -> None:
    directions = _altaz_to_directions(
        np.array([[0.0], [90.0]]),
        np.array([[0.0, 90.0, 180.0]]),
    )

    assert directions.shape == (6, 3)
    np.testing.assert_allclose(
        directions,
        np.array(
            [
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
            ]
        ),
        atol=1e-15,
    )
