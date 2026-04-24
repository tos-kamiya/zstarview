from __future__ import annotations

import math

from zstarview.astro import calculate_moon_render_data


def _signed_delta_deg(new_angle: float, old_angle: float) -> float:
    return (float(new_angle) - float(old_angle) + 180.0) % 360.0 - 180.0


def test_moon_screen_rotation_follows_zenith_view_azimuth() -> None:
    _sun_dir, base_rotation = calculate_moon_render_data(
        (20.0, 100.0),
        (45.0, 180.0),
        (90.0, 0.0),
    )
    _sun_dir, rotated_view_rotation = calculate_moon_render_data(
        (20.0, 100.0),
        (45.0, 180.0),
        (90.0, 45.0),
    )

    assert math.isclose(
        _signed_delta_deg(rotated_view_rotation, base_rotation),
        45.0,
        abs_tol=1.0e-6,
    )


def test_moon_screen_rotation_stays_stable_for_same_view() -> None:
    _sun_dir, rotation_a = calculate_moon_render_data(
        (15.0, 80.0),
        (40.0, 170.0),
        (60.0, 120.0),
    )
    _sun_dir, rotation_b = calculate_moon_render_data(
        (15.0, 80.0),
        (40.0, 170.0),
        (60.0, 120.0),
    )

    assert rotation_a == rotation_b
