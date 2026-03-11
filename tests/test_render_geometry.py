from __future__ import annotations

import math

from zstarview.render.geometry import (
    altaz_to_cylindrical_normalized_xy,
    shortest_azimuth_delta_deg,
)


def test_shortest_azimuth_delta_deg_wraps_across_north() -> None:
    assert math.isclose(shortest_azimuth_delta_deg(359.0, 1.0), -2.0, abs_tol=1e-9)
    assert math.isclose(shortest_azimuth_delta_deg(1.0, 359.0), 2.0, abs_tol=1e-9)


def test_altaz_to_cylindrical_normalized_xy_is_linear_in_azimuth_and_altitude() -> None:
    nx, ny = altaz_to_cylindrical_normalized_xy(40.0, 150.0, (30.0, 120.0))

    assert math.isclose(nx, 30.0 / 90.0, abs_tol=1e-9)
    assert math.isclose(ny, -(10.0 / 90.0), abs_tol=1e-9)
