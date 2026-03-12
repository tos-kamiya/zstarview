from __future__ import annotations

import math

from zstarview.render.draw import sample_constant_altitude_arc, split_urban_profile_into_flat_runs
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


def test_split_urban_profile_into_flat_runs_breaks_on_altitude_change() -> None:
    runs = split_urban_profile_into_flat_runs(
        [
            (-5.0, 10.0),
            (-5.0, 11.0),
            (-4.0, 12.0),
            (-4.0, 13.0),
            (-90.0, 14.0),
            (-4.0, 15.0),
            (-4.0, 16.0),
        ]
    )

    assert runs == [
        [(-5.0, 10.0), (-5.0, 11.0)],
        [(-4.0, 12.0), (-4.0, 13.0)],
        [(-4.0, 15.0), (-4.0, 16.0)],
    ]


def test_sample_constant_altitude_arc_wraps_shortest_path() -> None:
    arc = sample_constant_altitude_arc(
        altitude_deg=-5.0,
        start_azimuth_deg=359.0,
        end_azimuth_deg=1.0,
        step_deg=0.5,
    )

    assert arc[0] == (-5.0, 359.0)
    assert arc[-1] == (-5.0, 1.0)
    assert any(abs(az) < 1.0e-9 for _alt, az in arc)
