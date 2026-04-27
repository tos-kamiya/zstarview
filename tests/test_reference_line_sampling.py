from astropy.coordinates import EarthLocation
from astropy.time import Time

from zstarview.astro import (
    calculate_celestial_equator_points,
    calculate_ecliptic_points,
    calculate_horizon_points,
)


def test_celestial_reference_lines_use_dense_sampling() -> None:
    location = EarthLocation(lat=35.68, lon=139.76, height=0.0)
    time = Time("2026-04-27T00:00:00", scale="utc")

    equator_points = calculate_celestial_equator_points(location, time)
    ecliptic_points = calculate_ecliptic_points(location, time)
    horizon_points = calculate_horizon_points()

    assert len(equator_points) == 73
    assert len(ecliptic_points) == 73
    assert len(horizon_points) == 73
