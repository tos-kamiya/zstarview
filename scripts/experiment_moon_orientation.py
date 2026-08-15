"""Compare the current Moon rotation with a celestial-north projection."""

from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone
from typing import Any

import astropy.time
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from skyfield.api import Topos, load

from zstarview.astro import (
    altaz_to_normalized_xy,
    calculate_moon_north_up_screen_rotation,
    calculate_planets,
    load_ephemeris,
)


LAT = 35.65305006707789
LON = 139.67368395238296
EDGE_FOV_DEG = 90.0
START = datetime(2026, 8, 1, 12, tzinfo=timezone.utc)
SAMPLE_DAYS = 3
SAMPLES = 10
VIEW_CENTERS = (
    (0.0, 0.0),
    (30.0, 0.0),
    (60.0, 0.0),
    (30.0, 90.0),
    (30.0, 180.0),
    (30.0, 270.0),
    (70.0, 45.0),
)


def _angle(x: float, y: float) -> float:
    return math.degrees(math.atan2(y, x))


def _signed_delta(new_angle: float, old_angle: float) -> float:
    return (new_angle - old_angle + 180.0) % 360.0 - 180.0


def _celestial_north_screen_rotation(
    target: datetime,
    view_center: tuple[float, float],
    planets: Any,
) -> float:
    """Project the apparent celestial-north tangent direction onto the screen."""
    timescale = load.timescale()
    observer = planets["earth"] + Topos(latitude_degrees=LAT, longitude_degrees=LON)
    apparent = observer.at(timescale.from_datetime(target)).observe(planets["moon"]).apparent()
    ra, dec, _distance = apparent.radec()
    moon_coord = SkyCoord(ra=ra.hours * 15.0 * u.deg, dec=dec.degrees * u.deg, frame="icrs")
    north_coord = moon_coord.directional_offset_by(0.0 * u.deg, 0.01 * u.deg)
    astropy_target = astropy.time.Time(target)
    earth_location = EarthLocation(
        lat=LAT * u.deg,
        lon=LON * u.deg,
        height=0.0 * u.m,
    )
    frame = AltAz(obstime=astropy_target, location=earth_location)
    moon_altaz_coord = moon_coord.transform_to(frame)
    north_altaz_coord = north_coord.transform_to(frame)
    moon_xy = altaz_to_normalized_xy(
        moon_altaz_coord.alt.deg,
        moon_altaz_coord.az.deg,
        view_center,
        edge_fov_deg=EDGE_FOV_DEG,
    )
    north_xy = altaz_to_normalized_xy(
        north_altaz_coord.alt.deg,
        north_altaz_coord.az.deg,
        view_center,
        edge_fov_deg=EDGE_FOV_DEG,
    )
    screen_north_angle = _angle(north_xy[0] - moon_xy[0], north_xy[1] - moon_xy[1])
    return screen_north_angle + 90.0


def main() -> None:
    planets = load_ephemeris()
    errors: list[float] = []
    print("date,view_alt,view_az,formula_rotation,reference_rotation,error")
    for index in range(SAMPLES):
        target = START + timedelta(days=index * SAMPLE_DAYS)
        bodies = calculate_planets(
            LAT,
            LON,
            0.0,
            astropy.time.Time(target),
            VIEW_CENTERS[0],
            planets,
            content_fov_deg=EDGE_FOV_DEG,
        )
        moon = {body.name: body for body in bodies}["moon"]
        for view_center in VIEW_CENTERS:
            formula_rotation = calculate_moon_north_up_screen_rotation(
                (moon.alt, moon.az),
                view_center,
                observer_latitude_deg=LAT,
                edge_fov_deg=EDGE_FOV_DEG,
            )
            reference_rotation = _celestial_north_screen_rotation(
                target,
                view_center,
                planets,
            )
            error = abs(_signed_delta(formula_rotation, reference_rotation))
            errors.append(error)
            print(
                f"{target.date()},{view_center[0]:.0f},{view_center[1]:.0f},"
                f"{formula_rotation:.3f},{reference_rotation:.3f},{error:.3f}"
            )
    print(
        f"summary,n={len(errors)},mean_abs_error={sum(errors) / len(errors):.4f},"
        f"max_abs_error={max(errors):.4f},within_0_5={sum(value <= 0.5 for value in errors)}"
    )


if __name__ == "__main__":
    main()
