from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

import astropy.time
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord

from zstarview.meteors.projection import (
    _clip_segment_to_geometric_horizon,
    project_meteor_observations_to_celestial,
)
from zstarview.meteors.types import MeteorObservation


def _observation(**overrides) -> MeteorObservation:
    base = MeteorObservation(
        trajectory_id="test-meteor",
        beginning_utc=datetime(2026, 8, 10, 12, tzinfo=timezone.utc),
        begin_lat_deg=35.0,
        begin_lon_deg=139.0,
        begin_height_km=105.0,
        end_lat_deg=35.02,
        end_lon_deg=139.02,
        end_height_km=85.0,
        duration_s=0.5,
        peak_abs_magnitude=-1.0,
        initial_speed_km_s=50.0,
        shower_code="PER",
    )
    return replace(base, **overrides)


def test_same_location_trail_is_fixed_near_zenith_in_icrs() -> None:
    observation = _observation(
        end_lat_deg=35.0,
        end_lon_deg=139.0,
    )

    trails = project_meteor_observations_to_celestial(
        [observation],
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
    )

    assert len(trails) == 1
    trail = trails[0]
    observer = EarthLocation(lat=35.0 * u.deg, lon=139.0 * u.deg, height=0.0 * u.m)
    frame = AltAz(obstime=astropy.time.Time(observation.beginning_utc), location=observer)
    begin_altaz = SkyCoord(
        ra=trail.begin_ra_deg * u.deg,
        dec=trail.begin_dec_deg * u.deg,
        frame="icrs",
    ).transform_to(frame)
    assert begin_altaz.alt.deg > 89.9


def test_projection_rejects_observations_outside_candidate_radius() -> None:
    trails = project_meteor_observations_to_celestial(
        [_observation(begin_lat_deg=-35.0, end_lat_deg=-35.1)],
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
    )

    assert trails == ()


def test_horizon_clip_discards_fully_below_segment() -> None:
    observer = np.asarray([10.0, 0.0, 0.0])
    up = np.asarray([1.0, 0.0, 0.0])

    clipped = _clip_segment_to_geometric_horizon(
        np.asarray([9.0, 1.0, 0.0]),
        np.asarray([8.0, 2.0, 0.0]),
        observer_xyz=observer,
        up=up,
    )

    assert clipped is None


def test_horizon_clip_interpolates_crossing_endpoint() -> None:
    observer = np.asarray([10.0, 0.0, 0.0])
    up = np.asarray([1.0, 0.0, 0.0])

    clipped = _clip_segment_to_geometric_horizon(
        np.asarray([9.0, 1.0, 0.0]),
        np.asarray([11.0, 3.0, 0.0]),
        observer_xyz=observer,
        up=up,
    )

    assert clipped is not None
    begin, end = clipped
    assert begin == pytest.approx(np.asarray([10.0, 2.0, 0.0]))
    assert end == pytest.approx(np.asarray([11.0, 3.0, 0.0]))
