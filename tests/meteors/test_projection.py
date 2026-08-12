from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

from zstarview.meteors.projection import (
    project_meteor_observations_to_altaz,
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


def test_same_location_trail_is_fixed_near_zenith_in_event_altaz() -> None:
    observation = _observation(
        end_lat_deg=35.0,
        end_lon_deg=139.0,
    )

    trails = project_meteor_observations_to_altaz(
        [observation],
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
    )

    assert len(trails) == 1
    trail = trails[0]
    assert trail.begin_alt_deg > 89.9


def test_projection_rejects_observations_outside_candidate_radius() -> None:
    trails = project_meteor_observations_to_altaz(
        [_observation(begin_lat_deg=-35.0, end_lat_deg=-35.1)],
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=0.0,
    )

    assert trails == ()
