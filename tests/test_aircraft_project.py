import astropy.time

from zstarview.aircraft.project import project_aircraft_snapshots
from zstarview.aircraft.types import AircraftSnapshot


def test_project_aircraft_snapshots_places_same_location_overhead() -> None:
    time_obj = astropy.time.Time("2026-03-18T12:00:00Z")
    snapshot = AircraftSnapshot(
        icao24="abcd01",
        callsign="TEST123",
        latitude=35.47,
        longitude=133.05,
        baro_altitude_m=1000.0,
        velocity_mps=250.0,
        heading_deg=90.0,
        vertical_rate_mps=0.0,
        on_ground=False,
        last_contact_unix=1710000000,
    )

    points = project_aircraft_snapshots(
        [snapshot],
        observer_lat=35.47,
        observer_lon=133.05,
        observer_height_m=0.0,
        time_obj=time_obj,
    )

    assert len(points) == 1
    assert points[0].icao24 == "abcd01"
    assert points[0].alt_deg > 80.0
    assert points[0].distance_km > 0.9


def test_project_aircraft_snapshots_reports_farther_distance_for_farther_target() -> None:
    time_obj = astropy.time.Time("2026-03-18T12:00:00Z")
    near = AircraftSnapshot(
        icao24="near01",
        callsign="NEAR",
        latitude=35.47,
        longitude=133.05,
        baro_altitude_m=1000.0,
        velocity_mps=250.0,
        heading_deg=90.0,
        vertical_rate_mps=0.0,
        on_ground=False,
        last_contact_unix=1710000000,
    )
    far = AircraftSnapshot(
        icao24="far001",
        callsign="FAR",
        latitude=36.00,
        longitude=134.00,
        baro_altitude_m=10000.0,
        velocity_mps=250.0,
        heading_deg=90.0,
        vertical_rate_mps=0.0,
        on_ground=False,
        last_contact_unix=1710000000,
    )

    points = project_aircraft_snapshots(
        [near, far],
        observer_lat=35.47,
        observer_lon=133.05,
        observer_height_m=0.0,
        time_obj=time_obj,
    )

    assert len(points) == 2
    points_by_id = {point.icao24: point for point in points}
    assert points_by_id["far001"].distance_km > points_by_id["near01"].distance_km
