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
        last_contact_unix=int(time_obj.unix),
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
    assert points[0].age_seconds == 0.0
    assert points[0].alpha_scale == 1.0
    assert len(points[0].trail_alt_az_points) == 5
    assert len(points[0].trail_geodetic_points) == 5
    assert points[0].trail_alt_az_points[0][1] != points[0].trail_alt_az_points[-1][1]


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
        last_contact_unix=int(time_obj.unix),
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
        last_contact_unix=int(time_obj.unix),
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


def test_project_aircraft_snapshots_fades_after_ninety_seconds() -> None:
    base_time = astropy.time.Time("2026-03-18T12:00:00Z")
    current_time = astropy.time.Time("2026-03-18T12:02:00Z")
    snapshot = AircraftSnapshot(
        icao24="fade01",
        callsign="FADE",
        latitude=35.47,
        longitude=133.05,
        baro_altitude_m=1000.0,
        velocity_mps=250.0,
        heading_deg=90.0,
        vertical_rate_mps=0.0,
        on_ground=False,
        last_contact_unix=int(base_time.unix),
    )

    points = project_aircraft_snapshots(
        [snapshot],
        observer_lat=35.47,
        observer_lon=133.05,
        observer_height_m=0.0,
        time_obj=current_time,
    )

    assert len(points) == 1
    assert 119.0 <= points[0].age_seconds <= 121.0
    assert 0.3 < points[0].alpha_scale < 1.0


def test_project_aircraft_snapshots_builds_eight_second_motion_segment() -> None:
    base_time = astropy.time.Time("2026-03-18T12:00:00Z")
    current_time = astropy.time.Time("2026-03-18T12:00:20Z")
    snapshot = AircraftSnapshot(
        icao24="line01",
        callsign="LINE",
        latitude=35.47,
        longitude=133.05,
        baro_altitude_m=5000.0,
        velocity_mps=250.0,
        heading_deg=90.0,
        vertical_rate_mps=0.0,
        on_ground=False,
        last_contact_unix=int(base_time.unix),
    )

    points = project_aircraft_snapshots(
        [snapshot],
        observer_lat=35.47,
        observer_lon=133.05,
        observer_height_m=0.0,
        time_obj=current_time,
    )

    assert len(points) == 1
    point = points[0]
    assert len(point.trail_alt_az_points) == 5
    assert len(point.trail_geodetic_points) == 5
    assert point.trail_alt_az_points[2] == (point.alt_deg, point.az_deg)
    assert len(set(point.trail_alt_az_points)) == 5
