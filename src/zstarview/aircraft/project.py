from __future__ import annotations

import math
from typing import Iterable

import astropy.time
from astropy import units as u
from astropy.coordinates import EarthLocation

from .types import AircraftOverlayPoint, AircraftSnapshot


def project_aircraft_snapshots(
    snapshots: Iterable[AircraftSnapshot],
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
) -> list[AircraftOverlayPoint]:
    observer_location = EarthLocation(
        lat=float(observer_lat) * u.deg,
        lon=float(observer_lon) * u.deg,
        height=float(observer_height_m) * u.m,
    )
    observer_xyz = observer_location.to_geocentric()
    obs_x = float(observer_xyz[0].to_value(u.m))
    obs_y = float(observer_xyz[1].to_value(u.m))
    obs_z = float(observer_xyz[2].to_value(u.m))
    lat_rad = math.radians(float(observer_lat))
    lon_rad = math.radians(float(observer_lon))
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)

    overlay_points: list[AircraftOverlayPoint] = []
    for snapshot in snapshots:
        target_location = EarthLocation(
            lat=float(snapshot.latitude) * u.deg,
            lon=float(snapshot.longitude) * u.deg,
            height=float(snapshot.baro_altitude_m or 0.0) * u.m,
        )
        target_xyz = target_location.to_geocentric()
        dx = float(target_xyz[0].to_value(u.m)) - obs_x
        dy = float(target_xyz[1].to_value(u.m)) - obs_y
        dz = float(target_xyz[2].to_value(u.m)) - obs_z

        east_m = (-sin_lon * dx) + (cos_lon * dy)
        north_m = (-sin_lat * cos_lon * dx) - (sin_lat * sin_lon * dy) + (cos_lat * dz)
        up_m = (cos_lat * cos_lon * dx) + (cos_lat * sin_lon * dy) + (sin_lat * dz)

        horizontal_m = math.hypot(east_m, north_m)
        distance_km = math.sqrt((horizontal_m * horizontal_m) + (up_m * up_m)) / 1000.0
        alt_deg = math.degrees(math.atan2(up_m, horizontal_m))
        if alt_deg <= 0.0:
            continue
        az_deg = math.degrees(math.atan2(east_m, north_m)) % 360.0
        overlay_points.append(
            AircraftOverlayPoint(
                icao24=snapshot.icao24,
                callsign=snapshot.callsign,
                alt_deg=alt_deg,
                az_deg=az_deg,
                distance_km=distance_km,
            )
        )
    return overlay_points
