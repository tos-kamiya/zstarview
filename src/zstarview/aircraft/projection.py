from __future__ import annotations

import math
from dataclasses import dataclass

from astropy import units as u
from astropy.coordinates import EarthLocation


@dataclass(frozen=True, slots=True)
class ObserverProjectionState:
    obs_x: float
    obs_y: float
    obs_z: float
    sin_lat: float
    cos_lat: float
    sin_lon: float
    cos_lon: float


def make_observer_projection_state(
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
) -> ObserverProjectionState:
    observer_location = EarthLocation(
        lat=float(observer_lat) * u.deg,
        lon=float(observer_lon) * u.deg,
        height=float(observer_height_m) * u.m,
    )
    observer_xyz = observer_location.to_geocentric()
    lat_rad = math.radians(float(observer_lat))
    lon_rad = math.radians(float(observer_lon))
    return ObserverProjectionState(
        obs_x=float(observer_xyz[0].to_value(u.m)),
        obs_y=float(observer_xyz[1].to_value(u.m)),
        obs_z=float(observer_xyz[2].to_value(u.m)),
        sin_lat=math.sin(lat_rad),
        cos_lat=math.cos(lat_rad),
        sin_lon=math.sin(lon_rad),
        cos_lon=math.cos(lon_rad),
    )


def project_geodetic_to_altaz(
    target_lat: float,
    target_lon: float,
    target_alt_m: float,
    *,
    observer_state: ObserverProjectionState,
) -> tuple[float, float, float]:
    target_location = EarthLocation(
        lat=float(target_lat) * u.deg,
        lon=float(target_lon) * u.deg,
        height=float(target_alt_m) * u.m,
    )
    target_xyz = target_location.to_geocentric()
    dx = float(target_xyz[0].to_value(u.m)) - float(observer_state.obs_x)
    dy = float(target_xyz[1].to_value(u.m)) - float(observer_state.obs_y)
    dz = float(target_xyz[2].to_value(u.m)) - float(observer_state.obs_z)

    east_m = (-float(observer_state.sin_lon) * dx) + (float(observer_state.cos_lon) * dy)
    north_m = (
        (-float(observer_state.sin_lat) * float(observer_state.cos_lon) * dx)
        - (float(observer_state.sin_lat) * float(observer_state.sin_lon) * dy)
        + (float(observer_state.cos_lat) * dz)
    )
    up_m = (
        (float(observer_state.cos_lat) * float(observer_state.cos_lon) * dx)
        + (float(observer_state.cos_lat) * float(observer_state.sin_lon) * dy)
        + (float(observer_state.sin_lat) * dz)
    )

    horizontal_m = math.hypot(east_m, north_m)
    distance_km = math.sqrt((horizontal_m * horizontal_m) + (up_m * up_m)) / 1000.0
    alt_deg = math.degrees(math.atan2(up_m, horizontal_m))
    az_deg = math.degrees(math.atan2(east_m, north_m)) % 360.0
    return alt_deg, az_deg, distance_km
