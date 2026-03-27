from __future__ import annotations

import math
from dataclasses import dataclass

from astropy import units as u
from astropy.coordinates import EarthLocation


@dataclass(frozen=True)
class PlaceTargetProjection:
    alt_deg: float
    az_deg: float
    distance_km: float
    target_latitude_deg: float
    target_longitude_deg: float
    target_height_m: float


def project_place_target_to_altaz(
    *,
    observer_latitude_deg: float,
    observer_longitude_deg: float,
    observer_height_m: float,
    target_latitude_deg: float,
    target_longitude_deg: float,
    target_height_m: float = 0.0,
) -> PlaceTargetProjection:
    observer_location = EarthLocation(
        lat=float(observer_latitude_deg) * u.deg,
        lon=float(observer_longitude_deg) * u.deg,
        height=float(observer_height_m) * u.m,
    )
    observer_xyz = observer_location.to_geocentric()
    obs_x = float(observer_xyz[0].to_value(u.m))
    obs_y = float(observer_xyz[1].to_value(u.m))
    obs_z = float(observer_xyz[2].to_value(u.m))

    observer_lat_rad = math.radians(float(observer_latitude_deg))
    observer_lon_rad = math.radians(float(observer_longitude_deg))
    sin_lat = math.sin(observer_lat_rad)
    cos_lat = math.cos(observer_lat_rad)
    sin_lon = math.sin(observer_lon_rad)
    cos_lon = math.cos(observer_lon_rad)

    target_location = EarthLocation(
        lat=float(target_latitude_deg) * u.deg,
        lon=float(target_longitude_deg) * u.deg,
        height=float(target_height_m) * u.m,
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
    if horizontal_m <= 1e-9:
        az_deg = 0.0
    else:
        az_deg = math.degrees(math.atan2(east_m, north_m)) % 360.0
    return PlaceTargetProjection(
        alt_deg=alt_deg,
        az_deg=az_deg,
        distance_km=distance_km,
        target_latitude_deg=float(target_latitude_deg),
        target_longitude_deg=float(target_longitude_deg),
        target_height_m=float(target_height_m),
    )


__all__ = ["PlaceTargetProjection", "project_place_target_to_altaz"]
