from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class PlaceTargetProjection:
    alt_deg: float
    az_deg: float
    distance_km: float
    target_latitude_deg: float
    target_longitude_deg: float
    target_height_m: float


WGS84_SEMI_MAJOR_AXIS_M = 6378137.0
WGS84_FLATTENING = 1.0 / 298.257223563
WGS84_ECCENTRICITY_SQUARED = WGS84_FLATTENING * (2.0 - WGS84_FLATTENING)


def _geodetic_to_ecef_numpy(
    latitude_deg: np.ndarray | float,
    longitude_deg: np.ndarray | float,
    height_m: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lat_deg, lon_deg, height = np.broadcast_arrays(
        np.asarray(latitude_deg, dtype=np.float64),
        np.asarray(longitude_deg, dtype=np.float64),
        np.asarray(height_m, dtype=np.float64),
    )
    lat_rad = np.deg2rad(lat_deg)
    lon_rad = np.deg2rad(lon_deg)
    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)
    sin_lon = np.sin(lon_rad)
    cos_lon = np.cos(lon_rad)
    prime_vertical = WGS84_SEMI_MAJOR_AXIS_M / np.sqrt(
        1.0 - (WGS84_ECCENTRICITY_SQUARED * sin_lat * sin_lat)
    )
    x = (prime_vertical + height) * cos_lat * cos_lon
    y = (prime_vertical + height) * cos_lat * sin_lon
    z = (prime_vertical * (1.0 - WGS84_ECCENTRICITY_SQUARED) + height) * sin_lat
    return x, y, z


def _project_place_target_arrays_to_altaz(
    *,
    observer_latitude_deg: float,
    observer_longitude_deg: float,
    observer_height_m: float,
    target_latitude_deg: np.ndarray | float,
    target_longitude_deg: np.ndarray | float,
    target_height_m: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    obs_x, obs_y, obs_z = _geodetic_to_ecef_numpy(
        float(observer_latitude_deg),
        float(observer_longitude_deg),
        float(observer_height_m),
    )
    observer_lat_rad = math.radians(float(observer_latitude_deg))
    observer_lon_rad = math.radians(float(observer_longitude_deg))
    sin_lat = math.sin(observer_lat_rad)
    cos_lat = math.cos(observer_lat_rad)
    sin_lon = math.sin(observer_lon_rad)
    cos_lon = math.cos(observer_lon_rad)

    tgt_x, tgt_y, tgt_z = _geodetic_to_ecef_numpy(
        target_latitude_deg,
        target_longitude_deg,
        target_height_m,
    )
    dx = tgt_x - float(obs_x)
    dy = tgt_y - float(obs_y)
    dz = tgt_z - float(obs_z)

    east_m = (-sin_lon * dx) + (cos_lon * dy)
    north_m = (-sin_lat * cos_lon * dx) - (sin_lat * sin_lon * dy) + (cos_lat * dz)
    up_m = (cos_lat * cos_lon * dx) + (cos_lat * sin_lon * dy) + (sin_lat * dz)

    horizontal_m = np.hypot(east_m, north_m)
    distance_km = np.sqrt((horizontal_m * horizontal_m) + (up_m * up_m)) / 1000.0
    alt_deg = np.degrees(np.arctan2(up_m, horizontal_m))
    az_deg = np.degrees(np.arctan2(east_m, north_m)) % 360.0
    az_deg = np.where(horizontal_m <= 1.0e-9, 0.0, az_deg)
    return alt_deg, az_deg, distance_km


def project_place_target_to_altaz(
    *,
    observer_latitude_deg: float,
    observer_longitude_deg: float,
    observer_height_m: float,
    target_latitude_deg: float,
    target_longitude_deg: float,
    target_height_m: float = 0.0,
) -> PlaceTargetProjection:
    alt_deg, az_deg, distance_km = _project_place_target_arrays_to_altaz(
        observer_latitude_deg=observer_latitude_deg,
        observer_longitude_deg=observer_longitude_deg,
        observer_height_m=observer_height_m,
        target_latitude_deg=float(target_latitude_deg),
        target_longitude_deg=float(target_longitude_deg),
        target_height_m=float(target_height_m),
    )
    return PlaceTargetProjection(
        alt_deg=float(np.asarray(alt_deg)),
        az_deg=float(np.asarray(az_deg)),
        distance_km=float(np.asarray(distance_km)),
        target_latitude_deg=float(target_latitude_deg),
        target_longitude_deg=float(target_longitude_deg),
        target_height_m=float(target_height_m),
    )


def project_place_targets_to_altaz(
    *,
    observer_latitude_deg: float,
    observer_longitude_deg: float,
    observer_height_m: float,
    target_latitude_deg: Sequence[float] | np.ndarray,
    target_longitude_deg: Sequence[float] | np.ndarray,
    target_height_m: Sequence[float] | np.ndarray | float = 0.0,
) -> tuple[PlaceTargetProjection, ...]:
    alt_deg, az_deg, distance_km = _project_place_target_arrays_to_altaz(
        observer_latitude_deg=observer_latitude_deg,
        observer_longitude_deg=observer_longitude_deg,
        observer_height_m=observer_height_m,
        target_latitude_deg=target_latitude_deg,
        target_longitude_deg=target_longitude_deg,
        target_height_m=target_height_m,
    )
    target_lat_arr = np.asarray(target_latitude_deg, dtype=np.float64)
    target_lon_arr = np.asarray(target_longitude_deg, dtype=np.float64)
    target_height_arr = np.asarray(target_height_m, dtype=np.float64)
    target_lat_arr, target_lon_arr, target_height_arr, alt_deg, az_deg, distance_km = np.broadcast_arrays(
        target_lat_arr,
        target_lon_arr,
        target_height_arr,
        np.asarray(alt_deg, dtype=np.float64),
        np.asarray(az_deg, dtype=np.float64),
        np.asarray(distance_km, dtype=np.float64),
    )
    projections: list[PlaceTargetProjection] = []
    for target_lat, target_lon, target_height, alt, az, distance in zip(
        target_lat_arr.flat,
        target_lon_arr.flat,
        target_height_arr.flat,
        alt_deg.flat,
        az_deg.flat,
        distance_km.flat,
        strict=False,
    ):
        projections.append(
            PlaceTargetProjection(
                alt_deg=float(alt),
                az_deg=float(az),
                distance_km=float(distance),
                target_latitude_deg=float(target_lat),
                target_longitude_deg=float(target_lon),
                target_height_m=float(target_height),
            )
        )
    return tuple(projections)


__all__ = [
    "PlaceTargetProjection",
    "project_place_target_to_altaz",
    "project_place_targets_to_altaz",
]
