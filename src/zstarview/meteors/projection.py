"""Project observed GMN atmospheric paths into event-time Alt/Az coordinates."""

from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation

from .constants import GMN_CANDIDATE_RADIUS_KM
from .types import MeteorObservation, MeteorTrail


def project_meteor_observations_to_altaz(
    observations: Iterable[MeteorObservation],
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    candidate_radius_km: float = GMN_CANDIDATE_RADIUS_KM,
) -> tuple[MeteorTrail, ...]:
    observer = EarthLocation(
        lat=float(observer_lat) * u.deg,
        lon=float(observer_lon) * u.deg,
        height=float(observer_height_m) * u.m,
    )
    observer_xyz = _earth_location_xyz_m(observer)
    trails: list[MeteorTrail] = []
    for observation in observations:
        if not _within_candidate_radius(
            observation,
            observer_lat=float(observer_lat),
            observer_lon=float(observer_lon),
            radius_km=float(candidate_radius_km),
        ):
            continue
        begin_xyz = _geodetic_xyz_m(
            observation.begin_lat_deg,
            observation.begin_lon_deg,
            observation.begin_height_km,
        )
        end_xyz = _geodetic_xyz_m(
            observation.end_lat_deg,
            observation.end_lon_deg,
            observation.end_height_km,
        )
        begin_alt, begin_az = _line_of_sight_to_altaz(
            begin_xyz - observer_xyz,
            observer_lat_deg=float(observer_lat),
            observer_lon_deg=float(observer_lon),
        )
        end_alt, end_az = _line_of_sight_to_altaz(
            end_xyz - observer_xyz,
            observer_lat_deg=float(observer_lat),
            observer_lon_deg=float(observer_lon),
        )
        trails.append(
            MeteorTrail(
                trajectory_id=observation.trajectory_id,
                beginning_utc=observation.beginning_utc,
                begin_alt_deg=begin_alt,
                begin_az_deg=begin_az,
                end_alt_deg=end_alt,
                end_az_deg=end_az,
                duration_s=observation.duration_s,
                peak_abs_magnitude=observation.peak_abs_magnitude,
                shower_code=observation.shower_code,
            )
        )
    return tuple(trails)


def project_meteor_observations_to_celestial(
    observations: Iterable[MeteorObservation],
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    candidate_radius_km: float = GMN_CANDIDATE_RADIUS_KM,
) -> tuple[MeteorTrail, ...]:
    """Compatibility wrapper for the former projection function name."""
    return project_meteor_observations_to_altaz(
        observations,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        observer_height_m=observer_height_m,
        candidate_radius_km=candidate_radius_km,
    )


def _line_of_sight_to_altaz(
    vector_m: np.ndarray,
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> tuple[float, float]:
    east, north, up = _enu_basis(
        observation_lat_deg=observer_lat_deg,
        observation_lon_deg=observer_lon_deg,
    )
    east_m = float(np.dot(vector_m, east))
    north_m = float(np.dot(vector_m, north))
    up_m = float(np.dot(vector_m, up))
    horizontal_m = math.hypot(east_m, north_m)
    alt_deg = math.degrees(math.atan2(up_m, horizontal_m))
    az_deg = math.degrees(math.atan2(east_m, north_m)) % 360.0
    return alt_deg, az_deg


def _within_candidate_radius(
    observation: MeteorObservation,
    *,
    observer_lat: float,
    observer_lon: float,
    radius_km: float,
) -> bool:
    midpoint = _spherical_midpoint(
        observation.begin_lat_deg,
        observation.begin_lon_deg,
        observation.end_lat_deg,
        observation.end_lon_deg,
    )
    locations = (
        (observation.begin_lat_deg, observation.begin_lon_deg),
        (observation.end_lat_deg, observation.end_lon_deg),
        midpoint,
    )
    return any(
        _great_circle_distance_km(observer_lat, observer_lon, lat, lon) <= radius_km
        for lat, lon in locations
    )


def _great_circle_distance_km(
    lat1_deg: float,
    lon1_deg: float,
    lat2_deg: float,
    lon2_deg: float,
) -> float:
    lat1 = math.radians(lat1_deg)
    lat2 = math.radians(lat2_deg)
    dlat = lat2 - lat1
    dlon = math.radians(lon2_deg - lon1_deg)
    a = math.sin(dlat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * math.asin(min(1.0, math.sqrt(a)))


def _spherical_midpoint(
    begin_lat_deg: float,
    begin_lon_deg: float,
    end_lat_deg: float,
    end_lon_deg: float,
) -> tuple[float, float]:
    first = _unit_sphere_xyz(begin_lat_deg, begin_lon_deg)
    second = _unit_sphere_xyz(end_lat_deg, end_lon_deg)
    midpoint = first + second
    norm = float(np.linalg.norm(midpoint))
    if norm <= 1e-12:
        return begin_lat_deg, begin_lon_deg
    midpoint /= norm
    return (
        math.degrees(math.asin(float(midpoint[2]))),
        math.degrees(math.atan2(float(midpoint[1]), float(midpoint[0]))),
    )


def _unit_sphere_xyz(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    return np.asarray(
        [math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)],
        dtype=float,
    )


def _geodetic_xyz_m(lat_deg: float, lon_deg: float, height_km: float) -> np.ndarray:
    return _earth_location_xyz_m(
        EarthLocation(
            lat=float(lat_deg) * u.deg,
            lon=float(lon_deg) * u.deg,
            height=float(height_km) * u.km,
        )
    )


def _earth_location_xyz_m(location: EarthLocation) -> np.ndarray:
    xyz = location.to_geocentric()
    return np.asarray([component.to_value(u.m) for component in xyz], dtype=float)


def _enu_basis(
    *,
    observation_lat_deg: float,
    observation_lon_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lat = math.radians(observation_lat_deg)
    lon = math.radians(observation_lon_deg)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)
    east = np.asarray([-sin_lon, cos_lon, 0.0])
    north = np.asarray([-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat])
    up = np.asarray([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])
    return east, north, up
