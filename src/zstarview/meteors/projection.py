"""Project observed GMN atmospheric paths into fixed celestial coordinates."""

from __future__ import annotations

import math
from collections.abc import Iterable

import astropy.time
import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS, AltAz, EarthLocation, SkyCoord

from .constants import GMN_CANDIDATE_RADIUS_KM
from .types import CelestialMeteorTrail, MeteorObservation


def project_meteor_observations_to_celestial(
    observations: Iterable[MeteorObservation],
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    candidate_radius_km: float = GMN_CANDIDATE_RADIUS_KM,
) -> tuple[CelestialMeteorTrail, ...]:
    observer = EarthLocation(
        lat=float(observer_lat) * u.deg,
        lon=float(observer_lon) * u.deg,
        height=float(observer_height_m) * u.m,
    )
    observer_xyz = _earth_location_xyz_m(observer)
    _, _, up = _enu_basis(
        observation_lat_deg=float(observer_lat),
        observation_lon_deg=float(observer_lon),
    )
    trails: list[CelestialMeteorTrail] = []
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
        clipped = _clip_segment_to_geometric_horizon(
            begin_xyz,
            end_xyz,
            observer_xyz=observer_xyz,
            up=up,
        )
        if clipped is None:
            continue
        begin_visible, end_visible = clipped
        begin_ra, begin_dec = _line_of_sight_to_icrs(
            begin_visible - observer_xyz,
            observer=observer,
            observation=observation,
        )
        end_ra, end_dec = _line_of_sight_to_icrs(
            end_visible - observer_xyz,
            observer=observer,
            observation=observation,
        )
        trails.append(
            CelestialMeteorTrail(
                trajectory_id=observation.trajectory_id,
                beginning_utc=observation.beginning_utc,
                begin_ra_deg=begin_ra,
                begin_dec_deg=begin_dec,
                end_ra_deg=end_ra,
                end_dec_deg=end_dec,
                duration_s=observation.duration_s,
                peak_abs_magnitude=observation.peak_abs_magnitude,
                shower_code=observation.shower_code,
            )
        )
    return tuple(trails)


def _clip_segment_to_geometric_horizon(
    begin_xyz: np.ndarray,
    end_xyz: np.ndarray,
    *,
    observer_xyz: np.ndarray,
    up: np.ndarray,
) -> tuple[np.ndarray, np.ndarray] | None:
    begin_up = float(np.dot(begin_xyz - observer_xyz, up))
    end_up = float(np.dot(end_xyz - observer_xyz, up))
    if begin_up < 0.0 and end_up < 0.0:
        return None
    clipped_begin = begin_xyz
    clipped_end = end_xyz
    if begin_up < 0.0 <= end_up:
        fraction = -begin_up / (end_up - begin_up)
        clipped_begin = begin_xyz + fraction * (end_xyz - begin_xyz)
    elif end_up < 0.0 <= begin_up:
        fraction = -begin_up / (end_up - begin_up)
        clipped_end = begin_xyz + fraction * (end_xyz - begin_xyz)
    return clipped_begin, clipped_end


def _line_of_sight_to_icrs(
    vector_m: np.ndarray,
    *,
    observer: EarthLocation,
    observation: MeteorObservation,
) -> tuple[float, float]:
    east, north, up = _enu_basis(
        observation_lat_deg=float(observer.lat.to_value(u.deg)),
        observation_lon_deg=float(observer.lon.to_value(u.deg)),
    )
    east_m = float(np.dot(vector_m, east))
    north_m = float(np.dot(vector_m, north))
    up_m = float(np.dot(vector_m, up))
    horizontal_m = math.hypot(east_m, north_m)
    alt_deg = math.degrees(math.atan2(up_m, horizontal_m))
    az_deg = math.degrees(math.atan2(east_m, north_m)) % 360.0
    frame = AltAz(
        obstime=astropy.time.Time(observation.beginning_utc),
        location=observer,
    )
    icrs = SkyCoord(az=az_deg * u.deg, alt=alt_deg * u.deg, frame=frame).transform_to(ICRS())
    return float(icrs.ra.deg) % 360.0, float(icrs.dec.deg)


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
