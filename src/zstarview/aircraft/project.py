from __future__ import annotations

import math
from collections.abc import Iterable

import astropy.time

from ..aircraft_constants import (
    AIRCRAFT_FADE_END_SECONDS,
    AIRCRAFT_FADE_START_SECONDS,
    AIRCRAFT_MIN_ALPHA_SCALE,
    AIRCRAFT_TRAIL_HALF_SPAN_SECONDS,
)
from .projection import (
    ObserverProjectionState,
    make_observer_projection_state,
    project_geodetic_to_altaz,
)
from .types import AircraftOverlayPoint, AircraftSnapshot


def project_aircraft_snapshots(
    snapshots: Iterable[AircraftSnapshot],
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    time_obj: astropy.time.Time,
) -> list[AircraftOverlayPoint]:
    now_unix = float(time_obj.unix)
    observer_state = make_observer_projection_state(
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        observer_height_m=observer_height_m,
    )

    overlay_points: list[AircraftOverlayPoint] = []
    for snapshot in snapshots:
        predicted_lat, predicted_lon, predicted_alt_m, age_seconds = _predict_snapshot_geodetic(snapshot, now_unix=now_unix)
        alt_deg, az_deg, distance_km = project_geodetic_to_altaz(
            predicted_lat,
            predicted_lon,
            predicted_alt_m,
            observer_state=observer_state,
        )
        if alt_deg <= 0.0:
            continue
        trail_alt_az_points, trail_geodetic_points = _project_trail_points(
            snapshot,
            age_seconds=age_seconds,
            observer_state=observer_state,
        )
        overlay_points.append(
            AircraftOverlayPoint(
                icao24=snapshot.icao24,
                callsign=snapshot.callsign,
                alt_deg=alt_deg,
                az_deg=az_deg,
                trail_alt_az_points=trail_alt_az_points,
                distance_km=distance_km,
                age_seconds=age_seconds,
                alpha_scale=_aircraft_alpha_scale(age_seconds),
                trail_geodetic_points=trail_geodetic_points,
            )
        )
    return overlay_points


def _project_trail_points(
    snapshot: AircraftSnapshot,
    *,
    age_seconds: float,
    observer_state: ObserverProjectionState,
) -> tuple[tuple[tuple[float, float], ...], tuple[tuple[float, float, float], ...]]:
    half_span = AIRCRAFT_TRAIL_HALF_SPAN_SECONDS
    sample_offsets = (-half_span, -(half_span / 2.0), 0.0, half_span / 2.0, half_span)
    alt_az_points: list[tuple[float, float]] = []
    geodetic_points: list[tuple[float, float, float]] = []
    for offset_seconds in sample_offsets:
        sample_lat, sample_lon, sample_alt_m = _predict_snapshot_at_age(
            snapshot,
            age_seconds=max(0.0, age_seconds + offset_seconds),
        )
        sample_alt_deg, sample_az_deg, _ = project_geodetic_to_altaz(
            sample_lat,
            sample_lon,
            sample_alt_m,
            observer_state=observer_state,
        )
        alt_az_points.append((sample_alt_deg, sample_az_deg))
        geodetic_points.append((float(sample_lat), float(sample_lon), float(sample_alt_m)))
    return tuple(alt_az_points), tuple(geodetic_points)


def _predict_snapshot_geodetic(
    snapshot: AircraftSnapshot,
    *,
    now_unix: float,
) -> tuple[float, float, float, float]:
    age_seconds = max(0.0, now_unix - float(snapshot.last_contact_unix or now_unix))
    predicted_lat, predicted_lon, predicted_alt_m = _predict_snapshot_at_age(snapshot, age_seconds=age_seconds)
    return predicted_lat, predicted_lon, predicted_alt_m, age_seconds


def _predict_snapshot_at_age(
    snapshot: AircraftSnapshot,
    *,
    age_seconds: float,
) -> tuple[float, float, float]:
    predicted_lat = float(snapshot.latitude)
    predicted_lon = float(snapshot.longitude)
    predicted_alt_m = float(snapshot.baro_altitude_m or 0.0)

    velocity_mps = float(snapshot.velocity_mps or 0.0)
    heading_deg = float(snapshot.heading_deg or 0.0)
    vertical_rate_mps = float(snapshot.vertical_rate_mps or 0.0)

    if velocity_mps > 0.0:
        heading_rad = math.radians(heading_deg)
        north_m = math.cos(heading_rad) * velocity_mps * age_seconds
        east_m = math.sin(heading_rad) * velocity_mps * age_seconds
        lat_rad = math.radians(predicted_lat)
        predicted_lat += north_m / 111320.0
        cos_lat = max(0.01, math.cos(lat_rad))
        predicted_lon += east_m / (111320.0 * cos_lat)

    if vertical_rate_mps != 0.0:
        predicted_alt_m += vertical_rate_mps * age_seconds
        predicted_alt_m = max(0.0, predicted_alt_m)

    return predicted_lat, predicted_lon, predicted_alt_m


def _aircraft_alpha_scale(age_seconds: float) -> float:
    age = max(0.0, float(age_seconds))
    if age <= AIRCRAFT_FADE_START_SECONDS:
        return 1.0
    if age >= AIRCRAFT_FADE_END_SECONDS:
        return AIRCRAFT_MIN_ALPHA_SCALE
    span = AIRCRAFT_FADE_END_SECONDS - AIRCRAFT_FADE_START_SECONDS
    progress = (age - AIRCRAFT_FADE_START_SECONDS) / span
    return 1.0 + (AIRCRAFT_MIN_ALPHA_SCALE - 1.0) * progress
