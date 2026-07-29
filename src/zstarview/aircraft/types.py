from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AircraftSnapshot:
    icao24: str
    callsign: str | None
    latitude: float
    longitude: float
    baro_altitude_m: float | None
    velocity_mps: float | None
    heading_deg: float | None
    vertical_rate_mps: float | None
    on_ground: bool
    last_contact_unix: int | None


@dataclass(frozen=True)
class AircraftOverlayPoint:
    icao24: str
    callsign: str | None
    alt_deg: float
    az_deg: float
    trail_alt_az_points: tuple[tuple[float, float], ...]
    distance_km: float
    age_seconds: float
    alpha_scale: float
    trail_geodetic_points: tuple[tuple[float, float, float], ...] = ()
