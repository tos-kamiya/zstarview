from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class AircraftSnapshot:
    icao24: str
    callsign: Optional[str]
    latitude: float
    longitude: float
    baro_altitude_m: Optional[float]
    velocity_mps: Optional[float]
    heading_deg: Optional[float]
    vertical_rate_mps: Optional[float]
    on_ground: bool
    last_contact_unix: Optional[int]


@dataclass(frozen=True)
class AircraftOverlayPoint:
    icao24: str
    callsign: Optional[str]
    alt_deg: float
    az_deg: float
    trail_alt_az_points: tuple[tuple[float, float], ...]
    distance_km: float
    age_seconds: float
    alpha_scale: float
