from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Iterable, Optional
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from ..aircraft_constants import AIRCRAFT_BBOX_DELTA_DEG
from .types import AircraftSnapshot

OPENSKY_STATES_ALL_URL = "https://opensky-network.org/api/states/all"
MIN_AIRBORNE_SPEED_MPS = 30.0


@dataclass(frozen=True)
class AircraftBoundingBox:
    min_lat: float
    max_lat: float
    min_lon: float
    max_lon: float

    @property
    def area_sq_deg(self) -> float:
        return max(0.0, self.max_lat - self.min_lat) * max(0.0, self.max_lon - self.min_lon)

    def to_query_params(self) -> dict[str, str]:
        return {
            "lamin": f"{self.min_lat:.6f}",
            "lamax": f"{self.max_lat:.6f}",
            "lomin": f"{self.min_lon:.6f}",
            "lomax": f"{self.max_lon:.6f}",
        }

    def to_url(self, base_url: str = OPENSKY_STATES_ALL_URL) -> str:
        return f"{base_url}?{urlencode(self.to_query_params())}"


def build_observer_bbox(
    observer_lat: float,
    observer_lon: float,
    *,
    delta_deg: float = AIRCRAFT_BBOX_DELTA_DEG,
) -> AircraftBoundingBox:
    delta = max(0.0, float(delta_deg))
    lat = float(observer_lat)
    lon = float(observer_lon)
    min_lat = max(-90.0, lat - delta)
    max_lat = min(90.0, lat + delta)
    min_lon = max(-180.0, lon - delta)
    max_lon = min(180.0, lon + delta)
    return AircraftBoundingBox(
        min_lat=min_lat,
        max_lat=max_lat,
        min_lon=min_lon,
        max_lon=max_lon,
    )


def normalize_opensky_state_vectors(
    payload: dict[str, Any],
    *,
    min_airborne_speed_mps: float = MIN_AIRBORNE_SPEED_MPS,
) -> list[AircraftSnapshot]:
    states = payload.get("states")
    if not isinstance(states, list):
        return []

    normalized: list[AircraftSnapshot] = []
    for row in states:
        snapshot = _parse_state_vector(row, min_airborne_speed_mps=min_airborne_speed_mps)
        if snapshot is not None:
            normalized.append(snapshot)
    return normalized


def fetch_opensky_states(
    bbox: AircraftBoundingBox,
    *,
    auth_header: Optional[str] = None,
    timeout_s: float = 20.0,
    base_url: str = OPENSKY_STATES_ALL_URL,
) -> list[AircraftSnapshot]:
    headers = {"Accept": "application/json"}
    if auth_header:
        headers["Authorization"] = auth_header
    request = Request(bbox.to_url(base_url), headers=headers)
    with urlopen(request, timeout=float(timeout_s)) as response:
        payload = json.load(response)
    if not isinstance(payload, dict):
        return []
    return normalize_opensky_state_vectors(payload)


def _parse_state_vector(
    row: Any,
    *,
    min_airborne_speed_mps: float,
) -> Optional[AircraftSnapshot]:
    if not isinstance(row, Iterable) or isinstance(row, (str, bytes, dict)):
        return None
    values = list(row)
    if len(values) < 12:
        return None

    icao24 = _clean_text(values[0])
    latitude = _to_float(values[6])
    longitude = _to_float(values[5])
    if not icao24 or latitude is None or longitude is None:
        return None

    on_ground = bool(values[8])
    velocity_mps = _to_float(values[9])
    if on_ground:
        return None
    if velocity_mps is not None and velocity_mps < float(min_airborne_speed_mps):
        return None

    last_contact = _to_int(values[4])
    return AircraftSnapshot(
        icao24=icao24,
        callsign=_clean_text(values[1]),
        latitude=latitude,
        longitude=longitude,
        baro_altitude_m=_to_float(values[7]),
        velocity_mps=velocity_mps,
        heading_deg=_to_float(values[10]),
        vertical_rate_mps=_to_float(values[11]),
        on_ground=on_ground,
        last_contact_unix=last_contact,
    )


def _clean_text(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _to_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None
