# -*- coding: utf-8 -*-
"""
Utilities for selecting the best geostationary satellite for a given location.
"""

import math
from typing import Dict, Tuple, List

# Defines the sub-satellite longitude (the longitude directly below the satellite)
# for each available geostationary satellite.
LEGACY_SATELLITE_ALIASES: Dict[str, str] = {
    "G16": "G19",
}
GOES_SATELLITES: Tuple[str, ...] = ("G19", "G18")
SUPPORTED_SATELLITES: Tuple[str, ...] = GOES_SATELLITES + ("HIMAWARI",)
SAT_LON: Dict[str, float] = {
    "G19": -75.2,
    "G18": -137.0,
    "HIMAWARI": 140.7,
}
MAX_VISIBLE_CENTRAL_ANGLE_DEG = 81.3


def _sat_lon_map() -> Dict[str, float]:
    return SAT_LON


def normalize_satellite_name(sat_name: str) -> str:
    return LEGACY_SATELLITE_ALIASES.get(sat_name, sat_name)


def central_angle_deg(lat_deg: float, lon_deg: float, sub_lon_deg: float) -> float:
    """
    Calculates the central angle between an observer and a satellite's nadir.

    This angle is the great-circle distance on a spherical Earth. It is used to
    determine if the satellite is above the horizon for the observer.

    Args:
        lat_deg: Latitude of the observer in degrees.
        lon_deg: Longitude of the observer in degrees.
        sub_lon_deg: Sub-satellite longitude in degrees.

    Returns:
        The central angle in degrees.
    """
    lat_rad = math.radians(lat_deg)

    # Normalize the longitude difference to the range [-180, 180].
    lon_diff_rad = math.radians((lon_deg - sub_lon_deg + 540) % 360 - 180)

    # Use the spherical law of cosines for the great-circle distance.
    # For a geostationary satellite, the sub-satellite latitude is 0.
    cos_angle = math.cos(lat_rad) * math.cos(lon_diff_rad)

    # Clamp the value to handle potential floating-point inaccuracies near the poles.
    cos_angle = max(-1.0, min(1.0, cos_angle))

    return math.degrees(math.acos(cos_angle))


def is_satellite_visible(
    lat: float,
    lon: float,
    sat_name: str,
    max_angle_deg: float = MAX_VISIBLE_CENTRAL_ANGLE_DEG,
) -> bool:
    """Return True if a geostationary satellite is above the visibility angle."""
    sat_lon_map = _sat_lon_map()
    sat_name = normalize_satellite_name(sat_name)
    if sat_name not in sat_lon_map:
        return False
    return central_angle_deg(lat, lon, sat_lon_map[sat_name]) <= max_angle_deg


def visible_satellites(
    lat: float,
    lon: float,
    sat_names: Tuple[str, ...],
    max_angle_deg: float = MAX_VISIBLE_CENTRAL_ANGLE_DEG,
) -> List[str]:
    """Return visible satellites from sat_names ordered by smaller central angle first."""
    sat_lon_map = _sat_lon_map()
    normalized_names: List[str] = []
    for sat in sat_names:
        normalized = normalize_satellite_name(sat)
        if normalized not in normalized_names:
            normalized_names.append(normalized)
    visible: List[Tuple[float, str]] = []
    for sat in normalized_names:
        if sat not in sat_lon_map:
            continue
        angle = central_angle_deg(lat, lon, sat_lon_map[sat])
        if angle <= max_angle_deg:
            visible.append((angle, sat))
    visible.sort()
    return [sat for _, sat in visible]


def pick_satellite(
    lat: float,
    lon: float,
    priority: Tuple[str, ...] = ("AUTO",),
) -> str:
    """
    Selects the best satellite for a given location based on visibility and priority.

    If "AUTO" is in the priority list, it selects the satellite with the smallest
    central angle (i.e., highest in the sky) that is within the visibility range.
    The visibility limit for a geostationary satellite is approximately 81.3 degrees.

    If priority is an explicit list (e.g., ("HIMAWARI", "G19")), it returns the
    first satellite from that list that is a valid satellite name.

    Args:
        lat: Observer's latitude in degrees.
        lon: Observer's longitude in degrees.
        priority: A tuple defining the satellite selection strategy.

    Returns:
        The name of the selected satellite (e.g., "HIMAWARI").
    """
    # --- Automatic Selection Mode ---
    if "AUTO" in priority:
        sat_lon_map = _sat_lon_map()
        candidates: List[Tuple[float, str]] = []
        for sat, sub_lon in sat_lon_map.items():
            angle = central_angle_deg(lat, lon, sub_lon)
            if angle <= MAX_VISIBLE_CENTRAL_ANGLE_DEG:
                candidates.append((angle, sat))

        if not candidates:
            # If no satellite is technically visible, find the one that is closest.
            candidates = [(central_angle_deg(lat, lon, slon), s) for s, slon in sat_lon_map.items()]

        # Sort by angle (smallest first) to find the best candidate.
        candidates.sort()
        return candidates[0][1]

    # --- Manual Priority Mode ---
    sat_lon_map = _sat_lon_map()
    # Iterate through the user-provided list and return the first valid satellite.
    for sat_name in priority:
        sat_name = normalize_satellite_name(sat_name)
        if sat_name in sat_lon_map:
            return sat_name

    # As a final fallback, if the priority list is misconfigured or empty.
    return "HIMAWARI"
