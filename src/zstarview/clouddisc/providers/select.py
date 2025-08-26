"""
Utilities for selecting the best satellite.
"""

import math
from typing import Tuple, List

# Sub-satellite longitudes for available geostationary satellites.
SAT_LON = {"G16": -75.2, "G18": -137.0, "HIMAWARI": 140.7}


def central_angle_deg(lat_deg: float, lon_deg: float, sub_lon_deg: float) -> float:
    """
    Calculates the central angle between a point on Earth and a satellite's nadir.
    This angle determines if the satellite is visible from the given location.

    Args:
        lat_deg: Latitude of the observer in degrees.
        lon_deg: Longitude of the observer in degrees.
        sub_lon_deg: Sub-satellite longitude in degrees.

    Returns:
        The central angle in degrees.
    """
    lat_rad = math.radians(lat_deg)

    # Normalize the longitude difference to the range [-180, 180]
    lon_diff_rad = math.radians((lon_deg - sub_lon_deg + 540) % 360 - 180)

    # Great-circle distance formula for a geostationary satellite
    cos_angle = math.cos(lat_rad) * math.cos(lon_diff_rad)
    # Clamp the value to handle potential floating point inaccuracies
    cos_angle = max(-1.0, min(1.0, cos_angle))

    return math.degrees(math.acos(cos_angle))


def pick_satellite(lat: float, lon: float, priority: Tuple[str, ...] = ("AUTO",)) -> str:
    """
    Selects the best satellite for a given location based on visibility and priority.

    If priority contains "AUTO", it selects the satellite with the smallest central
    angle that is within the visibility range (<= 81.3 degrees). If no satellite
    is strictly visible, it picks the one with the smallest angle anyway (the closest).

    If priority is an explicit list (e.g., ("HIMAWARI", "G18")), it returns the
    first satellite from that list that is a valid satellite name.

    Args:
        lat: Observer's latitude.
        lon: Observer's longitude.
        priority: A tuple defining the satellite selection strategy.

    Returns:
        The name of the selected satellite.
    """
    if "AUTO" in priority:
        candidates: List[Tuple[float, str]] = []
        for sat, sub_lon in SAT_LON.items():
            angle = central_angle_deg(lat, lon, sub_lon)
            # 81.3 degrees is the approximate limit for visibility of a geostationary satellite
            if angle <= 81.3:
                candidates.append((angle, sat))

        if not candidates:
            # If no satellite is visible, find the one with the minimum angle anyway.
            candidates = [(central_angle_deg(lat, lon, slon), s) for s, slon in SAT_LON.items()]

        # Sort by angle to find the best candidate
        candidates.sort()
        return candidates[0][1]

    # Use the explicit priority list
    for sat_name in priority:
        if sat_name == "AUTO":
            continue
        if sat_name in SAT_LON:
            return sat_name

    # Default fallback if priority list is misconfigured or empty
    return "HIMAWARI"
