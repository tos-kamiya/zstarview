# -*- coding: utf-8 -*-
"""
Utilities for selecting the best geostationary satellite for a given location.
"""

import math
from typing import Tuple, List

# Defines the sub-satellite longitude (the longitude directly below the satellite)
# for each available geostationary satellite.
SAT_LON = {"G16": -75.2, "G18": -137.0, "HIMAWARI": 140.7}


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


def pick_satellite(lat: float, lon: float, priority: Tuple[str, ...] = ("AUTO",)) -> str:
    """
    Selects the best satellite for a given location based on visibility and priority.

    If "AUTO" is in the priority list, it selects the satellite with the smallest
    central angle (i.e., highest in the sky) that is within the visibility range.
    The visibility limit for a geostationary satellite is approximately 81.3 degrees.

    If priority is an explicit list (e.g., ("HIMAWARI", "G18")), it returns the
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
        candidates: List[Tuple[float, str]] = []
        for sat, sub_lon in SAT_LON.items():
            angle = central_angle_deg(lat, lon, sub_lon)
            # 81.3 degrees is the approximate maximum central angle for a geostationary
            # satellite to be visible above the horizon.
            if angle <= 81.3:
                candidates.append((angle, sat))

        if not candidates:
            # If no satellite is technically visible, find the one that is closest.
            candidates = [(central_angle_deg(lat, lon, slon), s) for s, slon in SAT_LON.items()]

        # Sort by angle (smallest first) to find the best candidate.
        candidates.sort()
        return candidates[0][1]

    # --- Manual Priority Mode ---
    # Iterate through the user-provided list and return the first valid satellite.
    for sat_name in priority:
        if sat_name in SAT_LON:
            return sat_name

    # As a final fallback, if the priority list is misconfigured or empty.
    return "HIMAWARI"
