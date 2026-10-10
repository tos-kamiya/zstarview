"""Display approximations for lunar and ground cloud illumination."""

from __future__ import annotations

import math
from collections.abc import Iterable

NIGHT_AMBIENT_STRENGTH = 0.08
MOON_LIGHT_STRENGTH = 0.25
MOON_LIGHT_RGB = (0.96, 0.975, 1.0)
GROUND_LIGHT_STRENGTH = 0.25
GROUND_LIGHT_RGB = (1.0, 0.90, 0.72)
GROUND_RADIANCE_SCALE = 20.0
NIGHT_EXPOSURE = 2.0
LIGHT_EARTH_RADIUS_KM = 6371.0
DIRECT_LIGHT_HORIZON_FADE_DEG = 0.25


def cloud_horizon_depression_deg(height_km: float) -> float:
    """Spherical, unrefracted horizon depression at a cloud cell's height."""
    if not math.isfinite(height_km) or height_km < 0.0:
        raise ValueError("cloud height must be finite and nonnegative")
    return math.degrees(math.acos(LIGHT_EARTH_RADIUS_KM / (LIGHT_EARTH_RADIUS_KM + height_km)))


def cloud_light_visibility(alt_deg: float, height_km: float) -> float:
    """Fade direct light just above the height-dependent geometric horizon."""
    if not math.isfinite(alt_deg):
        return 0.0
    clearance = alt_deg + cloud_horizon_depression_deg(height_km)
    fraction = max(0.0, min(clearance / DIRECT_LIGHT_HORIZON_FADE_DEG, 1.0))
    return fraction * fraction * (3.0 - 2.0 * fraction)


def moon_cloud_strength(
    alt_deg: float, elongation_deg: float, height_km: float = 0.0
) -> float:
    """Display strength above the cloud horizon; 0=new, 180=full elongation."""
    if not math.isfinite(alt_deg) or not math.isfinite(elongation_deg):
        return 0.0
    effective_alt = alt_deg + cloud_horizon_depression_deg(height_km)
    if effective_alt <= 0.0:
        return 0.0
    phase = (1.0 - math.cos(math.radians(elongation_deg))) * 0.5
    horizon = min(effective_alt / 5.0, 1.0)
    horizon = horizon * horizon * (3.0 - 2.0 * horizon)
    altitude = math.sqrt(max(0.0, math.sin(math.radians(min(effective_alt, 90.0)))))
    return MOON_LIGHT_STRENGTH * phase * phase * horizon * altitude


def cloud_moon_input(bodies: Iterable) -> tuple[float, float, float] | None:
    """Extract altitude, azimuth and solar elongation, including offscreen Moon."""
    for body in bodies:
        if body.name == "moon" and body.phase_angle is not None:
            return float(body.alt), float(body.az), float(body.phase_angle)
    return None


def cloud_illumination_key(bodies: Iterable) -> tuple:
    """Quantize light changes to 0.1 degree for asynchronous cloud refresh."""
    bodies = tuple(bodies)
    sun = (-90.0, 0.0)
    for body in bodies:
        if body.name == "sun":
            sun = (float(body.alt), float(body.az))
    return cloud_light_input_key(sun, cloud_moon_input(bodies))


def cloud_light_input_key(sun_altaz, moon_altaz_phase) -> tuple:
    """Identify the display lighting of a completed cloud request."""
    sun_altaz = sun_altaz or (-90.0, 0.0)
    return (
        round(float(sun_altaz[0]), 1), round(float(sun_altaz[1]) % 360, 1),
        None if moon_altaz_phase is None else (
            round(float(moon_altaz_phase[0]), 1),
            round(float(moon_altaz_phase[1]) % 360, 1),
            round(float(moon_altaz_phase[2]), 1),
        ),
    )
