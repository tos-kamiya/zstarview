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


def moon_cloud_strength(alt_deg: float, elongation_deg: float) -> float:
    """Display strength; elongation is 0 at new moon and 180 at full moon."""
    if not math.isfinite(alt_deg) or not math.isfinite(elongation_deg):
        return 0.0
    if alt_deg <= 0.0:
        return 0.0
    phase = (1.0 - math.cos(math.radians(elongation_deg))) * 0.5
    horizon = min(alt_deg / 5.0, 1.0)
    horizon = horizon * horizon * (3.0 - 2.0 * horizon)
    altitude = math.sqrt(max(0.0, math.sin(math.radians(min(alt_deg, 90.0)))))
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
