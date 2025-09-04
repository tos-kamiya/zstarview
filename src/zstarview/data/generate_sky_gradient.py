#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates a chart of the sky color gradient for different viewing angles.

This script produces an image showing how the sky color model changes based on
the viewing altitude (from -90 to +90 degrees) and azimuth relative to the
sun (directly towards, 90 degrees from, and opposite).
"""

import math
import argparse
from typing import Tuple

import numpy as np
from PIL import Image

# --- Script Parameters ---
IMG_WIDTH = 300
IMG_HEIGHT = 181
SUN_ALTITUDE_DEG = 60.0  # Sun's altitude (-5 for sunset, 90 for zenith)
EXPOSURE = 1.0  # Overall brightness adjustment
SATURATION = 1.2  # Color saturation
GROUND_COLOR = np.array([0.1, 0.1, 0.1])  # Color below the horizon
OUTPUT_FILENAME = "sky_gradient_chart.png"


# --- Helper Functions ---
def _deg2rad(d: float) -> float:
    """Converts degrees to radians."""
    return d * math.pi / 180.0


def _clamp01(x: float) -> float:
    """Clamps a value to the [0, 1] range."""
    return max(0.0, min(1.0, x))


def _angle_between(alt1_deg: float, az1_deg: float, alt2_deg: float, az2_deg: float) -> float:
    """Calculates the angular separation between two points in alt-az coordinates."""
    a1, z1 = _deg2rad(alt1_deg), _deg2rad(az1_deg)
    a2, z2 = _deg2rad(alt2_deg), _deg2rad(az2_deg)
    cos_g = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(z2 - z1)
    cos_g = max(-1.0, min(1.0, cos_g))  # Clamp to handle potential floating point errors
    return math.acos(cos_g)


def _lerp_color(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    """Linearly interpolates between two colors represented as NumPy arrays."""
    t = _clamp01(t)
    return a + (b - a) * t


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """Performs a smooth Hermite interpolation between 0 and 1."""
    t = _clamp01((x - edge0) / (edge1 - edge0))
    return t * t * (3.0 - 2.0 * t)


# --- Phenomenological Sky and Sun Color Model ---


def get_sun_color(sun_alt_deg: float) -> np.ndarray:
    """
    Determines the color of sunlight based on the sun's altitude.

    Args:
        sun_alt_deg: The sun's altitude in degrees (-90 to +90).

    Returns:
        A NumPy array representing the (r, g, b) color in the range [0, 1].
    """
    # Define key colors for different sun altitudes.
    zenith_color = np.array([0.3, 0.55, 0.98])  # Deep blue for when the sun is at zenith.
    horizon_color = np.array([0.95, 0.50, 0.30])  # Orange for when the sun is at the horizon.
    night_color = np.array([0.01, 0.02, 0.05])  # Dark blue for night.

    # --- Calculate daytime color ---
    # Normalize sun altitude from -7° (below horizon) to 90° (zenith) to a [0,1] range.
    t = float(np.clip((sun_alt_deg + 7.0) / 97.0, 0.0, 1.0))
    # Apply a power function to create a more natural, non-linear transition.
    t = math.pow(t, 0.35)
    day_color = _lerp_color(horizon_color, zenith_color, t)

    # --- Fade to night color ---
    # Create a rapid fade to darkness as the sun sets (from +1° to -5°).
    fade = float(np.clip((-sun_alt_deg + 1.0) / 6.0, 0.0, 1.0))
    return _lerp_color(day_color, night_color, fade)


def get_sky_color(view_alt_deg: float, view_az_deg: float, sun_alt_deg: float, sun_az_deg: float) -> np.ndarray:
    """
    Step 2: Determine the final sky color based on viewing angle and sun position.
    """
    # 1. Get the base sunlight color.
    sun_color = get_sun_color(sun_alt_deg)

    # For deep night, the sky is simply black.
    if sun_alt_deg < -10:
        return np.array([0.0, 0.0, 0.0])

    # 2. Determine brightness based on the angle between the view direction and the sun.
    sun_angle = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)  # [0..pi]
    cosg = math.cos(sun_angle)
    brightness = (1.0 + cosg) * 0.5  # 0..1
    brightness = brightness**2.0  # Emphasize the sun-facing direction

    # 3. Adjust tone based on viewing altitude (darker at zenith, whitish at horizon).
    t = _clamp01(view_alt_deg / 90.0)  # 0 (horizon) -> 1 (zenith)
    zenith_dim = 1.0 - 0.35 * t  # Becomes darker (0.65x) towards the zenith.
    horizon_mix = (1.0 - t) * 0.2  # Becomes whiter (20%) towards the horizon.

    # 4. Apply a twilight factor to smoothly transition from day to night.
    twilight = _smoothstep(-10.0, 0.0, sun_alt_deg) if sun_alt_deg < 0.0 else 1.0

    # 5. Composite the final color.
    base_sky = sun_color * brightness
    horizon_color = np.array([1.0, 1.0, 1.0])
    rgb = _lerp_color(base_sky, horizon_color, horizon_mix * twilight)
    rgb *= zenith_dim

    return np.clip(rgb, 0.0, 1.0)


# --- Main Execution ---
def main() -> None:
    """Generates and saves the sky gradient image."""
    print("Generating sky gradient image...")
    print(f"Settings: Sun Altitude = {SUN_ALTITUDE_DEG}°")

    img = Image.new("RGB", (IMG_WIDTH, IMG_HEIGHT))
    sun_az = 0.0
    view_azimuths = {"sun": sun_az, "side": 90.0, "opposite": 180.0}

    for y in range(IMG_HEIGHT):
        view_alt = 90.0 - y  # Top of image is zenith (90), bottom is nadir (-90)

        for i, (name, view_az) in enumerate(view_azimuths.items()):
            # Get color from the phenomenological model.
            rgb = get_sky_color(view_alt, view_az, SUN_ALTITUDE_DEG, sun_az)

            # Apply post-processing for saturation and exposure.
            gray = np.dot(rgb, [0.299, 0.587, 0.114])
            saturated_rgb = _lerp_color(np.array([gray, gray, gray]), rgb, SATURATION)
            final_rgb = np.clip(saturated_rgb * EXPOSURE, 0.0, 1.0)

            # Fade to ground color below the horizon.
            if view_alt < 0:
                t = _clamp01((view_alt + 5.0) / 5.0)  # Fade between -5 and 0 degrees
                final_rgb = _lerp_color(GROUND_COLOR, final_rgb, t)

            # Convert to integer color for PIL.
            color_int = tuple((final_rgb * 255).astype(int))

            # Draw a vertical strip for this azimuth.
            strip_width = IMG_WIDTH // 3
            for x_offset in range(strip_width):
                x = i * strip_width + x_offset
                img.putpixel((x, y), color_int)

    img.save(OUTPUT_FILENAME)
    print(f"Image saved as '{OUTPUT_FILENAME}'")


if __name__ == "__main__":
    main()
