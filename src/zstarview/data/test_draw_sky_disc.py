#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Minimal tester for the sky disc rendering logic.

This script renders a sky-color disc to a QImage and saves it as a PNG file.
It is useful for visually debugging the projection and color models without
running the full application.

Example:
    python test_draw_sky_disc.py --width 256 --height 256 --fov 200 \
      --sun-alt 15 --sun-az 270 --center-alt 90 --center-az 0 \
      --outfile sky_test.png
"""

import argparse
import math
from dataclasses import dataclass
from typing import Callable, Optional, Tuple

import numpy as np
from PySide6.QtCore import QRect
from PySide6.QtGui import QColor, QImage, QPainter, QPixmap
from PySide6.QtWidgets import QApplication

# Default field of view for the projection if not specified.
FIELD_OF_VIEW_DEG = 100.0


@dataclass
class ScreenGeometry:
    """Simple container for screen space geometry."""

    center: Tuple[int, int]
    radius: int


# -----------------------------
# Math and Color Helpers
# -----------------------------


def _clamp01(x: float) -> float:
    """Clamps a value to the [0, 1] range."""
    return 0.0 if x < 0.0 else 1.0 if x > 1.0 else x


def _lerp_color(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    """Linearly interpolates between two colors represented as NumPy arrays."""
    return a + (b - a) * t


def _angle_between(alt1_deg: float, az1_deg: float, alt2_deg: float, az2_deg: float) -> float:
    """Calculates the angular separation between two points in alt-az coordinates."""
    a1, z1 = math.radians(alt1_deg), math.radians(az1_deg)
    a2, z2 = math.radians(alt2_deg), math.radians(az2_deg)
    d_az = z2 - z1

    # Ensure d_az is in [-pi, pi] for numerical stability
    d_az = (d_az + math.pi) % (2 * math.pi) - math.pi

    cos_g = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(d_az)
    cos_g = max(-1.0, min(1.0, cos_g))  # Clamp to handle potential floating point errors
    return math.acos(cos_g)


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """Performs a smooth Hermite interpolation between 0 and 1."""
    t = _clamp01((x - edge0) / (edge1 - edge0))
    return t * t * (3.0 - 2.0 * t)


# -----------------------------
# Sky and Sun Color Model
# -----------------------------


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
    Calculates a simple, stable heuristic sky color.

    The color is returned as linear RGB in the [0, 1] range. This model is
    designed to be numerically stable, especially when viewing the zenith.

    Returns:
        A NumPy array for the (r, g, b) color.
    """
    if sun_alt_deg <= -10.0:
        return np.array([0.0, 0.0, 0.0])

    # 1. Get the base sunlight color.
    sun_color = get_sun_color(sun_alt_deg)

    # 2. Determine brightness based on the angle between the view direction and the sun.
    sun_angle = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)
    cosg = math.cos(sun_angle)
    brightness = (1.0 + cosg) * 0.5  # Range [0, 1]
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


# A small epsilon to avoid singularities at the zenith.
EPS = 1e-3


def clamp_alt(alt: float, alt_min: float = -60.0, alt_max: float = 89.5) -> float:
    """Clamps altitude to a safe range, avoiding the exact zenith."""
    if alt < alt_min:
        return alt_min
    if alt > alt_max:
        return alt_max
    # Avoid the singularity at exactly 90 degrees.
    if abs(alt - 90.0) < EPS:
        return 90.0 - EPS
    return alt


# -----------------------------
# Projection: screen -> (alt, az)
# -----------------------------
def screen_to_altaz_equidistant(x: int, y: int, geometry: ScreenGeometry, view_center: Tuple[float, float], fov_deg: float) -> Tuple[float, float]:
    """
    Inverse azimuthal equidistant projection.

    Maps a screen coordinate (x, y) to a celestial coordinate (altitude, azimuth).
    - The disc radius maps to half the field of view.
    - Altitude: 0° is the horizon, 90° is the zenith.
    - Azimuth: 0° is North, 90° is East (clockwise).
    """
    cx, cy = geometry.center
    R = geometry.radius
    dx, dy = x - cx, y - cy
    rho = math.hypot(dx, dy)

    if rho == 0 or R <= 0:
        return view_center

    # Angle on screen (0 is up/North, clockwise is positive)
    psi = math.atan2(dx, -dy)

    # Angular distance from the center of the view
    rho_ratio = min(rho / R, 1.0)
    theta = math.radians(rho_ratio * fov_deg)

    # Spherical trigonometry to find the new alt/az
    alt_c, az_c = view_center
    alt_c = clamp_alt(alt_c)
    lat1, lon1 = math.radians(alt_c), math.radians(az_c)

    sin_lat1, cos_lat1 = math.sin(lat1), math.cos(lat1)
    sin_theta, cos_theta = math.sin(theta), math.cos(theta)
    sin_psi, cos_psi = math.sin(psi), math.cos(psi)

    sin_lat2 = sin_lat1 * cos_theta + cos_lat1 * sin_theta * cos_psi
    sin_lat2 = max(-1.0, min(1.0, sin_lat2))
    lat2 = math.asin(sin_lat2)
    cos_lat2 = math.cos(lat2)

    if abs(cos_lat2) < 1e-6:
        lon2 = lon1  # Azimuth is undefined at zenith/nadir
    else:
        y_ = sin_psi * sin_theta * cos_lat1
        x_ = cos_theta - sin_lat1 * sin_lat2
        lon2 = lon1 + math.atan2(y_, x_)

    alt = math.degrees(lat2)
    az = (math.degrees(lon2) + 360.0) % 360.0
    return alt, az


# -----------------------------
# The Disc Renderer Under Test
# -----------------------------
def draw_sky_color_disc(
    painter: QPainter,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    sun_alt: float,
    sun_az: float,
    *,
    fov_deg: float,
    pixel_step: int = 1,
    scale: float = 1.0,
    exposure: float = 1.0,
    saturation: float = 1.2,
    ground_color: Tuple[float, float, float] = (0.1, 0.1, 0.1),
    to_altaz: Optional[Callable] = None,
) -> None:
    """Renders the sky color disc onto a QPainter."""
    if to_altaz is None:
        to_altaz = screen_to_altaz_equidistant

    R = int(geometry.radius)
    if R <= 0 or scale <= 0.0 or pixel_step < 1:
        return

    cx, cy = geometry.center
    buf_R = int(round(R * scale))
    if buf_R < 2:
        return

    buf_w = buf_h = buf_R * 2
    inv_scale = 1.0 / scale
    x0, y0 = cx - R, cy - R
    r2_buf = buf_R * buf_R

    img = QImage(buf_w, buf_h, QImage.Format.Format_RGB32)
    img.fill(QColor(0, 0, 0))

    # Iterate over each pixel in the buffer, calculate its color, and set it.
    for by in range(0, buf_h, pixel_step):
        y = int(round(y0 + by * inv_scale))
        dy_b = by - buf_R
        for bx in range(0, buf_w, pixel_step):
            dx_b = bx - buf_R
            if (dx_b * dx_b + dy_b * dy_b) > r2_buf:
                continue  # Skip pixels outside the disc

            x = int(round(x0 + bx * inv_scale))

            # 1. Project screen pixel to sky coordinate
            alt, az = to_altaz(x, y, geometry, view_center, fov_deg)

            # 2. Get base color from the model
            rgb = get_sky_color(alt, az, sun_alt, sun_az)

            # 3. Apply saturation and exposure adjustments
            gray = np.dot(rgb, [0.299, 0.587, 0.114])
            rgb = _lerp_color(np.array([gray, gray, gray]), rgb, saturation)
            rgb *= exposure

            # 4. Fade to ground color below the horizon
            if alt < 0.0:
                t = _clamp01((alt + 5.0) / 5.0)  # Fade between -5 and 0 degrees
                rgb = _lerp_color(np.array(ground_color), rgb, t)

            # 5. Convert to final color and set pixel
            rgb = np.clip(rgb, 0.0, 1.0)
            img.setPixel(bx, by, QColor.fromRgbF(rgb[0], rgb[1], rgb[2]).rgb())

    # Draw the generated image to the painter.
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    target_rect = QRect(cx - R, cy - R, 2 * R, 2 * R)
    painter.drawImage(target_rect, img)
    painter.restore()


# -----------------------------
# Main: Generate a Test PNG
# -----------------------------
def main() -> None:
    """Parses arguments and runs the test rendering."""
    ap = argparse.ArgumentParser(description="Test utility for sky disc rendering.")
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=320)
    ap.add_argument("--fov", type=float, default=FIELD_OF_VIEW_DEG)
    ap.add_argument("--sun-alt", type=float, default=20.0)
    ap.add_argument("--sun-az", type=float, default=270.0, help="(0=N, 90=E)")
    ap.add_argument("--center-alt", type=float, default=90.0)
    ap.add_argument("--center-az", type=float, default=0.0)
    ap.add_argument("--scale", type=float, default=1.0, help="Internal render scale.")
    ap.add_argument("--pixel-step", type=int, default=1, help="Render every Nth pixel.")
    ap.add_argument("--exposure", type=float, default=1.0)
    ap.add_argument("--saturation", type=float, default=1.2)
    ap.add_argument("--outfile", type=str, default="sky_test.png")
    args = ap.parse_args()

    # A QApplication is needed to handle QImage and QPainter.
    app = QApplication([])

    W, H = args.width, args.height
    # A simplified geometry calculation for this test script.
    margin_x, margin_y = 10, 10
    radius = max(2, (W - margin_x * 2) // 2)
    center = (W // 2, H // 2)
    geometry = ScreenGeometry(center=center, radius=radius)

    # Prepare an image canvas.
    canvas = QImage(W, H, QImage.Format.Format_RGB32)
    canvas.fill(QColor(0, 0, 0))

    painter = QPainter(canvas)

    # Draw the sky-color disc.
    draw_sky_color_disc(
        painter=painter,
        geometry=geometry,
        view_center=(args.center_alt, args.center_az),
        sun_alt=args.sun_alt,
        sun_az=args.sun_az,
        fov_deg=args.fov,
        pixel_step=args.pixel_step,
        scale=args.scale,
        exposure=args.exposure,
        saturation=args.saturation,
        ground_color=(0.1, 0.1, 0.1),
    )

    painter.end()

    # Save the final image.
    ok = canvas.save(args.outfile, "PNG")
    print(f"Saved: {args.outfile}, Success: {ok}")


if __name__ == "__main__":
    main()
