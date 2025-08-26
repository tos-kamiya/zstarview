import math
from typing import Optional, Tuple

import numpy as np
from PySide6.QtGui import QColor, QImage, QPainter, Qt

from ..astro import altaz_to_normalized_xy
from ..types import ScreenGeometry, SolarEclipseInfo
from .draw import normalized_to_screen_xy


TURBIDITY = 5   # 2 (clear blue sky) to 10 (hazy white sky)


def _ground_cutoff(alt: float, alpha: float, fade_hi: float = 0.0, fade_lo: float = -2.0) -> float:
    """
    Returns an alpha value based on altitude, for fading near the horizon.
       alt >= fade_hi  -> alpha (equivalent to 1)
       fade_lo < alt < fade_hi -> linear fade
       alt <= fade_lo -> 0 (not drawn)

    Args:
        alt: The altitude in degrees.
        alpha: The base alpha value.
        fade_hi: The altitude above which alpha is at its maximum.
        fade_lo: The altitude below which alpha is zero.

    Returns:
        The calculated alpha value.
    """
    if alt <= fade_lo:
        return 0.0
    if alt >= fade_hi:
        return alpha
    t = (alt - fade_lo) / (fade_hi - fade_lo)  # [fade_hi, fade_lo] -> [0,1]
    return alpha * t


def _fwd_offset_altaz(alt_c: float, az_c: float, theta_deg: float, psi_deg: float) -> Tuple[float, float]:
    """
    Calculates the (alt, az) of a point that is at an angular distance 'theta'
    and direction 'psi' from a center point (alt_c, az_c).

    Args:
        alt_c: Center altitude in degrees.
        az_c: Center azimuth in degrees.
        theta_deg: Angular distance from the center in degrees.
        psi_deg: Direction from the center in degrees (0° North, 90° East, clockwise).

    Returns:
        A tuple of (altitude, azimuth) in degrees for the new point.
    """
    phi1 = math.radians(alt_c)
    lambda1 = math.radians(az_c)
    theta = math.radians(theta_deg)
    psi = math.radians(psi_deg)

    sin_phi1, cos_phi1 = math.sin(phi1), math.cos(phi1)
    sin_theta, cos_theta = math.sin(theta), math.cos(theta)

    sin_phi2 = sin_phi1 * cos_theta + cos_phi1 * sin_theta * math.cos(psi)
    sin_phi2 = max(-1.0, min(1.0, sin_phi2))
    phi2 = math.asin(sin_phi2)

    y = math.sin(psi) * sin_theta * cos_phi1
    x = cos_theta - sin_phi1 * sin_phi2
    lambda2 = lambda1 + math.atan2(y, x)

    alt = math.degrees(phi2)
    az = (math.degrees(lambda2) + 360.0) % 360.0
    return alt, az


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """
    Performs a smooth Hermite interpolation between 0 and 1 when edge0 < x < edge1.
    """
    # Hermite smoothstep
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _angle_between(alt1_deg: float, az1_deg: float, alt2_deg: float, az2_deg: float) -> float:
    """
    Calculates the angular distance between two points on the celestial sphere.

    Args:
        alt1_deg: Altitude of the first point in degrees.
        az1_deg: Azimuth of the first point in degrees.
        alt2_deg: Altitude of the second point in degrees.
        az2_deg: Azimuth of the second point in degrees.

    Returns:
        The angular distance in radians.
    """
    a1, z1 = math.radians(alt1_deg), math.radians(az1_deg)
    a2, z2 = math.radians(alt2_deg), math.radians(az2_deg)
    cos_g = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(z2 - z1)
    cos_g = max(-1.0, min(1.0, cos_g))
    return math.acos(cos_g)


def _lerp_color(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    """
    Linearly interpolates between two colors.
    """
    return a + (b - a) * t


def get_sun_color(sun_alt_deg: float) -> np.ndarray:
    """
    Step 1: Determine the color of sunlight based on the sun's altitude.

    Args:
        sun_alt_deg: The sun's altitude in degrees.

    Returns:
        A numpy array of (r, g, b) float values representing the sun's color.
    """
    # Define colors
    zenith_color = np.array([0.3, 0.48, 0.96])  # Color at zenith (blue)
    horizon_color = np.array([1.0, 0.61, 0.32])  # Color at horizon (orange)
    night_color = np.array([0.01, 0.02, 0.05])  # Night color (dark blue)

    # Normalize sun altitude from -7 degrees (sunset) to 90 degrees (zenith) to a 0-1 range
    t = np.clip((sun_alt_deg + 7.0) / 97.0, 0.0, 1.0)
    t = math.pow(t, 0.4)

    # Daytime color (horizon to zenith)
    day_color = _lerp_color(horizon_color, zenith_color, t)

    # Mix day and night colors (to represent the rapid darkening near the horizon)
    fade = np.clip((sun_alt_deg + 2.0) / 12.0, 0.0, 1.0)  # Fade between -2 and -10 degrees

    return _lerp_color(night_color, day_color, fade)


def get_sky_color(view_altaz: Tuple[float, float], sun_altaz: Tuple[float, float]) -> np.ndarray:
    """
    Calculate the sky color for a given viewing direction and sun position.

    This function models the sky color based on:
    - Sun altitude and azimuth (position of the sun).
    - Viewing direction (altitude, azimuth).
    - Atmospheric turbidity, which controls how clear (deep blue) or hazy
      (whitish and desaturated) the sky appears.

    Args:
        view_altaz: Tuple (altitude_deg, azimuth_deg) of the viewing direction.
        sun_altaz: Tuple (altitude_deg, azimuth_deg) of the sun's position.

    Returns:
        np.ndarray: RGB color (float values 0..1) representing the sky color.
    """
    sun_alt_deg, sun_az_deg = sun_altaz
    view_alt_deg, view_az_deg = view_altaz

    if sun_alt_deg <= -10.0:
        # Completely dark if the sun is more than 10° below the horizon
        return np.array([0.0, 0.0, 0.0])

    sun_color = get_sun_color(sun_alt_deg)

    # --- Turbidity effect (0..1 normalized) ---
    tau = np.clip((TURBIDITY - 2.0) / (10.0 - 2.0), 0.0, 1.0)  # 2.0..10.0 -> 0..1

    # 1) Brightness in the sun-facing direction
    #    Higher turbidity reduces directional contrast (lower exponent).
    sun_angle = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)  # [0..pi]
    cosg = math.cos(sun_angle)
    f = max(0.0, cosg)  # set to 0 on the antisolar side
    exp = 2.0 - 0.6 * tau  # tau=0 -> 2.0, tau=1 -> 1.4
    brightness = f ** exp

    # 2) Altitude tone factor
    t = np.clip(view_alt_deg / 90.0, 0.0, 1.0)  # 0(horizon) -> 1(zenith)

    # Darkening toward zenith
    # Higher turbidity reduces the darkening (sky looks brighter/whiter at the zenith).
    zenith_dim_coef = 0.30 * (1.0 - 0.5 * tau)   # halved at tau=1
    zenith_dim = 1.0 - zenith_dim_coef * t

    # White mixing near the horizon
    # Higher turbidity increases the amount of white mixed in.
    horizon_k = 0.25 * (0.8 + 0.6 * tau)
    horizon_mix = (1.0 - t) * horizon_k

    # Twilight correction (smooth transition -10..0°)
    if sun_alt_deg < 0.0:
        twilight = _smoothstep(-10.0, 0.0, sun_alt_deg)
    else:
        twilight = 1.0

    # Composite color
    lit = sun_color * brightness
    color = _lerp_color(lit, np.array([1.0, 1.0, 1.0]), horizon_mix * twilight)
    color *= zenith_dim

    return np.clip(color, 0.0, 1.0)


def grade_color(
    color: np.ndarray,
    *,
    saturation: float = 1.0,
    exposure: float = 1.0,
    gamma: float = 1.0,
) -> np.ndarray:
    """
    Sky-disc-only color grading:
      1) luma-based saturation (BT.601/709)
      2) exposure scaling
      3) simple power-law gamma (on current space)
    Assumes color is a numpy array in [0,1] and *not* premultiplied.
    """
    # --- Luma ---
    luma = np.dot(color, np.array([0.299, 0.587, 0.114]))  # BT.601

    # --- Saturation: lerp(gray, original, s) ---
    color = luma + (color - luma) * saturation

    # --- Exposure ---
    color *= exposure

    # --- Gamma (simple power-law) ---
    if gamma > 0.0 and gamma != 1.0:
        inv = 1.0 / gamma
        color = np.power(np.maximum(0.0, color), inv)

    return np.clip(color, 0.0, 1.0)


def draw_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],  # (alt, az)
    sun_altaz: Tuple[float, float],  # (alt, az)
    *,
    exposure: float = 1.0,
    saturation: float = 1.2,
    alpha: float = 1.0,
    eclipse_factor: float = 1.0,
    mask_fov_deg: float = 90.0,
    # --- Sampling density (knobs for quality vs. speed) ---
    sample_step_px: int = 10,  # Target pixel interval (basis for determining Δθ)
    min_ang_samples: int = 8,  # Minimum number of samples for each ring
    # --- Parameters for stabilizing Δθ estimation ---
    deriv_probe_deg: float = 0.25,  # Small angle (degrees) for finite difference of dr/dθ
    min_theta_step_deg: float = 0.2,  # Lower limit for Δθ (degrees)
    max_theta_step_deg: float = 6.0,  # Upper limit for Δθ (degrees)
) -> QImage:
    """
    Draws the sky color disc with dynamic sampling and returns it as a QImage.

    The radial step (Δθ) is dynamically determined so that the on-screen
    pixel distance is around `sample_step_px`. The number of samples in the
    circumferential direction is determined by the "measured radius" of the
    ring and the sample pitch.

    Args:
        geometry: The screen geometry.
        view_center: The (alt, az) of the view center.
        sun_altaz: The (alt, az) of the sun.
        exposure: Exposure adjustment for the final color.
        saturation: Saturation adjustment for the final color.
        alpha: Overall alpha transparency.
        eclipse_factor: Multiplicative factor to darken the sky during a solar eclipse (set < 1.0 to simulate dimming).
        mask_fov_deg: Full field-of-view angle for the visibility mask.
        sample_step_px: Target pixel distance steps.
        min_ang_samples: Minimum number of circumferential samples.
        deriv_probe_deg: Probe angle for derivative estimation.
        min_theta_step_deg: Minimum angular step.
        max_theta_step_deg: Maximum angular step.

    Returns:
        A QImage of the rendered sky disc.
    """
    assert altaz_to_normalized_xy and normalized_to_screen_xy and get_sky_color

    R = int(geometry.radius)
    if R < 2:
        return QImage(2 * R, 2 * R, QImage.Format.Format_ARGB32_Premultiplied)

    sample_step_px = max(2, min(R // 128, sample_step_px))
    tile_size = int(sample_step_px * 1.5 + 1)

    cx, cy = geometry.center

    img_w = img_h = R * 2
    img = QImage(img_w, img_h, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(QColor(0, 0, 0, 0))

    ip = QPainter(img)
    ip.setPen(Qt.PenStyle.NoPen)

    # Pseudo clip to a circle
    ip.setCompositionMode(QPainter.CompositionMode_Source)
    ip.setBrush(QColor(0, 0, 0, 255))
    ip.drawEllipse(0, 0, 2 * R, 2 * R)
    ip.setCompositionMode(QPainter.CompositionMode_SourceAtop)  # no alpha drawing hereafter

    # Draw the background disc
    ip.fillRect(0, 0, 2 * R, 2 * R, QColor(0, 0, 0, 255))

    # Clamp to avoid zenith singularity
    EPS = 1e-3
    alt_c, az_c = view_center
    alt_c = max(-(90.0 - EPS), min(90.0 - EPS, alt_c))

    # The outer circumference is always 90° to match the normalization spec
    theta_max = mask_fov_deg

    # Function to "measure" the local radius r_px(θ)
    def screen_radius_px(theta_deg: float) -> float:
        alt_p, az_p = _fwd_offset_altaz(alt_c, az_c, theta_deg, 0.0)
        nx, ny = altaz_to_normalized_xy(alt_p, az_p, (alt_c, az_c))
        # Convert to image coordinates relative to image center (R, R)
        sx, sy = normalized_to_screen_xy(nx, ny, geometry)
        return max(0.0, math.hypot(sx - cx, sy - cy))

    gamma = (1.0 - alpha) * 0.2 + 1.0 if alpha < 1.0 else 1.0

    # Advance angle from 0 to 90° (Δθ is dynamic)
    theta = 0.0
    half = max(1, int(tile_size // 2))
    while True:
        # Current ring radius (px)
        r_px = screen_radius_px(theta)

        # Number of circumferential samples
        circumference = max(1.0, 2.0 * math.pi * r_px)
        n_ang = max(min_ang_samples, int(round(circumference / max(1.0, float(sample_step_px)))))

        # Fill the ring
        for i in range(n_ang):
            psi_deg = (360.0 * i) / n_ang  # 0° North, 90° East, clockwise

            # (θ,ψ) -> (alt,az)
            alt, az = _fwd_offset_altaz(alt_c, az_c, theta, psi_deg)

            cutoff = _ground_cutoff(alt, alpha, fade_hi=0.5, fade_lo=-2.0)
            if cutoff <= 0.0:
                continue

            # (alt,az) -> screen -> image coordinates
            nx, ny = altaz_to_normalized_xy(alt, az, (alt_c, az_c))
            sx, sy = normalized_to_screen_xy(nx, ny, geometry)
            # Convert to image coordinates (origin at top-left of the QImage)
            xi = int(round(sx - (cx - R)))
            yi = int(round(sy - (cy - R)))

            if xi < 0 or xi >= img_w or yi < 0 or yi >= img_h:
                continue

            color = get_sky_color((alt, az), sun_altaz)
            color = grade_color(
                color,
                saturation=saturation,
                exposure=exposure,
                gamma=gamma,
            )

            color *= cutoff * eclipse_factor
            ip.fillRect(xi - half, yi - half, 2 * half, 2 * half, QColor.fromRgbF(*color, 1.0))

        # Termination condition (ensure the 90° ring is always drawn)
        if theta >= theta_max - 1e-6:
            break

        # ---- Determine Δθ from screen distance: Δθ ≈ sample_step_px / (dr_px/dθ) ----
        # Estimate dr/dθ using finite differences
        probe = min(deriv_probe_deg, theta_max - theta)
        r_next = screen_radius_px(theta + probe)
        dr_dtheta = (r_next - r_px) / max(1e-6, probe)

        if dr_dtheta <= 1e-6:
            dtheta = max_theta_step_deg  # Safe side (almost no change / numerical error at edges)
        else:
            dtheta = float(sample_step_px) / dr_dtheta

        # Limit the angle step
        dtheta = max(min_theta_step_deg, min(max_theta_step_deg, dtheta))

        theta += dtheta

    ip.end()

    return img
