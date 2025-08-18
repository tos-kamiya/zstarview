import math
from typing import Tuple

from PySide6.QtGui import QColor, QImage, QPainter, QPainterPath

from ..astro import altaz_to_normalized_xy
from ..types import ScreenGeometry
from .draw import (
    normalized_to_screen_xy
)


def _alpha_from_alt(alt: float, alpha: float, fade_hi: float = 0.0, fade_lo: float = -2.0) -> float:
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
    t = (alt - fade_lo) / (fade_hi - fade_lo)   # [fade_hi, fade_lo] -> [0,1]
    return alpha * t


def _fwd_offset_altaz(alt_c: float, az_c: float,
                      theta_deg: float, psi_deg: float) -> Tuple[float, float]:
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
    phi1  = math.radians(alt_c)
    lambda1  = math.radians(az_c)
    theta   = math.radians(theta_deg)
    psi   = math.radians(psi_deg)

    sin_phi1, cos_phi1 = math.sin(phi1), math.cos(phi1)
    sin_theta, cos_theta   = math.sin(theta), math.cos(theta)

    sin_phi2 = sin_phi1 * cos_theta + cos_phi1 * sin_theta * math.cos(psi)
    sin_phi2 = max(-1.0, min(1.0, sin_phi2))
    phi2 = math.asin(sin_phi2)

    y = math.sin(psi) * sin_theta * cos_phi1
    x = cos_theta - sin_phi1 * sin_phi2
    lambda2 = lambda1 + math.atan2(y, x)

    alt = math.degrees(phi2)
    az  = (math.degrees(lambda2) + 360.0) % 360.0
    return alt, az


def _clamp01(x: float) -> float:
    """Clamps a value to the range [0, 1]."""
    return max(0.0, min(1.0, x))


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """
    Performs a smooth Hermite interpolation between 0 and 1 when edge0 < x < edge1.
    """
    # Hermite smoothstep
    t = _clamp01((x - edge0) / (edge1 - edge0))
    return t * t * (3.0 - 2.0 * t)


def _deg2rad(d: float) -> float:
    """Converts degrees to radians."""
    return d * math.pi / 180.0


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
    a1, z1 = _deg2rad(alt1_deg), _deg2rad(az1_deg)
    a2, z2 = _deg2rad(alt2_deg), _deg2rad(az2_deg)
    cos_g = (math.sin(a1)*math.sin(a2) + math.cos(a1)*math.cos(a2)*math.cos(z2 - z1))
    cos_g = max(-1.0, min(1.0, cos_g))
    return math.acos(cos_g)


def _lerp_color(a: Tuple[float,float,float],
                b: Tuple[float,float,float], t: float) -> Tuple[float,float,float]:
    """
    Linearly interpolates between two colors.
    """
    return (a[0] + (b[0]-a[0]) * t,
            a[1] + (b[1]-a[1]) * t,
            a[2] + (b[2]-a[2]) * t)


def get_sun_color(sun_alt_deg: float) -> Tuple[float, float, float]:
    """
    Step 1: Determine the color of sunlight based on the sun's altitude.

    Args:
        sun_alt_deg: The sun's altitude in degrees.

    Returns:
        A tuple of (r, g, b) float values representing the sun's color.
    """
    # Define colors
    zenith_color = (0.3, 0.5, 1.0)  # Color at zenith (blue)
    horizon_color = (1.0, 0.8, 0.4) # Color at horizon (orange)
    night_color = (0.01, 0.02, 0.05) # Night color (dark blue)

    # Normalize sun altitude from -5 degrees (sunset) to 90 degrees (zenith) to a 0-1 range
    t = _clamp01((sun_alt_deg + 5.0) / 95.0)

    # Daytime color (horizon to zenith)
    day_color = _lerp_color(horizon_color, zenith_color, t)

    # Mix day and night colors (to represent the rapid darkening near the horizon)
    fade = _clamp01(sun_alt_deg / 10.0) # Fade between 0 and 10 degrees

    return _lerp_color(night_color, day_color, fade)


def get_sky_color(view_altaz: Tuple[float, float], sun_altaz: Tuple[float, float]) -> Tuple[float, float, float]:
    """
    Calculates the sky color for a given viewing direction and sun position.

    Args:
        view_altaz: A tuple of (altitude, azimuth) for the viewing direction.
        sun_altaz: A tuple of (altitude, azimuth) for the sun's position.

    Returns:
        A tuple of (r, g, b) float values for the sky color.
    """
    # After sunset, it gradually darkens (completely dark at -10°, towards day at 0°)
    sun_alt_deg, sun_az_deg = sun_altaz
    view_alt_deg, view_az_deg = view_altaz

    if sun_alt_deg <= -10.0:
        return (0.0, 0.0, 0.0)

    # Basic sun color (assumed 0..1: preferably linear RGB)
    sun_color = get_sun_color(sun_alt_deg)

    # 1) Brightness based on angle to the sun (stable even at zenith)
    gamma = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)  # [0..pi]
    cosg = math.cos(gamma)
    brightness = (1.0 + cosg) * 0.5          # 0..1
    brightness = brightness ** 2.0            # Emphasize the sun-facing direction

    # 2) Tone based on altitude (darker at zenith, whitish at horizon)
    t = _clamp01(view_alt_deg / 90.0)         # 0(horizon) -> 1(zenith)
    zenith_darkness = 0.5 + 0.5 * t           # 0.5..1.0 (darker towards zenith)
    horizon_whiteness = (1.0 - t) * 0.3       # 0.3..0.0 (whiter towards horizon)

    # 3) Twilight correction (interpolate -10..0° to 0..1)
    if sun_alt_deg < 0.0:
        twilight = _smoothstep(-10.0, 0.0, sun_alt_deg)  # 0..1
    else:
        twilight = 1.0

    # Composite (simple: mix of white + sun contribution)
    r = sun_color[0] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight
    g = sun_color[1] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight
    b = sun_color[2] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight

    # 4) Clip (final safety)
    return (_clamp01(r), _clamp01(g), _clamp01(b))


def draw_sky_color_disc(
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],   # (alt, az)
    sun_altaz: Tuple[float, float],     # (alt, az)
    *,
    exposure: float = 1.0,
    saturation: float = 1.2,
    alpha: float = 1.0,
    # --- Sampling density (knobs for quality vs. speed) ---
    ring_step_px: int = 6,            # Target pixel interval in the radial direction (basis for determining Δθ)
    sample_pitch_px: float = 6.0,     # Target sample interval on each ring (px)
    min_ang_samples: int = 8,         # Minimum number of samples for each ring
    dot_size: int = 12,               # Side length of the tile (square) in px: roughly 2*ring_step_px is a guideline
    # --- Parameters for stabilizing Δθ estimation ---
    deriv_probe_deg: float = 0.25,    # Small angle (degrees) for finite difference of dr/dθ
    min_theta_step_deg: float = 0.2,  # Lower limit for Δθ (degrees)
    max_theta_step_deg: float = 6.0,  # Upper limit for Δθ (degrees)
) -> QImage:
    """
    Draws the sky color disc with dynamic sampling and returns it as a QImage.

    The radial step (Δθ) is dynamically determined so that the on-screen
    pixel distance is around `ring_step_px`. The number of samples in the
    circumferential direction is determined by the "measured radius" of the
    ring and the sample pitch.

    Args:
        geometry: The screen geometry.
        view_center: The (alt, az) of the view center.
        sun_altaz: The (alt, az) of the sun.
        exposure: Exposure adjustment for the final color.
        saturation: Saturation adjustment for the final color.
        alpha: Overall alpha transparency.
        ring_step_px: Target pixel distance for radial steps.
        sample_pitch_px: Target pixel distance for circumferential samples.
        min_ang_samples: Minimum number of circumferential samples.
        dot_size: The size of the colored dots to draw.
        deriv_probe_deg: Probe angle for derivative estimation.
        min_theta_step_deg: Minimum angular step.
        max_theta_step_deg: Maximum angular step.

    Returns:
        A QImage of the rendered sky disc.
    """
    assert altaz_to_normalized_xy and normalized_to_screen_xy and get_sky_color

    R = int(geometry.radius)
    if R < 2:
        return QImage(2*R, 2*R, QImage.Format.Format_ARGB32_Premultiplied)

    cx, cy = geometry.center

    img_w = img_h = R * 2
    img = QImage(img_w, img_h, QImage.Format.Format_ARGB32)
    img.fill(QColor(0, 0, 0, 0))  # Transparent

    ip = QPainter(img)
    ip.setRenderHint(QPainter.RenderHint.Antialiasing, False)
    ip.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, False)

    # Clip to a circle
    path = QPainterPath()
    path.addEllipse(0, 0, 2*R, 2*R)  # Circle at (0,0) with radius R in image coords
    ip.setClipPath(path)

    # Clamp to avoid zenith singularity
    EPS = 0.01
    alt_c, az_c = view_center
    if alt_c <= -90.0: alt_c = -(90.0 - EPS)
    if alt_c >=  90.0: alt_c =  (90.0 - EPS)

    # The outer circumference is always 90° to match the normalization spec
    theta_max = 90.0

    # Function to "measure" the local radius r_px(θ)
    def screen_radius_px(theta_deg: float) -> float:
        alt_p, az_p = _fwd_offset_altaz(alt_c, az_c, theta_deg, 0.0)
        nx, ny = altaz_to_normalized_xy(alt_p, az_p, (alt_c, az_c))
        # Convert to image coordinates relative to image center (R, R)
        sx, sy = normalized_to_screen_xy(nx, ny, geometry)
        return max(0.0, math.hypot(sx - cx, sy - cy))

    # Advance angle from 0 to 90° (Δθ is dynamic)
    theta = 0.0
    half = max(1, int(dot_size // 2))
    while True:
        # Current ring radius (px)
        r_px = screen_radius_px(theta)

        # Number of circumferential samples
        circumference = max(1.0, 2.0 * math.pi * r_px)
        n_ang = max(min_ang_samples, int(round(circumference / max(1.0, float(sample_pitch_px)))))

        # Fill the ring
        for i in range(n_ang):
            psi_deg = (360.0 * i) / n_ang  # 0° North, 90° East, clockwise

            # (θ,ψ) -> (alt,az)
            alt, az = _fwd_offset_altaz(alt_c, az_c, theta, psi_deg)

            # (alt,az) -> screen -> image coordinates
            nx, ny = altaz_to_normalized_xy(alt, az, (alt_c, az_c))
            sx, sy = normalized_to_screen_xy(nx, ny, geometry)
            # Convert to image coordinates (origin at top-left of the QImage)
            xi = int(round(sx - (cx - R)))
            yi = int(round(sy - (cy - R)))

            if xi < 0 or xi >= img_w or yi < 0 or yi >= img_h:
                continue

            aa = _alpha_from_alt(alt, alpha, fade_hi=1.0, fade_lo=-1.0)
            if aa <= 0.0:
                continue

            rr, gg, bb = get_sky_color((alt, az), sun_altaz)
            gray = rr*0.299 + gg*0.587 + bb*0.114
            rr, gg, bb = _lerp_color((gray, gray, gray), (rr, gg, bb), saturation)
            rr *= exposure; gg *= exposure; bb *= exposure
            rr = _clamp01(rr); gg = _clamp01(gg); bb = _clamp01(bb)

            ip.fillRect(xi - half, yi - half, 2*half, 2*half, QColor.fromRgbF(rr, gg, bb, aa))

        # Termination condition (ensure the 90° ring is always drawn)
        if theta >= theta_max - 1e-6:
            break

        # ---- Determine Δθ from screen distance: Δθ ≈ ring_step_px / (dr_px/dθ) ----
        # Estimate dr/dθ using finite differences
        probe = min(deriv_probe_deg, theta_max - theta)
        r_next = screen_radius_px(theta + probe)
        dr_dtheta = (r_next - r_px) / max(1e-6, probe)

        if dr_dtheta <= 1e-6:
            dtheta = max_theta_step_deg  # Safe side (almost no change / numerical error at edges)
        else:
            dtheta = float(ring_step_px) / dr_dtheta

        # Limit the angle step
        dtheta = max(min_theta_step_deg, min(max_theta_step_deg, dtheta))

        theta += dtheta

    ip.end()
    return img