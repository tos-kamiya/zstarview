#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Minimal tester for draw_sky_color_disc.
- Renders a sky-color disc to a QImage and saves it as PNG.
- Center of view is zenith (alt=90°) by default to reproduce/check the N–S seam.
- Requires: PySide6
"""

"""
python test_draw_sky_disc.py --width 256 --height 256 --fov 200 \
  --sun-alt 15 --sun-az 270 --center-alt 90 --center-az 0 \
  --outfile sky_test.png
"""

import math
import argparse
from dataclasses import dataclass
from typing import Tuple

from PySide6.QtGui import QImage, QColor, QPainter
from PySide6.QtCore import QRect
from PySide6.QtWidgets import QApplication

# -----------------------------
# Config & helpers
# -----------------------------
FIELD_OF_VIEW_DEG = 200.0  # adjust as needed


@dataclass
class ScreenGeometry:
    center: Tuple[int, int]
    radius: int


def _clamp01(x: float) -> float:
    return 0.0 if x < 0.0 else 1.0 if x > 1.0 else x


def _lerp_color(a: Tuple[float, float, float], b: Tuple[float, float, float], t: float) -> Tuple[float, float, float]:
    # Simple linear interpolation; if you want "saturation boost" behavior,
    # pass t>1.0 (works as extrapolation)
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t)


def _angle_between(alt1_deg, az1_deg, alt2_deg, az2_deg) -> float:
    a1, z1 = math.radians(alt1_deg), math.radians(az1_deg)
    a2, z2 = math.radians(alt2_deg), math.radians(az2_deg)
    d_az = z2 - z1
    d_az = (d_az + math.pi) % (2 * math.pi) - math.pi
    cos_g = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(d_az)
    cos_g = max(-1.0, min(1.0, cos_g))
    return math.acos(cos_g)


# -----------------------------
# Sky color model (simple)
# -----------------------------
def get_sun_color(sun_alt_deg: float) -> Tuple[float, float, float]:
    """
    Very crude daylight tint by solar altitude.
    Return linear RGB in [0,1].
    - High sun: bluer sky light
    - Low sun: warmer tint
    """
    t = _clamp01((sun_alt_deg + 6.0) / 36.0)  # ~[-6,30] → [0,1]
    # interpolate between warm (low) and cool (high)
    warm = (1.0, 0.85, 0.70)  # warm tint
    cool = (0.75, 0.85, 1.0)  # cool tint
    return _lerp_color(warm, cool, t)


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    t = _clamp01((x - edge0) / (edge1 - edge0))
    return t * t * (3.0 - 2.0 * t)


def get_sky_color(
    view_alt_deg: float, view_az_deg: float, sun_alt_deg: float, sun_az_deg: float
) -> Tuple[float, float, float]:
    """
    Simple, stable heuristic sky color. Linear RGB in [0,1].
    Designed to be numerically stable even when view_alt_deg == 90.
    """
    if sun_alt_deg <= -10.0:
        return (0.0, 0.0, 0.0)

    sun_color = get_sun_color(sun_alt_deg)

    gamma = _angle_between(view_alt_deg, view_az_deg, sun_alt_deg, sun_az_deg)
    cosg = math.cos(gamma)
    brightness = ((1.0 + cosg) * 0.5) ** 2.0  # emphasize near the sun

    t_alt = _clamp01(view_alt_deg / 90.0)  # 0 at horizon, 1 at zenith
    zenith_darkness = 0.5 + 0.5 * t_alt  # 0.5..1.0
    horizon_whiteness = (1.0 - t_alt) * 0.3  # 0.3..0.0

    twilight = 1.0 if sun_alt_deg >= 0.0 else _smoothstep(-10.0, 0.0, sun_alt_deg)

    r = sun_color[0] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight
    g = sun_color[1] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight
    b = sun_color[2] * brightness * zenith_darkness * twilight + horizon_whiteness * twilight

    return (_clamp01(r), _clamp01(g), _clamp01(b))


EPS = 1e-3  # 天頂の特異を避けるための余白


def clamp_alt(alt: float, alt_min: float = -60.0, alt_max: float = 89.5) -> float:
    # alt_max は 90 より少し小さく
    if alt < alt_min:
        return alt_min
    if alt > alt_max:
        return alt_max
    # ぴったり90°は避ける
    if abs(alt - 90.0) < EPS:
        return 90.0 - EPS
    return alt


# -----------------------------
# Projection: screen → (alt, az)
# -----------------------------
def screen_to_altaz_equidistant(
    x: int, y: int, geometry: "ScreenGeometry", view_center: Tuple[float, float], fov_deg: float  # (alt_c, az_c)
) -> Tuple[float, float]:
    """
    Azimuthal equidistant projection inverse:
    - Disk radius maps to FOV/2.
    - Alt: 0=horizon, 90=zenith
    - Az:  0=N, 90=E (clockwise)
    """
    cx, cy = geometry.center
    R = geometry.radius
    dx = x - cx
    dy = y - cy
    rho = math.hypot(dx, dy)
    if rho == 0 or R <= 0:
        return view_center

    psi = math.atan2(dx, -dy)  # compass bearing: up/N=0, CW positive

    rho_ratio = min(rho / R, 1.0)
    theta = math.radians(rho_ratio * (fov_deg * 0.5))  # angular distance from center

    alt_c, az_c = view_center
    alt_c = clamp_alt(alt_c)
    lat1 = math.radians(alt_c)
    lon1 = math.radians(az_c)

    sin_lat1 = math.sin(lat1)
    cos_lat1 = math.cos(lat1)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    sin_psi = math.sin(psi)
    cos_psi = math.cos(psi)

    sin_lat2 = sin_lat1 * cos_theta + cos_lat1 * sin_theta * cos_psi
    sin_lat2 = max(-1.0, min(1.0, sin_lat2))
    lat2 = math.asin(sin_lat2)
    cos_lat2 = math.cos(lat2)

    if abs(cos_lat2) < 1e-6:
        lon2 = lon1  # az undefined at zenith/nadir -> keep center az for continuity
    else:
        y_ = sin_psi * sin_theta * cos_lat1
        x_ = cos_theta - sin_lat1 * sin_lat2
        lon2 = lon1 + math.atan2(y_, x_)

    alt = math.degrees(lat2)
    az = (math.degrees(lon2) + 360.0) % 360.0
    return alt, az


# -----------------------------
# The disc renderer under test
# -----------------------------
def draw_sky_color_disc(
    painter: QPainter,
    geometry: "ScreenGeometry",
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
    to_altaz=None,
) -> None:
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
    x0 = cx - R
    y0 = cy - R
    r2 = buf_R * buf_R

    img = QImage(buf_w, buf_h, QImage.Format.Format_RGB32)
    img.fill(QColor(0, 0, 0))

    # raster
    for by in range(0, buf_h, pixel_step):
        y = int(round(y0 + by * inv_scale))
        dy_b = by - buf_R
        dy2 = dy_b * dy_b
        for bx in range(0, buf_w, pixel_step):
            dx_b = bx - buf_R
            if (dx_b * dx_b + dy2) > r2:
                continue  # outside disc

            x = int(round(x0 + bx * inv_scale))

            # 1) to alt/az
            alt, az = to_altaz(x, y, geometry, view_center, fov_deg)

            # 2) linear RGB base
            r, g, b = get_sky_color(alt, az, sun_alt, sun_az)

            # 3) saturation & exposure in linear space
            gray = r * 0.299 + g * 0.587 + b * 0.114
            r, g, b = _lerp_color((gray, gray, gray), (r, g, b), saturation)
            r *= exposure
            g *= exposure
            b *= exposure

            # ground fade: -5..0° → 0..1
            if alt < 0.0:
                t = _clamp01((alt + 5.0) / 5.0)
                r, g, b = _lerp_color(ground_color, (r, g, b), t)

            r = _clamp01(r)
            g = _clamp01(g)
            b = _clamp01(b)
            img.setPixel(bx, by, QColor.fromRgbF(r, g, b).rgb())

    painter.save()
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    target_rect = QRect(cx - R, cy - R, 2 * R, 2 * R)
    painter.drawImage(target_rect, img)
    painter.restore()


# -----------------------------
# Main: generate a small test PNG
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=320)
    ap.add_argument("--fov", type=float, default=FIELD_OF_VIEW_DEG)
    ap.add_argument("--sun-alt", type=float, default=20.0)
    ap.add_argument("--sun-az", type=float, default=270.0)  # W
    ap.add_argument("--center-alt", type=float, default=90.0)  # zenith
    ap.add_argument("--center-az", type=float, default=0.0)  # north
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--pixel-step", type=int, default=1)
    ap.add_argument("--exposure", type=float, default=1.0)
    ap.add_argument("--saturation", type=float, default=1.2)
    ap.add_argument("--outfile", type=str, default="sky_test.png")
    args = ap.parse_args()

    app = QApplication([])

    W, H = args.width, args.height
    # Compute geometry: center and radius similar to your get_screen_geometry
    margin_x, margin_y = 10, 10
    radius = max(2, (W - margin_x * 2) // 2)
    # Place center vertically at middle for testing (simplify)
    center = (W // 2, H // 2)
    geometry = ScreenGeometry(center=center, radius=radius)

    # Prepare an image canvas
    canvas = QImage(W, H, QImage.Format.Format_RGB32)
    canvas.fill(QColor(0, 0, 0))

    painter = QPainter(canvas)

    # Optional: draw your dark radial background first (simplified here)
    painter.fillRect(0, 0, W, H, QColor(0, 0, 0))

    # Draw the sky-color disc
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

    # Save
    ok = canvas.save(args.outfile, "PNG")
    print(f"Saved: {args.outfile}, ok={ok}")


if __name__ == "__main__":
    main()
