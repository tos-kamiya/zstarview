# -*- coding: utf-8 -*-
"""
Image generation utilities, such as creating celestial body sprites.
"""

import math
from typing import Optional, Tuple

import numpy as np
from PIL import Image


def generate_moon_phase_image(
    size: int,
    sun_dir_3d: np.ndarray,
    view_dir_3d: np.ndarray,
    moon_color: Tuple[int, int, int] = (230, 230, 230),
    dark_color: Tuple[int, int, int] = (30, 30, 30),
    earthshine_factor: float = 0.15,
    tint_color: Optional[Tuple[int, int, int, int]] = None,  # RGBA
    edge_soft_px: float = 1.0,
) -> Image.Image:
    """
    Generates a spherical moon-phase image as an RGBA bitmap.

    This function renders a 2D image of the Moon, accurately simulating its phase
    based on the directions of the Sun and the viewer. It includes several
    visual features for realism and aesthetic quality.

    Features:
        - Pixels outside the lunar disc are fully transparent (alpha = 0).
        - Only the hemisphere visible from the viewer's direction is rendered.
        - Lit areas use Lambertian shading with `moon_color`.
        - Unlit but visible areas are faintly illuminated by `dark_color`,
          simulating earthshine.
        - An optional `tint_color` (e.g., for eclipses) can be blended over the surface.
        - A soft feathered edge is applied to the rim for anti-aliasing.

    Args:
        size: The output square image size in pixels (width = height = size).
        sun_dir_3d: A unit 3D vector pointing towards the Sun in the Moon's local frame.
        view_dir_3d: A unit 3D vector pointing towards the viewer in the same frame.
        moon_color: The RGB color for the sunlit surface.
        dark_color: The base RGB color for the unlit side (used for earthshine).
        earthshine_factor: A scalar multiplier for `dark_color` on the unlit side.
        tint_color: An optional RGBA color to overlay on the final image.
        edge_soft_px: The width of the feathering (in pixels) for the rim's alpha.

    Returns:
        A Pillow `Image` object in "RGBA" mode.
    """
    cx = cy = size // 2
    r = size // 2
    img = np.zeros((size, size, 4), dtype=np.uint8)  # RGBA buffer

    # Pre-cast colors and parameters for performance
    moon_rgb = np.array(moon_color, dtype=np.float32)
    dark_rgb = np.array(dark_color, dtype=np.float32)
    use_tint = tint_color is not None
    if use_tint:
        tr, tg, tb, ta = tint_color
        tint_rgb = np.array([tr, tg, tb], dtype=np.float32)
        tint_a = float(ta) / 255.0

    # Iterate over each pixel in the output image
    for y in range(size):
        dy_px = y - cy
        for x in range(size):
            dx_px = x - cx

            # Convert pixel coordinates to normalized disc coordinates (-1.0 to 1.0)
            dx = dx_px / r
            dy = dy_px / r
            d2 = dx * dx + dy * dy

            # Skip pixels outside the circular disc of the moon
            if d2 > 1.0:
                continue

            # Project the 2D disc coordinate to a 3D surface normal on a unit sphere.
            # The z-component points towards the viewer.
            dz = math.sqrt(max(0.0, 1.0 - d2))
            surf_normal = np.array([dx, dy, dz], dtype=np.float32)

            # Check if this part of the surface is visible to the viewer.
            # This culls the back side of the moon.
            view_dot = float(np.dot(surf_normal, view_dir_3d))
            if view_dot <= 0.0:
                continue

            # --- Shading Calculation ---
            # Calculate Lambertian shading based on the angle to the sun.
            light = float(np.dot(surf_normal, sun_dir_3d))
            if light > 0.0:
                # Sunlit side
                c = moon_rgb * light
            else:
                # Dark side (with earthshine)
                c = dark_rgb * earthshine_factor

            # --- Tinting ---
            # Apply an optional color tint (e.g., for a lunar eclipse).
            if use_tint and tint_a > 0.0:
                c = (1.0 - tint_a) * c + tint_a * tint_rgb

            # --- Anti-aliasing ---
            # Feather the alpha channel at the moon's rim to reduce jagged edges.
            if edge_soft_px > 0.0:
                dist_to_rim_px = (1.0 - math.sqrt(d2)) * r
                alpha_factor = np.clip(dist_to_rim_px / edge_soft_px, 0.0, 1.0)
                a = 255.0 * alpha_factor
            else:
                a = 255.0

            # Write the final RGBA value to the image buffer, clamping values.
            img[y, x, :3] = np.clip(c, 0.0, 255.0).astype(np.uint8)
            img[y, x, 3] = int(a)

    return Image.fromarray(img, mode="RGBA")
