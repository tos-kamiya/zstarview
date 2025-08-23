from typing import Optional, Tuple

import math
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
    Generate a spherical moon-phase image as an RGBA bitmap.

    - Pixels outside the lunar disc are fully transparent (alpha = 0).
    - The visible hemisphere only (dot(view, normal) >= 0) is shaded.
    - Lit areas use Lambertian shading with `moon_color`;
      unlit but visible areas use `dark_color * earthshine_factor` (simple earthshine).
    - Optional `tint_color` (RGBA) is alpha-blended on top of the shaded color.
    - A small feather (edge_soft_px in pixels) is applied near the rim to reduce jaggies.

    Args:
        size: Output square image size in pixels (width = height = size).
        sun_dir_3d: Unit 3D vector of sun direction in the moon-local frame.
        view_dir_3d: Unit 3D vector pointing toward the viewer (camera) in the same frame.
        moon_color: RGB color for the lit surface.
        dark_color: RGB base for the unlit side (earthshine).
        earthshine_factor: Scalar multiplier for `dark_color` on the unlit side.
        tint_color: Optional RGBA to overlay via alpha blending per pixel.
        edge_soft_px: Feather width (in pixels) for the rim anti-alias alpha.

    Returns:
        A Pillow `Image` in mode "RGBA".
    """
    cx = cy = size // 2
    r = size // 2
    img = np.zeros((size, size, 4), dtype=np.uint8)  # RGBA

    # Pre-cast colors to arrays
    moon_rgb = np.array(moon_color, dtype=np.float32)
    dark_rgb = np.array(dark_color, dtype=np.float32)
    use_tint = tint_color is not None
    if use_tint:
        tr, tg, tb, ta = tint_color
        tint_rgb = np.array([tr, tg, tb], dtype=np.float32)
        tint_a = float(ta) / 255.0

    # Iterate pixels; only fill those inside the unit disc
    for y in range(size):
        dy_px = y - cy
        for x in range(size):
            dx_px = x - cx

            # Normalized disc coordinates (-1..1)
            dx = dx_px / r
            dy = dy_px / r
            d2 = dx * dx + dy * dy
            if d2 > 1.0:
                # Outside the lunar disc: fully transparent
                continue

            # Surface normal z > 0 hemisphere (sphere projection)
            dz = math.sqrt(max(0.0, 1.0 - d2))
            surf = np.array([dx, dy, dz], dtype=np.float32)
            # view-facing? (if not, keep alpha 0)
            view_dot = float(np.dot(surf, view_dir_3d))
            if view_dot <= 0.0:
                continue

            # Lambertian lighting from sun
            light = float(np.dot(surf, sun_dir_3d))
            if light > 0.0:
                c = moon_rgb * light
            else:
                c = dark_rgb * earthshine_factor

            # Optional tint (simple over operator on RGB)
            if use_tint and tint_a > 0.0:
                c = (1.0 - tint_a) * c + tint_a * tint_rgb

            # Soft alpha near the rim (anti-aliasing)
            if edge_soft_px > 0.0:
                # Distance to rim in pixels
                dist_to_rim_px = (1.0 - math.sqrt(d2)) * r
                a = 255.0 * np.clip(dist_to_rim_px / edge_soft_px, 0.0, 1.0)
            else:
                a = 255.0

            # Write RGBA (clamp)
            img[y, x, :3] = np.clip(c, 0.0, 255.0).astype(np.uint8)
            img[y, x, 3] = int(a)

    return Image.fromarray(img, mode="RGBA")
