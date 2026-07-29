"""
Image generation utilities, such as creating celestial body sprites.
"""


import numpy as np


def generate_moon_phase_rgba(
    size: int,
    sun_dir_3d: np.ndarray,
    view_dir_3d: np.ndarray,
    moon_color: tuple[int, int, int] = (230, 230, 230),
    dark_color: tuple[int, int, int] = (30, 30, 30),
    earthshine_factor: float = 0.15,
    tint_color: tuple[int, int, int, int] | None = None,  # RGBA
    edge_soft_px: float = 1.0,
) -> np.ndarray:
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
        A NumPy array in RGBA format with shape ``(size, size, 4)``.
    """
    cx = cy = size // 2
    r = size // 2
    img = np.zeros((size, size, 4), dtype=np.uint8)  # RGBA buffer
    if r <= 0:
        return img

    # Pre-cast colors and vectors for vectorized shading.
    moon_rgb = np.asarray(moon_color, dtype=np.float32)
    dark_rgb = np.asarray(dark_color, dtype=np.float32)
    sun_dir = np.asarray(sun_dir_3d, dtype=np.float32)
    view_dir = np.asarray(view_dir_3d, dtype=np.float32)

    xs = (np.arange(size, dtype=np.float32) - float(cx)) / float(r)
    ys = (np.arange(size, dtype=np.float32) - float(cy)) / float(r)
    dx, dy = np.meshgrid(xs, ys)
    d2 = dx * dx + dy * dy
    disc_mask = d2 <= 1.0

    dz = np.zeros_like(dx, dtype=np.float32)
    dz[disc_mask] = np.sqrt(np.maximum(0.0, 1.0 - d2[disc_mask]))
    normals = np.stack((dx, dy, dz), axis=-1)

    # Culls the back side of the moon.
    view_dot = normals @ view_dir
    visible_mask = disc_mask & (view_dot > 0.0)
    if not np.any(visible_mask):
        return img

    light = normals @ sun_dir
    sunlit_mask = visible_mask & (light > 0.0)
    dark_mask = visible_mask & ~sunlit_mask

    rgb = np.zeros((size, size, 3), dtype=np.float32)
    rgb[sunlit_mask] = moon_rgb[None, :] * light[sunlit_mask, None]
    rgb[dark_mask] = dark_rgb[None, :] * float(earthshine_factor)

    if tint_color is not None:
        tr, tg, tb, ta = tint_color
        tint_a = float(ta) / 255.0
        if tint_a > 0.0:
            tint_rgb = np.asarray([tr, tg, tb], dtype=np.float32)
            rgb[visible_mask] = (1.0 - tint_a) * rgb[visible_mask] + tint_a * tint_rgb[None, :]

    alpha = np.zeros((size, size), dtype=np.float32)
    if edge_soft_px > 0.0:
        dist_to_rim_px = (1.0 - np.sqrt(np.clip(d2, 0.0, 1.0))) * float(r)
        alpha_factor = np.clip(dist_to_rim_px / float(edge_soft_px), 0.0, 1.0)
        alpha[visible_mask] = 255.0 * alpha_factor[visible_mask]
    else:
        alpha[visible_mask] = 255.0

    img[..., :3] = np.clip(rgb, 0.0, 255.0).astype(np.uint8)
    img[..., 3] = np.clip(alpha, 0.0, 255.0).astype(np.uint8)
    return img
