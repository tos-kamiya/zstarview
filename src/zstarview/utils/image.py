from __future__ import annotations

import numpy as np
from PIL import Image
from typing import Tuple


def generate_moon_phase_image(
    size: int,
    sun_dir_3d: np.ndarray,
    view_dir_3d: np.ndarray,
    moon_color: Tuple[int, int, int] = (230, 230, 230),
    dark_color: Tuple[int, int, int] = (30, 30, 30),
    earthshine_factor: float = 0.15,
) -> Image.Image:
    """Generate a spherical moon phase image with simple earthshine.

    Returns a PIL Image (RGB) of dimensions (size, size).
    """
    cx, cy = size // 2, size // 2
    r = size // 2
    img_array = np.zeros((size, size, 3), dtype=np.uint8)
    for y in range(size):
        for x in range(size):
            dx = (x - cx) / r
            dy = (y - cy) / r
            if dx**2 + dy**2 > 1:
                continue
            dz = np.sqrt(1 - dx**2 - dy**2)
            surf = np.array([dx, dy, dz])
            surf /= np.linalg.norm(surf)
            view = np.dot(surf, view_dir_3d)
            if view < 0:
                continue  # back side

            light = np.dot(surf, sun_dir_3d)
            if light > 0:
                c = np.array(moon_color) * light
            else:
                c = np.array(dark_color) * earthshine_factor
            c = np.clip(c, 0, 255)
            img_array[y, x] = c
    return Image.fromarray(img_array)
