"""Application-level appearance settings for voxel-rendered clouds."""

from __future__ import annotations

import numpy as np

# Ambient night cloud color at full cloud amount, in normalized RGB.
CLOUD_VOXEL_NIGHT_COLOR_RGB = (0.36, 0.385, 0.44)

# Voxel cloud images use a 513x513 source surface (256px radius).
CLOUD_VOXEL_CUTOUT_REFERENCE_SIZE = 513
CLOUD_VOXEL_CUTOUT_WIDTH_PX = 1
CLOUD_VOXEL_CUTOUT_COUNT = 10


def apply_cloud_voxel_cutout(
    cloud_rgba: np.ndarray,
    *,
    count: int = CLOUD_VOXEL_CUTOUT_COUNT,
) -> np.ndarray:
    """Return a voxel cloud raster with sparse, 45-degree cutout lines."""
    image = np.asarray(cloud_rgba)
    if image.ndim != 3 or image.shape[2] != 4 or int(count) <= 0:
        return image

    height, width = image.shape[:2]
    if height < 3 or width < 1:
        return image

    line_count = min(int(count), height + width - 2)
    line_width = max(
        1,
        int(
            round(
                CLOUD_VOXEL_CUTOUT_WIDTH_PX
                * min(height, width)
                / CLOUD_VOXEL_CUTOUT_REFERENCE_SIZE
            )
        ),
    )
    center_x = (width - 1) * 0.5
    center_y = (height - 1) * 0.5
    x = np.arange(width, dtype=np.float32)[None, :] - center_x
    y = np.arange(height, dtype=np.float32)[:, None] - center_y
    diagonal_coord = x - y
    max_offset = 0.5 * (min(height, width) - 1) * np.sqrt(2.0)
    line_offsets = np.rint(
        np.linspace(-max_offset, max_offset, line_count + 2)[1:-1]
    )
    half_width_in_diagonal = line_width / np.sqrt(2.0)
    cutout_mask = np.zeros((height, width), dtype=bool)
    for offset in line_offsets:
        cutout_mask |= np.abs(diagonal_coord - offset) <= half_width_in_diagonal

    output = image.copy()
    output[cutout_mask, 3] = 0
    return output
