from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

from ..paths import GEOSATELLITE_GRAY_COMMON_MASK_FILE

DEFAULT_GRAY_SPREAD = 0.33
DEFAULT_WHITE_BRIGHTNESS = 0.58
DEFAULT_MASK_PATH = Path(GEOSATELLITE_GRAY_COMMON_MASK_FILE)
DEFAULT_MASK_THRESHOLD = 127.0
DEFAULT_RADIUS = 1
DEFAULT_MAX_ITERATIONS = 256


def _rgb_features(rgb: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rgb = np.asarray(rgb, dtype=np.float32) / 255.0
    max_rgb = np.max(rgb, axis=2)
    min_rgb = np.min(rgb, axis=2)
    spread = max_rgb - min_rgb
    brightness = 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]
    return spread.astype(np.float32), brightness.astype(np.float32)


def build_common_mask(
    images: list[Image.Image],
    *,
    gray_spread: float = DEFAULT_GRAY_SPREAD,
    white_brightness: float = DEFAULT_WHITE_BRIGHTNESS,
) -> Image.Image:
    arrays = [np.asarray(image.convert("RGB"), dtype=np.uint8) for image in images]
    shape = arrays[0].shape
    for array in arrays[1:]:
        if array.shape != shape:
            raise ValueError("All input images must have the same dimensions.")

    common = np.ones(shape[:2], dtype=bool)
    for array in arrays:
        spread, brightness = _rgb_features(array)
        keep = (spread <= float(gray_spread)) & (brightness <= float(white_brightness))
        common &= keep

    mask = np.where(common, 255, 0).astype(np.uint8)
    return Image.fromarray(mask, mode="L")


def load_default_mask() -> Image.Image:
    with Image.open(DEFAULT_MASK_PATH) as image:
        return image.convert("L")


def _accumulate_neighborhood(values: np.ndarray, valid: np.ndarray, radius: int) -> tuple[np.ndarray, np.ndarray]:
    if radius < 1:
        raise ValueError("radius must be >= 1")

    padded_values = np.pad(values, radius, mode="constant", constant_values=0.0)
    padded_valid = np.pad(valid.astype(np.float32), radius, mode="constant", constant_values=0.0)
    height, width = values.shape
    sum_acc = np.zeros((height, width), dtype=np.float32)
    count_acc = np.zeros((height, width), dtype=np.float32)

    for dy in range(-radius, radius + 1):
        y0 = radius + dy
        y1 = y0 + height
        for dx in range(-radius, radius + 1):
            x0 = radius + dx
            x1 = x0 + width
            sum_acc += padded_values[y0:y1, x0:x1] * padded_valid[y0:y1, x0:x1]
            count_acc += padded_valid[y0:y1, x0:x1]

    return sum_acc, count_acc


def fill_masked_regions(
    image: np.ndarray,
    mask: np.ndarray,
    *,
    radius: int = DEFAULT_RADIUS,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
) -> np.ndarray:
    if image.shape != mask.shape:
        raise ValueError("Image and mask must have the same dimensions.")

    result = image.astype(np.float32, copy=True)
    working_mask = mask.astype(bool, copy=True)
    known = ~working_mask

    for _ in range(max_iterations):
        sum_acc, count_acc = _accumulate_neighborhood(result, known, radius)
        fillable = working_mask & (count_acc > 0.0)
        if not np.any(fillable):
            break
        result[fillable] = sum_acc[fillable] / count_acc[fillable]
        known[fillable] = True
        working_mask = working_mask & ~fillable
        if not np.any(working_mask):
            break

    if np.any(working_mask):
        fallback = float(np.mean(result[known])) if np.any(known) else 0.0
        result[working_mask] = fallback

    return np.clip(result, 0.0, 255.0)
