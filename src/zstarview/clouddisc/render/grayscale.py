# -*- coding: utf-8 -*-
"""
Functions for rendering brightness temperature data into grayscale images.

This module handles the final conversion from a numerical array of brightness
temperatures (in Kelvin) into a visual, 8-bit grayscale image suitable for display.
"""

import warnings

import numpy as np
from PIL import Image

# Suppress a specific, harmless warning that can occur deep inside the Satpy/Dask
# execution path when dealing with NaN values in logarithms.
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in log",
    category=RuntimeWarning,
    module="dask._task_spec",
)


def _bt_to_weight(bt: np.ndarray, bt_warm: float, bt_cold: float) -> np.ndarray:
    """
    Converts brightness temperature (BT) to a normalized weight from 0.0 to 1.0.

    The weight represents cloudiness, where 1.0 corresponds to the coldest
    (and therefore highest/thickest) clouds, and 0.0 corresponds to the warmest
    (clear sky or low clouds).

    Args:
        bt: A NumPy array of brightness temperatures in Kelvin.
        bt_warm: The temperature (K) to map to weight 0.0.
        bt_cold: The temperature (K) to map to weight 1.0.

    Returns:
        A NumPy array of normalized weights, clipped to the [0, 1] range.
    """
    # Suppress "invalid value" warnings for this calculation, as we handle NaNs later.
    with np.errstate(invalid="ignore"):
        # Linearly map the temperature range [bt_warm, bt_cold] to [0, 1].
        weight = (bt_warm - bt) / max(1e-6, (bt_warm - bt_cold))

    # Clean up the result by clipping to the valid range and replacing any NaNs.
    weight = np.clip(weight, 0.0, 1.0)
    weight = np.nan_to_num(weight, nan=0.0, posinf=1.0, neginf=0.0)
    return weight


def _smoothstep(edge0: float, edge1: float, x: np.ndarray) -> np.ndarray:
    """Vectorized smoothstep for gently ramping values in [edge0, edge1]."""
    t = np.clip((x - edge0) / max(1e-6, (edge1 - edge0)), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _suppress_low_cloud_weight(
    weight: np.ndarray,
    *,
    floor: float = 0.03,
    knee: float = 0.08,
) -> np.ndarray:
    """Suppress low cloud amounts so faint clear-sky stripes are less visible."""
    floor = float(np.clip(floor, 0.0, 0.95))
    knee = float(np.clip(max(knee, floor + 1e-6), floor + 1e-6, 0.99))

    # 1) Hard floor to remove tiny cloud/noise values.
    base = np.clip((weight - floor) / max(1e-6, (1.0 - floor)), 0.0, 1.0)
    # 2) Smooth gain near the knee so faint clouds do not pop abruptly.
    gain = _smoothstep(floor, knee, weight)
    return np.clip(base * gain, 0.0, 1.0)


def convert_bt_to_la_image(
    bt: np.ndarray,
    mask_inside: np.ndarray,
    bt_warm: float,
    bt_cold: float,
) -> Image.Image:
    """
    Converts a brightness temperature (BT) array to a grayscale-alpha (LA) image.

    The output is an "LA" mode Pillow Image, which has two 8-bit channels:
    - L: Luminance (grayscale value), fixed to white for visible pixels.
    - A: Alpha (transparency), representing cloud amount (0=clear, 255=thick cloud).

    Args:
        bt: A NumPy array of brightness temperatures in Kelvin.
        mask_inside: A boolean NumPy array where True indicates pixels inside the
                     visible disc (i.e., above the horizon and within the FOV).
        bt_warm: The temperature (K) mapped to alpha 0.
        bt_cold: The temperature (K) mapped to alpha 255.

    Returns:
        A Pillow Image object in "LA" mode.
    """
    # Ensure there is a valid temperature range to prevent division by zero.
    if bt_warm <= bt_cold + 0.5:
        mid = 0.5 * (bt_warm + bt_cold)
        bt_cold, bt_warm = mid - 0.25, mid + 0.25

    # 1. Convert BT to cloud amount in [0,1].
    weight = _bt_to_weight(bt, bt_warm, bt_cold)
    weight = _suppress_low_cloud_weight(weight)

    # 2. Keep cloud color white and encode cloud amount in alpha.
    l_channel = np.zeros_like(weight, dtype=np.uint8)
    l_channel[mask_inside] = 255
    a_channel = np.zeros_like(weight, dtype=np.uint8)
    a_channel[mask_inside] = (weight[mask_inside] * 255.0).astype(np.uint8)

    # 3. Stack the L and A channels together to create the final LA image data.
    la_data = np.dstack([l_channel, a_channel])

    return Image.fromarray(la_data, mode="LA")
