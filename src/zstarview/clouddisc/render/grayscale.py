"""
Renders brightness temperature data into grayscale images.
"""

import warnings

import numpy as np
from PIL import Image

# Suppress a specific, harmless warning that can occur deep inside the Satpy/Dask execution path.
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in log",
    category=RuntimeWarning,
    module="dask._task_spec",
)


def _bt_to_weight(bt: np.ndarray, bt_warm: float, bt_cold: float) -> np.ndarray:
    """
    Converts brightness temperature (BT) to a normalized weight (0.0 to 1.0).

    The weight represents cloudiness, where 1.0 is coldest (thick clouds)
    and 0.0 is warmest (clear sky).

    Args:
        bt: Numpy array of brightness temperatures in Kelvin.
        bt_warm: The temperature (K) to map to weight 0.0.
        bt_cold: The temperature (K) to map to weight 1.0.

    Returns:
        A numpy array of normalized weights (0.0 to 1.0).
    """
    # Suppress division-by-zero warnings for this calculation
    with np.errstate(invalid="ignore"):
        # w = 0..1 (warm -> 0, cold -> 1)
        weight = (bt_warm - bt) / max(1e-6, (bt_warm - bt_cold))

    # Clean up the result
    weight = np.clip(weight, 0.0, 1.0)
    weight = np.nan_to_num(weight, nan=0.0, posinf=1.0, neginf=0.0)
    return weight


def convert_bt_to_la_image(
    bt: np.ndarray,
    mask_inside: np.ndarray,
    bt_warm: float,
    bt_cold: float,
) -> Image.Image:
    """
    Converts a brightness temperature (BT) array to a grayscale-alpha (LA) image.

    This function can operate in two modes:
    1. Standard mode: The grayscale value (L) represents cloud brightness, and the
       alpha channel (A) is a hard mask for the horizon/FOV.
    2. `brightness_as_alpha` mode: The grayscale value (L) is solid white, and the
       alpha channel (A) represents the cloud brightness. This is useful for
       compositing the cloud layer onto other images.

    Args:
        bt: Numpy array of brightness temperatures in Kelvin.
        mask_inside: A boolean numpy array where True indicates pixels inside the
                     visible disc (above horizon and within FOV).
        bt_warm: The temperature (K) to map to black (0).
        bt_cold: The temperature (K) to map to white (255).

    Returns:
        A PIL Image object in "LA" mode.
    """

    if bt_warm <= bt_cold + 0.5:
        mid = 0.5 * (bt_warm + bt_cold)
        bt_cold, bt_warm = mid - 0.25, mid + 0.25

    # 1) Calculate base brightness (0-255) from temperature
    weight = _bt_to_weight(bt, bt_warm, bt_cold)
    l_gray = (weight * 255.0).astype(np.uint8)

    # 2) Create a hard mask for the alpha channel (inside=255, outside=0)
    a_mask = np.zeros_like(l_gray, dtype=np.uint8)
    a_mask[mask_inside] = 255

    l_channel = l_gray
    a_channel = a_mask
    la_data = np.dstack([l_channel, a_channel])

    return Image.fromarray(la_data, mode="LA")
