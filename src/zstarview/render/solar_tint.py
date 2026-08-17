"""Atmospheric RGB colorization for the HMI Sun hover image."""

from __future__ import annotations

from collections import OrderedDict
from functools import lru_cache

import numpy as np
from PySide6.QtGui import QImage

from .atmosphere import direct_atmospheric_transmission_rgb
from .atmospheric_tint import multiply_qimage_linear_rgb, srgb_to_linear
from .photometry import planet_marker_color

SOLAR_HOVER_REFERENCE_RGB = np.asarray((1.0, 0.94, 0.80), dtype=np.float32)
SOLAR_HOVER_MARKER_FALLBACK_STRENGTH = 0.45
_COLOR_CACHE_MAXSIZE = 8
_COLOR_ALTITUDE_STEP_DEG = 0.25
_COLOR_AOD_STEP = 0.01
_COLOR_HEIGHT_STEP_M = 10.0
_SOLAR_COLOR_CACHE: OrderedDict[tuple[int, int, int, int], QImage] = OrderedDict()


def _quantize(value: float, step: float) -> int:
    return int(round(float(value) / step))


def _subdued_solar_marker_multiplier() -> np.ndarray:
    marker = planet_marker_color("sun")
    marker_rgb = np.asarray(
        (marker.red(), marker.green(), marker.blue()), dtype=np.float32
    ) / 255.0
    return srgb_to_linear(marker_rgb) * SOLAR_HOVER_MARKER_FALLBACK_STRENGTH


@lru_cache(maxsize=32)
def solar_hover_color_multiplier(
    sun_altitude_q: int,
    observer_height_q: int,
    aerosol_optical_depth_q: int,
) -> tuple[float, float, float]:
    """Return a cached linear-RGB multiplier for quantized inputs."""
    transmission = direct_atmospheric_transmission_rgb(
        sun_altitude_q * _COLOR_ALTITUDE_STEP_DEG,
        observer_height_m=observer_height_q * _COLOR_HEIGHT_STEP_M,
        aerosol_optical_depth=aerosol_optical_depth_q * _COLOR_AOD_STEP,
    )
    atmospheric_multiplier = SOLAR_HOVER_REFERENCE_RGB * transmission
    multiplier = np.maximum(
        atmospheric_multiplier,
        _subdued_solar_marker_multiplier(),
    )
    return (float(multiplier[0]), float(multiplier[1]), float(multiplier[2]))


def colorize_solar_hover_image(
    image: QImage,
    *,
    image_id: int,
    sun_alt_deg: float,
    observer_height_m: float,
    aerosol_optical_depth: float,
) -> QImage:
    """Apply the atmospheric solar tint while preserving the input alpha."""
    key = (
        int(image_id),
        _quantize(sun_alt_deg, _COLOR_ALTITUDE_STEP_DEG),
        _quantize(observer_height_m, _COLOR_HEIGHT_STEP_M),
        _quantize(aerosol_optical_depth, _COLOR_AOD_STEP),
    )
    cached = _SOLAR_COLOR_CACHE.get(key)
    if cached is not None:
        _SOLAR_COLOR_CACHE.move_to_end(key)
        return cached

    multiplier = np.asarray(
        solar_hover_color_multiplier(key[1], key[2], key[3]),
        dtype=np.float32,
    )
    tinted_image = multiply_qimage_linear_rgb(image, multiplier)
    _SOLAR_COLOR_CACHE[key] = tinted_image
    _SOLAR_COLOR_CACHE.move_to_end(key)
    while len(_SOLAR_COLOR_CACHE) > _COLOR_CACHE_MAXSIZE:
        _SOLAR_COLOR_CACHE.popitem(last=False)
    return tinted_image


def clear_solar_hover_color_cache() -> None:
    """Clear colorized hover frames, primarily for tests and profile changes."""
    _SOLAR_COLOR_CACHE.clear()
    solar_hover_color_multiplier.cache_clear()


__all__ = [
    "SOLAR_HOVER_REFERENCE_RGB",
    "SOLAR_HOVER_MARKER_FALLBACK_STRENGTH",
    "clear_solar_hover_color_cache",
    "colorize_solar_hover_image",
    "solar_hover_color_multiplier",
]
