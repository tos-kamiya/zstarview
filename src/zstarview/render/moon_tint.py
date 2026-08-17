"""Atmospheric RGB colorization for NASA Moon hover images."""

from __future__ import annotations

from collections import OrderedDict
from functools import lru_cache

import numpy as np
from PySide6.QtGui import QImage

from .atmosphere import direct_atmospheric_transmission_rgb
from .atmospheric_tint import multiply_qimage_linear_rgb, srgb_to_linear
from .photometry import planet_marker_color

MOON_HOVER_MARKER_FALLBACK_STRENGTH = 0.35
_COLOR_CACHE_MAXSIZE = 8
_COLOR_ALTITUDE_STEP_DEG = 0.25
_COLOR_AOD_STEP = 0.01
_COLOR_HEIGHT_STEP_M = 10.0
_MOON_COLOR_CACHE: OrderedDict[tuple[int, int, int, int], QImage] = OrderedDict()


def _quantize(value: float, step: float) -> int:
    return int(round(float(value) / step))


def _subdued_moon_marker_multiplier() -> np.ndarray:
    marker = planet_marker_color("moon")
    marker_rgb = np.asarray(
        (marker.red(), marker.green(), marker.blue()), dtype=np.float32
    ) / 255.0
    return srgb_to_linear(marker_rgb) * MOON_HOVER_MARKER_FALLBACK_STRENGTH


@lru_cache(maxsize=32)
def moon_hover_color_multiplier(
    moon_altitude_q: int,
    observer_height_q: int,
    aerosol_optical_depth_q: int,
) -> tuple[float, float, float]:
    """Return a cached linear-RGB multiplier for quantized lunar conditions."""
    transmission = direct_atmospheric_transmission_rgb(
        moon_altitude_q * _COLOR_ALTITUDE_STEP_DEG,
        observer_height_m=observer_height_q * _COLOR_HEIGHT_STEP_M,
        aerosol_optical_depth=aerosol_optical_depth_q * _COLOR_AOD_STEP,
    )
    multiplier = np.maximum(transmission, _subdued_moon_marker_multiplier())
    return (float(multiplier[0]), float(multiplier[1]), float(multiplier[2]))


def colorize_moon_hover_image(
    image: QImage,
    *,
    moon_alt_deg: float,
    observer_height_m: float,
    aerosol_optical_depth: float,
) -> QImage:
    """Apply observer-side atmospheric transmission to a Moon hover image."""
    key = (
        int(image.cacheKey()),
        _quantize(moon_alt_deg, _COLOR_ALTITUDE_STEP_DEG),
        _quantize(observer_height_m, _COLOR_HEIGHT_STEP_M),
        _quantize(aerosol_optical_depth, _COLOR_AOD_STEP),
    )
    cached = _MOON_COLOR_CACHE.get(key)
    if cached is not None:
        _MOON_COLOR_CACHE.move_to_end(key)
        return cached

    multiplier = np.asarray(
        moon_hover_color_multiplier(key[1], key[2], key[3]),
        dtype=np.float32,
    )
    tinted_image = multiply_qimage_linear_rgb(image, multiplier)
    _MOON_COLOR_CACHE[key] = tinted_image
    _MOON_COLOR_CACHE.move_to_end(key)
    while len(_MOON_COLOR_CACHE) > _COLOR_CACHE_MAXSIZE:
        _MOON_COLOR_CACHE.popitem(last=False)
    return tinted_image


def clear_moon_hover_color_cache() -> None:
    """Clear colorized Moon hover frames, primarily for tests."""
    _MOON_COLOR_CACHE.clear()
    moon_hover_color_multiplier.cache_clear()


__all__ = [
    "MOON_HOVER_MARKER_FALLBACK_STRENGTH",
    "clear_moon_hover_color_cache",
    "colorize_moon_hover_image",
    "moon_hover_color_multiplier",
]
