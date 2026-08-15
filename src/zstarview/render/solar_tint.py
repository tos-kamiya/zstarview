"""Atmospheric RGB colorization for the HMI Sun hover image."""

from __future__ import annotations

from collections import OrderedDict
from functools import lru_cache

import numpy as np
from PySide6.QtGui import QImage

from .atmosphere import direct_solar_transmission_rgb
from .photometry import planet_marker_color
from .qt_image import np_rgba_to_qimage, qimage_to_np_rgba

SOLAR_HOVER_REFERENCE_RGB = np.asarray((1.0, 0.94, 0.80), dtype=np.float32)
SOLAR_HOVER_MARKER_FALLBACK_STRENGTH = 0.45
_COLOR_CACHE_MAXSIZE = 8
_COLOR_ALTITUDE_STEP_DEG = 0.25
_COLOR_AOD_STEP = 0.01
_COLOR_HEIGHT_STEP_M = 10.0
_SOLAR_COLOR_CACHE: OrderedDict[tuple[int, int, int, int], QImage] = OrderedDict()


def _srgb_to_linear(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    return np.where(
        values <= 0.04045,
        values / 12.92,
        np.power((values + 0.055) / 1.055, 2.4),
    )


def _linear_to_srgb(values: np.ndarray) -> np.ndarray:
    values = np.clip(np.asarray(values, dtype=np.float32), 0.0, 1.0)
    return np.where(
        values <= 0.0031308,
        values * 12.92,
        1.055 * np.power(values, 1.0 / 2.4) - 0.055,
    )


def _quantize(value: float, step: float) -> int:
    return int(round(float(value) / step))


def _subdued_solar_marker_multiplier() -> np.ndarray:
    marker = planet_marker_color("sun")
    marker_rgb = np.asarray(marker.getRgb()[:3], dtype=np.float32) / 255.0
    return _srgb_to_linear(marker_rgb) * SOLAR_HOVER_MARKER_FALLBACK_STRENGTH


@lru_cache(maxsize=32)
def solar_hover_color_multiplier(
    sun_altitude_q: int,
    observer_height_q: int,
    aerosol_optical_depth_q: int,
) -> tuple[float, float, float]:
    """Return a cached linear-RGB multiplier for quantized inputs."""
    transmission = direct_solar_transmission_rgb(
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
    rgba = qimage_to_np_rgba(image)
    alpha = rgba[..., 3].copy()
    linear_rgb = _srgb_to_linear(rgba[..., :3].astype(np.float32) / 255.0)
    tinted = _linear_to_srgb(linear_rgb * multiplier[None, None, :])
    output = np.empty_like(rgba)
    output[..., :3] = np.clip(np.rint(tinted * 255.0), 0.0, 255.0).astype(np.uint8)
    output[..., 3] = alpha
    tinted_image = np_rgba_to_qimage(output).convertToFormat(
        QImage.Format.Format_ARGB32_Premultiplied
    )
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
