"""Shared linear-RGB helpers for atmospheric body-image tinting."""

from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage

from .qt_image import np_rgba_to_qimage, qimage_to_np_rgba


def srgb_to_linear(values: np.ndarray) -> np.ndarray:
    """Convert normalized sRGB values to linear RGB."""
    values = np.asarray(values, dtype=np.float32)
    return np.where(
        values <= 0.04045,
        values / 12.92,
        np.power((values + 0.055) / 1.055, 2.4),
    )


def linear_to_srgb(values: np.ndarray) -> np.ndarray:
    """Convert normalized linear RGB values to display sRGB."""
    values = np.clip(np.asarray(values, dtype=np.float32), 0.0, 1.0)
    return np.where(
        values <= 0.0031308,
        values * 12.92,
        1.055 * np.power(values, 1.0 / 2.4) - 0.055,
    )


def multiply_qimage_linear_rgb(image: QImage, multiplier: np.ndarray) -> QImage:
    """Multiply a QImage in linear RGB while preserving its alpha channel."""
    rgba = qimage_to_np_rgba(image)
    alpha = rgba[..., 3].copy()
    linear_rgb = srgb_to_linear(rgba[..., :3].astype(np.float32) / 255.0)
    tinted = linear_to_srgb(
        linear_rgb * np.asarray(multiplier, dtype=np.float32)[None, None, :]
    )
    output = np.empty_like(rgba)
    output[..., :3] = np.clip(np.rint(tinted * 255.0), 0.0, 255.0).astype(np.uint8)
    output[..., 3] = alpha
    return np_rgba_to_qimage(output).convertToFormat(
        QImage.Format.Format_ARGB32_Premultiplied
    )


__all__ = ["linear_to_srgb", "multiply_qimage_linear_rgb", "srgb_to_linear"]
