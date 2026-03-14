# -*- coding: utf-8 -*-
"""
Utilities for converting between Qt and other common image formats.

This module provides helper functions to seamlessly convert between QImage/QPixmap,
NumPy arrays, and Pillow (PIL) Images. This is crucial for integrating data from
image processing libraries with the Qt-based user interface.
"""

import numpy as np
from PIL import Image
from PySide6.QtGui import QImage, QPixmap


def pil_to_qimage(img: Image.Image, premultiplied: bool = True) -> QImage:
    """
    Converts a Pillow (PIL) Image to a QImage.

    Args:
        img: The input Pillow image.
        premultiplied: If True, the output QImage will be in
            `ARGB32_Premultiplied` format. This is recommended for better
            alpha blending performance and quality in Qt. If False, the
            output will be `RGBA8888` (straight alpha).

    Returns:
        The converted QImage instance.
    """
    # Ensure the image is in RGBA format for consistency
    if img.mode != "RGBA":
        img = img.convert("RGBA")
    w, h = img.size
    buf = img.tobytes("raw", "RGBA")

    # Create a QImage from the buffer. The format is RGBA8888 (non-premultiplied).
    # Note: At this stage, QImage does not own the buffer `buf`.
    qimg = QImage(buf, w, h, w * 4, QImage.Format.Format_RGBA8888)

    if premultiplied:
        # This conversion creates a new QImage with premultiplied alpha, and Qt
        # takes ownership of the new memory.
        return qimg.convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
    else:
        # If not premultiplying, we must explicitly copy the data to ensure
        # the new QImage owns its memory, otherwise it could point to a garbage
        # collected buffer.
        return qimg.copy()


def _qimage_to_pil(qimg: QImage) -> Image.Image:
    """
    Converts a QImage to a Pillow (PIL) Image in RGBA format.

    Args:
        qimg: The input QImage.

    Returns:
        A Pillow Image in "RGBA" mode.
    """
    # Ensure the QImage is in a standard, non-premultiplied format.
    qimg = qimg.convertToFormat(QImage.Format.Format_RGBA8888)
    w, h = qimg.width(), qimg.height()

    # Get a pointer to the image data and create a NumPy array view.
    ptr = qimg.bits()
    arr = np.frombuffer(ptr, np.uint8).reshape((h, w, 4))

    # Create a Pillow Image from the NumPy array. A copy is made.
    return Image.fromarray(arr, "RGBA")


def pil2qpixmap(img: Image.Image, premultiplied: bool = True) -> QPixmap:
    """
    Convenience function to convert a Pillow (PIL) Image directly to a QPixmap.

    Args:
        img: The input Pillow image.
        premultiplied: If True, the intermediate QImage is created with
            premultiplied alpha for better rendering quality.

    Returns:
        The converted QPixmap object.
    """
    return QPixmap.fromImage(pil_to_qimage(img, premultiplied=premultiplied))


def qimage_to_np_rgba(qimg: QImage) -> np.ndarray:
    """
    Converts a QImage to a NumPy array (H, W, 4) in RGBA8888 format.

    The returned array is a deep copy, ensuring it is independent of the
    lifetime of the source QImage.

    Args:
        qimg: The input QImage.

    Returns:
        A NumPy array of shape (height, width, 4) with dtype=uint8,
        representing the RGBA pixel data.
    """
    # Ensure a consistent, non-premultiplied input format.
    qimg = qimg.convertToFormat(QImage.Format.Format_RGBA8888)
    w, h = qimg.width(), qimg.height()
    bpl = qimg.bytesPerLine()

    # In PySide6, bits() returns a memoryview. Convert to bytes for frombuffer.
    buf = bytes(qimg.constBits())
    arr = np.frombuffer(buf, dtype=np.uint8).reshape(h, bpl)

    # Remove any extra padding at the end of each scanline and ensure the
    # final array has the correct shape. A copy is made here.
    return arr[:, : w * 4].reshape(h, w, 4).copy()


def np_rgba_to_qimage(arr: np.ndarray) -> QImage:
    """
    Converts a NumPy RGBA array (H, W, 4) to a QImage.

    Args:
        arr: An input NumPy array of shape (height, width, 4) and dtype=uint8.

    Returns:
        A QImage in RGBA8888 format. The buffer is copied, so the QImage
        owns its own memory.
    """
    h, w, c = arr.shape
    assert c == 4 and arr.dtype == np.uint8, "Input must be a (H, W, 4) uint8 array"

    # Create a QImage that views the NumPy array's data.
    qimg = QImage(arr.data, w, h, w * 4, QImage.Format_RGBA8888)

    # Return a deep copy to ensure the QImage has its own memory, independent
    # of the NumPy array's lifetime.
    return qimg.copy()
