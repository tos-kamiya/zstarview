# -*- coding: utf-8 -*-
"""
Pillow-backed conversion helpers for Qt images and pixmaps.

Pillow-backed helpers are kept here for tooling compatibility. NumPy/QImage
helpers live in ``zstarview.render.qt_image`` so the main runtime path no
longer depends on this module.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtGui import QImage, QPixmap

if TYPE_CHECKING:
    from PIL import Image


def pil_to_qimage(img: "Image.Image", premultiplied: bool = True) -> QImage:
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


def _qimage_to_pil(qimg: QImage) -> "Image.Image":
    """
    Converts a QImage to a Pillow (PIL) Image in RGBA format.

    Args:
        qimg: The input QImage.

    Returns:
        A Pillow Image in "RGBA" mode.
    """
    from PIL import Image

    # Ensure the QImage is in a standard, non-premultiplied format.
    qimg = qimg.convertToFormat(QImage.Format.Format_RGBA8888)
    w, h = qimg.width(), qimg.height()

    import numpy as np

    ptr = qimg.bits()
    arr = np.frombuffer(ptr, np.uint8).reshape((h, w, 4))

    # Create a Pillow Image from the NumPy array. A copy is made.
    return Image.fromarray(arr, "RGBA")


def pil2qpixmap(img: "Image.Image", premultiplied: bool = True) -> QPixmap:
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
