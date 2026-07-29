"""
Qt image conversion helpers used by the main render/runtime paths.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage


def qimage_to_np_rgba(qimg: QImage) -> np.ndarray:
    """
    Converts a QImage to a NumPy array (H, W, 4) in RGBA8888 format.

    The returned array is a deep copy, ensuring it is independent of the
    lifetime of the source QImage.
    """
    qimg = qimg.convertToFormat(QImage.Format.Format_RGBA8888)
    w, h = qimg.width(), qimg.height()
    bpl = qimg.bytesPerLine()
    buf = bytes(qimg.constBits())
    arr = np.frombuffer(buf, dtype=np.uint8).reshape(h, bpl)
    return arr[:, : w * 4].reshape(h, w, 4).copy()


def np_rgba_to_qimage(arr: np.ndarray) -> QImage:
    """
    Converts a NumPy RGBA array (H, W, 4) to a QImage.
    """
    arr = np.asarray(arr)
    if arr.ndim != 3 or arr.shape[2] != 4:
        raise ValueError("Input must be a (H, W, 4) array")
    if arr.dtype != np.uint8:
        raise ValueError("Input must have dtype uint8")
    if not arr.flags.c_contiguous:
        arr = np.ascontiguousarray(arr)
    h, w, _ = arr.shape
    qimg = QImage(arr.data, w, h, w * 4, QImage.Format_RGBA8888)
    return qimg.copy()
