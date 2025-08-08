from __future__ import annotations

import numpy as np
from PIL import Image
from PyQt5.QtGui import QImage, QPixmap


def pil2qpixmap(img: Image.Image) -> QPixmap:
    """Convert a PIL Image to a PyQt5 QPixmap (RGBA)."""
    arr = np.array(img.convert("RGBA"))
    h, w, ch = arr.shape
    bytes_per_line = ch * w
    qimg = QImage(arr.data.tobytes(), w, h, bytes_per_line, QImage.Format.Format_RGBA8888)
    return QPixmap.fromImage(qimg)
