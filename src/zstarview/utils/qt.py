import numpy as np
from PIL import Image
from PySide6.QtGui import QImage, QPixmap


def pil_to_qimage(img: Image.Image, premultiplied: bool = True) -> QImage:
    """Convert a PIL Image to QImage. If premultiplied=True, return ARGB32_Premultiplied."""
    if img.mode != "RGBA":
        img = img.convert("RGBA")
    w, h = img.size
    buf = img.tobytes("raw", "RGBA")  # RGBA (straight alpha)

    # まずは RGBA8888（非プレマルチ）として QImage を作る（この時点ではバッファ非所有）
    qimg = QImage(buf, w, h, w * 4, QImage.Format.Format_RGBA8888)

    if premultiplied:
        # 並び変換 + プレマルチ化 + 所有権確保 を一括でやってくれる
        qimg = qimg.convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
    else:
        # 非プレマルチのまま使う場合は deep copy して所有権を Qt 側に渡す
        qimg = qimg.copy()

    return qimg

def pil2qpixmap(img: Image.Image, premultiplied: bool = True) -> QPixmap:
    """Convert a PIL Image to QPixmap. premultiplied=True is recommended for smooth edges."""
    return QPixmap.fromImage(pil_to_qimage(img, premultiplied=premultiplied))


def qimage_to_pil(qimg: QImage) -> Image.Image:
    """QImage → Pillow Image (RGBA)"""
    qimg = qimg.convertToFormat(QImage.Format.Format_RGBA8888)
    w, h = qimg.width(), qimg.height()
    ptr = qimg.bits()
    arr = np.frombuffer(ptr, np.uint8).reshape((h, w, 4)).copy()
    return Image.fromarray(arr, "RGBA")
