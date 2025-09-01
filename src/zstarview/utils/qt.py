import numpy as np
from PIL import Image
from PySide6.QtGui import QImage, QPixmap


def pil_to_qimage(img: Image.Image, premultiplied: bool = True) -> QImage:
    """
    Convert a Pillow (PIL) Image to a QImage.

    Args:
        img (Image.Image): Input Pillow image.
        premultiplied (bool, optional): If True, the returned QImage will be
            in ARGB32_Premultiplied format (recommended for better blending
            and performance in Qt). If False, the QImage will remain
            in RGBA8888 format.

    Returns:
        QImage: The converted QImage instance.
    """
    if img.mode != "RGBA":
        img = img.convert("RGBA")
    w, h = img.size
    buf = img.tobytes("raw", "RGBA")  # RGBA (straight alpha)

    # Create QImage in RGBA8888 (non-premultiplied).
    # At this point, the buffer is not owned by QImage.
    qimg = QImage(buf, w, h, w * 4, QImage.Format.Format_RGBA8888)

    if premultiplied:
        # Convert to ARGB32_Premultiplied.
        # This ensures proper alpha blending and transfers buffer ownership.
        qimg = qimg.convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
    else:
        # Make a deep copy to ensure QImage owns the buffer.
        qimg = qimg.copy()

    return qimg


def qimage_to_pil(qimg: QImage) -> Image.Image:
    """
    Convert a QImage to a Pillow (PIL) Image in RGBA format.

    Args:
        qimg (QImage): Input QImage.

    Returns:
        Image.Image: Converted Pillow image in RGBA mode.
    """
    qimg = qimg.convertToFormat(QImage.Format.Format_RGBA8888)
    w, h = qimg.width(), qimg.height()
    ptr = qimg.bits()
    arr = np.frombuffer(ptr, np.uint8).reshape((h, w, 4)).copy()
    return Image.fromarray(arr, "RGBA")


def pil2qpixmap(img: Image.Image, premultiplied: bool = True) -> QPixmap:
    """
    Convert a Pillow (PIL) Image to a QPixmap.

    Args:
        img (Image.Image): Input Pillow image.
        premultiplied (bool, optional): If True, conversion is done with
            ARGB32_Premultiplied format, which is recommended for smoother
            edges and proper blending.

    Returns:
        QPixmap: Converted QPixmap object.
    """
    return QPixmap.fromImage(pil_to_qimage(img, premultiplied=premultiplied))


def qimage_to_np_rgba(qimg: QImage) -> np.ndarray:
    """
    Convert a QImage to a NumPy array (H, W, 4) in RGBA8888 format.

    The returned array is a deep copy and independent from the QImage.

    Args:
        qimg (QImage): Input QImage.

    Returns:
        np.ndarray: A NumPy array of shape (height, width, 4), dtype=uint8,
            representing RGBA8888 pixel data.
    """
    qimg = qimg.convertToFormat(QImage.Format_RGBA8888)
    w, h = qimg.width(), qimg.height()
    bpl = qimg.bytesPerLine()
    # In PySide6, bits() returns a memoryview. Convert explicitly to bytes.
    buf = qimg.constBits().tobytes()
    arr = np.frombuffer(buf, dtype=np.uint8).reshape(h, bpl)
    arr = arr[:, : w * 4].reshape(h, w, 4).copy()  # Ensure independence from QImage
    return arr


def np_rgba_to_qimage(arr: np.ndarray) -> QImage:
    """
    Convert a NumPy RGBA array (H, W, 4) to a QImage in RGBA8888 format.

    Args:
        arr (np.ndarray): Input NumPy array of shape (height, width, 4) and dtype=uint8.

    Returns:
        QImage: Converted QImage in RGBA8888 format. The buffer is copied
            so that QImage owns the memory.
    """
    h, w, c = arr.shape
    assert c == 4 and arr.dtype == np.uint8, "Input must be (H,W,4) uint8 array"
    qimg = QImage(arr.data, w, h, w * 4, QImage.Format_RGBA8888)
    return qimg.copy()
