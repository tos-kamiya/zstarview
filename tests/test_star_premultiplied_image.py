import numpy as np
from PySide6.QtGui import QImage

from zstarview.render.stars import (
    _canvas_to_argb32_premultiplied_image,
    _rgba_to_argb32_premultiplied_image,
)


def test_rgba_image_uses_premultiplied_argb_and_preserves_thin_color() -> None:
    rgba = np.array([[[255, 200, 100, 10], [30, 20, 10, 0]]], dtype=np.uint8)

    image = _rgba_to_argb32_premultiplied_image(rgba)

    assert image.format() == QImage.Format.Format_ARGB32_Premultiplied
    assert image.pixelColor(0, 0).getRgb()[:4] == (255, 204, 102, 10)
    assert image.pixelColor(1, 0).getRgb()[:4] == (0, 0, 0, 0)


def test_canvas_values_are_stored_as_premultiplied_rgb() -> None:
    canvas = np.array([[[192, 128, 64], [0, 0, 0]]], dtype=np.uint8)

    image = _canvas_to_argb32_premultiplied_image(canvas)

    assert image.format() == QImage.Format.Format_ARGB32_Premultiplied
    assert image.pixelColor(0, 0).getRgb()[:4] == (255, 170, 85, 192)
    assert image.pixelColor(1, 0).getRgb()[:4] == (0, 0, 0, 0)
