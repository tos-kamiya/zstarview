from __future__ import annotations

import numpy as np
from PySide6.QtGui import QColor, QImage

from zstarview.render.moon_tint import (
    MOON_HOVER_MARKER_FALLBACK_STRENGTH,
    clear_moon_hover_color_cache,
    colorize_moon_hover_image,
    moon_hover_color_multiplier,
)
from zstarview.render.qt_image import qimage_to_np_rgba


def _source_image() -> QImage:
    image = QImage(3, 1, QImage.Format.Format_RGBA8888)
    image.setPixelColor(0, 0, QColor(30, 20, 10, 255))
    image.setPixelColor(1, 0, QColor(128, 128, 128, 190))
    image.setPixelColor(2, 0, QColor(255, 255, 255, 0))
    return image


def _colorize(*, altitude: float, aod: float = 0.15) -> QImage:
    return colorize_moon_hover_image(
        _source_image(),
        moon_alt_deg=altitude,
        observer_height_m=0.0,
        aerosol_optical_depth=aod,
    )


def test_colorize_moon_hover_image_preserves_alpha_and_source_color() -> None:
    clear_moon_hover_color_cache()
    rgba = qimage_to_np_rgba(_colorize(altitude=60.0))

    assert np.array_equal(rgba[:, :, 3], np.asarray([[255, 190, 0]], dtype=np.uint8))
    assert int(rgba[0, 0, 0]) > int(rgba[0, 0, 1]) > int(rgba[0, 0, 2])


def test_low_moon_hover_image_is_redder_than_high_moon_image() -> None:
    clear_moon_hover_color_cache()
    high = qimage_to_np_rgba(_colorize(altitude=60.0, aod=0.3))[0, 1, :3]
    low = qimage_to_np_rgba(_colorize(altitude=2.0, aod=0.3))[0, 1, :3]

    assert float(low[0]) / max(1.0, float(low[2])) > float(high[0]) / float(high[2])


def test_colorize_moon_hover_image_reuses_cached_qimage() -> None:
    clear_moon_hover_color_cache()
    source = _source_image()
    kwargs = {
        "moon_alt_deg": 30.0,
        "observer_height_m": 100.0,
        "aerosol_optical_depth": 0.2,
    }

    first = colorize_moon_hover_image(source, **kwargs)
    second = colorize_moon_hover_image(source, **kwargs)

    assert first.cacheKey() == second.cacheKey()


def test_below_horizon_moon_uses_neutral_marker_fallback() -> None:
    clear_moon_hover_color_cache()
    multiplier = np.asarray(moon_hover_color_multiplier(-4, 0, 15))

    assert np.all(multiplier > 0.0)
    assert np.allclose(multiplier, multiplier[0])
    assert float(multiplier[0]) < MOON_HOVER_MARKER_FALLBACK_STRENGTH
