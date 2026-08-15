from __future__ import annotations

import numpy as np
from PySide6.QtGui import QColor, QImage

from zstarview.render.qt_image import qimage_to_np_rgba
from zstarview.render.solar_tint import (
    SOLAR_HOVER_MARKER_FALLBACK_STRENGTH,
    clear_solar_hover_color_cache,
    colorize_solar_hover_image,
    solar_hover_color_multiplier,
)


def _source_image() -> QImage:
    image = QImage(3, 1, QImage.Format.Format_RGBA8888)
    image.setPixelColor(0, 0, QColor(0, 0, 0, 255))
    image.setPixelColor(1, 0, QColor(128, 128, 128, 190))
    image.setPixelColor(2, 0, QColor(255, 255, 255, 0))
    return image


def test_colorize_solar_hover_image_preserves_alpha_and_adds_warm_tint() -> None:
    clear_solar_hover_color_cache()
    tinted = colorize_solar_hover_image(
        _source_image(),
        image_id=1,
        sun_alt_deg=45.0,
        observer_height_m=0.0,
        aerosol_optical_depth=0.15,
    )
    rgba = qimage_to_np_rgba(tinted)

    assert np.array_equal(rgba[:, :, 3], np.asarray([[255, 190, 0]], dtype=np.uint8))
    assert tuple(rgba[0, 0, :3]) == (0, 0, 0)
    assert int(rgba[0, 1, 0]) > int(rgba[0, 1, 2])
    assert int(rgba[0, 1, 1]) > int(rgba[0, 1, 2])


def test_colorize_solar_hover_image_reuses_cached_qimage() -> None:
    clear_solar_hover_color_cache()
    kwargs = {
        "image_id": 2,
        "sun_alt_deg": 30.0,
        "observer_height_m": 100.0,
        "aerosol_optical_depth": 0.2,
    }

    first = colorize_solar_hover_image(_source_image(), **kwargs)
    second = colorize_solar_hover_image(_source_image(), **kwargs)

    assert first.cacheKey() == second.cacheKey()


def test_below_horizon_tint_falls_back_to_subdued_solar_marker_color() -> None:
    clear_solar_hover_color_cache()
    multiplier = np.asarray(
        solar_hover_color_multiplier(-4, 0, 15),
        dtype=np.float32,
    )
    marker = np.asarray((248, 196, 64), dtype=np.float32) / 255.0
    marker_linear = np.where(
        marker <= 0.04045,
        marker / 12.92,
        np.power((marker + 0.055) / 1.055, 2.4),
    )

    np.testing.assert_allclose(
        multiplier,
        marker_linear * SOLAR_HOVER_MARKER_FALLBACK_STRENGTH,
        rtol=1e-6,
        atol=1e-6,
    )
    assert float(multiplier.max()) > 0.0
