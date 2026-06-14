from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage, QPainter
from PySide6.QtWidgets import QApplication

import zstarview.gui.composite as render_composite
from zstarview import night_lights
from zstarview.types import ScreenGeometry

app = QApplication.instance() or QApplication([])


def test_blur_glow_mask_alpha_smooths_peak() -> None:
    alpha = np.zeros((5, 5), dtype=np.float32)
    alpha[2, 2] = 1.0

    blurred = render_composite._blur_glow_mask_alpha(alpha, passes=1)

    assert blurred.shape == alpha.shape
    assert blurred.dtype == np.float32
    assert 0.0 < float(blurred[2, 2]) < 1.0
    assert float(blurred[2, 2]) >= float(blurred[1, 2])
    assert np.all((blurred >= 0.0) & (blurred <= 1.0))


def test_glow_mask_to_qimage_uses_bright_base_color() -> None:
    mask = render_composite.GlowMask(
        alpha=np.asarray([[0.0, 0.5], [1.0, 0.25]], dtype=np.float32),
        scale=0.25,
    )

    image = render_composite._glow_mask_to_qimage(mask, (200, 100, 50))

    assert image.format() == QImage.Format.Format_ARGB32_Premultiplied
    center = image.pixelColor(0, 1)
    assert center.alpha() > 0
    assert center.red() >= center.green() >= center.blue()
    assert center.red() == 255
    assert center.green() == 128
    assert center.blue() == 64


def test_build_glow_mask_rasterizes_low_res_layers(monkeypatch) -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
        band_profiles=(
            night_lights.NightLightDistanceBandProfile(
                min_distance_km=0.5,
                max_distance_km=1.0,
                samples=(
                    night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
                ),
            ),
        ),
    )

    def _paint_dot(painter: QPainter, **_kwargs) -> None:
        painter.setPen(render_composite.Qt.PenStyle.NoPen)
        painter.setBrush(render_composite.Qt.GlobalColor.white)
        center = painter.viewport().center()
        painter.drawEllipse(center, 2, 2)

    monkeypatch.setattr(render_composite, "draw_night_light_glow_normal", _paint_dot)
    monkeypatch.setattr(render_composite, "draw_ridge_glow_normal", _paint_dot)

    mask = render_composite._build_glow_mask(
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        view_center=(0.0, 180.0),
        terrain_profile_altaz=[(0.0, 180.0)],
        terrain_secondary_ridges_altaz_layers=[[(0.0, 180.0)]],
        night_light_glow_profile=profile,
        night_light_opacity=0.5,
        ridge_glow_opacity=0.5,
        night_light_sun_alt_deg=-5.0,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
    )

    assert mask is not None
    assert mask.alpha.shape == (20, 20)
    assert mask.alpha.dtype == np.float32
    assert np.any(mask.alpha > 0.0)
    assert np.all((mask.alpha >= 0.0) & (mask.alpha <= 1.0))
