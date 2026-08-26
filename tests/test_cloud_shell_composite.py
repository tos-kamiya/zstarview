from __future__ import annotations

import datetime as dt

import numpy as np
from PySide6.QtGui import QImage, QPainter
from PySide6.QtWidgets import QApplication

from zstarview.clouddisc.altaz_grid import CloudAltAzGrid
from zstarview.gui import composite as composite_module
from zstarview.gui.composite import _combine_cloud_shell_rgba
from zstarview.types import ScreenGeometry


def _float_reference(
    layers: list[tuple[np.ndarray, tuple[int, int, int]]],
) -> np.ndarray:
    shape = layers[0][0].shape
    premultiplied = np.zeros(shape[:2] + (3,), dtype=np.float32)
    alpha = np.zeros(shape[:2], dtype=np.float32)
    for image, tint_rgb in reversed(layers):
        source_alpha = image[..., 3].astype(np.float32) / 255.0
        source_rgb = (
            image[..., :3].astype(np.float32)
            * np.asarray(tint_rgb, dtype=np.float32)[None, None, :]
            / 255.0
        )
        premultiplied = (
            source_rgb * source_alpha[..., None]
            + premultiplied * (1.0 - source_alpha[..., None])
        )
        alpha = source_alpha + alpha * (1.0 - source_alpha)
    rgb = np.zeros_like(premultiplied)
    np.divide(
        premultiplied,
        np.maximum(alpha[..., None], 1.0e-6),
        out=rgb,
        where=alpha[..., None] > 0.0,
    )
    result = np.zeros(shape, dtype=np.uint8)
    result[..., :3] = np.clip(np.round(rgb), 0, 255).astype(np.uint8)
    result[..., 3] = np.clip(np.round(alpha * 255.0), 0, 255).astype(np.uint8)
    return result


def test_combine_cloud_shell_rgba_matches_float_source_over() -> None:
    rng = np.random.default_rng(20260826)
    layers = [
        (
            rng.integers(0, 256, size=(13, 17, 4), dtype=np.uint8),
            tint,
        )
        for tint in ((148, 162, 171), (204, 214, 220), (255, 255, 255))
    ]

    actual = _combine_cloud_shell_rgba(layers)
    expected = _float_reference(layers)

    assert actual is not None
    np.testing.assert_allclose(actual, expected, atol=2, rtol=0)


def test_combine_cloud_shell_rgba_returns_none_without_layers() -> None:
    assert _combine_cloud_shell_rgba([]) is None


def test_draw_cloud_overlay_reuses_matching_raster(monkeypatch) -> None:
    _app = QApplication.instance() or QApplication([])
    amount = np.ones((4, 8), dtype=np.float32)
    grid = CloudAltAzGrid(
        amount=amount,
        missing_mask=np.zeros_like(amount, dtype=np.uint8),
        alt_min_deg=-5.0,
        alt_max_deg=90.0,
        az_min_deg=0.0,
        az_max_deg=360.0,
        observer_lat=35.0,
        observer_lon=139.0,
        satellite="test",
        product="test",
        time_utc=dt.datetime(2026, 8, 26, tzinfo=dt.UTC),
        shells_km=(6374.0,),
        source_key="test",
        coverage_ratio=1.0,
    )
    calls = 0

    def fake_render(*_args, **_kwargs) -> np.ndarray:
        nonlocal calls
        calls += 1
        result = np.zeros((24, 32, 4), dtype=np.uint8)
        result[8:16, 10:22] = (255, 255, 255, 255)
        return result

    monkeypatch.setattr(composite_module, "_render_cloud_grid_rgba", fake_render)
    compositor = composite_module.SkyCompositorCache()
    target = QImage(64, 48, QImage.Format.Format_ARGB32_Premultiplied)
    target.fill(0)
    painter = QPainter(target)
    kwargs = dict(
        geometry=ScreenGeometry((32, 24), 24),
        cloud_alpha=0.09,
        render_size=(32, 24),
        view_center=(30.0, 180.0),
        cloud_altaz_grid=grid,
        missing_mask=None,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        sun_alt_deg=20.0,
        theme=None,
    )
    compositor.draw_cloud_overlay(painter, **kwargs)
    compositor.draw_cloud_overlay(painter, **kwargs)
    painter.end()

    assert calls == 1
