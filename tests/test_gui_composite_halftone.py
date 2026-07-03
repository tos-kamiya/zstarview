# -*- coding: utf-8 -*-
"""Tests for halftone cloud sizing."""

from __future__ import annotations

import math
from types import SimpleNamespace

import numpy as np

from zstarview.gui import composite
from zstarview.gui.composite import _halftone_grid_delta, _halftone_level_diameters
from zstarview.types import ScreenGeometry, ViewProjection


def test_halftone_grid_delta_has_minimum_spacing() -> None:
    assert _halftone_grid_delta(600.0, 24) == 25.0
    assert _halftone_grid_delta(100.0, 24) == 22.0
    assert math.isclose(_halftone_grid_delta(1200.0, 24), 1200.0 / 24.0)


def test_halftone_level_diameters_scale_with_grid_spacing() -> None:
    diameters = _halftone_level_diameters(8.0, 1.0)
    assert diameters[0] == 0.0
    assert diameters[-1] > diameters[1]
    small = _halftone_level_diameters(8.0, 1.0)
    large = _halftone_level_diameters(120.0, 1.0)
    assert large[-1] > small[-1]


def test_render_halftone_cloud_uses_expanded_content_fov(monkeypatch) -> None:
    seen: dict[str, float] = {}

    class DummyImage:
        Format_ARGB32_Premultiplied = object()

        def __init__(self, width: int, height: int, _format: object) -> None:
            self.width = width
            self.height = height

        def fill(self, _value: object) -> None:
            return None

    class DummyPainter:
        RenderHint = SimpleNamespace(Antialiasing=object())

        def __init__(self, _image: DummyImage) -> None:
            self.calls: list[tuple[object, ...]] = []

        def setRenderHint(self, *args: object) -> None:
            self.calls.append(("setRenderHint", *args))

        def setClipPath(self, path: object) -> None:
            self.calls.append(("setClipPath", path))

        def setBrush(self, brush: object) -> None:
            self.calls.append(("setBrush", brush))

        def setPen(self, pen: object) -> None:
            self.calls.append(("setPen", pen))

        def drawEllipse(self, center: object, rx: float, ry: float) -> None:
            self.calls.append(("drawEllipse", center, rx, ry))

        def drawPoint(self, point: object) -> None:
            self.calls.append(("drawPoint", point))

        def end(self) -> None:
            self.calls.append(("end",))

    def fake_inverse_project_points(xs, ys, cx, cy, rr, view_center, edge_fov_deg, content_fov_deg):
        seen["content_fov_deg"] = float(content_fov_deg)
        count = len(xs)
        return (
            np.full(count, 45.0, dtype=np.float64),
            np.full(count, 180.0, dtype=np.float64),
            np.ones(count, dtype=bool),
        )

    def fake_sample_altaz_grid_amount(_grid, alts, _azs):
        return np.full(len(alts), 0.9, dtype=np.float32)

    monkeypatch.setattr(composite, "QImage", DummyImage)
    monkeypatch.setattr(composite, "QPainter", DummyPainter)
    monkeypatch.setattr(
        composite,
        "qimage_to_np_rgba",
        lambda image: np.zeros((image.height, image.width, 4), dtype=np.uint8),
    )
    monkeypatch.setattr(composite, "_halftone_grid_delta", lambda *_args, **_kwargs: 10.0)
    monkeypatch.setattr(composite, "_inverse_project_points", fake_inverse_project_points)
    monkeypatch.setattr(composite, "_sample_altaz_grid_amount", fake_sample_altaz_grid_amount)

    grid = SimpleNamespace(
        amount=np.ones((2, 2), dtype=np.float32),
        alt_min_deg=0.0,
        alt_max_deg=90.0,
        az_min_deg=0.0,
        az_max_deg=360.0,
    )
    geometry = ScreenGeometry(center=(16, 16), radius=10)
    projection = ViewProjection(view_center=(45.0, 180.0), edge_fov_deg=90.0, content_fov_deg=90.0)

    composite._render_halftone_cloud_rgba_from_altaz_grid(
        grid,
        32,
        32,
        composite.CLOUD_HATCH_DEFAULT,
        geometry=geometry,
        projection=projection,
    )

    assert seen["content_fov_deg"] == 102.0
