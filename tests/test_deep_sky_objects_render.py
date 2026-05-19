from __future__ import annotations

from types import SimpleNamespace

import astropy.time
import numpy as np
from PySide6.QtGui import QColor

from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.render import deep_sky_objects as render_deep_sky_objects
from zstarview.render.deep_sky_objects import DSO_LABEL_RGB
from zstarview.types import CelestialData, ScreenGeometry, ViewerData


class _DummyPainter:
    def __init__(self) -> None:
        self.brush_rgbs: list[tuple[int, int, int, int]] = []
        self.pen_rgbs: list[tuple[int, int, int, int]] = []

    def save(self) -> None:
        pass

    def restore(self) -> None:
        pass

    def setPen(self, *_args, **_kwargs) -> None:
        if _args and hasattr(_args[0], "color"):
            pen = _args[0]
            color = QColor(pen.color())
            self.pen_rgbs.append((color.red(), color.green(), color.blue(), color.alpha()))

    def setBrush(self, brush) -> None:
        if not hasattr(brush, "red"):
            return
        color = QColor(brush)
        self.brush_rgbs.append((color.red(), color.green(), color.blue(), color.alpha()))

    def drawPolygon(self, *_args, **_kwargs) -> None:
        pass


def _make_celestial_data() -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
        planets=[],
        stars={
            "name": [],
            "source_id": [],
            "alt": [],
            "az": [],
            "vmag": [],
            "bv": [],
            "size_factor": [],
            "color_factor_base": [],
        },
        deep_sky_objects={
            "id": np.array(["M31"], dtype=object),
            "name": np.array(["Andromeda Galaxy"], dtype=object),
            "type": np.array(["G"], dtype=object),
            "alt": np.array([45.0], dtype=float),
            "az": np.array([180.0], dtype=float),
            "vmag": np.array([3.4], dtype=float),
            "major_arcmin": np.array([30.0], dtype=float),
            "minor_arcmin": np.array([20.0], dtype=float),
            "pa_deg": np.array([0.0], dtype=float),
        },
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
        star_catalog_meta=None,
    )


def test_draw_deep_sky_shapes_uses_the_same_fill_color_in_all_themes(monkeypatch) -> None:
    painter = _DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    celestial_data = _make_celestial_data()

    monkeypatch.setattr(
        render_deep_sky_objects,
        "_dso_ellipse_polygon",
        lambda **_kwargs: SimpleNamespace(),
    )

    for preset in ("night", "black", "transparent", "white", "day"):
        theme = THEME_STYLES_BY_PRESET[preset]
        painter.brush_rgbs.clear()
        render_deep_sky_objects.draw_deep_sky_shapes(
            painter,
            geometry,
            celestial_data,
            viewer,
            theme=theme,
        )
        assert painter.brush_rgbs
        assert painter.brush_rgbs[0][:3] == (110, 185, 255)


def test_draw_dso_hover_info_uses_dso_label_color_in_all_themes() -> None:
    painter = _DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    highlighted_dso = (
        {
            "name": "Andromeda Galaxy",
            "alt": 45.0,
            "az": 180.0,
            "major_arcmin": 30.0,
            "minor_arcmin": 20.0,
            "pa_deg": 0.0,
        },
        SimpleNamespace(),
    )

    for preset in ("night", "black", "transparent", "white", "day"):
        theme = THEME_STYLES_BY_PRESET[preset]
        painter.pen_rgbs.clear()
        render_deep_sky_objects.draw_dso_hover_info(
            painter,
            geometry,
            viewer,
            highlighted_dso,
            text_font=object(),
            theme=theme,
        )
        assert painter.pen_rgbs
        assert painter.pen_rgbs[-1][:3] == DSO_LABEL_RGB
