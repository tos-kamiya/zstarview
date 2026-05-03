from __future__ import annotations

import astropy.time
import numpy as np
import pytest
from PySide6.QtCore import QPoint, QPointF
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QApplication

from zstarview.asterisms import Asterism
from zstarview.paths import ASTERISM_CLIP_FIELD_OF_VIEW_DEG, PALETTE_ASTERISM_RGB, THEME_STYLES_BY_PRESET
from zstarview.render import asterisms as render_asterisms
from zstarview.render import stars as render_stars
from zstarview.types import CelestialData, ScreenGeometry, StarCatalogMeta, ViewerData

app = QApplication.instance() or QApplication([])


class DummyPainter:
    def __init__(self) -> None:
        self.polyline_count = 0
        self.polylines = []
        self.pen_widths: list[float] = []
        self.pen_alphas: list[float] = []

    def save(self) -> None:
        pass

    def restore(self) -> None:
        pass

    def setPen(self, pen, *_args, **_kwargs) -> None:
        self.pen_widths.append(float(pen.widthF()))
        self.pen_alphas.append(float(pen.color().alphaF()))

    def drawPolyline(self, _poly) -> None:
        self.polyline_count += 1
        self.polylines.append([(point.x(), point.y()) for point in _poly])


def _celestial_data_with_asterism_star_positions() -> CelestialData:
    star_catalog_meta = StarCatalogMeta(
        name_indices=np.array([0, 1], dtype=np.int32),
        names=np.array(["Star A", "Star B"], dtype=object),
        source_id_indices=np.array([0, 1], dtype=np.int32),
        source_ids=np.array(["HIP1", "HIP2"], dtype=object),
    )
    return CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=[],
        stars={
            "star_index": np.array([0, 1], dtype=np.int32),
            "alt": np.array([45.0, 45.0], dtype=float),
            "az": np.array([170.0, 190.0], dtype=float),
            "vmag": np.array([1.0, 2.0], dtype=float),
            "bv": np.zeros(2, dtype=float),
            "size_factor": np.ones(2, dtype=float),
            "color_factor_base": np.ones(2, dtype=float),
        },
        deep_sky_objects={
            "id": np.array([], dtype=object),
            "name": np.array([], dtype=object),
            "type": np.array([], dtype=object),
            "alt": np.array([], dtype=float),
            "az": np.array([], dtype=float),
            "vmag": np.array([], dtype=float),
            "major_arcmin": np.array([], dtype=float),
            "minor_arcmin": np.array([], dtype=float),
            "pa_deg": np.array([], dtype=float),
        },
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
        star_catalog_meta=star_catalog_meta,
    )


def test_draw_asterisms_draws_dim_overlay_without_hover(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))

    monkeypatch.setattr(render_asterisms, "ASTERISMS", (asterism,))

    render_asterisms.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=None,
        text_font=QFont(),
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert painter.polyline_count == 3


def test_draw_asterisms_keeps_dim_overlay_base_widths_fixed(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))

    monkeypatch.setattr(render_asterisms, "ASTERISMS", (asterism,))

    render_asterisms.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=None,
        text_font=QFont(),
        theme=THEME_STYLES_BY_PRESET["night"],
        line_width_scale=2.0,
    )

    assert painter.polyline_count == 3
    assert painter.pen_widths == [8.0, 5.2, 2.8]


def test_draw_asterisms_scales_dim_overlay_alpha_with_visibility_boost(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))

    monkeypatch.setattr(render_asterisms, "ASTERISMS", (asterism,))

    render_asterisms.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=None,
        text_font=QFont(),
        theme=THEME_STYLES_BY_PRESET["night"],
        base_line_width_scale=2.0,
        base_line_alpha_scale=2.0,
    )

    assert painter.pen_alphas[:3] == [
        pytest.approx(5 / 255.0),
        pytest.approx(10 / 255.0),
        pytest.approx(32 / 255.0),
    ]
    assert painter.pen_widths[:3] == [4.0, 2.6, 2.8]


def test_draw_asterisms_dim_overlay_uses_softer_alpha(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))

    monkeypatch.setattr(render_asterisms, "ASTERISMS", (asterism,))

    render_asterisms.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=None,
        text_font=QFont(),
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert painter.pen_alphas == sorted(painter.pen_alphas)
    assert painter.pen_alphas[-1] < 0.1


def test_draw_asterisms_hover_adds_bright_overlay_and_label(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))
    label_candidates: list[dict[str, object]] = []

    monkeypatch.setattr(render_asterisms, "ASTERISMS", (asterism,))
    monkeypatch.setattr(render_asterisms, "pick_rotating_asterism", lambda *_args, **_kwargs: asterism)

    render_asterisms.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=({"source_id": "HIP1", "name": "Star A"}, QPointF(120.0, 90.0)),
        text_font=QFont(),
        label_candidates=label_candidates,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert painter.polyline_count == 6
    assert painter.pen_widths[-1] == 1.0
    assert [c["text"] for c in label_candidates] == ["Test Asterism"]
    label_style = label_candidates[0]["style"]
    assert (label_style.text_color.red(), label_style.text_color.green(), label_style.text_color.blue()) == PALETTE_ASTERISM_RGB


def test_draw_asterisms_deduplicates_shared_dim_segments(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    first = Asterism("first", "First", (("HIP1", "HIP2"),))
    second = Asterism("second", "Second", (("HIP2", "HIP1"),))

    monkeypatch.setattr(render_asterisms, "ASTERISMS", (first, second))

    render_asterisms.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=None,
        text_font=QFont(),
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert painter.polyline_count == 3


def test_draw_asterisms_clips_with_asterism_specific_wide_fov(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=[],
        stars={
            "star_index": np.array([0, 1], dtype=np.int32),
            "alt": np.array([45.0, -58.0], dtype=float),
            "az": np.array([180.0, 180.0], dtype=float),
            "vmag": np.array([1.0, 2.0], dtype=float),
            "bv": np.zeros(2, dtype=float),
            "size_factor": np.ones(2, dtype=float),
            "color_factor_base": np.ones(2, dtype=float),
        },
        deep_sky_objects={
            "id": np.array([], dtype=object),
            "name": np.array([], dtype=object),
            "type": np.array([], dtype=object),
            "alt": np.array([], dtype=float),
            "az": np.array([], dtype=float),
            "vmag": np.array([], dtype=float),
            "major_arcmin": np.array([], dtype=float),
            "minor_arcmin": np.array([], dtype=float),
            "pa_deg": np.array([], dtype=float),
        },
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
        star_catalog_meta=StarCatalogMeta(
            name_indices=np.array([0, 1], dtype=np.int32),
            names=np.array(["Star A", "Star B"], dtype=object),
            source_id_indices=np.array([0, 1], dtype=np.int32),
            source_ids=np.array(["HIP1", "HIP2"], dtype=object),
        ),
    )
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))

    monkeypatch.setattr(render_asterisms, "ASTERISMS", (asterism,))

    render_asterisms.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=None,
        text_font=QFont(),
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert painter.polyline_count == 3
    ys = [point[1] for polyline in painter.polylines for point in polyline]
    assert ys
    assert max(ys) <= geometry.center[1] + (geometry.radius * (ASTERISM_CLIP_FIELD_OF_VIEW_DEG / 90.0)) + 1.0e-6


def test_draw_asterisms_scales_line_widths_with_star_upscale(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))
    label_candidates: list[dict[str, object]] = []

    monkeypatch.setattr(render_asterisms, "ASTERISMS", (asterism,))
    monkeypatch.setattr(render_asterisms, "pick_rotating_asterism", lambda *_args, **_kwargs: asterism)

    render_asterisms.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=({"source_id": "HIP1", "name": "Star A"}, QPointF(120.0, 90.0)),
        text_font=QFont(),
        label_candidates=label_candidates,
        line_width_scale=2.0,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert painter.pen_widths[:6] == [8.0, 5.2, 2.8, 10.0, 6.4, 2.0]
    assert [c["text"] for c in label_candidates] == ["Test Asterism"]


def test_find_highlighted_object_accepts_unnamed_asterism_member(monkeypatch) -> None:
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=[],
        stars={
            "star_index": np.array([0], dtype=np.int32),
            "alt": np.array([45.0], dtype=float),
            "az": np.array([180.0], dtype=float),
            "vmag": np.array([3.0], dtype=float),
            "bv": np.zeros(1, dtype=float),
            "size_factor": np.ones(1, dtype=float),
            "color_factor_base": np.ones(1, dtype=float),
        },
        deep_sky_objects={
            "id": np.array([], dtype=object),
            "name": np.array([], dtype=object),
            "type": np.array([], dtype=object),
            "alt": np.array([], dtype=float),
            "az": np.array([], dtype=float),
            "vmag": np.array([], dtype=float),
            "major_arcmin": np.array([], dtype=float),
            "minor_arcmin": np.array([], dtype=float),
            "pa_deg": np.array([], dtype=float),
        },
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
        star_catalog_meta=StarCatalogMeta(
            name_indices=np.array([], dtype=np.int32),
            names=np.array([], dtype=object),
            source_id_indices=np.array([0], dtype=np.int32),
            source_ids=np.array(["HIP1"], dtype=object),
        ),
    )

    monkeypatch.setattr(render_stars, "ASTERISM_REQUIRED_SOURCE_IDS", frozenset({"HIP1"}))

    highlighted = render_stars.find_highlighted_object(
        celestial_data=celestial_data,
        viewer_data=viewer,
        mouse_pos=QPoint(120, 90),
        geometry=geometry,
    )

    assert highlighted is not None
    obj, _ = highlighted
    assert isinstance(obj, dict)
    assert obj["source_id"] == "HIP1"
