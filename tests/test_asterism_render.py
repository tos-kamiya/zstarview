from __future__ import annotations

import astropy.time
import numpy as np
from PySide6.QtCore import QPoint, QPointF
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QApplication

from zstarview.asterisms import Asterism
from zstarview.render import draw as render_draw
from zstarview.types import CelestialData, ScreenGeometry, ViewerData

app = QApplication.instance() or QApplication([])


class DummyPainter:
    def __init__(self) -> None:
        self.polyline_count = 0

    def save(self) -> None:
        pass

    def restore(self) -> None:
        pass

    def setPen(self, *_args, **_kwargs) -> None:
        pass

    def drawPolyline(self, _poly) -> None:
        self.polyline_count += 1


def _celestial_data_with_asterism_star_positions() -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=[],
        stars={
            "name": np.array(["Star A", "Star B"], dtype=object),
            "source_id": np.array(["HIP1", "HIP2"], dtype=object),
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
    )


def test_draw_asterisms_draws_dim_overlay_without_hover(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))

    monkeypatch.setattr(render_draw, "ASTERISMS", (asterism,))

    render_draw.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=None,
        text_font=QFont(),
    )

    assert painter.polyline_count == 2


def test_draw_asterisms_hover_adds_bright_overlay_and_label(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    asterism = Asterism("test", "Test Asterism", (("HIP1", "HIP2"),))
    label_candidates: list[dict[str, object]] = []

    monkeypatch.setattr(render_draw, "ASTERISMS", (asterism,))
    monkeypatch.setattr(render_draw, "pick_rotating_asterism", lambda *_args, **_kwargs: asterism)

    render_draw.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=({"source_id": "HIP1", "name": "Star A"}, QPointF(120.0, 90.0)),
        text_font=QFont(),
        label_candidates=label_candidates,
    )

    assert painter.polyline_count == 4
    assert [c["text"] for c in label_candidates] == ["Test Asterism"]


def test_draw_asterisms_deduplicates_shared_dim_segments(monkeypatch) -> None:
    painter = DummyPainter()
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = _celestial_data_with_asterism_star_positions()
    first = Asterism("first", "First", (("HIP1", "HIP2"),))
    second = Asterism("second", "Second", (("HIP2", "HIP1"),))

    monkeypatch.setattr(render_draw, "ASTERISMS", (first, second))

    render_draw.draw_asterisms(
        painter=painter,
        geometry=geometry,
        celestial_data=celestial_data,
        viewer_data=viewer,
        highlighted_object=None,
        text_font=QFont(),
    )

    assert painter.polyline_count == 2


def test_find_highlighted_object_accepts_unnamed_asterism_member(monkeypatch) -> None:
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    celestial_data = CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=[],
        stars={
            "name": np.array([""], dtype=object),
            "source_id": np.array(["HIP1"], dtype=object),
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
    )

    monkeypatch.setattr(render_draw, "ASTERISM_REQUIRED_SOURCE_IDS", frozenset({"HIP1"}))

    highlighted = render_draw.find_highlighted_object(
        celestial_data=celestial_data,
        viewer_data=viewer,
        mouse_pos=QPoint(120, 90),
        geometry=geometry,
    )

    assert highlighted is not None
    obj, _ = highlighted
    assert isinstance(obj, dict)
    assert obj["source_id"] == "HIP1"
