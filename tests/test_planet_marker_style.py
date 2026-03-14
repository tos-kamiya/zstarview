from __future__ import annotations

import astropy.time
import numpy as np
from PySide6.QtCore import QPoint, QPointF
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QApplication

from zstarview.render import draw as render_draw
from zstarview.types import CelestialData, PlanetBody, ScreenGeometry, ViewerData


def _empty_celestial_data(planets: list[PlanetBody]) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=planets,
        stars={
            "name": np.array([], dtype=object),
            "alt": np.array([], dtype=float),
            "az": np.array([], dtype=float),
            "vmag": np.array([], dtype=float),
            "bv": np.array([], dtype=float),
            "size_factor": np.array([], dtype=float),
            "color_factor_base": np.array([], dtype=float),
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


def _celestial_data_with_stars(stars: dict[str, np.ndarray]) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=[],
        stars=stars,
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


def _star_table(names: list[str], alt: float = 45.0, az: float = 180.0) -> dict[str, np.ndarray]:
    count = len(names)
    return {
        "name": np.array(names, dtype=object),
        "alt": np.full(count, alt, dtype=float),
        "az": np.full(count, az, dtype=float),
        "vmag": np.linspace(1.0, 5.0, count, dtype=float),
        "bv": np.zeros(count, dtype=float),
        "size_factor": np.ones(count, dtype=float),
        "color_factor_base": np.ones(count, dtype=float),
    }


def test_planets_are_drawn_with_disc_and_cross_markers(monkeypatch) -> None:
    disc_calls: list[tuple[float, int, tuple[int, int, int, int]]] = []
    bloom_calls: list[tuple[float, float | None]] = []
    cross_calls: list[tuple[float, float]] = []
    label_calls: list[str] = []

    def fake_draw_planet_disc(_painter, _pos, color, *, radius_px=1.0, alpha=255) -> None:
        disc_calls.append((radius_px, alpha, color.getRgb()))

    def fake_draw_planet_bloom(_painter, _pos, _color, *, core_radius_px=1.0, vmag=None) -> None:
        bloom_calls.append((core_radius_px, vmag))

    def fake_draw_gauge_cross(_painter, _color, _center, *, scale=1.0, pen_width=1.0) -> None:
        cross_calls.append((scale, pen_width))

    def fake_draw_outlined_text(_painter, text, _pos, _font, *_args, **_kwargs) -> None:
        label_calls.append(text)

    monkeypatch.setattr(render_draw, "draw_planet_disc", fake_draw_planet_disc)
    monkeypatch.setattr(render_draw, "draw_planet_bloom", fake_draw_planet_bloom)
    monkeypatch.setattr(render_draw, "draw_gauge_cross", fake_draw_gauge_cross)
    monkeypatch.setattr(render_draw, "draw_outlined_text", fake_draw_outlined_text)

    mars = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=True)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    geometry = ScreenGeometry(center=(100, 100), radius=80)

    render_draw.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=_empty_celestial_data([mars]),
        viewer_data=viewer,
        enlarge_moon=False,
    )

    assert len(disc_calls) == 1
    assert disc_calls[0][0] > 0.0
    assert disc_calls[0][1] > 0
    r, g, b, _ = disc_calls[0][2]
    assert r > g and r > b
    assert len(bloom_calls) == 1
    assert bloom_calls[0][0] > 0.0
    assert len(cross_calls) == 1
    assert cross_calls[0][0] < 1.0
    assert label_calls == ["Mars"]


def test_hover_can_identify_planet_name() -> None:
    mars = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=True)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    mouse_pos = QPoint(120, 90)

    highlighted = render_draw.find_highlighted_object(
        _empty_celestial_data([mars]),
        viewer,
        mouse_pos,
        geometry,
    )

    assert highlighted is not None
    obj, _ = highlighted
    assert getattr(obj, "name", "") == "mars"


def test_planet_draw_and_hover_ignore_horizon_visibility_flag(monkeypatch) -> None:
    disc_calls: list[tuple[float, int, tuple[int, int, int, int]]] = []

    def fake_draw_planet_disc(_painter, _pos, color, *, radius_px=1.0, alpha=255) -> None:
        disc_calls.append((radius_px, alpha, color.getRgb()))

    monkeypatch.setattr(render_draw, "draw_planet_disc", fake_draw_planet_disc)
    monkeypatch.setattr(render_draw, "draw_planet_bloom", lambda *_a, **_k: None)
    monkeypatch.setattr(render_draw, "draw_gauge_cross", lambda *_a, **_k: None)
    monkeypatch.setattr(render_draw, "draw_outlined_text", lambda *_a, **_k: None)

    mars = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=False)
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    geometry = ScreenGeometry(center=(120, 90), radius=70)

    render_draw.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=_empty_celestial_data([mars]),
        viewer_data=viewer,
        enlarge_moon=False,
    )
    assert len(disc_calls) == 1

    highlighted = render_draw.find_highlighted_object(
        _empty_celestial_data([mars]),
        viewer,
        QPoint(120, 90),
        geometry,
    )
    assert highlighted is not None
app = QApplication.instance() or QApplication([])


def test_hover_ignores_unnamed_stars() -> None:
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    mouse_pos = QPoint(120, 90)

    stars = _star_table(names=["", "Sirius"])
    highlighted = render_draw.find_highlighted_object(
        _celestial_data_with_stars(stars),
        viewer,
        mouse_pos,
        geometry,
    )

    assert highlighted is not None
    obj, _ = highlighted
    assert obj.get("name") == "Sirius"


def test_hover_returns_none_without_named_star() -> None:
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    mouse_pos = QPoint(120, 90)

    stars = _star_table(names=["", ""])
    highlighted = render_draw.find_highlighted_object(
        _celestial_data_with_stars(stars),
        viewer,
        mouse_pos,
        geometry,
    )

    assert highlighted is None


def test_overlay_skips_label_for_planet(monkeypatch) -> None:
    class DummyPainter:
        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawEllipse(self, *_args, **_kwargs) -> None:
            pass

    painter = DummyPainter()
    viewer = ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo", view_center=(45.0, 180.0))
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    planet = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=True)
    highlighted_object = (planet, QPointF(120.0, 90.0))

    label_calls: list[str] = []

    def fake_draw_outlined_text(_painter, text, *_args, **_kwargs) -> None:
        label_calls.append(text)

    monkeypatch.setattr(render_draw, "draw_outlined_text", fake_draw_outlined_text)

    render_draw.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=highlighted_object,
        text_font=QFont(),
    )

    assert "mars" not in label_calls


def test_overlay_info_includes_location_height_and_explicit_observer_height(monkeypatch) -> None:
    class DummyPainter:
        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawEllipse(self, *_args, **_kwargs) -> None:
            pass

    painter = DummyPainter()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="t/Tokyo Skytree",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
        location_height_label="Tower height",
        location_height_m=634.0,
        show_observer_height=True,
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)

    label_calls: list[str] = []

    def fake_draw_outlined_text(_painter, text, *_args, **_kwargs) -> None:
        label_calls.append(text)

    monkeypatch.setattr(render_draw, "draw_outlined_text", fake_draw_outlined_text)

    render_draw.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=None,
        text_font=QFont(),
    )

    assert "t/Tokyo Skytree" in label_calls
    assert "Tower height 634 m" in label_calls
    assert "Observer height 12 m" in label_calls
    assert label_calls.index("t/Tokyo Skytree") < label_calls.index("2026-02-27 00:00:00 UTC")
    assert label_calls.index("Tower height 634 m") < label_calls.index("2026-02-27 00:00:00 UTC")
    assert label_calls.index("Observer height 12 m") < label_calls.index("2026-02-27 00:00:00 UTC")
    assert label_calls.index("2026-02-27 00:00:00 UTC") < label_calls.index("Alt 45°  Az 180° (S)")
    assert label_calls.index("Alt 45°  Az 180° (S)") < label_calls.index("Vmag limit 6.0")
