from __future__ import annotations

import astropy.time
import numpy as np
from PySide6.QtCore import QPoint

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
        },
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
    )


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
    assert label_calls == ["mars"]


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
