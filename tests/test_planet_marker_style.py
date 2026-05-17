from __future__ import annotations

import astropy.time
import math
import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRectF
from PySide6.QtGui import QColor, QFont, QFontMetrics, QImage, QPainter
from PySide6.QtWidgets import QApplication

from zstarview.render import background as render_background
from zstarview.render import aircraft as render_aircraft
from zstarview.render import guides as render_guides
from zstarview.render import overlay_info as render_overlay_info
from zstarview.render import satellites as render_satellites
from zstarview.render import solar_system as render_solar_system
from zstarview.render import stars as render_stars
from zstarview.render import text as render_text
from zstarview.paths import PALETTE_AIRCRAFT_AND_SATELLITE_RGB
from zstarview.render.deep_sky_objects import DSO_LABEL_RGB
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.aircraft.types import AircraftOverlayPoint
from zstarview.satellites.types import SatelliteOverlayPoint
from zstarview.types import (
    CelestialData,
    PlanetBody,
    ScreenGeometry,
    StarCatalogMeta,
    ViewerData,
)


def _empty_star_catalog_meta() -> StarCatalogMeta:
    return StarCatalogMeta(
        name_indices=np.array([], dtype=np.int32),
        names=np.array([], dtype=object),
        source_id_indices=np.array([], dtype=np.int32),
        source_ids=np.array([], dtype=object),
    )


def _empty_celestial_data(planets: list[PlanetBody]) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        planets=planets,
        stars={
            "star_index": np.array([], dtype=np.int32),
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
        star_catalog_meta=_empty_star_catalog_meta(),
    )


def _celestial_data_with_stars(
    stars: dict[str, np.ndarray], star_catalog_meta: StarCatalogMeta | None = None
) -> CelestialData:
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
        star_catalog_meta=star_catalog_meta,
    )


def _star_table(
    names: list[str],
    *,
    source_ids: list[str] | None = None,
    alt: float = 45.0,
    az: float = 180.0,
) -> tuple[dict[str, np.ndarray], StarCatalogMeta]:
    count = len(names)
    stars = {
        "star_index": np.arange(count, dtype=np.int32),
        "alt": np.full(count, alt, dtype=float),
        "az": np.full(count, az, dtype=float),
        "vmag": np.linspace(1.0, 5.0, count, dtype=float),
        "bv": np.zeros(count, dtype=float),
        "size_factor": np.ones(count, dtype=float),
        "color_factor_base": np.ones(count, dtype=float),
    }
    source_values = [""] * count if source_ids is None else source_ids
    name_indices = np.array(
        [idx for idx, name in enumerate(names) if str(name).strip()], dtype=np.int32
    )
    source_id_indices = np.array(
        [idx for idx, value in enumerate(source_values) if str(value).strip()],
        dtype=np.int32,
    )
    meta = StarCatalogMeta(
        name_indices=name_indices,
        names=np.array([names[idx] for idx in name_indices], dtype=object),
        source_id_indices=source_id_indices,
        source_ids=np.array(
            [source_values[idx] for idx in source_id_indices], dtype=object
        ),
    )
    return stars, meta


def test_planets_are_drawn_with_disc_and_cross_markers(monkeypatch) -> None:
    disc_calls: list[tuple[float, int, tuple[int, int, int, int]]] = []
    bloom_calls: list[tuple[float, float | None]] = []
    cross_calls: list[tuple[float, float]] = []
    label_calls: list[str] = []

    def fake_draw_planet_disc(
        _painter, _pos, color, *, radius_px=1.0, alpha=255
    ) -> None:
        disc_calls.append((radius_px, alpha, color.getRgb()))

    def fake_draw_planet_bloom(
        _painter, _pos, _color, *, core_radius_px=1.0, vmag=None
    ) -> None:
        bloom_calls.append((core_radius_px, vmag))

    def fake_draw_gauge_cross(
        _painter, _color, _center, *, scale=1.0, pen_width=1.0
    ) -> None:
        cross_calls.append((scale, pen_width))

    def fake_draw_outlined_text(_painter, text, _pos, _font, *_args, **_kwargs) -> None:
        label_calls.append(text)

    monkeypatch.setattr(render_solar_system, "draw_planet_disc", fake_draw_planet_disc)
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", fake_draw_planet_bloom
    )
    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", fake_draw_gauge_cross)
    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", fake_draw_outlined_text
    )

    mars = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=True)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(100, 100), radius=80)

    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=_empty_celestial_data([mars]),
        viewer_data=viewer,
        enlarge_moon=False,
        theme=THEME_STYLES_BY_PRESET["night"],
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


def test_planet_label_is_skipped_when_body_marker_is_outside_viewport(
    monkeypatch,
) -> None:
    monkeypatch.setattr(render_solar_system, "draw_planet_disc", lambda *_a, **_k: None)
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", lambda *_a, **_k: None
    )
    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", lambda *_a, **_k: None)

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        mars = PlanetBody(name="mars", alt=0.0, az=180.0, symbol="♂", is_visible=True)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
        )
        geometry = ScreenGeometry(center=(20, 20), radius=80)
        label_candidates: list[dict[str, object]] = []

        render_solar_system.draw_solar_system_bodies(
            painter=painter,
            geometry=geometry,
            celestial_data=_empty_celestial_data([mars]),
            viewer_data=viewer,
            enlarge_moon=False,
            label_candidates=label_candidates,
            theme=THEME_STYLES_BY_PRESET["night"],
        )
    finally:
        painter.end()

    assert label_candidates == []


def test_hover_can_identify_planet_name() -> None:
    mars = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=True)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    mouse_pos = QPoint(120, 90)

    highlighted = render_stars.find_highlighted_object(
        _empty_celestial_data([mars]),
        viewer,
        mouse_pos,
        geometry,
    )

    assert highlighted is not None
    obj, _ = highlighted
    assert obj.name == "mars"


def test_enlarge_moon_scales_display_radius_by_five(monkeypatch) -> None:
    moon_draw_radii: list[float] = []

    def fake_draw_moon(_painter, _center, radius_px, *_args, **_kwargs) -> None:
        moon_draw_radii.append(float(radius_px))

    monkeypatch.setattr(render_solar_system, "draw_moon", fake_draw_moon)
    monkeypatch.setattr(
        render_solar_system, "draw_gauge_cross", lambda *_args, **_kwargs: None
    )

    sun = PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True)
    moon = PlanetBody(name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(100, 100), radius=80)
    celestial = _empty_celestial_data([sun, moon])

    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=celestial,
        viewer_data=viewer,
        enlarge_moon=False,
        theme=THEME_STYLES_BY_PRESET["night"],
        label_candidates=[],
    )
    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=celestial,
        viewer_data=viewer,
        enlarge_moon=True,
        label_candidates=[],
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert len(moon_draw_radii) == 2
    assert moon_draw_radii[0] == 2.5
    assert moon_draw_radii[1] == 12.5


def test_outline_bright_bodies_keeps_enlarged_moon_filled(monkeypatch) -> None:
    moon_outline_radii: list[float] = []
    moon_draw_radii: list[float] = []
    planet_outline_radii: list[float] = []
    cross_calls: list[tuple[float, float]] = []

    def fake_draw_moon_outline(_painter, _center, radius_px, _color) -> None:
        moon_outline_radii.append(float(radius_px))

    def fake_draw_moon(_painter, _center, radius_px, *_args, **_kwargs) -> None:
        moon_draw_radii.append(float(radius_px))

    def fake_draw_planet_outline(_painter, _pos, _color, *, radius_px=1.0) -> None:
        planet_outline_radii.append(float(radius_px))

    def fake_draw_gauge_cross(_painter, _color, _center, *, scale=1.0, pen_width=1.0) -> None:
        cross_calls.append((float(scale), float(pen_width)))

    monkeypatch.setattr(render_solar_system, "draw_moon_outline", fake_draw_moon_outline)
    monkeypatch.setattr(render_solar_system, "draw_planet_outline", fake_draw_planet_outline)
    monkeypatch.setattr(render_solar_system, "draw_planet_disc", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(render_solar_system, "draw_planet_bloom", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(render_solar_system, "draw_moon", fake_draw_moon)
    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", fake_draw_gauge_cross)
    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", lambda *_args, **_kwargs: None
    )

    sun = PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True)
    moon = PlanetBody(name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True)
    mars = PlanetBody(
        name="mars",
        alt=45.0,
        az=180.0,
        symbol="♂",
        is_visible=True,
        vmag=0.0,
    )
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(100, 100), radius=80)
    celestial = _empty_celestial_data([sun, moon, mars])

    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=celestial,
        viewer_data=viewer,
        enlarge_moon=True,
        outline_bright_bodies=True,
        label_candidates=[],
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert moon_outline_radii == []
    assert moon_draw_radii == [12.5]
    assert len(planet_outline_radii) == 1
    assert all(radius > 0.0 for radius in planet_outline_radii)
    assert len(cross_calls) == 3


def test_hovered_moon_is_filled_even_in_outline_mode(monkeypatch) -> None:
    moon_outline_radii: list[float] = []
    moon_draw_radii: list[float] = []

    def fake_draw_moon_outline(_painter, _center, radius_px, _color) -> None:
        moon_outline_radii.append(float(radius_px))

    def fake_draw_moon(_painter, _center, radius_px, *_args, **_kwargs) -> None:
        moon_draw_radii.append(float(radius_px))

    def fake_draw_gauge_cross(_painter, _color, _center, **_kwargs) -> None:
        return None

    monkeypatch.setattr(render_solar_system, "draw_moon_outline", fake_draw_moon_outline)
    monkeypatch.setattr(render_solar_system, "draw_moon", fake_draw_moon)
    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", fake_draw_gauge_cross)

    sun = PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True)
    moon = PlanetBody(name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(100, 100), radius=80)
    celestial = _empty_celestial_data([sun, moon])

    render_solar_system.draw_hovered_moon_overlay(
        painter=object(),
        geometry=geometry,
        celestial_data=celestial,
        viewer_data=viewer,
        highlighted_object=(moon, QPointF(100.0, 100.0)),
        outline_bright_bodies=True,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert moon_outline_radii == []
    assert moon_draw_radii == [12.5]


def test_marker_scale_applies_to_planets_and_moon(monkeypatch) -> None:
    planet_draw_radii: list[float] = []
    moon_draw_radii: list[float] = []

    def fake_draw_planet_disc(_painter, _pos, _color, *, radius_px=1.0, alpha=255) -> None:
        planet_draw_radii.append(float(radius_px))

    def fake_draw_planet_bloom(_painter, _pos, _color, *, core_radius_px=1.0, vmag=None) -> None:
        return None

    def fake_draw_moon(_painter, _center, radius_px, *_args, **_kwargs) -> None:
        moon_draw_radii.append(float(radius_px))

    monkeypatch.setattr(render_solar_system, "draw_planet_disc", fake_draw_planet_disc)
    monkeypatch.setattr(render_solar_system, "draw_planet_bloom", fake_draw_planet_bloom)
    monkeypatch.setattr(render_solar_system, "draw_moon", fake_draw_moon)
    monkeypatch.setattr(
        render_solar_system, "draw_gauge_cross", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", lambda *_args, **_kwargs: None
    )

    sun = PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True)
    moon = PlanetBody(name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True)
    mars = PlanetBody(
        name="mars",
        alt=45.0,
        az=180.0,
        symbol="♂",
        is_visible=True,
        vmag=0.0,
    )
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(100, 100), radius=80)
    celestial = _empty_celestial_data([sun, moon, mars])

    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=celestial,
        viewer_data=viewer,
        enlarge_moon=False,
        marker_scale=1.0,
        theme=THEME_STYLES_BY_PRESET["night"],
    )
    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=celestial,
        viewer_data=viewer,
        enlarge_moon=False,
        marker_scale=2.0,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert len(planet_draw_radii) == 2
    assert planet_draw_radii[1] == planet_draw_radii[0] * 2.0
    assert len(moon_draw_radii) == 2
    assert moon_draw_radii[1] == moon_draw_radii[0] * 2.0


def test_planet_draw_and_hover_ignore_horizon_visibility_flag(monkeypatch) -> None:
    disc_calls: list[tuple[float, int, tuple[int, int, int, int]]] = []

    def fake_draw_planet_disc(
        _painter, _pos, color, *, radius_px=1.0, alpha=255
    ) -> None:
        disc_calls.append((radius_px, alpha, color.getRgb()))

    monkeypatch.setattr(render_solar_system, "draw_planet_disc", fake_draw_planet_disc)
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", lambda *_a, **_k: None
    )
    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", lambda *_a, **_k: None)
    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", lambda *_a, **_k: None
    )

    mars = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=False)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)

    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=_empty_celestial_data([mars]),
        viewer_data=viewer,
        enlarge_moon=False,
        theme=THEME_STYLES_BY_PRESET["night"],
    )
    assert len(disc_calls) == 1

    highlighted = render_stars.find_highlighted_object(
        _empty_celestial_data([mars]),
        viewer,
        QPoint(120, 90),
        geometry,
    )
    assert highlighted is not None


def test_draw_gauge_cross_respects_small_scale() -> None:
    image = QImage(64, 64, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        render_guides.draw_gauge_cross(
            painter,
            QColor(255, 255, 255, 255),
            QPointF(32.0, 32.0),
            scale=0.13,
            pen_width=1.0,
        )
    finally:
        painter.end()

    alpha = np.frombuffer(image.bits(), dtype=np.uint8).reshape((64, 64, 4))[..., 3]
    ys, xs = np.nonzero(alpha)
    assert len(xs) > 0
    half_span_x = max(abs(int(x) - 32) for x in xs)
    half_span_y = max(abs(int(y) - 32) for y in ys)
    assert half_span_x <= 2
    assert half_span_y <= 2


def test_satellite_overlay_draws_below_horizon_marker_when_in_fov(monkeypatch) -> None:
    cross_calls: list[tuple[float, float]] = []

    monkeypatch.setattr(
        render_satellites,
        "draw_gauge_cross",
        lambda _painter, _color, _center, *, scale=1.0, pen_width=1.0: (
            cross_calls.append((scale, pen_width))
        ),
    )

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_points=[
                SatelliteOverlayPoint(
                    group_key="iss",
                    satellite_name="ISS (ZARYA)",
                    alt_deg=-40.0,
                    az_deg=151.0,
                    marker_scale=0.42,
                )
            ],
            view_center=(0.0, 151.0),
            opacity=1.0,
            label_candidates=[],
            theme=THEME_STYLES_BY_PRESET["night"],
        )
    finally:
        painter.end()

    assert cross_calls == [(0.42, 1.0)]


def test_satellite_overlay_scales_marker_with_window_scale(monkeypatch) -> None:
    cross_calls: list[tuple[float, float]] = []

    monkeypatch.setattr(
        render_satellites,
        "draw_gauge_cross",
        lambda _painter, _color, _center, *, scale=1.0, pen_width=1.0: (
            cross_calls.append((scale, pen_width))
        ),
    )

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_points=[
                SatelliteOverlayPoint(
                    group_key="iss",
                    satellite_name="ISS (ZARYA)",
                    alt_deg=-40.0,
                    az_deg=151.0,
                    marker_scale=0.42,
                )
            ],
            view_center=(0.0, 151.0),
            opacity=1.0,
            label_candidates=[],
            marker_scale=2.5,
            theme=THEME_STYLES_BY_PRESET["night"],
        )
    finally:
        painter.end()

    assert cross_calls == [(1.05, 1.0)]


def test_satellite_overlay_keeps_overscan_position_beyond_90_deg(monkeypatch) -> None:
    positions: list[tuple[float, float]] = []

    def _record_cross(_painter, _color, center, *, scale=1.0, pen_width=1.0) -> None:
        positions.append((float(center.x()), float(center.y())))

    monkeypatch.setattr(render_satellites, "draw_gauge_cross", _record_cross)

    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(120, 120), radius=80),
            satellite_points=[
                SatelliteOverlayPoint(
                    group_key="iss",
                    satellite_name="ISS (ZARYA)",
                    alt_deg=-50.0,
                    az_deg=180.0,
                    marker_scale=0.42,
                )
            ],
            view_center=(45.0, 180.0),
            opacity=1.0,
            label_candidates=[],
            edge_fov_deg=110.0,
            content_fov_deg=110.0,
            theme=THEME_STYLES_BY_PRESET["night"],
        )
    finally:
        painter.end()

    assert len(positions) == 1
    x, y = positions[0]
    assert math.isclose(x, 120.0, abs_tol=0.2)
    assert math.isclose(y, 189.1, abs_tol=0.3)


app = QApplication.instance() or QApplication([])


def test_satellite_overlay_does_not_add_labels(monkeypatch) -> None:
    monkeypatch.setattr(
        render_satellites, "draw_gauge_cross", lambda *_args, **_kwargs: None
    )

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        label_candidates: list[dict[str, object]] = []
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_points=[
                SatelliteOverlayPoint(
                    group_key="iss",
                    satellite_name="ISS (ZARYA)",
                    alt_deg=10.0,
                    az_deg=151.0,
                    marker_scale=0.42,
                )
            ],
            view_center=(0.0, 151.0),
            opacity=1.0,
            label_candidates=label_candidates,
            theme=THEME_STYLES_BY_PRESET["night"],
        )
    finally:
        painter.end()

    assert label_candidates == []


def test_satellite_overlay_info_shows_hover_name(monkeypatch) -> None:
    class DummyPainter:
        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawEllipse(self, *_args, **_kwargs) -> None:
            raise AssertionError("satellite hover should not draw a circle")

    painter = DummyPainter()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    hovered_satellite = (
        SatelliteOverlayPoint(
            group_key="horizons",
            satellite_name="JWST",
            alt_deg=12.0,
            az_deg=220.0,
            marker_scale=0.3,
        ),
        QPointF(120.0, 90.0),
    )
    label_calls: list[tuple[str, QColor]] = []

    def fake_draw_outlined_text(_painter, text, *_args, **kwargs) -> None:
        style = kwargs["style"]
        label_calls.append((str(text), style.text_color))

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=None,
        highlighted_satellite=hovered_satellite,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert any(
        text == "JWST"
        and (color.red(), color.green(), color.blue())
        == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
        for text, color in label_calls
    )


def test_overlay_info_shows_star_and_satellite_labels_independently() -> None:
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
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    highlighted_object = ({"name": "Vega"}, QPointF(100.0, 80.0))
    highlighted_satellite = (
        SatelliteOverlayPoint(
            group_key="iss",
            satellite_name="ISS",
            alt_deg=12.0,
            az_deg=220.0,
            marker_scale=0.42,
        ),
        QPointF(130.0, 95.0),
    )
    label_calls: list[tuple[str, QColor]] = []

    def fake_draw_outlined_text(_painter, text, *_args, **kwargs) -> None:
        style = kwargs["style"]
        label_calls.append((str(text), style.text_color))

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=highlighted_object,
        highlighted_satellite=highlighted_satellite,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    label_names = {text for text, _color in label_calls}
    assert "Vega" in label_names
    assert "ISS" in label_names
    assert any(
        text == "ISS"
        and (color.red(), color.green(), color.blue())
        == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
        for text, color in label_calls
    )


def test_overlay_info_colors_dso_labels_like_the_dso_marker() -> None:
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
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    highlighted_dso = (
        {
            "name": "Andromeda Galaxy",
            "alt": 45.0,
            "az": 180.0,
            "major_arcmin": 30.0,
            "minor_arcmin": 20.0,
            "pa_deg": 0.0,
        },
        QPointF(120.0, 90.0),
    )
    label_calls: list[tuple[str, QColor]] = []

    def fake_draw_outlined_text(_painter, text, *_args, **kwargs) -> None:
        style = kwargs["style"]
        label_calls.append((str(text), style.text_color))

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=highlighted_dso,
        highlighted_object=None,
        highlighted_satellite=None,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert any(
        text == "Andromeda Galaxy"
        and (color.red(), color.green(), color.blue()) == DSO_LABEL_RGB
        for text, color in label_calls
    )


def test_find_highlighted_satellite_prefers_nearest_marker(monkeypatch) -> None:
    monkeypatch.setattr(render_satellites, "is_in_fov", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(
        render_satellites,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(alt), float(az)),
    )
    monkeypatch.setattr(
        render_satellites,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    geometry = ScreenGeometry(center=(120, 90), radius=70)
    view_center = (45.0, 180.0)
    satellite_points = [
        SatelliteOverlayPoint(
            group_key="iss",
            satellite_name="ISS",
            alt_deg=120.0,
            az_deg=90.0,
            marker_scale=0.42,
        ),
        SatelliteOverlayPoint(
            group_key="horizons",
            satellite_name="JWST",
            alt_deg=123.0,
            az_deg=90.0,
            marker_scale=0.30,
        ),
    ]

    highlighted = render_satellites.find_highlighted_satellite(
        satellite_points,
        QPointF(121.0, 90.0),
        geometry,
        view_center,
    )

    assert highlighted is not None
    satellite, _ = highlighted
    assert satellite.satellite_name == "ISS"


def test_aircraft_label_uses_black_theme_style_in_day_theme() -> None:
    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def font(self) -> QFont:
            return QFont()

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

    label_candidates: list[dict[str, object]] = []
    render_aircraft.draw_aircraft_overlay(
        painter=_Painter(),
        geometry=ScreenGeometry(center=(20, 20), radius=20),
        aircraft_points=[
            AircraftOverlayPoint(
                icao24="abc123",
                callsign="TEST123",
                alt_deg=10.0,
                az_deg=151.0,
                trail_alt_az_points=((10.0, 151.0), (10.2, 151.3)),
                distance_km=5.0,
                age_seconds=10.0,
                alpha_scale=1.0,
            )
        ],
        view_center=(0.0, 151.0),
        opacity=1.0,
        label_candidates=label_candidates,
        theme=THEME_STYLES_BY_PRESET["day"],
        content_fov_deg=180.0,
    )

    assert len(label_candidates) == 1
    style = label_candidates[0]["style"]
    expected_rgb = PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert (
        style.text_color.red(),
        style.text_color.green(),
        style.text_color.blue(),
    ) == expected_rgb
    assert (
        style.outline_width
        == render_text.resolve_text_style(THEME_STYLES_BY_PRESET["night"], QFont()).outline_width
    )


def test_hover_ignores_unnamed_stars() -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    mouse_pos = QPoint(120, 90)

    stars, meta = _star_table(names=["", "Sirius"])
    highlighted = render_stars.find_highlighted_object(
        _celestial_data_with_stars(stars, meta),
        viewer,
        mouse_pos,
        geometry,
    )

    assert highlighted is not None
    obj, _ = highlighted
    assert obj.get("name") == "Sirius"


def test_hover_returns_none_without_named_star() -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    mouse_pos = QPoint(120, 90)

    stars, meta = _star_table(names=["", ""])
    highlighted = render_stars.find_highlighted_object(
        _celestial_data_with_stars(stars, meta),
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
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    planet = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=True)
    highlighted_object = (planet, QPointF(120.0, 90.0))

    label_calls: list[str] = []

    def fake_draw_outlined_text(_painter, text, *_args, **_kwargs) -> None:
        label_calls.append(text)

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=highlighted_object,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert "mars" not in label_calls


def test_overlay_info_includes_location_height_and_explicit_observer_height(
    monkeypatch,
) -> None:
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
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
        show_observer_height=True,
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)

    label_calls: list[str] = []

    def fake_draw_outlined_text(_painter, text, *_args, **_kwargs) -> None:
        label_calls.append(text)

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=None,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert label_calls == [
        "t/Tokyo Skytree",
        "Lat: 35.00000, Lon: 139.00000 | Height: ground 35 m, building 634 m, add 12 m",
        "2026-02-27 00:00:00 UTC",
        "Alt 45°  Az 180° (S)",
    ]


def test_overlay_info_first_line_top_margin_matches_left_margin_when_cursor_is_lower_half(
    monkeypatch,
) -> None:
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
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
        show_observer_height=True,
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    text_font = QFont()

    first_label_pos = None

    def fake_draw_outlined_text(_painter, text, pos, *_args, **_kwargs) -> None:
        nonlocal first_label_pos
        if text == "t/Tokyo Skytree" and first_label_pos is None:
            first_label_pos = pos

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=None,
        text_font=text_font,
        viewport_rect=QRectF(0.0, 0.0, 240.0, 180.0),
        mouse_pos=QPoint(10, 170),
        bottom_left=False,
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert first_label_pos is not None
    fm = QFontMetrics(text_font)
    bounds = fm.tightBoundingRect("Ag")
    top_margin = float(first_label_pos.y()) + float(bounds.top())
    left_margin = float(fm.lineSpacing())
    assert abs(top_margin - left_margin) <= 1.0


def test_overlay_info_moves_to_bottom_when_cursor_is_in_upper_half(monkeypatch) -> None:
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
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
        show_observer_height=True,
    )

    first_label_pos = None

    def fake_draw_outlined_text(_painter, text, pos, *_args, **_kwargs) -> None:
        nonlocal first_label_pos
        if text == "t/Tokyo Skytree" and first_label_pos is None:
            first_label_pos = pos

    render_overlay_info.draw_overlay_info(
        painter,
        ScreenGeometry(center=(120, 90), radius=70),
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=None,
        text_font=QFont(),
        viewport_rect=QRectF(0.0, 0.0, 240.0, 180.0),
        mouse_pos=QPoint(10, 20),
        bottom_left=True,
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert first_label_pos is not None
    fm = QFontMetrics(QFont())
    bounds = fm.tightBoundingRect("Ag")
    line_height = int(fm.lineSpacing() * 1.2)
    bottom_margin = 180.0 - (
        float(first_label_pos.y()) + float(bounds.bottom()) + 3.0 * line_height
    )
    left_margin = float(fm.lineSpacing())
    assert abs(bottom_margin - left_margin) <= 2.0


def test_format_overlay_info_lines_matches_static_overlay_order() -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="t/Tokyo Skytree",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
        show_observer_height=True,
    )

    assert render_background.format_overlay_info_lines(
        _empty_celestial_data([]), viewer, 6.0
    ) == [
        "t/Tokyo Skytree",
        "Lat: 35.00000, Lon: 139.00000 | Height: ground 35 m, building 634 m, add 12 m",
        "2026-02-27 00:00:00 UTC",
        "Alt 45°  Az 180° (S)",
    ]
