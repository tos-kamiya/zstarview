
from tests._planet_marker_support import *

from datetime import datetime, timezone

from zstarview.solar_hover import SolarHoverImage


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


def test_sun_label_uses_desaturated_moon_tone(monkeypatch) -> None:
    label_calls: list[
        tuple[str, tuple[int, int, int, int], tuple[int, int, int, int], float]
    ] = []

    def fake_draw_outlined_text(
        _painter, text, _pos, _font, *_args, style=None, **_kwargs
    ):
        assert style is not None
        label_calls.append(
            (
                str(text),
                style.text_color.getRgb(),
                style.outline_color.getRgb(),
                float(style.outline_width),
            )
        )

    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", fake_draw_outlined_text
    )
    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", lambda *_a, **_k: None)
    monkeypatch.setattr(render_solar_system, "draw_planet_disc", lambda *_a, **_k: None)
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", lambda *_a, **_k: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_outline", lambda *_a, **_k: None
    )

    sun = PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True)
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
        celestial_data=_empty_celestial_data([sun]),
        viewer_data=viewer,
        enlarge_moon=False,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert len(label_calls) == 1
    text, text_rgb, _outline_rgb, _outline_width = label_calls[0]
    assert text == "Sun"
    assert (
        text_rgb[:3]
        == render_text.blend_color_toward_white(
            render_solar_system.planet_marker_color("moon"),
            amount=render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
        ).getRgb()[:3]
    )


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
    draw_order: list[str] = []

    def fake_draw_moon(_painter, _center, radius_px, *_args, **_kwargs) -> None:
        moon_draw_radii.append(float(radius_px))
        draw_order.append("moon")

    monkeypatch.setattr(render_solar_system, "draw_moon", fake_draw_moon)
    monkeypatch.setattr(
        render_solar_system,
        "draw_moon_outline",
        lambda *_args, **_kwargs: draw_order.append("limb"),
    )
    monkeypatch.setattr(
        render_solar_system,
        "draw_gauge_cross",
        lambda *_args, **_kwargs: draw_order.append("cross"),
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
    assert draw_order[-3:] == ["cross", "limb", "moon"]


def test_suppress_moon_marker_skips_base_moon_rendering(monkeypatch) -> None:
    moon_draw_calls: list[float] = []

    monkeypatch.setattr(
        render_solar_system,
        "draw_moon",
        lambda _painter, _center, radius_px, *_args, **_kwargs: moon_draw_calls.append(
            float(radius_px)
        ),
    )
    monkeypatch.setattr(
        render_solar_system, "draw_gauge_cross", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_disc", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", lambda *_args, **_kwargs: None
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

    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=geometry,
        celestial_data=_empty_celestial_data([sun, moon]),
        viewer_data=viewer,
        enlarge_moon=False,
        suppress_moon_marker=True,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert moon_draw_calls == []


def test_outline_bright_bodies_keeps_enlarged_moon_filled(monkeypatch) -> None:
    moon_outline_radii: list[float] = []
    moon_draw_radii: list[float] = []
    planet_outline_radii: list[float] = []
    cross_calls: list[tuple[float, float]] = []

    def fake_draw_moon_outline(
        _painter, _center, radius_px, _color, **_kwargs
    ) -> None:
        moon_outline_radii.append(float(radius_px))

    def fake_draw_moon(_painter, _center, radius_px, *_args, **_kwargs) -> None:
        moon_draw_radii.append(float(radius_px))

    def fake_draw_planet_outline(_painter, _pos, _color, *, radius_px=1.0) -> None:
        planet_outline_radii.append(float(radius_px))

    def fake_draw_gauge_cross(
        _painter, _color, _center, *, scale=1.0, pen_width=1.0
    ) -> None:
        cross_calls.append((float(scale), float(pen_width)))

    monkeypatch.setattr(
        render_solar_system, "draw_moon_outline", fake_draw_moon_outline
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_outline", fake_draw_planet_outline
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_disc", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", lambda *_args, **_kwargs: None
    )
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

    assert moon_outline_radii == [12.5]
    assert moon_draw_radii == [12.5]
    assert len(planet_outline_radii) == 1
    assert all(radius > 0.0 for radius in planet_outline_radii)
    assert len(cross_calls) == 3


def test_outline_bright_bodies_draws_moon_phase_outline(monkeypatch) -> None:
    phase_outline_calls: list[tuple[float, float]] = []

    def fake_draw_moon_phase_outline(
        _painter,
        _center,
        radius_px,
        *,
        sun_dir_in_moon_frame,
        screen_rotation_deg,
        color,
    ) -> None:
        del color
        phase_outline_calls.append(
            (float(radius_px), float(np.linalg.norm(sun_dir_in_moon_frame)))
        )
        assert math.isfinite(float(screen_rotation_deg))

    monkeypatch.setattr(
        render_solar_system,
        "draw_moon_phase_outline",
        fake_draw_moon_phase_outline,
    )
    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", lambda *_a, **_k: None)
    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", lambda *_a, **_k: None
    )

    sun = PlanetBody(name="sun", alt=20.0, az=180.0, symbol="☉", is_visible=True)
    moon = PlanetBody(name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )

    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        celestial_data=_empty_celestial_data([sun, moon]),
        viewer_data=viewer,
        enlarge_moon=False,
        outline_bright_bodies=True,
        label_candidates=[],
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert len(phase_outline_calls) == 1
    assert phase_outline_calls[0][0] == 2.5
    assert phase_outline_calls[0][1] > 0.0


def test_instrument_presentation_draws_planets_with_outlines_and_sun_as_cross(
    monkeypatch,
) -> None:
    disc_alphas: list[int] = []
    moon_draw_calls: list[float] = []
    outline_calls: list[float] = []
    cross_scales: list[float] = []

    monkeypatch.setattr(
        render_solar_system,
        "draw_planet_disc",
        lambda *_args, alpha=0, **_kwargs: disc_alphas.append(int(alpha)),
    )
    monkeypatch.setattr(
        render_solar_system,
        "draw_moon",
        lambda _painter, _center, radius_px, *_args, **_kwargs: moon_draw_calls.append(
            float(radius_px)
        ),
    )
    monkeypatch.setattr(
        render_solar_system,
        "draw_planet_outline",
        lambda _painter, _pos, _color, *, radius_px=1.0: outline_calls.append(
            float(radius_px)
        ),
    )
    monkeypatch.setattr(
        render_solar_system, "draw_moon_outline", lambda *_a, **_k: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", lambda *_a, **_k: None
    )
    monkeypatch.setattr(
        render_solar_system,
        "draw_gauge_cross",
        lambda _painter, _color, _center, *, scale=1.0, pen_width=1.0: (
            cross_scales.append(float(scale))
        ),
    )

    sun = PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True)
    moon = PlanetBody(name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True)
    mars = PlanetBody(
        name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=True, vmag=0.0
    )
    labels: list[dict[str, object]] = []
    render_solar_system.draw_solar_system_bodies(
        painter=object(),
        geometry=ScreenGeometry(center=(100, 100), radius=80),
        celestial_data=_empty_celestial_data([sun, moon, mars]),
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
        ),
        enlarge_moon=False,
        outline_bright_bodies=True,
        label_candidates=labels,
        theme=THEME_STYLES_BY_PRESET["atlas-white"],
        instrument_presentation=True,
    )

    assert disc_alphas == [255]
    assert len(moon_draw_calls) == 1
    assert len(outline_calls) == 1
    assert cross_scales == [1.0, 1.0, 0.55]
    assert len(labels) == 3
    assert all(
        candidate["style"].text_color.getRgb() == (24, 24, 24, 255)
        for candidate in labels
    )


def test_hovered_moon_is_filled_even_in_outline_mode(monkeypatch) -> None:
    moon_outline_radii: list[float] = []
    moon_draw_radii: list[float] = []
    cross_calls: list[bool] = []

    def fake_draw_moon_outline(
        _painter, _center, radius_px, _color, **_kwargs
    ) -> None:
        moon_outline_radii.append(float(radius_px))

    def fake_draw_moon(_painter, _center, radius_px, *_args, **_kwargs) -> None:
        moon_draw_radii.append(float(radius_px))

    def fake_draw_gauge_cross(_painter, _color, _center, **_kwargs) -> None:
        cross_calls.append(True)

    monkeypatch.setattr(
        render_solar_system, "draw_moon_outline", fake_draw_moon_outline
    )
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

    assert moon_outline_radii == [12.5]
    assert moon_draw_radii == [12.5]
    assert cross_calls == [True]


def test_hovered_sun_keeps_external_image_above_cross(monkeypatch) -> None:
    draw_order: list[str] = []

    def fake_draw_hmi_solar_image(*_args, **_kwargs) -> float:
        draw_order.append("image")
        return 20.0

    def fake_draw_gauge_cross(_painter, _color, _center, **_kwargs) -> None:
        draw_order.append("cross")

    monkeypatch.setattr(
        render_solar_system, "draw_hmi_solar_image", fake_draw_hmi_solar_image
    )
    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", fake_draw_gauge_cross)
    monkeypatch.setattr(
        render_solar_system,
        "calculate_solar_north_up_screen_rotation",
        lambda *_args, **_kwargs: 0.0,
    )
    monkeypatch.setattr(
        render_solar_system, "bundled_aod550_or_default", lambda *_args: 0.15
    )
    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", lambda *_args, **_kwargs: None
    )

    sun = PlanetBody(name="sun", alt=45.0, az=180.0, symbol="sun", is_visible=True)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(100, 100), radius=80)
    image = SolarHoverImage(
        image=QImage(8, 8, QImage.Format.Format_ARGB32),
        time_utc=datetime(2026, 2, 27, tzinfo=timezone.utc),
        source_radius_px=4.0,
        image_id=1,
    )

    render_solar_system.draw_hovered_sun_overlay(
        painter=object(),
        geometry=geometry,
        celestial_data=_empty_celestial_data([sun]),
        viewer_data=viewer,
        highlighted_object=(sun, QPointF(100.0, 100.0)),
        time_obj=None,
        marker_scale=2.0,
        text_font=QFont(),
        theme=THEME_STYLES_BY_PRESET["night"],
        external_solar_image=image,
    )

    assert draw_order == ["cross", "image"]


def test_marker_scale_applies_to_planets_and_moon(monkeypatch) -> None:
    planet_draw_radii: list[float] = []
    moon_draw_radii: list[float] = []

    def fake_draw_planet_disc(
        _painter, _pos, _color, *, radius_px=1.0, alpha=255
    ) -> None:
        planet_draw_radii.append(float(radius_px))

    def fake_draw_planet_bloom(
        _painter, _pos, _color, *, core_radius_px=1.0, vmag=None
    ) -> None:
        return None

    def fake_draw_moon(_painter, _center, radius_px, *_args, **_kwargs) -> None:
        moon_draw_radii.append(float(radius_px))

    monkeypatch.setattr(render_solar_system, "draw_planet_disc", fake_draw_planet_disc)
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", fake_draw_planet_bloom
    )
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


def test_day_and_white_themes_use_desaturated_planet_colors_for_solar_system_labels(
    monkeypatch,
) -> None:
    cross_colors: list[tuple[int, int, int]] = []
    label_colors: dict[str, tuple[int, int, int]] = {}
    label_outlines: list[tuple[int, int, int, int, float]] = []

    def fake_draw_gauge_cross(_painter, color, _center, *, scale=1.0, pen_width=1.0):
        cross_colors.append((color.red(), color.green(), color.blue()))

    def fake_draw_outlined_text(
        _painter, text, _pos, _font, *_args, style=None, **_kwargs
    ):
        if style is not None:
            label_colors[str(text)] = (
                style.text_color.red(),
                style.text_color.green(),
                style.text_color.blue(),
            )
            label_outlines.append((*style.outline_color.getRgb(), style.outline_width))

    monkeypatch.setattr(render_solar_system, "draw_gauge_cross", fake_draw_gauge_cross)
    monkeypatch.setattr(
        render_solar_system, "draw_outlined_text", fake_draw_outlined_text
    )
    monkeypatch.setattr(
        render_solar_system, "draw_moon", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_disc", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_bloom", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_solar_system, "draw_planet_outline", lambda *_args, **_kwargs: None
    )

    sun = PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True)
    moon = PlanetBody(name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True)
    mercury = PlanetBody(
        name="mercury",
        alt=45.0,
        az=180.0,
        symbol="☿",
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
    celestial = _empty_celestial_data([sun, moon, mercury])

    for preset in ("day", "white"):
        cross_colors.clear()
        label_colors.clear()
        render_solar_system.draw_solar_system_bodies(
            painter=object(),
            geometry=geometry,
            celestial_data=celestial,
            viewer_data=viewer,
            enlarge_moon=False,
            theme=THEME_STYLES_BY_PRESET[preset],
        )

        assert cross_colors
        assert all(color == (180, 180, 180) for color in cross_colors)
        expected_sun_rgb = render_text.blend_color_toward_white(
            render_solar_system.planet_marker_color("moon"),
            amount=render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
        ).getRgb()[:3]
        expected_moon_rgb = render_text.blend_color_toward_white(
            render_solar_system.planet_marker_color("moon"),
            amount=render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
        ).getRgb()[:3]
        expected_mercury_rgb = render_text.blend_color_toward_white(
            render_solar_system.planet_marker_color("mercury"),
            amount=render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
        ).getRgb()[:3]
        assert label_colors["Sun"] == expected_sun_rgb
        assert label_colors["Moon"] == expected_moon_rgb
        assert label_colors["Mercury"] == expected_mercury_rgb
        assert label_outlines
        assert all(
            r == 0 and g == 0 and b == 0 and a == 76 and width == 3.0
            for r, g, b, a, width in label_outlines
        )


def test_day_and_white_label_outlines_match_night_theme() -> None:
    night_style = render_text.resolve_label_text_style(
        THEME_STYLES_BY_PRESET["night"],
        QFont(),
    )
    night_outline = night_style.outline_color.getRgb()

    for preset in ("day", "white"):
        style = render_text.resolve_label_text_style(
            THEME_STYLES_BY_PRESET[preset],
            QFont(),
        )
        assert style.outline_color.getRgb() == night_outline
        assert style.outline_width == night_style.outline_width


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
