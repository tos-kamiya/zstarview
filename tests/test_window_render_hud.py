
from tests._window_render_support import *

import zstarview.render.meteors as render_meteors_module


def test_render_hud_overlay_draws_persistent_search_label(monkeypatch) -> None:
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        pipeline_module,
        "_draw_hover_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_static_observation_overlay",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_status_line_text",
        lambda *_args, **_kwargs: None,
    )

    def fake_draw_search_target_overlay(*_args, **kwargs) -> None:
        captured.update(kwargs)
        label_candidates = kwargs.get("label_candidates")
        if isinstance(label_candidates, list):
            label_candidates.append({"text": "Dubhe", "priority": 15})

    monkeypatch.setattr(
        pipeline_module.render_search_overlay,
        "draw_search_target_overlay",
        fake_draw_search_target_overlay,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args: captured.update({"labels": _args[1]}),
    )

    scene = replace(
        _make_scene(),
        viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        ),
    )
    pipeline_module.render_hud_overlay_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(100, 100), radius=80),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(star_render_expected_width=600),
        hud=_make_hud(),
        highlighted_object=None,
        highlighted_dso=None,
        label_candidates=[],
        search_overlay_target=SearchJumpTarget(
            label="Dubhe",
            kind="star",
            sort_key=(0.0, "dubhe"),
            alt_deg=50.0,
            az_deg=20.0,
            persistent_keep_marker=True,
        ),
    )

    assert captured["draw_label"] is True
    assert captured["labels"] == [{"text": "Dubhe", "priority": 15}]


def test_collect_visible_named_star_labels_returns_only_named_visible_stars() -> None:
    scene = _make_scene(
        celestial_data=CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[],
            stars={
                "star_index": np.array([11, 12], dtype=np.int32),
                "alt": np.array([45.0, 10.0]),
                "az": np.array([180.0, 10.0]),
                "vmag": np.array([5.96, 1.0]),
                "bv": np.array([0.869, 0.0]),
                "size_factor": np.array([1.0, 1.0]),
                "color_factor_base": np.array([1.0, 1.0]),
            },
            deep_sky_objects={},
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
            star_catalog_meta=StarCatalogMeta(
                name_indices=np.array([11], dtype=np.int32),
                names=np.array(["Dubhe"], dtype=object),
                source_id_indices=np.array([], dtype=np.int32),
                source_ids=np.array([], dtype=object),
            ),
        ),
    )

    labels = pipeline_module.render_stars.collect_visible_named_star_labels(
        scene.celestial_data,
        scene.viewer,
        SimpleNamespace(center=(200, 200), radius=200),
        star_base_radius=4.0,
        draw_vmag_limit=6.0,
        viewport_size=(400, 400),
    )

    assert len(labels) == 1
    assert labels[0][0] == "Dubhe"
    assert labels[0][1].x() == pytest.approx(200.0)
    assert labels[0][1].y() == pytest.approx(200.0)
    assert labels[0][2] == (255, 210, 161)


def test_render_hud_overlay_draws_simplified_named_star_labels_at_fixed_offset(
    monkeypatch,
) -> None:
    captured: list[tuple[str, float, float, int]] = []

    monkeypatch.setattr(
        pipeline_module,
        "_draw_hover_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_static_observation_overlay",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_status_line_text",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "collect_visible_named_star_labels",
        lambda *_args, **_kwargs: [
            ("Dubhe", QPointF(120.0, 80.0), (10, 20, 30)),
            ("Merak", QPointF(150.0, 100.0), (40, 50, 60)),
        ],
    )

    def fake_draw_outlined_text(*_args, **kwargs) -> None:
        text = _args[1]
        pos = _args[2]
        captured.append(
            (text, float(pos.x()), float(pos.y()), kwargs["style"].text_color.alpha())
        )
        style = kwargs["style"]
        assert style.font is not None
        assert style.outline_width == 0.0
        assert style.outline_color.alpha() == 0
        if text == "Dubhe":
            expected = pipeline_module.render_text.blend_color_toward_white(
                QColor(10, 20, 30),
                amount=pipeline_module.render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
            )
            assert style.text_color.red() == expected.red()
            assert style.text_color.green() == expected.green()
            assert style.text_color.blue() == expected.blue()
            assert style.text_color.alpha() == round(255 * 0.4)
        if text == "Merak":
            expected = pipeline_module.render_text.blend_color_toward_white(
                QColor(40, 50, 60),
                amount=pipeline_module.render_text.LABEL_COLOR_WHITE_BLEND_AMOUNT,
            )
            assert style.text_color.red() == expected.red()
            assert style.text_color.green() == expected.green()
            assert style.text_color.blue() == expected.blue()
            assert style.text_color.alpha() == round(255 * 0.4)

    monkeypatch.setattr(
        pipeline_module.render_text,
        "draw_outlined_text",
        fake_draw_outlined_text,
    )

    scene = _make_scene()
    style = _make_style(text_font=QFont(), status_line_font=QFont())
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
    pipeline_module.render_hud_overlay_into_painter(
        painter=painter,
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(200, 200), radius=200),
            QRect(0, 0, 400, 400),
        ),
        scene=scene,
        style=style,
        hud=_make_hud(simplified_view_enabled=True),
        highlighted_object=None,
        highlighted_dso=None,
        label_candidates=[],
    )
    painter.end()

    assert len(captured) == 2
    expected_positions = {
        "Dubhe": QPointF(120.0, 80.0),
        "Merak": QPointF(150.0, 100.0),
    }
    for text, x, y, alpha in captured:
        assert alpha == round(255 * 0.4)
        bounds = pipeline_module.render_text._text_bounds_at_baseline(
            text,
            QFont(),
            QPointF(0.0, 0.0),
        )
        expected = expected_positions[text]
        assert x + float(bounds.left()) == pytest.approx(float(expected.x()))
        assert y + float(bounds.bottom()) == pytest.approx(float(expected.y()))


def test_render_hud_overlay_skips_simplified_labels_when_disabled(monkeypatch) -> None:
    labels_drawn: list[str] = []

    monkeypatch.setattr(
        pipeline_module,
        "_draw_hover_overlay_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_static_observation_overlay",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_status_line_text",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "collect_visible_named_star_labels",
        lambda *_args, **_kwargs: [
            ("Dubhe", QPointF(120.0, 80.0), (10, 20, 30)),
        ],
    )

    def fake_draw_outlined_text(*_args, **_kwargs) -> None:
        labels_drawn.append(str(_args[1]))

    monkeypatch.setattr(
        pipeline_module.render_text,
        "draw_outlined_text",
        fake_draw_outlined_text,
    )

    scene = _make_scene()
    style = _make_style(text_font=QFont(), status_line_font=QFont())
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
    pipeline_module.render_hud_overlay_into_painter(
        painter=painter,
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(200, 200), radius=200),
            QRect(0, 0, 400, 400),
        ),
        scene=scene,
        style=style,
        hud=_make_hud(
            simplified_view_enabled=True, simplified_view_labels_enabled=False
        ),
        highlighted_object=None,
        highlighted_dso=None,
        label_candidates=[],
    )
    painter.end()

    assert labels_drawn == []


def test_render_fast_overlay_layers_skip_meteors_in_fast_mode(monkeypatch) -> None:
    draw_meteors = Mock()
    monkeypatch.setattr(render_meteors_module, "draw_meteor_trails", draw_meteors)
    monkeypatch.setattr(
        pipeline_module, "_draw_satellite_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_aircraft_layer", lambda *_args, **_kwargs: None
    )

    scene = _make_scene()
    style = _make_style(
        satellite_opacity=0.0,
        aircraft_opacity=0.0,
        meteor_opacity=0.4,
        tropical_cyclone_opacity=0.0,
        text_font=QFont(),
    )
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
    try:
        pipeline_module.render_fast_overlay_layers_into_painter(
            painter=painter,
            frame=_make_frame(
                scene,
                SimpleNamespace(center=(200, 200), radius=200),
                QRect(0, 0, 400, 400),
            ),
            scene=scene,
            style=style,
            draw_labels=False,
            fast_mode=True,
        )
    finally:
        painter.end()

    draw_meteors.assert_not_called()


def test_render_fast_overlay_layers_passes_simplified_satellite_labels(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        pipeline_module.render_satellites,
        "draw_satellite_overlay",
        lambda *_args, **kwargs: captured.update(
            {"draw_simplified_labels": kwargs.get("draw_simplified_labels")}
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_aircraft,
        "draw_aircraft_overlay",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_tropical_cyclones,
        "draw_tropical_cyclone_overlay",
        lambda *_args, **_kwargs: None,
    )

    scene = _make_scene()
    style = _make_style(
        satellite_opacity=1.0,
        aircraft_opacity=0.0,
        tropical_cyclone_opacity=0.0,
        text_font=QFont(),
    )
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
    try:
        pipeline_module.render_fast_overlay_layers_into_painter(
            painter=painter,
            frame=_make_frame(
                scene,
                SimpleNamespace(center=(200, 200), radius=200),
                QRect(0, 0, 400, 400),
            ),
            scene=scene,
            style=style,
            draw_labels=False,
            draw_simplified_satellite_labels=True,
        )
    finally:
        painter.end()

    assert captured == {"draw_simplified_labels": True}


def test_render_hud_overlay_forwards_simplified_satellite_label_flag(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        pipeline_module.render_overlay_info,
        "draw_overlay_info",
        lambda *_args, **kwargs: captured.update(
            {
                "draw_simplified_satellite_labels": kwargs[
                    "draw_simplified_satellite_labels"
                ],
                "highlighted_satellite": _args[9],
            }
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_simplified_named_star_labels",
        lambda *_args, **_kwargs: None,
    )

    scene = _make_scene()
    style = _make_style(text_font=QFont())
    img = QImage(400, 400, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(img)
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
    try:
        pipeline_module.render_hud_overlay_into_painter(
            painter=painter,
            frame=_make_frame(
                scene,
                SimpleNamespace(center=(200, 200), radius=200),
                QRect(0, 0, 400, 400),
            ),
            scene=scene,
            style=style,
            hud=_make_hud(
                simplified_view_enabled=True,
                simplified_view_labels_enabled=True,
            ),
            highlighted_object=None,
            highlighted_dso=None,
            highlighted_satellite=highlighted_satellite,
            label_candidates=[],
        )
    finally:
        painter.end()

    assert captured["draw_simplified_satellite_labels"] is True
    assert captured["highlighted_satellite"] == highlighted_satellite


def test_render_scene_hides_sky_bitmap_during_viewport_interaction(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}
    scene = replace(_make_scene(), sky_disc_image=object())

    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_clear_background_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_background_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_guide_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **kwargs: captured.update(
            {
                "cloud_disc_alpha": kwargs["style"].cloud_disc_alpha,
                "sky_disc_image": kwargs["scene"].sky_disc_image,
                "draw_sky_disc": kwargs["draw_sky_disc"],
            }
        ),
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_viewport_interaction_layers",
        lambda *_args, **_kwargs: None,
    )

    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(100, 100), radius=80),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(cloud_disc_alpha=0.2),
        hud=_make_hud(viewport_interaction_mode=True),
        compositor=object(),
    )

    assert captured == {
        "cloud_disc_alpha": 0.0,
        "sky_disc_image": scene.sky_disc_image,
        "draw_sky_disc": False,
    }


def test_draw_guide_layer_draws_zenith_marker(monkeypatch) -> None:
    calls: list[str] = []

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: calls.append("direction"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: calls.append("zenith"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_celestial_pole_markers",
        lambda *_args, **_kwargs: calls.append("poles"),
    )

    pipeline_module._draw_guide_layer(
        painter=object(),
        geometry=SimpleNamespace(center=(100, 100), radius=80),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(),
        style=_make_style(show_guidelines=True),
    )

    assert calls == ["direction", "zenith", "poles"]


def test_draw_instrument_guide_layer_draws_reference_lines(monkeypatch) -> None:
    calls: list[str] = []

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: calls.append("reference"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_grid_overlay",
        lambda *_args, **_kwargs: calls.append("grid"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_guide_layer",
        lambda *_args, **_kwargs: calls.append("annotations"),
    )

    atlas_pipeline_module._draw_instrument_guide_layer(
        painter=object(),
        geometry=SimpleNamespace(center=(100, 100), radius=80),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(),
        style=_make_style(show_guidelines=True),
    )

    assert calls == ["grid", "reference", "annotations"]


def test_render_base_scene_can_skip_fast_overlays(monkeypatch) -> None:
    calls: list[str] = []
    planet_kwargs: dict[str, object] = {}
    scene = _make_scene()
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    viewport_rect = SimpleNamespace(width=lambda: 200, height=lambda: 200)

    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_clear_background_layer",
        lambda *_args, **_kwargs: calls.append("clear"),
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_background_layer",
        lambda *_args, **_kwargs: calls.append("background"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_guide_layer",
        lambda *_args, **_kwargs: calls.append("guide"),
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **_kwargs: calls.append("sky-cloud"),
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_terrain_layers",
        lambda *_args, **_kwargs: calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: calls.append("stars"),
    )

    def fake_draw_planet_layer(*_args, **kwargs) -> None:
        planet_kwargs.update(kwargs)
        calls.append("planets")

    monkeypatch.setattr(pipeline_module, "_draw_planet_layer", fake_draw_planet_layer)
    monkeypatch.setattr(
        pipeline_module,
        "_draw_satellite_layer",
        lambda *_args, **_kwargs: calls.append("satellites"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_aircraft_layer",
        lambda *_args, **_kwargs: calls.append("aircraft"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "render_text",
        pipeline_module.render_text,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: calls.append("labels"),
    )

    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(scene, geometry, viewport_rect),
        scene=scene,
        style=_make_style(),
        hud=_make_hud(),
        compositor=object(),
        draw_fast_overlays=False,
    )

    assert calls == [
        "clear",
        "background",
        "sky-cloud",
        "guide",
        "terrain",
        "stars",
        "planets",
        "labels",
    ]
    assert planet_kwargs.get("draw_markers", True) is True


def test_draw_planet_layer_passes_marker_scale(monkeypatch) -> None:
    seen_marker_scales: list[float] = []

    monkeypatch.setattr(
        pipeline_module.render_solar_system,
        "draw_solar_system_bodies",
        lambda *_args, **kwargs: seen_marker_scales.append(
            float(kwargs.get("marker_scale", 1.0))
        ),
    )

    pipeline_module._draw_planet_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=_make_scene(),
        style=_make_style(star_render_expected_width=600),
        enlarge_moon=False,
        label_candidates=[],
    )

    assert seen_marker_scales == [
        pipeline_module.compute_star_render_upscale_factor(1200, 600)
    ]


def test_draw_hover_overlay_passes_marker_scale_to_moon(monkeypatch) -> None:
    seen_marker_scales: list[float] = []

    monkeypatch.setattr(
        render_solar_system_module,
        "draw_hovered_moon_overlay",
        lambda *_args, **kwargs: seen_marker_scales.append(
            float(kwargs.get("marker_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        render_overlay_info_module,
        "draw_overlay_info",
        lambda *_args, **_kwargs: None,
    )

    scene = _make_scene(
        celestial_data=CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[
                PlanetBody(name="sun", alt=45.0, az=180.0, symbol="☉", is_visible=True),
                PlanetBody(
                    name="moon", alt=45.0, az=180.0, symbol="☾", is_visible=True
                ),
            ],
            stars={
                "star_index": np.array([], dtype=np.int64),
                "alt": np.array([], dtype=np.float64),
                "az": np.array([], dtype=np.float64),
                "vmag": np.array([], dtype=np.float64),
                "bv": np.array([], dtype=np.float64),
                "size_factor": np.array([], dtype=np.float64),
                "color_factor_base": np.array([], dtype=np.float64),
            },
            deep_sky_objects={
                "id": np.array([], dtype=np.int64),
                "name": np.array([], dtype=object),
                "type": np.array([], dtype=object),
                "alt": np.array([], dtype=np.float64),
                "az": np.array([], dtype=np.float64),
                "vmag": np.array([], dtype=np.float64),
                "major_arcmin": np.array([], dtype=np.float64),
                "minor_arcmin": np.array([], dtype=np.float64),
                "pa_deg": np.array([], dtype=np.float64),
            },
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
        )
    )

    pipeline_module._draw_hover_overlay_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=QRect(0, 0, 1200, 1200),
        scene=scene,
        style=_make_style(
            star_render_expected_width=600,
            show_asterisms=False,
            show_observation_info=False,
        ),
        highlighted_object=({"name": "moon"}, QPointF(10.0, 10.0)),
        highlighted_dso=None,
    )

    assert seen_marker_scales == [
        pipeline_module.compute_star_render_upscale_factor(1200, 600)
    ]


def test_draw_hover_overlay_requires_direction_marker_hover(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "resolve_direction_marker_hover",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_altitude_ring_overlay",
        lambda *_args, **_kwargs: calls.append("background-rings"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_grid_overlay",
        lambda *_args, **_kwargs: calls.append("grid"),
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_solar_system,
        "draw_hovered_moon_overlay",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_dso_hover_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        render_overlay_info_module,
        "draw_overlay_info",
        lambda *_args, **_kwargs: None,
    )

    pipeline_module._draw_hover_overlay_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=QRect(0, 0, 1200, 1200),
        scene=_make_scene(celestial_data=object()),
        style=_make_style(
            show_guidelines=True,
            sky_disc_altaz_rings_hover="altaz",
        ),
        highlighted_object=None,
        highlighted_dso=None,
        mouse_pos=QPoint(600, 600),
    )

    assert calls == []


def test_draw_static_observation_overlay_skips_static_info_when_disabled(
    monkeypatch,
) -> None:
    draw_overlay_info = Mock()
    monkeypatch.setattr(
        pipeline_module.render_overlay_info, "draw_overlay_info", draw_overlay_info
    )

    pipeline_module._draw_static_observation_overlay(
        painter=object(),
        geometry=SimpleNamespace(radius=100),
        viewport_rect=SimpleNamespace(height=lambda: 200),
        scene=_make_scene(),
        style=_make_style(show_observation_info=False),
        mouse_pos=None,
        overlay_info_bottom_left=False,
        highlighted_object=None,
        highlighted_dso=None,
        enlarge_moon=False,
        label_reservations=[],
        label_candidates=[],
    )

    draw_overlay_info.assert_not_called()


def test_draw_background_layer_skips_gradient_when_disabled(monkeypatch) -> None:
    draw_radial_background = Mock()
    draw_window_border = Mock()
    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_radial_background",
        draw_radial_background,
    )
    monkeypatch.setattr(
        pipeline_module.render_background, "draw_window_border", draw_window_border
    )

    zstarview_pipeline_module._draw_background_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=100),
        viewport_rect=SimpleNamespace(),
        scene=_make_scene(),
        style=_make_style(show_background_gradient=False),
    )

    draw_radial_background.assert_not_called()
    draw_window_border.assert_not_called()


def test_draw_background_layer_skips_custom_frame_when_disabled(monkeypatch) -> None:
    draw_radial_background = Mock()
    draw_window_border = Mock()
    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_radial_background",
        draw_radial_background,
    )
    monkeypatch.setattr(
        pipeline_module.render_background, "draw_window_border", draw_window_border
    )

    zstarview_pipeline_module._draw_background_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=100),
        viewport_rect=QRect(0, 0, 200, 200),
        scene=_make_scene(),
        style=_make_style(show_custom_window_frame=False),
    )

    draw_radial_background.assert_called_once()
    draw_window_border.assert_not_called()


def test_draw_hover_overlay_layer_enlarges_hovered_moon_by_name(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_solar_system,
        "draw_hovered_moon_overlay",
        lambda *_args, **_kwargs: calls.append("moon-hover"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_dso_hover_layer",
        lambda *_args, **_kwargs: calls.append("dso-hover"),
    )
    monkeypatch.setattr(
        render_overlay_info_module,
        "draw_overlay_info",
        lambda *_args, **_kwargs: calls.append("overlay-info"),
    )

    pipeline_module._draw_hover_overlay_layer(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=QRect(0, 0, 1200, 1200),
        scene=_make_scene(celestial_data=object()),
        style=_make_style(),
        highlighted_object=({"name": "moon"}, object()),
        highlighted_dso=None,
    )

    assert calls == ["moon-hover", "dso-hover", "overlay-info"]


def test_draw_sky_reference_lines_uses_render_view_center_projection() -> None:
    calls: list[tuple[float, float]] = []

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

    render_guides_module.draw_sky_reference_lines(
        _Painter(),
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        celestial_data=SimpleNamespace(
            celestial_equator_points=[(1.0, 2.0)],
            ecliptic_points=[(3.0, 4.0)],
            horizon_points=[(5.0, 6.0)],
        ),
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(55.0, 200.0),
            observer_height_m=10.0,
        ),
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (
            calls.append(_view_center) or (alt, az)
        ),
    )

    assert len(calls) == 6
    assert calls == [(55.0, 200.0)] * 6


def test_draw_sky_reference_lines_uses_wider_dash_patterns(monkeypatch) -> None:
    dash_patterns: list[list[int]] = []
    pen_styles: list[object] = []
    pen_alpha_values: list[int] = []
    pen_widths: list[float] = []

    class _FakePen:
        def __init__(self, color, _width, style=None) -> None:
            self._dash_pattern: list[int] = []
            self._style = style
            self._alpha = color.alpha() if hasattr(color, "alpha") else None
            self._width = float(_width)

        def setCosmetic(self, *_args, **_kwargs) -> None:
            pass

        def setDashPattern(self, pattern) -> None:
            self._dash_pattern = list(pattern)

        def setCapStyle(self, *_args, **_kwargs) -> None:
            pass

        def setJoinStyle(self, *_args, **_kwargs) -> None:
            pass

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            dash_patterns.append(list(getattr(pen, "_dash_pattern", [])))
            pen_styles.append(getattr(pen, "_style", None))
            pen_alpha_values.append(getattr(pen, "_alpha", None))
            pen_widths.append(getattr(pen, "_width", None))

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

    monkeypatch.setattr(render_guides_module, "QPen", _FakePen)

    render_guides_module.draw_sky_reference_lines(
        _Painter(),
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        celestial_data=SimpleNamespace(
            celestial_equator_points=[(0.10, 0.20), (0.12, 0.22)],
            ecliptic_points=[(0.30, 0.40), (0.32, 0.42)],
            horizon_points=[(0.50, 0.60), (0.52, 0.62)],
        ),
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(55.0, 200.0),
            observer_height_m=10.0,
        ),
        altaz_to_normalized_xy_func=lambda alt, az, _view_center, **_kwargs: (alt, az),
    )

    assert dash_patterns[0::3] == [[], [], []]
    assert dash_patterns[1::3] == [[], [], []]
    assert dash_patterns[2::3] == [[16, 6], [4, 6], [10, 1]]
    assert pen_styles[0::3] == [
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
    ]
    assert pen_styles[1::3] == [
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
    ]
    assert pen_styles[2::3] == [None, None, None]
    assert pen_alpha_values[0::3] == [18, 18, 18]
    assert pen_alpha_values[1::3] == [30, 30, 30]
    assert pen_alpha_values[2::3] == [255, 255, 255]
    assert [round(width, 3) for width in pen_widths[0::3]] == [1.100, 1.254, 1.100]
    assert [round(width, 3) for width in pen_widths[1::3]] == [0.750, 0.855, 0.750]
    assert [round(width, 3) for width in pen_widths[2::3]] == [0.510, 0.627, 0.550]


def test_draw_direction_labels_uses_horizon_line_color(monkeypatch) -> None:
    seen_colors: list[tuple[int, int, int]] = []

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setFont(self, _font) -> None:
            pass

        def setPen(self, _pen) -> None:
            pass

        def drawLine(self, *_args, **_kwargs) -> None:
            pass

        def viewport(self):
            return QRect(0, 0, 200, 200)

    def _capture(*_args, style=None, **_kwargs):
        seen_colors.append(
            (style.text_color.red(), style.text_color.green(), style.text_color.blue())
        )
        assert style.outline_width == 0.0
        assert style.outline_color.alpha() == 0

    monkeypatch.setattr(render_guides_module, "draw_outlined_text", _capture)
    monkeypatch.setattr(
        render_guides_module, "is_in_fov", lambda *_args, **_kwargs: True
    )
    monkeypatch.setattr(
        render_guides_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_guides_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    render_guides_module.draw_direction_labels(
        _Painter(),
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(30.0, 40.0),
            edge_fov_deg=95.0,
            content_fov_deg=180.0,
        ),
        text_font=QFont(),
        theme=window_module.THEME_STYLES_BY_PRESET["white"],
    )

    assert seen_colors
    assert all(
        color == render_guides_module.HORIZON_LINE_COLOR for color in seen_colors
    )


def test_draw_zenith_marker_uses_direction_grid_color_for_all_themes(
    monkeypatch,
) -> None:
    seen_colors: list[tuple[int, int, int]] = []

    class _FakePen:
        def __init__(self, color, _width) -> None:
            seen_colors.append((color.red(), color.green(), color.blue()))

    class _Painter:
        def setPen(self, _pen) -> None:
            pass

        def drawLine(self, *_args, **_kwargs) -> None:
            pass

    monkeypatch.setattr(render_guides_module, "QPen", _FakePen)
    monkeypatch.setattr(
        render_guides_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_guides_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        render_guides_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    for theme in THEME_STYLES_BY_PRESET.values():
        render_guides_module.draw_zenith_marker(
            _Painter(),
            geometry=SimpleNamespace(center=(0, 0), radius=1),
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(30.0, 40.0),
                edge_fov_deg=95.0,
                content_fov_deg=180.0,
            ),
            theme=theme,
        )

    assert seen_colors
    assert seen_colors == [
        tuple(theme.guide_style.grid_rgb or theme.guide_style.horizon_rgb)
        for theme in THEME_STYLES_BY_PRESET.values()
        for _ in range(2)
    ]


def test_draw_celestial_pole_markers_uses_celestial_equator_color_for_all_themes(
    monkeypatch,
) -> None:
    seen_colors: list[tuple[int, int, int]] = []
    seen_positions: list[tuple[float, float]] = []

    class _FakePen:
        def __init__(self, color, _width) -> None:
            seen_colors.append((color.red(), color.green(), color.blue()))

    class _Painter:
        def setPen(self, _pen) -> None:
            pass

        def drawLine(self, *_args, **_kwargs) -> None:
            pass

    monkeypatch.setattr(render_guides_module, "QPen", _FakePen)
    monkeypatch.setattr(
        render_guides_module,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        render_guides_module,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (
            seen_positions.append((float(alt), float(az))) or (float(az), float(alt))
        ),
    )
    monkeypatch.setattr(
        render_guides_module,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    for theme in THEME_STYLES_BY_PRESET.values():
        render_guides_module.draw_celestial_pole_markers(
            _Painter(),
            geometry=SimpleNamespace(center=(0, 0), radius=1),
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(30.0, 40.0),
                edge_fov_deg=95.0,
                content_fov_deg=180.0,
            ),
        )

    assert seen_colors
    assert all(
        color == render_guides_module.CELESTIAL_EQUATOR_COLOR for color in seen_colors
    )
    assert seen_positions[0::2] == [(35.0, 0.0)] * len(THEME_STYLES_BY_PRESET)
    assert seen_positions[1::2] == [(-35.0, 180.0)] * len(THEME_STYLES_BY_PRESET)
