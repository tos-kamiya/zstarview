
from tests._window_render_support import *


def test_instrument_presentation_uses_stable_context_layers(monkeypatch) -> None:
    calls: list[str] = []

    monkeypatch.setattr(
        render_instrument_background_module,
        "draw_instrument_background",
        lambda *_args, **_kwargs: calls.append("instrument-background"),
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_background_layer",
        lambda *_args, **_kwargs: calls.append("background"),
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **_kwargs: calls.append("sky-cloud"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_guide_layer",
        lambda *_args, **_kwargs: calls.append("guide"),
    )
    monkeypatch.setattr(
        atlas_pipeline_module,
        "_draw_instrument_guide_layer",
        lambda *_args, **_kwargs: calls.append("guide"),
    )
    monkeypatch.setattr(
        atlas_pipeline_module,
        "_draw_instrument_context_layers",
        lambda *_args, **_kwargs: calls.append("instrument-context"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: calls.append("stars"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: calls.append("planets"),
    )
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
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: calls.append("labels"),
    )

    scene = _make_scene()
    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(radius=600),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(presentation_id="instrument"),
        hud=_make_hud(viewport_interaction_mode=True),
        compositor=object(),
    )

    assert calls == [
        "instrument-background",
        "guide",
        "instrument-context",
        "stars",
        "planets",
        "satellites",
        "aircraft",
        "labels",
    ]


def test_instrument_presentation_does_not_use_shared_background(monkeypatch) -> None:
    called = {"instrument": 0, "background": 0}

    def _bump(key: str) -> None:
        called[key] += 1

    monkeypatch.setattr(
        render_instrument_background_module,
        "draw_instrument_background",
        lambda *_args, **_kwargs: _bump("instrument"),
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_background_layer",
        lambda *_args, **_kwargs: _bump("background"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_guide_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        atlas_pipeline_module,
        "_draw_instrument_guide_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        atlas_pipeline_module,
        "_draw_instrument_context_layers",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_satellite_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_aircraft_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: None,
    )

    scene = _make_scene()
    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(radius=600),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(presentation_id="instrument"),
        hud=_make_hud(),
        compositor=object(),
    )

    assert called == {"instrument": 1, "background": 0}


def test_instrument_simplified_view_hides_ground_tint_and_urban_outline(
    monkeypatch,
) -> None:
    calls: list[object] = []

    monkeypatch.setattr(
        render_instrument_background_module,
        "draw_instrument_background",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        render_terrain_module,
        "draw_ground_tint",
        lambda *_args, **_kwargs: calls.append("ground-tint"),
    )
    monkeypatch.setattr(
        atlas_pipeline_module,
        "_draw_instrument_guide_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        atlas_pipeline_module,
        "_draw_instrument_context_layers",
        lambda *_args, **kwargs: calls.append(
            ("context", kwargs["simplified_view_active"])
        ),
    )
    monkeypatch.setattr(
        atlas_pipeline_module,
        "_draw_instrument_cloud_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        render_instrument_background_module,
        "draw_instrument_time_of_day_marker",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_star_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_planet_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_satellite_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_aircraft_layer", lambda *_args, **_kwargs: None
    )

    scene = _make_scene()
    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(radius=600),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(presentation_id="instrument"),
        hud=_make_hud(
            simplified_view_enabled=True,
            simplified_view_labels_enabled=False,
        ),
        compositor=object(),
    )

    assert calls == [("context", True)]


def test_viewer_data_for_render_uses_render_view_center() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=123.0,
        height_add_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(60.0, 210.0))

    got = SkyWindow._viewer_data_for_render(dummy)
    assert got.view_center == (60.0, 210.0)
    assert got.location == (35.0, 139.0)
    assert got.observer_height_m == 123.0
    assert got.height_add_m == 1.7
    assert got.location_height_label is None
    assert got.location_height_m == 0.0


def test_render_style_uses_window_observation_info_toggle() -> None:
    dummy = _WindowStub()
    dummy.visual_preset = "night"
    dummy.text_font = object()
    dummy.status_line_font = object()
    dummy._frameless_window = False
    dummy.show_observation_info = False
    dummy.show_dso = True
    dummy.show_asterisms = False
    dummy.show_guidelines = True
    dummy.star_base_radius = 4.0
    dummy.star_visibility_boost = 1.0
    dummy.vmag_limit = 6.0
    dummy.cloud_disc_alpha = 0.0
    dummy.satellite_opacity = 0.0
    dummy.terrain_horizon_opacity = 0.25
    dummy.urban_outline_opacity = 0.2
    dummy.ridge_glow_opacity = 0.04
    dummy.show_urban_outline_layer = True
    dummy.aircraft_opacity = 0.0
    dummy._star_render_expected_width = 600

    style = SkyWindow._render_style(dummy)

    assert style.show_custom_window_frame is False
    assert style.show_observation_info is False


def test_render_hud_state_uses_upper_third_to_switch_overlay_to_bottom_left() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = QPoint(10, 20)
    dummy.state.overlay_info_bottom_left = False
    dummy.height = lambda: 300
    dummy._status_line_message = lambda: None

    hud = SkyWindow._render_hud_state(dummy)

    assert hud.overlay_info_bottom_left is True
    assert dummy.state.overlay_info_bottom_left is True


def test_render_hud_state_keeps_overlay_anchor_in_middle_third() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = QPoint(10, 150)
    dummy.state.overlay_info_bottom_left = True
    dummy.height = lambda: 300
    dummy._status_line_message = lambda: None

    hud = SkyWindow._render_hud_state(dummy)

    assert hud.overlay_info_bottom_left is True
    assert dummy.state.overlay_info_bottom_left is True


def test_render_hud_state_uses_lower_third_to_switch_overlay_to_top_left() -> None:
    dummy = _WindowStub()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = QPoint(10, 280)
    dummy.state.overlay_info_bottom_left = True
    dummy.height = lambda: 300
    dummy._status_line_message = lambda: None

    hud = SkyWindow._render_hud_state(dummy)

    assert hud.overlay_info_bottom_left is False
    assert dummy.state.overlay_info_bottom_left is False


def test_render_hud_state_preserves_pinned_overlay_anchor() -> None:
    dummy = _WindowStub(observation_info_pinned=True)
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = QPoint(10, 20)
    dummy.state.overlay_info_bottom_left = False
    dummy.height = lambda: 300
    dummy._status_line_message = lambda: None

    hud = SkyWindow._render_hud_state(dummy)

    assert hud.overlay_info_bottom_left is False
    assert dummy.state.overlay_info_bottom_left is False


def test_status_line_text_always_uses_night_style(monkeypatch) -> None:
    class DummyFontMetrics:
        def lineSpacing(self) -> int:
            return 12

    class DummyPainter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setFont(self, _font) -> None:
            pass

        def fontMetrics(self):
            return DummyFontMetrics()

    calls: list[tuple[str, bool]] = []

    def fake_get_text_style(theme, *, status_line: bool = False):
        calls.append((theme, status_line))
        return (object(), object())

    monkeypatch.setattr(render_text_module, "get_text_style", fake_get_text_style)
    monkeypatch.setattr(
        render_text_module, "draw_outlined_text", lambda *_args, **_kwargs: None
    )

    render_text_module._draw_status_line_text(
        painter=DummyPainter(),
        message="loading",
        status_line_font=QFont(),
        viewport_rect=SimpleNamespace(bottom=lambda: 100),
        theme=THEME_STYLES_BY_PRESET["day"],
    )

    assert calls == [(THEME_STYLES_BY_PRESET["day"], True)]


def test_wrap_text_lines_fits_each_hud_line_to_available_width() -> None:
    font = QFont()
    metrics = render_text_module.QFontMetrics(font)
    max_width = metrics.horizontalAdvance("Skytree")

    lines = render_text_module.wrap_text_lines("Tokyo Skytree", font, max_width)

    assert lines == ["Tokyo", "Skytree"]
    assert all(metrics.horizontalAdvance(line) <= max_width for line in lines)


def test_status_line_text_draws_multiple_lines_from_bottom(monkeypatch) -> None:
    class DummyFontMetrics:
        def lineSpacing(self) -> int:
            return 12

    class DummyPainter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setFont(self, _font) -> None:
            pass

        def fontMetrics(self):
            return DummyFontMetrics()

    drawn: list[tuple[str, float]] = []
    monkeypatch.setattr(
        render_text_module,
        "draw_outlined_text",
        lambda _painter, text, pos, _font, **_kwargs: drawn.append(
            (text, float(pos.y()))
        ),
    )

    render_text_module._draw_status_line_text(
        painter=DummyPainter(),
        message="fixed\ndynamic",
        status_line_font=QFont(),
        viewport_rect=SimpleNamespace(bottom=lambda: 100),
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert drawn == [("dynamic", 97.0), ("fixed", 85.0)]
