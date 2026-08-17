
from tests._window_render_support import *


def test_render_cached_frame_image_reuses_existing_image() -> None:
    render_calls: list[str] = []

    dummy = _WindowStub()
    dummy._present_frame_cache_key = None
    dummy._present_frame_cache_image = None
    dummy.size = lambda: window_render_module.QImage(
        24, 16, QImage.Format.Format_ARGB32_Premultiplied
    ).size()

    def render_fn(frame_painter) -> None:
        render_calls.append("render")
        frame_painter.fillRect(0, 0, 24, 16, window_render_module.Qt.GlobalColor.black)

    image_a = window_render_module.SkyWindowRenderMixin._render_cached_frame_image(
        dummy,
        frame_key=("same",),
        render_fn=render_fn,
        cache_key_attr="_present_frame_cache_key",
        cache_image_attr="_present_frame_cache_image",
    )
    image_b = window_render_module.SkyWindowRenderMixin._render_cached_frame_image(
        dummy,
        frame_key=("same",),
        render_fn=render_fn,
        cache_key_attr="_present_frame_cache_key",
        cache_image_attr="_present_frame_cache_image",
    )

    assert render_calls == ["render"]
    assert image_a.cacheKey() == image_b.cacheKey()


def test_render_fast_frame_image_downsamples_base_scene(monkeypatch) -> None:
    base_frame_sizes: list[tuple[int, int]] = []
    call_order: list[str] = []
    fast_overlay_modes: list[bool] = []

    def _capture_base_scene(*_args, **kwargs) -> None:
        call_order.append("base")
        frame = kwargs["frame"]
        base_frame_sizes.append(
            (
                int(frame.viewport_rect.width()),
                int(frame.viewport_rect.height()),
            )
        )

    def _capture_fast_overlays(*_args, **kwargs) -> None:
        call_order.append("fast-overlays")
        fast_overlay_modes.append(bool(kwargs.get("fast_mode", False)))

    monkeypatch.setattr(
        window_render_module,
        "render_base_scene_into_painter",
        _capture_base_scene,
    )
    monkeypatch.setattr(
        window_render_module,
        "render_fast_overlay_layers_into_painter",
        _capture_fast_overlays,
    )
    monkeypatch.setattr(
        window_render_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: call_order.append("labels"),
    )

    dummy = _WindowStub(
        state=SkyWindowState(
            render_view_center=(45.0, 180.0),
            viewport_interaction_mode=True,
        ),
    )
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy._compositor = _DummyCompositor()
    dummy._fast_frame_base_cache_key = None
    dummy._fast_frame_base_cache_image = None
    dummy._fast_frame_cache_key = None
    dummy._fast_frame_cache_image = None
    dummy.client_rect = lambda: QRect(0, 0, 1600, 900)
    dummy.client_size = lambda: QSize(1600, 900)

    scene = _make_scene(viewer=dummy.viewer_data)
    style = _make_style(show_custom_window_frame=True)
    hud = _make_hud(viewport_interaction_mode=True, status_message="fast")
    frame = window_render_module.FrameContext(
        viewer=scene.viewer,
        time_obj=scene.time_obj,
        geometry=render_geometry.get_screen_geometry(
            1600,
            900,
            scene.viewer.view_alt_deg,
        ),
        viewport_rect=QRect(0, 0, 1600, 900),
    )

    image = window_render_module.SkyWindowRenderMixin._render_fast_frame_image(
        dummy,
        base_frame_key=("frame",),
        frame=frame,
        render_inputs=window_render_module.RenderInputs(
            scene=scene,
            style=style,
            hud=hud,
        ),
    )

    assert base_frame_sizes == [(600, 338)]
    assert call_order == ["base", "fast-overlays", "labels"]
    assert fast_overlay_modes == [True]
    assert image.size() == QSize(1600, 900)


def test_render_fast_frame_image_disables_labels(monkeypatch) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        window_render_module,
        "render_base_scene_into_painter",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module,
        "render_fast_overlay_layers_into_painter",
        lambda *_args, **kwargs: captured.update(
            {"draw_labels": kwargs.get("draw_labels")}
        ),
    )
    monkeypatch.setattr(
        window_render_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )

    dummy = _WindowStub(
        state=SkyWindowState(
            render_view_center=(45.0, 180.0),
            viewport_interaction_mode=True,
        ),
    )
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy._compositor = _DummyCompositor()
    dummy._fast_frame_base_cache_key = None
    dummy._fast_frame_base_cache_image = None
    dummy._fast_frame_cache_key = None
    dummy._fast_frame_cache_image = None
    dummy.client_rect = lambda: QRect(0, 0, 1600, 900)
    dummy.client_size = lambda: QSize(1600, 900)

    scene = _make_scene(viewer=dummy.viewer_data)
    style = _make_style(show_custom_window_frame=True)
    hud = _make_hud(viewport_interaction_mode=True, status_message="fast")
    frame = window_render_module.FrameContext(
        viewer=scene.viewer,
        time_obj=scene.time_obj,
        geometry=render_geometry.get_screen_geometry(
            1600,
            900,
            scene.viewer.view_alt_deg,
        ),
        viewport_rect=QRect(0, 0, 1600, 900),
    )

    window_render_module.SkyWindowRenderMixin._render_fast_frame_image(
        dummy,
        base_frame_key=("frame",),
        frame=frame,
        render_inputs=window_render_module.RenderInputs(
            scene=scene,
            style=style,
            hud=hud,
        ),
    )

    assert captured == {"draw_labels": False}


def test_compose_periodic_debug_snapshot_image_includes_volatile_overlay(
    monkeypatch,
) -> None:
    dummy = _WindowStub()
    dummy._draw_volatile_overlay_layers = lambda painter, **_kwargs: painter.fillRect(
        0, 0, 1, 1, Qt.GlobalColor.red
    )

    present_frame = QImage(8, 8, QImage.Format.Format_ARGB32_Premultiplied)
    present_frame.fill(Qt.GlobalColor.black)
    scene = _make_scene()
    style = _make_style()
    hud = _make_hud()
    frame = _make_frame(
        scene,
        SimpleNamespace(center=(100, 100), radius=80),
        QRect(0, 0, 8, 8),
    )

    composed = window_render_module.SkyWindowRenderMixin._compose_periodic_debug_snapshot_image(
        dummy,
        present_frame,
        hover_targets=window_render_module.HoverTargets(
            object=(object(), QPointF(1, 1)),
        ),
        frame=frame,
        render_inputs=window_render_module.RenderInputs(
            scene=scene,
            style=style,
            hud=hud,
        ),
    )

    assert present_frame.pixelColor(0, 0) == QColor(Qt.GlobalColor.black)
    assert composed.pixelColor(0, 0) == QColor(Qt.GlobalColor.red)


def test_draw_background_layer_can_skip_menu_button(monkeypatch) -> None:
    border_menu_flags: list[bool] = []

    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_radial_background",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_background,
        "draw_window_border",
        lambda *_args, **kwargs: border_menu_flags.append(
            bool(kwargs.get("draw_menu_button", True))
        ),
    )

    zstarview_pipeline_module._draw_background_layer(
        painter=object(),
        geometry=SimpleNamespace(),
        viewport_rect=QRect(0, 0, 1600, 900),
        scene=_make_scene(),
        style=_make_style(show_custom_window_frame=True),
        draw_menu_button=False,
    )

    assert border_menu_flags == [False]


def test_viewport_interaction_hides_menu_button() -> None:
    class _MenuButton:
        def __init__(self) -> None:
            self.visible = True

        def setVisible(self, value: bool) -> None:
            self.visible = value

    dummy = _WindowStub()
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.menu_button = _MenuButton()

    SkyWindow._sync_viewport_interaction_chrome_visibility(dummy)
    assert dummy.menu_button.visible is True

    dummy.state.viewport_interaction_mode = True
    SkyWindow._sync_viewport_interaction_chrome_visibility(dummy)
    assert dummy.menu_button.visible is False


def test_paint_event_skips_rendering_while_startup_overlay_visible(
    monkeypatch,
) -> None:
    class _VisibleOverlay:
        def isVisible(self) -> bool:
            return True

    class _FailPainter:
        def __init__(self, *_args, **_kwargs) -> None:
            raise AssertionError(
                "QPainter should not be constructed while startup overlay is visible"
            )

    dummy = _WindowStub(
        _startup_log_overlay=_VisibleOverlay(),
        state=SkyWindowState(render_view_center=(45.0, 180.0)),
    )

    monkeypatch.setattr(window_render_module, "QPainter", _FailPainter)

    window_render_module.SkyWindowRenderMixin.paintEvent(dummy, SimpleNamespace())


def test_render_frame_cache_key_ignores_volatile_overlay_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.aircraft_opacity = 0.4
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._status_line_message = lambda: "initial"
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(), missing_mask=object(), cloud_amount_field=None
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.state.water_overlay_dots = [object()]
    dummy.state.mouse_pos = None
    dummy.state.jump_highlight_name = None
    dummy.state.jump_highlight_altaz = None
    dummy.state.jump_highlight_until_ms = 0.0
    dummy.state.twinkle_targets = ((0, 0.5),)

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    dummy.state.mouse_pos = SimpleNamespace(x=lambda: 10, y=lambda: 20)
    dummy.state.jump_highlight_name = "Vega"
    dummy.state.jump_highlight_altaz = (20.0, 30.0)
    dummy.state.jump_highlight_until_ms = 12345.0
    dummy.state.twinkle_targets = ((1, 0.25),)
    dummy._status_line_message = lambda: "changed"

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    assert key_a == key_b


def test_render_frame_cache_key_tracks_inverted_city_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.aircraft_opacity = 0.0
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy.inverted_city_enabled = False
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(), missing_mask=object(), cloud_amount_field=None
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.urban_outlines = [object()]
    dummy.state.water_overlay_dots = [object()]

    key_normal = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )
    dummy.inverted_city_enabled = True
    key_inverted = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    assert key_normal != key_inverted


def test_render_frame_cache_key_tracks_water_overlay_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.aircraft_opacity = 0.4
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(), missing_mask=object(), cloud_amount_field=None
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.state.water_overlay_dots = [object()]

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    dummy.state.water_overlay_dots = [object(), object()]

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
    )

    assert key_a != key_b


def test_render_frame_cache_key_ignores_projected_tropical_cyclone_state_for_base_cache() -> (
    None
):
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.aircraft_opacity = 0.4
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(), missing_mask=object(), cloud_amount_field=None
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.state.water_overlay_dots = [object()]
    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text=None,
    )
    dummy._current_time_obj = lambda: astropy.time.Time(
        "2026-04-18T12:00:00", scale="utc"
    )

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm",
    )
    dummy._current_time_obj = lambda: astropy.time.Time(
        "2026-04-18T12:00:03", scale="utc"
    )

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    assert key_a == key_b


def test_present_frame_cache_key_ignores_volatile_overlay_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.satellite_opacity = 0.4
    dummy.aircraft_opacity = 0.3
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=SimpleNamespace(source_cache_key=123),
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]

    base_key = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )
    key_a = SkyWindow._present_frame_cache_key(
        dummy,
        base_frame_key=base_key,
        hud=_make_hud(status_message="initial"),
    )

    dummy.state.mouse_pos = QPoint(10, 20)
    dummy.state.jump_highlight_name = "Vega"
    dummy.state.jump_highlight_altaz = (20.0, 30.0)
    dummy.state.jump_highlight_until_ms = 12345.0

    key_b = SkyWindow._present_frame_cache_key(
        dummy,
        base_frame_key=base_key,
        hud=_make_hud(mouse_pos=QPoint(10, 20), status_message="changed"),
    )

    assert key_a == key_b


def test_render_frame_cache_key_ignores_fast_overlay_state_for_base_cache() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.satellite_opacity = 0.4
    dummy.aircraft_opacity = 0.3
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=SimpleNamespace(source_cache_key=123),
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm-a",
    )

    key_a = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    dummy.satellite_opacity = 0.9
    dummy.aircraft_opacity = 0.8
    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm-b",
    )

    key_b = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )

    assert key_a == key_b


def test_present_frame_cache_key_tracks_projected_tropical_cyclone_state() -> None:
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    dummy = _WindowStub()
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.visual_preset = "night"
    dummy.show_dso = True
    dummy.show_asterisms = True
    dummy.show_guidelines = True
    dummy.vmag_limit = 6.0
    dummy.sky_disc_alpha = 1.0
    dummy.cloud_disc_alpha = 0.2
    dummy.satellite_opacity = 0.4
    dummy.aircraft_opacity = 0.3
    dummy.terrain_horizon_opacity = 0.5
    dummy.urban_outline_opacity = 0.2
    dummy.show_urban_outline_layer = True
    dummy._render_cache_stamp = lambda value: (
        window_render_module.SkyWindowRenderMixin._render_cache_stamp(dummy, value)
    )
    dummy.cloud_state = SimpleNamespace(
        image=object(),
        missing_mask=object(),
        cloud_amount_field=SimpleNamespace(source_cache_key=123),
    )
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.sky_disc_image = object()
    dummy.state.terrain_horizon_profile = [(1.0, 2.0)]
    dummy.state.urban_outlines = [object()]
    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm-a",
    )
    dummy._current_time_obj = lambda: astropy.time.Time(
        "2026-04-18T12:00:00", scale="utc"
    )

    base_key = SkyWindow._render_frame_cache_key(
        dummy,
        geometry=geometry,
        celestial_data=celestial_data,
        render_viewer=viewer,
        include_fast_overlays=False,
    )
    key_a = SkyWindow._present_frame_cache_key(
        dummy,
        base_frame_key=base_key,
        hud=_make_hud(status_message="initial"),
    )

    dummy.tropical_cyclone_state = SimpleNamespace(
        snapshots=(object(),),
        snapshot_collection=None,
        banner_text="storm-b",
    )
    key_b = SkyWindow._present_frame_cache_key(
        dummy,
        base_frame_key=base_key,
        hud=_make_hud(status_message="initial"),
    )

    assert key_a != key_b
