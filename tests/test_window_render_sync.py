from __future__ import annotations

import astropy.time
from types import SimpleNamespace
from PySide6.QtGui import QImage

import zstarview.render.pipeline as pipeline_module
import zstarview.ui.window_render as window_render_module
from zstarview.types import CelestialData, UrbanOutlinePolyline, ViewerData
from zstarview.ui.window import SkyWindow
from zstarview.ui.window_state import SkyWindowState


class _DummyTimer:
    def __init__(self, active: bool) -> None:
        self._active = active
        self.started_with: list[int] = []

    def isActive(self) -> bool:  # noqa: N802 - Qt naming
        return self._active

    def start(self, ms: int) -> None:
        self._active = True
        self.started_with.append(ms)


class _DummySignal:
    def __init__(self) -> None:
        self.calls = 0

    def emit(self) -> None:
        self.calls += 1


class _DummyCompositor:
    def __init__(self) -> None:
        self.invalidated = False

    def invalidate(self) -> None:
        self.invalidated = True


def _make_scene(
    *,
    viewer: ViewerData | None = None,
    celestial_data: CelestialData | object | None = None,
    terrain_horizon_profile: object | None = None,
    urban_outlines: object | None = None,
) -> pipeline_module.RenderSceneData:
    if viewer is None:
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        )
    if celestial_data is None:
        celestial_data = CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[],
            stars={"name": [], "source_id": [], "alt": [], "az": [], "vmag": [], "bv": [], "size_factor": [], "color_factor_base": []},
            deep_sky_objects={},
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
        )
    return pipeline_module.RenderSceneData(
        viewer=viewer,
        celestial_data=celestial_data,
        sky_disc_image=None,
        cloud_image=None,
        cloud_missing_mask=None,
        cloud_stripe_density=None,
        terrain_horizon_profile=terrain_horizon_profile,
        urban_outlines=urban_outlines,
        aircraft_overlay_points=None,
    )


def _make_style(**overrides) -> pipeline_module.RenderStyle:
    values = {
        "visual_preset": "night",
        "text_font": object(),
        "status_line_font": object(),
        "show_dso": False,
        "show_asterisms": False,
        "enlarge_moon": False,
        "star_base_radius": 4.0,
        "star_visibility_boost": 1.0,
        "vmag_limit": 6.0,
        "cloud_disc_alpha": 0.0,
        "terrain_horizon_opacity": 0.25,
        "urban_outline_opacity": 0.2,
        "show_urban_outline_layer": True,
        "aircraft_opacity": 0.0,
        "star_render_expected_width": 600,
    }
    values.update(overrides)
    return pipeline_module.RenderStyle(**values)


def _make_hud(**overrides) -> pipeline_module.RenderHudState:
    values = {
        "mouse_pos": None,
        "viewport_interaction_mode": False,
        "viewport_interaction_stars": None,
        "status_message": None,
    }
    values.update(overrides)
    return pipeline_module.RenderHudState(**values)


def test_viewer_data_for_render_uses_render_view_center() -> None:
    dummy = SimpleNamespace()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=123.0,
    )
    dummy.state = SkyWindowState(render_view_center=(60.0, 210.0))

    got = SkyWindow._viewer_data_for_render(dummy)
    assert got.view_center == (60.0, 210.0)
    assert got.location == (35.0, 139.0)
    assert got.observer_height_m == 123.0
    assert got.location_height_label is None
    assert got.location_height_m is None
    assert got.show_observer_height is False


def test_on_sky_data_calculated_updates_render_snapshot_once() -> None:
    dummy = SimpleNamespace()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 0
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.request_sky_data_update = lambda *_args, **_kwargs: None
    dummy._safe_request_cloud_repaint = lambda: None

    update_calls: list[str] = []
    dummy.update = lambda: update_calls.append("update")

    celestial = object()
    sky_disc = object()
    payload = {
        "celestial": celestial,
        "sky_disc": sky_disc,
        "view_center": (15.0, 120.0),
        "render_width_px": 640,
        "render_height_px": 480,
        "render_generation": 0,
    }
    SkyWindow._on_sky_data_calculated(dummy, payload)

    assert dummy.state.render_view_center == (15.0, 120.0)
    assert dummy.state.celestial_data is celestial
    assert dummy.state.sky_disc_image is sky_disc
    assert dummy._compositor.invalidated is True
    assert update_calls == ["update"]


def test_on_sky_data_calculated_preserves_render_center_during_viewport_interaction() -> None:
    dummy = SimpleNamespace()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(40.0, 150.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(40.0, 150.0),
        viewport_interaction_mode=True,
    )
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 0
    dummy.width = lambda: 640
    dummy.height = lambda: 480
    dummy.request_sky_data_update = lambda *_args, **_kwargs: None
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.update = lambda: None

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (15.0, 120.0),
            "render_width_px": 640,
            "render_height_px": 480,
            "render_generation": 0,
        },
    )

    assert dummy.state.render_view_center == (40.0, 150.0)


def test_on_sky_data_calculated_discards_stale_generation_and_requests_refresh() -> None:
    dummy = SimpleNamespace()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = None
    dummy.cloud_disc_alpha = 0.2
    dummy.sky_update_interval = 60
    dummy.initial_data_loaded = _DummySignal()
    dummy._is_shutting_down = False
    dummy._disc_generation = 3
    dummy.width = lambda: 800
    dummy.height = lambda: 500
    requests: list[str] = []
    dummy.request_sky_data_update = lambda *_args, **_kwargs: requests.append("refresh")
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.update = lambda: (_ for _ in ()).throw(AssertionError("should not repaint stale sky payload"))

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (15.0, 120.0),
            "render_width_px": 640,
            "render_height_px": 480,
            "render_generation": 2,
        },
    )

    assert dummy.state.sky_disc_image is None
    assert dummy._compositor.invalidated is False
    assert requests == ["refresh"]


def test_draw_viewport_interaction_layers_limits_stars_to_bright_subset(monkeypatch) -> None:
    calls: list[tuple[str, object]] = []

    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: calls.append(("reference", None)),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_terrain_horizon_line",
        lambda *_args, **_kwargs: calls.append(("terrain", None)),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: calls.append(("urban", None)),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_direction_labels",
        lambda *_args, **_kwargs: calls.append(("direction", None)),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: calls.append(("zenith", None)),
    )

    monkeypatch.setattr(
        pipeline_module,
        "draw_star_layer",
        lambda *_args, **kwargs: calls.append(("stars", kwargs.get("draw_vmag_limit"))),
    )
    pipeline_module.draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(celestial_data=object()),
        style=_make_style(),
        hud=_make_hud(),
    )

    assert calls == [
        ("reference", None),
        ("stars", 4.0),
        ("terrain", None),
    ]


def test_draw_cached_frame_reuses_existing_image() -> None:
    draws: list[tuple[int, int]] = []
    render_calls: list[str] = []

    class _Painter:
        def drawImage(self, x: int, y: int, image: QImage) -> None:  # noqa: N802 - Qt naming
            draws.append((x, y))
            assert not image.isNull()

    dummy = SimpleNamespace()
    dummy._frame_cache_key = None
    dummy._frame_cache_image = None
    dummy.size = lambda: window_render_module.QImage(32, 24, QImage.Format.Format_ARGB32_Premultiplied).size()

    def render_fn(frame_painter) -> None:
        render_calls.append("render")
        frame_painter.fillRect(0, 0, 32, 24, window_render_module.Qt.GlobalColor.black)

    painter = _Painter()

    window_render_module.SkyWindowRenderMixin._draw_cached_frame(dummy, painter, ("same",), render_fn)
    window_render_module.SkyWindowRenderMixin._draw_cached_frame(dummy, painter, ("same",), render_fn)

    assert render_calls == ["render"]
    assert draws == [(0, 0), (0, 0)]


def test_draw_viewport_interaction_layers_prefers_interaction_star_subset(monkeypatch) -> None:
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_terrain_horizon_line",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )

    full_stars = {"name": ["full"]}
    interaction_stars = {"name": ["bright-only"]}
    seen_stars: list[object] = []
    monkeypatch.setattr(
        pipeline_module,
        "draw_star_layer",
        lambda _p, **kwargs: seen_stars.append(kwargs["scene"].celestial_data.stars),
    )
    pipeline_module.draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(
            celestial_data=CelestialData(
                time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
                planets=[],
                stars=full_stars,
                deep_sky_objects={},
                celestial_equator_points=[],
                ecliptic_points=[],
                horizon_points=[],
            ),
        ),
        style=_make_style(),
        hud=_make_hud(viewport_interaction_stars=interaction_stars),
    )

    assert seen_stars == [interaction_stars]


def test_draw_urban_outline_layer_skips_when_hidden(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: calls.append("urban"),
    )

    pipeline_module.draw_urban_outline_layer(
        painter=object(),
        geometry=object(),
        scene=_make_scene(
            urban_outlines=[UrbanOutlinePolyline(points=[(-1.0, 10.0), (-2.0, 12.0)], height_m=50.0)]
        ),
        style=_make_style(show_urban_outline_layer=False),
    )

    assert calls == []


def test_draw_viewport_interaction_layers_draws_terrain_profile(monkeypatch) -> None:
    seen_profiles: list[object] = []
    seen_view_centers: list[object] = []
    seen_line_width_scales: list[float] = []
    expected_line_width_scale = pipeline_module.compute_star_render_upscale_factor(1200, 600)

    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_terrain_horizon_line",
        lambda _p, _g, profile, view_center, **kwargs: (
            seen_profiles.append(profile),
            seen_view_centers.append(view_center),
            seen_line_width_scales.append(float(kwargs.get("line_width_scale", 1.0))),
        ),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    terrain_profile = [(1.0, 10.0), (2.0, 20.0)]
    monkeypatch.setattr(pipeline_module, "draw_star_layer", lambda *_args, **_kwargs: None)
    pipeline_module.draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=terrain_profile,
        ),
        style=_make_style(),
        hud=_make_hud(),
    )

    assert seen_profiles == [terrain_profile]
    assert seen_view_centers == [(50.0, 210.0)]
    assert seen_line_width_scales == [expected_line_width_scale]


def test_draw_viewport_interaction_layers_skips_urban_outlines(monkeypatch) -> None:
    seen_profiles: list[object] = []
    seen_view_centers: list[object] = []

    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_terrain_horizon_line",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_urban_outlines",
        lambda _p, _g, profile, view_center, **_kwargs: (
            seen_profiles.append(profile),
            seen_view_centers.append(view_center),
        ),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    urban_outlines = [[(-1.0, 10.0), (-2.0, 20.0)]]
    monkeypatch.setattr(pipeline_module, "draw_star_layer", lambda *_args, **_kwargs: None)
    pipeline_module.draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            urban_outlines=urban_outlines,
        ),
        style=_make_style(),
        hud=_make_hud(),
    )

    assert seen_profiles == []
    assert seen_view_centers == []


def test_draw_terrain_layers_only_scales_asterisms_and_terrain(monkeypatch) -> None:
    calls: dict[str, list[float]] = {
        "asterisms": [],
        "terrain": [],
        "reference": [],
        "direction": [],
        "zenith": [],
        "urban": [],
    }
    expected_line_width_scale = pipeline_module.compute_star_render_upscale_factor(1200, 600)

    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_deep_sky_shapes",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_dso_hover_info",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_asterisms",
        lambda *_args, **kwargs: calls["asterisms"].append(float(kwargs.get("line_width_scale", 1.0))),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_terrain_horizon_line",
        lambda *_args, **kwargs: calls["terrain"].append(float(kwargs.get("line_width_scale", 1.0))),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_sky_reference_lines",
        lambda *_args, **kwargs: calls["reference"].append(float(kwargs.get("line_width_scale", 1.0))),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_direction_labels",
        lambda *_args, **kwargs: calls["direction"].append(float(kwargs.get("line_width_scale", 1.0))),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_zenith_marker",
        lambda *_args, **kwargs: calls["zenith"].append(float(kwargs.get("line_width_scale", 1.0))),
    )

    monkeypatch.setattr(
        pipeline_module,
        "draw_urban_outline_layer",
        lambda *_args, **kwargs: calls["urban"].append(float(kwargs.get("line_width_scale", 1.0))),
    )
    pipeline_module.draw_terrain_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=[(1.0, 10.0), (2.0, 20.0)],
        ),
        style=_make_style(show_asterisms=True),
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls["asterisms"] == [expected_line_width_scale]
    assert calls["terrain"] == [expected_line_width_scale]
    assert calls["reference"] == [1.0]
    assert calls["direction"] == []
    assert calls["zenith"] == []
    assert calls["urban"] == [1.0]


def test_draw_terrain_layers_does_not_draw_dso_hover_info(monkeypatch) -> None:
    dso_hover_calls: list[str] = []

    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_deep_sky_shapes",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_dso_hover_info",
        lambda *_args, **_kwargs: dso_hover_calls.append("dso-hover"),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_asterisms",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_terrain_horizon_line",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    monkeypatch.setattr(pipeline_module, "draw_urban_outline_layer", lambda *_args, **_kwargs: None)
    pipeline_module.draw_terrain_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=_make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=[(1.0, 10.0), (2.0, 20.0)],
        ),
        style=_make_style(show_dso=True),
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert dso_hover_calls == []


def test_render_scene_draws_dso_hover_immediately_before_overlay(monkeypatch) -> None:
    calls: list[str] = []

    monkeypatch.setattr(pipeline_module, "clear_background_layer", lambda *_args, **_kwargs: calls.append("clear"))
    monkeypatch.setattr(pipeline_module, "draw_background_layer", lambda *_args, **_kwargs: calls.append("background"))
    monkeypatch.setattr(pipeline_module, "draw_guide_layer", lambda *_args, **_kwargs: calls.append("guide"))
    monkeypatch.setattr(pipeline_module, "draw_sky_cloud_layers", lambda *_args, **_kwargs: calls.append("sky-cloud"))
    monkeypatch.setattr(pipeline_module, "draw_terrain_layers", lambda *_args, **_kwargs: calls.append("terrain"))
    monkeypatch.setattr(pipeline_module, "draw_star_layer", lambda *_args, **_kwargs: calls.append("stars"))
    monkeypatch.setattr(pipeline_module, "draw_aircraft_layer", lambda *_args, **_kwargs: calls.append("aircraft"))
    monkeypatch.setattr(pipeline_module, "draw_planet_layer", lambda *_args, **_kwargs: calls.append("planets"))
    monkeypatch.setattr(pipeline_module, "draw_dso_hover_layer", lambda *_args, **_kwargs: calls.append("dso-hover"))
    monkeypatch.setattr(pipeline_module, "draw_overlay_layer", lambda *_args, **_kwargs: calls.append("overlay"))
    monkeypatch.setattr(pipeline_module, "draw_label_layer", lambda *_args, **_kwargs: calls.append("labels"))
    monkeypatch.setattr(pipeline_module, "draw_status_line", lambda *_args, **_kwargs: calls.append("status"))

    geometry = SimpleNamespace(center=(100, 100), radius=80)
    viewport_rect = SimpleNamespace(width=lambda: 200, height=lambda: 200)
    scene = pipeline_module.RenderSceneData(
        viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        ),
        celestial_data=CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[],
            stars={"name": [], "source_id": [], "alt": [], "az": [], "vmag": [], "bv": [], "size_factor": [], "color_factor_base": []},
            deep_sky_objects={},
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
        ),
        sky_disc_image=None,
        cloud_image=None,
        cloud_missing_mask=None,
        cloud_stripe_density=None,
        terrain_horizon_profile=None,
        urban_outlines=None,
        aircraft_overlay_points=None,
    )
    style = pipeline_module.RenderStyle(
        visual_preset="night",
        text_font=object(),
        status_line_font=object(),
        show_dso=True,
        show_asterisms=False,
        enlarge_moon=False,
        star_base_radius=4.0,
        star_visibility_boost=1.0,
        vmag_limit=6.0,
        cloud_disc_alpha=0.0,
        terrain_horizon_opacity=0.0,
        urban_outline_opacity=0.0,
        show_urban_outline_layer=False,
        aircraft_opacity=0.0,
        star_render_expected_width=600,
    )
    hud = pipeline_module.RenderHudState(
        mouse_pos=None,
        viewport_interaction_mode=False,
        viewport_interaction_stars=None,
        status_message=None,
    )

    pipeline_module.render_scene_into_painter(
        painter=object(),
        geometry=geometry,
        viewport_rect=viewport_rect,
        scene=scene,
        style=style,
        hud=hud,
        compositor=object(),
        highlighted_object=({"name": "Vega"}, object()),
        highlighted_dso=({"name": "M31"}, object()),
    )

    assert calls == [
        "clear",
        "background",
        "sky-cloud",
        "guide",
        "terrain",
        "stars",
        "aircraft",
        "planets",
        "dso-hover",
        "overlay",
        "labels",
        "status",
    ]


def test_draw_sky_reference_lines_uses_render_view_center_projection(monkeypatch) -> None:
    calls: list[tuple[float, float]] = []

    monkeypatch.setattr(
        window_render_module.render_draw,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "altaz_to_normalized_xy",
        lambda alt, az, view_center: calls.append(view_center) or (alt, az),
    )

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

    window_render_module.render_draw.draw_sky_reference_lines(
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
    )

    assert calls == [(55.0, 200.0), (55.0, 200.0), (55.0, 200.0)]


def test_draw_urban_outlines_simplifies_narrow_outline_to_horizontal_segment(monkeypatch) -> None:
    monkeypatch.setattr(
        window_render_module.render_draw,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.polylines: list[list[tuple[float, float]]] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, poly) -> None:
            self.polylines.append([(poly.at(i).x(), poly.at(i).y()) for i in range(poly.count())])

    painter = _Painter()
    window_render_module.render_draw.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[UrbanOutlinePolyline(points=[(-10.0, 10.0), (-12.0, 10.3)], height_m=50.0)],
        view_center=(45.0, 180.0),
        opacity=0.38,
    )

    assert len(painter.polylines) == 1
    assert len(painter.polylines[0]) == 2
    assert painter.polylines[0][0][1] == painter.polylines[0][1][1]


def test_draw_urban_outlines_reduces_alpha_for_short_buildings(monkeypatch) -> None:
    monkeypatch.setattr(
        window_render_module.render_draw,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )

    class _Painter:
        def __init__(self) -> None:
            self.alpha_values: list[int] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.alpha_values.append(int(pen.color().alpha()))

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    window_render_module.render_draw.draw_urban_outlines(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        urban_outlines=[
            UrbanOutlinePolyline(points=[(-10.0, 10.0), (-12.0, 12.0)], height_m=0.0),
            UrbanOutlinePolyline(points=[(-10.0, 20.0), (-12.0, 22.0)], height_m=50.0),
        ],
        view_center=(45.0, 180.0),
        opacity=0.2,
    )

    assert painter.alpha_values == [13, 51]


def test_draw_terrain_horizon_line_scales_line_widths(monkeypatch) -> None:
    monkeypatch.setattr(
        window_render_module.render_draw,
        "is_in_fov",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center: (float(az), float(alt)),
    )
    monkeypatch.setattr(
        window_render_module.render_draw,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    class _Painter:
        def __init__(self) -> None:
            self.pen_widths: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, pen) -> None:
            self.pen_widths.append(float(pen.widthF()))

        def drawPolyline(self, _poly) -> None:
            pass

    painter = _Painter()
    window_render_module.render_draw.draw_terrain_horizon_line(
        painter,
        geometry=SimpleNamespace(center=(0, 0), radius=1),
        terrain_profile_altaz=[(0.0, 0.0), (0.1, 0.1)],
        view_center=(45.0, 180.0),
        opacity=0.38,
        line_width_scale=2.0,
    )

    assert painter.pen_widths[:2] == [6.0, 2.0]
