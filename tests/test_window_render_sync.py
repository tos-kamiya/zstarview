from __future__ import annotations

import astropy.time
from types import SimpleNamespace

import zstarview.ui.window_render as window_render_module
from zstarview.types import CelestialData, ViewerData
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
    dummy.request_sky_data_update = lambda *_args, **_kwargs: None
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.update = lambda: None

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (15.0, 120.0),
        },
    )

    assert dummy.state.render_view_center == (40.0, 150.0)


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

    dummy = SimpleNamespace()
    dummy.text_font = object()
    dummy.visual_preset = "night"
    dummy.terrain_horizon_opacity = 0.25
    dummy.show_urban_outline_layer = True
    dummy._star_render_expected_width = 600
    dummy.compute_star_render_upscale_factor = lambda *_args, **_kwargs: 2.0
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = None
    dummy._draw_star_layer = lambda *_args, **kwargs: calls.append(("stars", kwargs.get("draw_vmag_limit")))

    SkyWindow._draw_viewport_interaction_layers(
        dummy,
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        celestial_data=object(),
        render_viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        ),
    )

    assert calls == [
        ("reference", None),
        ("stars", 4.0),
        ("terrain", None),
        ("direction", None),
        ("zenith", None),
    ]


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
    dummy = SimpleNamespace()
    dummy.text_font = object()
    dummy.visual_preset = "night"
    dummy.terrain_horizon_opacity = 0.25
    dummy.show_urban_outline_layer = True
    dummy._star_render_expected_width = 600
    dummy.compute_star_render_upscale_factor = lambda *_args, **_kwargs: 2.0
    dummy.state = SkyWindowState(
        render_view_center=(45.0, 180.0),
        viewport_interaction_stars=interaction_stars,
    )
    dummy.state.mouse_pos = None
    seen_stars: list[object] = []
    dummy._draw_star_layer = lambda _p, _g, celestial_data, _rv, **_kwargs: seen_stars.append(celestial_data.stars)
    dummy._draw_urban_outline_layer = lambda *_args, **_kwargs: None

    SkyWindow._draw_viewport_interaction_layers(
        dummy,
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        celestial_data=CelestialData(
            time=astropy.time.Time("2026-03-09T00:00:00", scale="utc"),
            planets=[],
            stars=full_stars,
            deep_sky_objects={},
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
        ),
        render_viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        ),
    )

    assert seen_stars == [interaction_stars]


def test_draw_urban_outline_layer_skips_when_hidden(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        window_render_module.render_draw,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: calls.append("urban"),
    )

    dummy = SimpleNamespace()
    dummy.show_urban_outline_layer = False
    dummy.state = SkyWindowState(
        render_view_center=(45.0, 180.0),
        urban_outlines=[[(-1.0, 10.0), (-2.0, 12.0)]],
    )

    SkyWindow._draw_urban_outline_layer(
        dummy,
        painter=object(),
        geometry=object(),
        render_viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            observer_height_m=1.7,
        ),
    )

    assert calls == []


def test_draw_viewport_interaction_layers_draws_terrain_profile(monkeypatch) -> None:
    seen_profiles: list[object] = []
    seen_view_centers: list[object] = []
    seen_line_width_scales: list[float] = []

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
    dummy = SimpleNamespace()
    dummy.text_font = object()
    dummy.visual_preset = "night"
    dummy.terrain_horizon_opacity = 0.25
    dummy.show_urban_outline_layer = True
    dummy._star_render_expected_width = 600
    dummy.compute_star_render_upscale_factor = lambda *_args, **_kwargs: 2.0
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = None
    dummy.state.terrain_horizon_profile = terrain_profile
    dummy._draw_star_layer = lambda *_args, **_kwargs: None
    dummy._draw_urban_outline_layer = lambda *_args, **_kwargs: None

    SkyWindow._draw_viewport_interaction_layers(
        dummy,
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        celestial_data=object(),
        render_viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(50.0, 210.0),
            observer_height_m=1.7,
        ),
    )

    assert seen_profiles == [terrain_profile]
    assert seen_view_centers == [(50.0, 210.0)]
    assert seen_line_width_scales == [2.0]


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
    dummy = SimpleNamespace()
    dummy.text_font = object()
    dummy.visual_preset = "night"
    dummy.terrain_horizon_opacity = 0.25
    dummy.show_urban_outline_layer = True
    dummy._star_render_expected_width = 600
    dummy.compute_star_render_upscale_factor = lambda *_args, **_kwargs: 2.0
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = None
    dummy.state.urban_outlines = urban_outlines
    dummy._draw_star_layer = lambda *_args, **_kwargs: None
    dummy._draw_urban_outline_layer = (
        lambda painter, geometry, render_viewer: SkyWindow._draw_urban_outline_layer(
            dummy,
            painter=painter,
            geometry=geometry,
            render_viewer=render_viewer,
        )
    )

    SkyWindow._draw_viewport_interaction_layers(
        dummy,
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        celestial_data=object(),
        render_viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(50.0, 210.0),
            observer_height_m=1.7,
        ),
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

    dummy = SimpleNamespace()
    dummy.show_dso = False
    dummy.show_asterisms = True
    dummy.text_font = object()
    dummy.visual_preset = "night"
    dummy.terrain_horizon_opacity = 0.25
    dummy._star_render_expected_width = 600
    dummy.compute_star_render_upscale_factor = lambda *_args, **_kwargs: 2.0
    dummy.state = SkyWindowState(render_view_center=(45.0, 180.0))
    dummy.state.mouse_pos = None
    dummy.state.terrain_horizon_profile = [(1.0, 10.0), (2.0, 20.0)]
    dummy._draw_urban_outline_layer = lambda *_args, **kwargs: calls["urban"].append(
        float(kwargs.get("line_width_scale", 1.0))
    )

    SkyWindow._draw_terrain_layers(
        dummy,
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        celestial_data=object(),
        render_viewer=ViewerData(
            location=(35.0, 139.0),
            timezone_name="Asia/Tokyo",
            city_name="Tokyo",
            view_center=(50.0, 210.0),
            observer_height_m=1.7,
        ),
        highlighted_object=None,
        highlighted_dso=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls["asterisms"] == [2.0]
    assert calls["terrain"] == [2.0]
    assert calls["reference"] == [1.0]
    assert calls["direction"] == [1.0]
    assert calls["zenith"] == [1.0]
    assert calls["urban"] == [1.0]


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
        urban_outlines=[[(-10.0, 10.0), (-12.0, 10.3)]],
        view_center=(45.0, 180.0),
        opacity=0.38,
    )

    assert len(painter.polylines) == 1
    assert len(painter.polylines[0]) == 2
    assert painter.polylines[0][0][1] == painter.polylines[0][1][1]


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
