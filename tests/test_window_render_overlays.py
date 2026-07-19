# ruff: noqa: F403, F405

from tests._window_render_support import *


def test_jump_to_jpl_major_body_target_keeps_overlay_without_refresh() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=None,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "resolve_jpl_target_state_vector",
            lambda _target, **_kwargs: (
                datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
                (1.0, 2.0, 3.0),
                (0.1, 0.2, 0.3),
            ),
        )
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda _target, **_kwargs: (15.0, 123.0),
        )
        SkyWindow._jump_to_search_target(
            dummy,
            SearchJumpTarget(
                label="Mars",
                kind="jpl_body",
                sort_key=(0.0, "mars"),
                subtitle="major body",
                object_key="499",
                command="499",
                jpl_group="mb",
                persistent_keep_marker=True,
            ),
        )

    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.label == "Mars"
    assert dummy.state.persistent_search_next_refresh_utc == datetime(
        2026, 4, 18, 13, 0, tzinfo=timezone.utc
    )


def test_reproject_satellite_overlay_falls_back_to_disk_cache() -> None:
    dummy = _WindowStub()
    dummy.satellite_opacity = 1.0
    dummy.viewer_data = ViewerData(
        location=(40.7128, -74.0060),
        timezone_name="America/New_York",
        city_name="New York City",
        view_center=(0.0, 151.0),
        observer_height_m=10.0,
    )
    dummy.satellite_state = SimpleNamespace(
        records_by_group={"iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]},
        element_epoch_utc=datetime.now(timezone.utc),
    )
    dummy.state = SkyWindowState(
        render_view_center=(0.0, 151.0),
        satellite_projection_next_refresh_utc=None,
    )
    dummy._enabled_satellite_groups = ("iss",)
    dummy._satellite_validity_remaining_ms = lambda: 1000
    dummy._current_time_obj = lambda: astropy.time.Time("2026-03-23T12:13:24Z")
    dummy.update = Mock()

    SkyWindow.reproject_satellite_overlay(dummy)

    assert dummy.state.satellite_projection_next_refresh_utc is not None
    dummy.update.assert_called_once()


def test_on_sky_data_calculated_discards_stale_generation_and_requests_refresh() -> (
    None
):
    dummy = _WindowStub()
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
    dummy.update = lambda: (_ for _ in ()).throw(
        AssertionError("should not repaint stale sky payload")
    )

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (15.0, 120.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 2,
        },
    )

    assert dummy.state.sky_disc_image is None
    assert dummy._compositor.invalidated is False
    assert requests == ["refresh"]


def test_draw_viewport_interaction_layers_limits_stars_to_bright_subset(
    monkeypatch,
) -> None:
    calls: list[tuple[str, object]] = []

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: calls.append(("reference", None)),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: calls.append(("terrain", None)),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: calls.append(("urban", None)),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: calls.append(("direction", None)),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: calls.append(("zenith", None)),
    )

    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **kwargs: calls.append(("stars", kwargs.get("draw_vmag_limit"))),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **kwargs: calls.append(("planets", kwargs.get("draw_labels"))),
    )
    zstarview_pipeline_module._draw_viewport_interaction_layers(
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
        ("planets", False),
        ("terrain", None),
    ]


def test_resolve_hover_targets_keeps_star_and_satellite_candidates_independent(
    monkeypatch,
) -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=1.7,
    )
    geometry = SimpleNamespace(center=(100, 100), radius=80)
    celestial_data = SimpleNamespace(time=None)
    mouse_pos = QPoint(12, 34)

    star_hit = ({"name": "Vega"}, QPointF(10.0, 30.0))
    satellite_hit = (
        SimpleNamespace(satellite_name="ISS"),
        QPointF(14.0, 36.0),
    )

    monkeypatch.setattr(
        window_render_module.render_stars,
        "find_highlighted_object",
        lambda *_args, **_kwargs: star_hit,
    )
    monkeypatch.setattr(
        window_render_module.render_deep_sky_objects,
        "find_highlighted_dso",
        lambda *_args, **_kwargs: {"name": "DSO"},
    )
    monkeypatch.setattr(
        window_render_module.render_satellites,
        "find_highlighted_satellite",
        lambda *_args, **_kwargs: satellite_hit,
    )

    hover_targets = window_render_module._resolve_hover_targets(
        celestial_data=celestial_data,
        render_viewer=viewer,
        mouse_pos=mouse_pos,
        geometry=geometry,
        satellite_records_by_group={"iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]},
        show_dso=True,
    )

    assert hover_targets.object == star_hit
    assert hover_targets.dso == {"name": "DSO"}
    assert hover_targets.satellite == satellite_hit
    assert hover_targets.tropical_cyclone is None


def test_draw_viewport_interaction_layers_prefers_interaction_star_subset(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
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
        "_draw_star_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_fast",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_normal",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )

    full_stars = {"name": ["full"]}
    interaction_stars = {"name": ["bright-only"]}
    seen_stars: list[object] = []
    monkeypatch.setattr(
        pipeline_module,
        "_draw_star_layer",
        lambda _p, **kwargs: seen_stars.append(kwargs["scene"].celestial_data.stars),
    )
    zstarview_pipeline_module._draw_viewport_interaction_layers(
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


def test_draw_viewport_interaction_layers_skips_water_when_terrain_horizon_hidden(
    monkeypatch,
) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: calls.append("water"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_fast",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_normal",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: None,
    )

    scene = replace(
        _make_scene(terrain_horizon_profile=[(1.0, 10.0)]),
        water_overlay_dots=[object()],
    )
    zstarview_pipeline_module._draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=scene,
        style=_make_style(terrain_horizon_opacity=0.0, water_overlay_opacity=0.5),
        hud=_make_hud(),
    )

    assert calls == ["terrain"]


def test_draw_viewport_interaction_layers_prefers_scene_water_overlay_points(
    monkeypatch,
) -> None:
    terrain_calls: list[str] = []
    seen_water_points: list[object] = []
    sentinel_water_points = [object(), object()]

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: terrain_calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda _p, _g, _viewer, water_points, *_args, **_kwargs: (
            seen_water_points.append(water_points)
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_fast",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_normal",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_planet_layer",
        lambda *_args, **_kwargs: None,
    )

    scene = replace(
        _make_scene(terrain_horizon_profile=[(1.0, 10.0)]),
        water_overlay_dots=sentinel_water_points,
    )
    zstarview_pipeline_module._draw_viewport_interaction_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        viewport_rect=SimpleNamespace(width=lambda: 200, height=lambda: 200),
        scene=scene,
        style=_make_style(terrain_horizon_opacity=0.25, water_overlay_opacity=0.5),
        hud=_make_hud(),
    )

    assert terrain_calls == ["terrain"]
    assert seen_water_points == [sentinel_water_points]


def test_render_base_scene_skips_water_when_terrain_horizon_hidden(monkeypatch) -> None:
    calls: list[str] = []
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
        zstarview_pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_guide_layer", lambda *_args, **_kwargs: None
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
    monkeypatch.setattr(
        pipeline_module, "_draw_urban_outline_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_star_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_hover_overlay_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_label_candidates",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_status_line_text",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_fast",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_stars,
        "draw_stars_normal",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: calls.append("terrain"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: calls.append("water"),
    )

    scene = replace(
        _make_scene(terrain_horizon_profile=[(1.0, 10.0)]),
        water_overlay_dots=[object()],
    )

    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setCompositionMode(self, *_args, **_kwargs) -> None:
            pass

        def fillRect(self, *_args, **_kwargs) -> None:
            pass

    pipeline_module.render_base_scene_into_painter(
        painter=_Painter(),
        frame=_make_frame(
            scene,
            SimpleNamespace(radius=600),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(terrain_horizon_opacity=0.0, water_overlay_opacity=0.5),
        compositor=object(),
        hud=_make_hud(),
    )

    assert calls == []


def test_draw_urban_outline_layer_skips_when_hidden(monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda *_args, **_kwargs: calls.append("urban"),
    )

    pipeline_module._draw_urban_outline_layer(
        painter=object(),
        geometry=object(),
        scene=_make_scene(
            urban_outlines=[
                UrbanOutlinePolyline(points=[(-1.0, 10.0), (-2.0, 12.0)], height_m=50.0)
            ]
        ),
        style=_make_style(show_urban_outline_layer=False),
    )

    assert calls == []


def test_draw_viewport_interaction_layers_draws_terrain_profile(monkeypatch) -> None:
    seen_main_profiles: list[object] = []
    seen_view_centers: list[object] = []
    seen_line_width_scales: list[float] = []
    fast_calls: list[object] = []
    expected_line_width_scale = pipeline_module.compute_star_render_upscale_factor(
        1200, 600
    )

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda _p, _g, viewer, main_profile, _main_distances, **kwargs: (
            seen_main_profiles.append(main_profile),
            seen_view_centers.append(viewer.view_center),
            seen_line_width_scales.append(float(kwargs["spec"].line_width_scale)),
            fast_calls.append(True),
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    terrain_profile = [(1.0, 10.0), (2.0, 20.0)]
    monkeypatch.setattr(
        pipeline_module, "_draw_star_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_planet_layer", lambda *_args, **_kwargs: None
    )
    zstarview_pipeline_module._draw_viewport_interaction_layers(
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

    assert seen_main_profiles == [terrain_profile]
    assert seen_view_centers == [(50.0, 210.0)]
    assert seen_line_width_scales == [expected_line_width_scale]
    assert fast_calls == [True]


def test_draw_viewport_interaction_layers_skips_urban_outlines(monkeypatch) -> None:
    seen_profiles: list[object] = []
    seen_view_centers: list[object] = []

    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_urban_outlines",
        lambda _p, _g, profile, view_center, **_kwargs: (
            seen_profiles.append(profile),
            seen_view_centers.append(view_center),
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    urban_outlines = [[(-1.0, 10.0), (-2.0, 20.0)]]
    monkeypatch.setattr(
        pipeline_module, "_draw_star_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_planet_layer", lambda *_args, **_kwargs: None
    )
    zstarview_pipeline_module._draw_viewport_interaction_layers(
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


def test_draw_terrain_layers_scales_asterisms_but_keeps_urban_outline_widths_fixed(
    monkeypatch,
) -> None:
    calls: dict[str, list[float]] = {
        "asterisms": [],
        "dso": [],
        "terrain": [],
        "terrain_secondary": [],
        "reference": [],
        "direction": [],
        "zenith": [],
        "urban": [],
    }
    expected_line_width_scale = pipeline_module.compute_star_render_upscale_factor(
        1200, 600
    )

    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_deep_sky_shapes",
        lambda *_args, **kwargs: calls["dso"].append(
            float(kwargs.get("opacity_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_dso_hover_info",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **kwargs: calls["asterisms"].append(
            float(
                kwargs.get("base_line_width_scale", kwargs.get("line_width_scale", 1.0))
            )
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **kwargs: calls["terrain"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **kwargs: calls["terrain_secondary"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **kwargs: calls["reference"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **kwargs: calls["direction"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **kwargs: calls["zenith"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )

    monkeypatch.setattr(
        pipeline_module,
        "_draw_urban_outline_layer",
        lambda *_args, **kwargs: calls["urban"].append(
            float(kwargs.get("line_width_scale", 1.0))
        ),
    )
    zstarview_pipeline_module._draw_terrain_layers(
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
            terrain_secondary_ridges_altaz_layers=[[(1.0, 10.0), (2.0, 20.0)]],
            terrain_secondary_ridges_distances_m_layers=[[10_000.0, 12_000.0]],
        ),
        style=_make_style(
            show_dso=True, show_asterisms=True, asterism_visibility_boost=2.0
        ),
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls["dso"] == [1.0]
    assert calls["asterisms"] == [expected_line_width_scale * 2.0]
    assert calls["terrain"] == []
    assert calls["terrain_secondary"] == [expected_line_width_scale]
    assert calls["reference"] == [1.0]
    assert calls["direction"] == []
    assert calls["zenith"] == []
    assert calls["urban"] == [1.0]


def test_draw_terrain_layers_dims_dso_and_asterisms_in_simplified_view(
    monkeypatch,
) -> None:
    calls: dict[str, list[float]] = {
        "asterisms": [],
        "dso": [],
    }

    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_deep_sky_shapes",
        lambda *_args, **kwargs: calls["dso"].append(
            float(kwargs.get("opacity_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **kwargs: calls["asterisms"].append(
            float(kwargs.get("base_line_alpha_scale", 1.0))
        ),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_urban_outline_layer",
        lambda *_args, **_kwargs: None,
    )

    zstarview_pipeline_module._draw_terrain_layers(
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
            terrain_secondary_ridges_altaz_layers=[[(1.0, 10.0), (2.0, 20.0)]],
            terrain_secondary_ridges_distances_m_layers=[[10_000.0, 12_000.0]],
        ),
        style=_make_style(
            show_dso=True, show_asterisms=True, asterism_visibility_boost=2.0
        ),
        simplified_view_active=True,
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls["dso"] == [0.4]
    assert calls["asterisms"] == [0.8]


def test_draw_terrain_layers_does_not_draw_dso_hover_info(monkeypatch) -> None:
    dso_hover_calls: list[str] = []

    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_deep_sky_shapes",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_dso_hover_info",
        lambda *_args, **_kwargs: dso_hover_calls.append("dso-hover"),
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_direction_labels",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_zenith_marker",
        lambda *_args, **_kwargs: None,
    )

    monkeypatch.setattr(
        pipeline_module, "_draw_urban_outline_layer", lambda *_args, **_kwargs: None
    )
    zstarview_pipeline_module._draw_terrain_layers(
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


def test_draw_terrain_layers_skips_secondary_layers_while_simplified_view_active(
    monkeypatch,
) -> None:
    calls: list[str] = []

    monkeypatch.setattr(
        pipeline_module.render_deep_sky_objects,
        "draw_deep_sky_shapes",
        lambda *_args, **_kwargs: calls.append("dso"),
    )
    monkeypatch.setattr(
        pipeline_module.render_asterisms,
        "draw_asterisms",
        lambda *_args, **_kwargs: calls.append("asterisms"),
    )
    monkeypatch.setattr(
        pipeline_module.render_guides,
        "draw_sky_reference_lines",
        lambda *_args, **_kwargs: calls.append("guides"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "_draw_terrain_profile_layer",
        lambda *_args, **_kwargs: calls.append("main"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_terrain_secondary_ridges",
        lambda *_args, **_kwargs: calls.append("secondary"),
    )
    monkeypatch.setattr(
        pipeline_module.render_terrain,
        "draw_water_overlay_dots",
        lambda *_args, **_kwargs: calls.append("water"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_urban_outline_layer",
        lambda *_args, **_kwargs: calls.append("urban"),
    )

    scene = replace(
        _make_scene(
            viewer=ViewerData(
                location=(35.0, 139.0),
                timezone_name="Asia/Tokyo",
                city_name="Tokyo",
                view_center=(50.0, 210.0),
                observer_height_m=1.7,
            ),
            celestial_data=object(),
            terrain_horizon_profile=[(1.0, 10.0), (2.0, 20.0)],
            terrain_secondary_ridges_altaz_layers=[[(1.0, 10.0), (2.0, 20.0)]],
            terrain_secondary_ridges_distances_m_layers=[[10_000.0, 12_000.0]],
        ),
        water_overlay_dots=[object()],
    )

    zstarview_pipeline_module._draw_terrain_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=600),
        scene=scene,
        style=_make_style(
            show_dso=True,
            show_asterisms=True,
            show_guidelines=True,
            terrain_horizon_opacity=0.25,
            water_overlay_opacity=0.5,
            show_urban_outline_layer=True,
        ),
        simplified_view_active=True,
        highlighted_object=None,
        label_reservations=[],
        label_candidates=[],
    )

    assert calls == ["dso", "asterisms", "guides", "main"]


def test_render_scene_draws_dso_hover_immediately_before_overlay(monkeypatch) -> None:
    calls: list[str] = []

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
        pipeline_module,
        "_draw_static_observation_overlay",
        lambda *_args, **_kwargs: calls.append("overlay"),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_hover_overlay_layer",
        lambda *_args, **_kwargs: calls.append("hover"),
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
    monkeypatch.setattr(
        pipeline_module.render_text,
        "_draw_status_line_text",
        lambda *_args, **_kwargs: calls.append("status"),
    )

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
            stars={
                "name": [],
                "source_id": [],
                "alt": [],
                "az": [],
                "vmag": [],
                "bv": [],
                "size_factor": [],
                "color_factor_base": [],
            },
            deep_sky_objects={},
            celestial_equator_points=[],
            ecliptic_points=[],
            horizon_points=[],
        ),
        sky_disc_image=None,
        cloud_missing_mask=None,
        cloud_altaz_grid=None,
        terrain_horizon_profile=None,
        terrain_horizon_profile_distances_m=None,
        terrain_secondary_ridges_altaz_layers=None,
        terrain_secondary_ridges_distances_m_layers=None,
        urban_outlines=None,
        satellite_records_by_group=None,
        aircraft_snapshots=None,
    )
    style = pipeline_module.RenderStyle(
        theme=THEME_STYLES_BY_PRESET["night"],
        visual_preset="night",
        text_font=object(),
        status_line_font=object(),
        show_background_gradient=True,
        show_custom_window_frame=False,
        show_observation_info=True,
        show_dso=True,
        show_asterisms=False,
        show_guidelines=True,
        enlarge_moon=False,
        bright_bodies_mode="outline",
        star_base_radius=4.0,
        star_visibility_boost=1.0,
        asterism_visibility_boost=1.0,
        earth_guide_visibility_boost=1.0,
        vmag_limit=6.0,
        sky_disc_altaz_rings="off",
        sky_disc_altaz_rings_hover="altaz",
        cloud_disc_alpha=0.0,
        satellite_opacity=0.0,
        terrain_horizon_opacity=0.0,
        earth_guide_opacity=0.0,
        urban_outline_opacity=0.0,
        show_urban_outline_layer=False,
        aircraft_opacity=0.0,
        tropical_cyclone_opacity=0.4,
        show_tropical_cyclone_overlay=True,
        star_render_expected_width=600,
    )
    hud = pipeline_module.RenderHudState(
        mouse_pos=None,
        overlay_info_bottom_left=False,
        viewport_interaction_mode=False,
        viewport_interaction_stars=None,
        status_message=None,
    )

    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(scene, geometry, viewport_rect),
        scene=scene,
        style=style,
        hud=hud,
        compositor=object(),
    )
    pipeline_module.render_hud_overlay_into_painter(
        painter=object(),
        frame=_make_frame(scene, geometry, viewport_rect),
        scene=scene,
        style=style,
        hud=hud,
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
        "planets",
        "satellites",
        "aircraft",
        "labels",
        "hover",
        "overlay",
    ]


def test_render_scene_reduces_layers_during_simplified_view(monkeypatch) -> None:
    captured: dict[str, object] = {}

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
        zstarview_pipeline_module,
        "_draw_sky_cloud_layers",
        lambda *_args, **kwargs: captured.update(
            {
                "cloud_disc_alpha": kwargs["style"].cloud_disc_alpha,
                "earth_guide_opacity": kwargs["style"].earth_guide_opacity,
                "sky_disc_image": kwargs["scene"].sky_disc_image,
            }
        ),
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_guide_layer", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        zstarview_pipeline_module,
        "_draw_terrain_layers",
        lambda *_args, **kwargs: captured.update(
            {"simplified_view_active": kwargs["simplified_view_active"]}
        ),
    )
    monkeypatch.setattr(
        pipeline_module,
        "_draw_main_terrain_profile_layer",
        lambda *_args, **_kwargs: captured.update({"main_terrain_profile": True}),
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

    scene = replace(
        _make_scene(),
        sky_disc_image=object(),
        night_light_glow_profile=object(),
    )
    pipeline_module.render_base_scene_into_painter(
        painter=object(),
        frame=_make_frame(
            scene,
            SimpleNamespace(center=(100, 100), radius=80),
            SimpleNamespace(width=lambda: 200, height=lambda: 200),
        ),
        scene=scene,
        style=_make_style(cloud_disc_alpha=0.2, earth_guide_opacity=0.25),
        hud=_make_hud(simplified_view_enabled=True),
        compositor=object(),
    )

    assert captured == {
        "cloud_disc_alpha": 0.0,
        "earth_guide_opacity": 0.0,
        "sky_disc_image": scene.sky_disc_image,
        "simplified_view_active": True,
    }


def test_draw_sky_cloud_layers_skips_night_lights_while_simplified_view_active(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    class _Compositor:
        def draw(self, *_args, **kwargs) -> None:
            captured.update(
                {
                    "night_light_glow_profile": kwargs["night_light_glow_profile"],
                    "night_light_opacity": kwargs["night_light_opacity"],
                }
            )

    zstarview_pipeline_module._draw_sky_cloud_layers(
        painter=object(),
        geometry=SimpleNamespace(radius=80),
        scene=replace(_make_scene(), night_light_glow_profile=object()),
        style=_make_style(night_light_opacity=0.12, ridge_glow_opacity=0.34),
        compositor=_Compositor(),
        star_render_surface_size=(200, 200),
        simplified_view_active=True,
    )

    assert captured == {
        "night_light_glow_profile": None,
        "night_light_opacity": 0.0,
    }
