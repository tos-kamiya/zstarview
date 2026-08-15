
from tests._window_render_support import *


def test_on_sky_data_calculated_updates_render_snapshot_once() -> None:
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
        "view_center": (20.0, 30.0),
        "geometry": render_geometry.get_screen_geometry(
            640, 480, dummy.viewer_data.view_alt_deg
        ),
        "render_generation": 0,
    }
    SkyWindow._on_sky_data_calculated(dummy, payload)

    assert dummy.state.render_view_center == (20.0, 30.0)
    assert dummy.state.celestial_data is celestial
    assert dummy.state.sky_disc_image is sky_disc
    assert dummy._compositor.invalidated is True
    assert update_calls == ["update"]


def test_on_sky_data_calculated_preserves_render_center_during_viewport_interaction() -> (
    None
):
    dummy = _WindowStub()
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
    dummy._startup_initial_load_started = False
    dummy._startup_initial_data_loaded = True
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
            "view_center": (40.0, 150.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.render_view_center == (40.0, 150.0)


def test_on_sky_data_calculated_triggers_release_followup_updates() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(40.0, 150.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(40.0, 150.0),
        viewport_interaction_release_pending=True,
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
    dummy.request_client_update = lambda: None
    dummy.reproject_tropical_cyclone_overlay = Mock()

    class _MenuButton:
        def __init__(self) -> None:
            self.visible = False

        def setVisible(self, value: bool) -> None:
            self.visible = value

    dummy.menu_button = _MenuButton()
    cloud_calls: list[str] = []
    dummy.start_background_cloud_update = lambda **kwargs: cloud_calls.append(
        str(kwargs.get("reason"))
    )
    dummy.start_background_terrain_horizon_update = lambda **kwargs: cloud_calls.append(
        str(kwargs.get("reason"))
    )

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (40.0, 150.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.viewport_interaction_release_pending is False
    assert dummy.state.viewport_interaction_mode is False
    assert dummy.menu_button.visible is True
    assert cloud_calls == ["view-change-release", "view-change-release"]
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with()


def test_on_sky_data_calculated_uses_idle_completion_reason() -> None:
    dummy = _WindowStub()
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
        viewport_interaction_release_pending=True,
        viewport_interaction_completion_reason="view-change-idle",
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
    dummy.request_client_update = lambda: None
    dummy.reproject_tropical_cyclone_overlay = Mock()

    class _MenuButton:
        def __init__(self) -> None:
            self.visible = False

        def setVisible(self, value: bool) -> None:
            self.visible = value

    dummy.menu_button = _MenuButton()
    cloud_calls: list[str] = []
    dummy.start_background_cloud_update = lambda **kwargs: cloud_calls.append(
        str(kwargs.get("reason"))
    )
    dummy.start_background_terrain_horizon_update = lambda **kwargs: cloud_calls.append(
        str(kwargs.get("reason"))
    )

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (40.0, 150.0),
            "geometry": render_geometry.get_screen_geometry(
                640,
                480,
                dummy.viewer_data.view_alt_deg,
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.viewport_interaction_release_pending is False
    assert dummy.state.viewport_interaction_completion_reason is None
    assert dummy.state.viewport_interaction_mode is False
    assert dummy.menu_button.visible is True
    assert cloud_calls == ["view-change-idle", "view-change-idle"]
    dummy.reproject_tropical_cyclone_overlay.assert_called_once_with()


def test_on_sky_data_calculated_keeps_existing_cloud_refresh_deadline() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    cloud_deadline = datetime(2026, 5, 25, 0, 10, tzinfo=timezone.utc)
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        cloud_next_refresh_utc=cloud_deadline,
    )
    dummy._compositor = _DummyCompositor()
    dummy._sky_data_update_timer = _DummyTimer(active=True)
    dummy._cloud_update_timer = _DummyTimer(active=False)
    dummy._clouddisc = object()
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
            "view_center": (20.0, 30.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.cloud_next_refresh_utc == cloud_deadline


def test_apply_startup_delta_t_disables_tropical_cyclone_layer_when_time_shifted() -> (
    None
):
    dummy = _WindowStub(
        tropical_cyclone_opacity=0.4,
        show_tropical_cyclone_overlay=True,
        _tropical_cyclone_controller=object(),
        _tropical_cyclone_requested_enabled=True,
        _tropical_cyclone_opacity_when_enabled=0.4,
        _action_toggle_satellites=Mock(),
        _action_toggle_aircraft=Mock(),
        _action_toggle_tropical_cyclone=Mock(),
    )

    window_module.SkyWindow.apply_startup_delta_t(dummy, timedelta(hours=-10))

    assert dummy.tropical_cyclone_opacity == 0.0
    assert dummy.show_tropical_cyclone_overlay is True
    dummy._action_toggle_tropical_cyclone.setEnabled.assert_called_once_with(False)
    dummy._action_toggle_tropical_cyclone.setChecked.assert_not_called()


def test_apply_startup_delta_t_disables_precipitation_when_time_shifted() -> None:
    action = Mock()
    dummy = _WindowStub(
        precipitation_opacity=0.6,
        _precipitation_controller=object(),
        _precipitation_requested_enabled=True,
        _precipitation_opacity_when_enabled=0.6,
        _action_toggle_precipitation=action,
        _action_toggle_satellites=Mock(),
        _action_toggle_aircraft=Mock(),
        _action_toggle_tropical_cyclone=Mock(),
    )

    window_module.SkyWindow.apply_startup_delta_t(dummy, timedelta(hours=10))

    assert dummy.precipitation_opacity == 0.0
    action.setEnabled.assert_called_once_with(False)
    action.setChecked.assert_called_once_with(False)


def test_on_sky_data_calculated_discards_stale_view_center_after_jump() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(14.25, 87.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(
        render_view_center=(14.25, 87.0),
        viewport_interaction_mode=False,
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
    retry_calls: list[str] = []
    dummy.request_sky_data_update = lambda *_args, **kwargs: retry_calls.append(
        str(kwargs.get("reason"))
    )
    dummy._safe_request_cloud_repaint = lambda: None
    dummy.update = lambda: None

    SkyWindow._on_sky_data_calculated(
        dummy,
        {
            "celestial": object(),
            "sky_disc": object(),
            "view_center": (0.0, 180.0),
            "geometry": render_geometry.get_screen_geometry(
                640, 480, dummy.viewer_data.view_alt_deg
            ),
            "render_generation": 0,
        },
    )

    assert dummy.state.render_view_center == (14.25, 87.0)
    assert retry_calls == ["stale-view-center"]


def test_schedule_satellite_retry_after_failure_uses_two_hour_backoff() -> None:
    dummy = _WindowStub()
    dummy.satellite_opacity = 0.5
    dummy._is_shutting_down = False
    dummy.state = SimpleNamespace(satellite_next_refresh_utc=None)
    dummy._satellite_layer_enabled = lambda: True

    SkyWindow._schedule_satellite_retry_after_failure(dummy)

    assert dummy.state.satellite_next_refresh_utc is not None
    assert dummy.state.satellite_next_refresh_utc > datetime.now(timezone.utc)


def test_satellite_validity_remaining_ms_uses_refresh_time() -> None:
    dummy = _WindowStub()
    dummy.satellite_state = SimpleNamespace(
        refreshed_at_utc=datetime.now(timezone.utc),
        element_epoch_utc=datetime(2020, 1, 1, tzinfo=timezone.utc),
    )

    remaining_ms = SkyWindow._satellite_validity_remaining_ms(dummy)

    assert remaining_ms is not None
    assert remaining_ms > 0


def test_on_satellite_failed_schedules_failure_backoff() -> None:
    dummy = _WindowStub()
    dummy.satellite_opacity = 0.5
    dummy.satellite_state = SimpleNamespace(set_error_banner=Mock())
    retry_calls: list[str] = []
    dummy._schedule_satellite_retry_after_failure = lambda: retry_calls.append("retry")
    dummy.update = Mock()

    SkyWindow._on_satellite_failed(dummy, {"banner": "Satellites: timed out"})

    dummy.satellite_state.set_error_banner.assert_called_once_with(
        "Satellites: timed out"
    )
    assert retry_calls == ["retry"]
    dummy.update.assert_called_once()


def test_jump_to_satellite_target_uses_cached_satellite_records_below_horizon(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        window_module, "find_satellite_altaz", lambda *_args, **_kwargs: (-12.0, 123.0)
    )
    monkeypatch.setattr(
        window_module.QTimer,
        "singleShot",
        lambda _ms, func: func(),
    )

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy.satellite_state = SimpleNamespace(
        records_by_group={"iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]},
        set_banner=Mock(),
    )
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()
    dummy._current_time_obj = lambda: astropy.time.Time("2026-03-22T12:00:00Z")
    dummy._find_satellite_jump_altaz = lambda key: SkyWindow._find_satellite_jump_altaz(
        dummy, key
    )

    SkyWindow._jump_to_search_target(
        dummy,
        SearchJumpTarget(
            label="ISS",
            ra_hours=0.0,
            dec_deg=0.0,
            kind="satellite",
            sort_key=(99.0, "iss"),
            subtitle="Satellite",
            object_key="ISS",
        ),
    )

    assert dummy.viewer_data.view_center == (-12.0, 123.0)
    assert dummy.state.jump_highlight_name == "ISS"
    assert dummy.state.jump_highlight_altaz == (-12.0, 123.0)
    dummy.request_sky_data_update.assert_called_once()


def test_find_satellite_jump_altaz_falls_back_to_disk_cache(monkeypatch) -> None:
    monkeypatch.setattr(
        window_module,
        "find_satellite_altaz",
        lambda records_by_group, **_kwargs: (
            (-40.0, 151.0) if records_by_group.get("iss") else None
        ),
    )
    monkeypatch.setattr(
        window_module,
        "load_satellite_cache",
        lambda group_key: (
            SimpleNamespace(records=[{"OBJECT_NAME": "ISS (ZARYA)"}])
            if group_key == "iss"
            else None
        ),
    )

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(40.7128, -74.0060),
        timezone_name="America/New_York",
        city_name="New York City",
        view_center=(20.0, 30.0),
        observer_height_m=10.0,
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={})
    dummy._enabled_satellite_groups = ("iss",)
    dummy._current_time_obj = lambda: astropy.time.Time("2026-03-23T12:13:24Z")
    dummy._load_cached_satellite_records = lambda groups: (
        SkyWindow._load_cached_satellite_records(dummy, groups)
    )

    altaz = SkyWindow._find_satellite_jump_altaz(dummy, "ISS")

    assert altaz == (-40.0, 151.0)


def test_jump_to_satellite_target_sets_banner_when_not_available() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy.update = Mock()
    dummy._load_cached_satellite_records = lambda _groups: {}
    dummy._find_satellite_jump_altaz = lambda key: SkyWindow._find_satellite_jump_altaz(
        dummy, key
    )

    SkyWindow._jump_to_search_target(
        dummy,
        SearchJumpTarget(
            label="ISS",
            ra_hours=0.0,
            dec_deg=0.0,
            kind="satellite",
            sort_key=(99.0, "iss"),
            subtitle="Satellite",
            object_key="ISS",
        ),
    )

    dummy.satellite_state.set_banner.assert_called_once_with(
        "Satellites: ISS not available"
    )
    dummy.update.assert_called_once()


def test_jump_to_place_target_uses_projected_altaz(monkeypatch) -> None:
    monkeypatch.setattr(
        window_module,
        "_project_place_targets_to_altaz",
        lambda **kwargs: (
            PlaceTargetProjection(
                alt_deg=-3.5,
                az_deg=145.0,
                distance_km=12.0,
                target_latitude_deg=float(kwargs["target_latitude_deg"][0]),
                target_longitude_deg=float(kwargs["target_longitude_deg"][0]),
                target_height_m=float(kwargs["target_height_m"][0]),
            ),
        ),
    )
    monkeypatch.setattr(
        window_module.QTimer,
        "singleShot",
        lambda _ms, func: func(),
    )

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy.state = SkyWindowState(render_view_center=(20.0, 30.0))
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()

    SkyWindow._jump_to_search_target(
        dummy,
        SearchJumpTarget(
            label="Tokyo Station",
            kind="place",
            sort_key=(0.0, "tokyo station"),
            subtitle="Place / railway / station",
            latitude_deg=35.681236,
            longitude_deg=139.767125,
        ),
    )

    # Allow the view center to go below the horizon (>= -45.0°) per policy.
    assert dummy.viewer_data.view_center == (-3.5, 145.0)
    assert dummy.state.jump_highlight_name == "Tokyo Station"
    assert dummy.state.jump_highlight_altaz == (-3.5, 145.0)
    dummy.request_sky_data_update.assert_called_once()


def test_jump_to_jpl_small_body_target_can_set_persistent_overlay(caplog) -> None:
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
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)
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
            lambda _target, **_kwargs: (12.5, 220.0),
        )
        with caplog.at_level("INFO"):
            SkyWindow._jump_to_search_target(
                dummy,
                SearchJumpTarget(
                    label="Ceres",
                    kind="jpl_small_body",
                    sort_key=(0.0, "ceres"),
                    subtitle="Asteroid / 1 Ceres",
                    object_key="20000001",
                    command="DES=20000001;",
                    persistent_keep_marker=True,
                ),
            )

    assert dummy.viewer_data.view_center == (12.5, 220.0)
    assert dummy.state.jump_highlight_name == "Ceres"
    assert dummy.state.jump_highlight_altaz == (12.5, 220.0)
    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.label == "Ceres"
    assert dummy.state.persistent_search_target.persistent_keep_marker is True
    assert dummy.state.persistent_search_next_refresh_utc == datetime(
        2026, 4, 18, 13, 0, tzinfo=timezone.utc
    )
    assert (
        "JPL persistent target set: label=Ceres kind=jpl_small_body group=<none>"
        in caplog.text
    )
    assert "target_time_utc=2026-04-18T12:00:00+00:00" in caplog.text
    assert "alt=12.5 az=220.0 command=DES=20000001;" in caplog.text


def test_jump_to_jpl_small_body_target_uses_state_vector_when_present(
    monkeypatch,
) -> None:
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
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)
    dummy._sync_view_altitude_actions = Mock()
    dummy._begin_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.update = Mock()
    dummy._current_time_obj = lambda: astropy.time.Time("2026-04-18T12:00:00Z")

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
        lambda _target, **_kwargs: (11.5, 221.0),
    )

    SkyWindow._jump_to_search_target(
        dummy,
        SearchJumpTarget(
            label="Ceres",
            kind="jpl_small_body",
            sort_key=(0.0, "ceres"),
            subtitle="Asteroid / 1 Ceres",
            object_key="20000001",
            command="DES=20000001;",
            persistent_keep_marker=True,
        ),
    )

    assert dummy.viewer_data.view_center == (11.5, 221.0)
    assert dummy.state.jump_highlight_altaz == (11.5, 221.0)
    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.horizons_position_km == (1.0, 2.0, 3.0)
    assert dummy.state.persistent_search_target.horizons_velocity_km_s == (
        0.1,
        0.2,
        0.3,
    )


def test_jump_to_jpl_small_body_target_honors_fixed_search_axes() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(5.0, 210.0),
        observer_height_m=1.7,
    )
    dummy._search_view_center_base = (5.0, 210.0)
    dummy._search_view_center_alt_specified = True
    dummy._search_view_center_az_specified = False
    dummy.state = SkyWindowState(
        render_view_center=(5.0, 210.0),
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
            lambda _target, **_kwargs: (12.5, 220.0),
        )
        SkyWindow._jump_to_search_target(
            dummy,
            SearchJumpTarget(
                label="Ceres",
                kind="jpl_small_body",
                sort_key=(0.0, "ceres"),
                subtitle="Asteroid / 1 Ceres",
                object_key="20000001",
                command="DES=20000001;",
            ),
        )

    assert dummy.viewer_data.view_center == (5.0, 220.0)


def test_jump_to_jpl_small_body_target_can_disable_fixed_search_axes() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(5.0, 210.0),
        observer_height_m=1.7,
    )
    dummy._search_view_center_base = (5.0, 210.0)
    dummy._search_view_center_alt_specified = True
    dummy._search_view_center_az_specified = True
    dummy.state = SkyWindowState(
        render_view_center=(5.0, 210.0),
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
            lambda _target, **_kwargs: (12.5, 220.0),
        )
        SkyWindow._jump_to_search_target(
            dummy,
            SearchJumpTarget(
                label="Ceres",
                kind="jpl_small_body",
                sort_key=(0.0, "ceres"),
                subtitle="Asteroid / 1 Ceres",
                object_key="20000001",
                command="DES=20000001;",
                preserve_cli_view_center=False,
            ),
        )

    assert dummy.viewer_data.view_center == (12.5, 220.0)


def test_jump_to_jpl_small_body_target_without_keep_flags_clears_overlay() -> None:
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
        persistent_search_target=SearchJumpTarget(
            label="Old",
            kind="jpl_small_body",
            sort_key=(0.0, "old"),
            alt_deg=10.0,
            az_deg=30.0,
            persistent_keep_marker=True,
        ),
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
            lambda _target, **_kwargs: (12.5, 220.0),
        )
        SkyWindow._jump_to_search_target(
            dummy,
            SearchJumpTarget(
                label="Ceres",
                kind="jpl_small_body",
                sort_key=(0.0, "ceres"),
                subtitle="Asteroid / 1 Ceres",
                object_key="20000001",
                command="DES=20000001;",
            ),
        )

    assert dummy.state.persistent_search_target is None


def test_jpl_small_body_failure_reschedules_one_hour_later() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    current_target = SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "ceres"),
        subtitle="Asteroid / 1 Ceres",
        object_key="20000001",
        command="DES=20000001;",
        alt_deg=12.5,
        az_deg=220.0,
        target_time_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
        persistent_keep_marker=True,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=current_target,
        persistent_search_next_refresh_utc=datetime(
            2026, 4, 18, 13, 0, tzinfo=timezone.utc
        ),
        persistent_search_reference_time_utc=datetime(
            2026, 4, 18, 12, 0, tzinfo=timezone.utc
        ),
    )
    dummy.satellite_state = SimpleNamespace(records_by_group={}, set_banner=Mock())
    dummy.request_client_update = Mock()
    dummy._schedule_persistent_search_refresh = Mock()

    SkyWindow._on_jpl_failed(
        dummy,
        {
            "target": current_target,
            "target_time_utc": datetime(2026, 4, 18, 13, 0, tzinfo=timezone.utc),
            "refreshed_at_utc": datetime(2026, 4, 18, 13, 2, tzinfo=timezone.utc),
            "banner": "JPL: timed out",
            "error": "timed out",
            "reason": "timer",
        },
    )

    assert dummy.state.persistent_search_last_error == "timed out"
    assert dummy.state.persistent_search_last_refresh_utc == datetime(
        2026, 4, 18, 13, 2, tzinfo=timezone.utc
    )
    assert dummy.state.persistent_search_next_refresh_utc == datetime(
        2026, 4, 18, 14, 2, tzinfo=timezone.utc
    )
    dummy._schedule_persistent_search_refresh.assert_called_once()
    dummy.request_client_update.assert_called_once()


def test_on_jpl_ready_logs_refreshed_persistent_target(caplog) -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._current_time_obj = lambda: astropy.time.Time("2026-04-18T13:00:00Z")
    current_target = SearchJumpTarget(
        label="Voyager 1",
        kind="jpl_small_body",
        sort_key=(0.0, "voyager 1"),
        subtitle="spacecraft",
        object_key="-31",
        command="DES=-31;",
        alt_deg=48.6,
        az_deg=245.6,
        horizons_epoch_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
        horizons_position_km=(1.0, 2.0, 3.0),
        horizons_velocity_km_s=(0.1, 0.2, 0.3),
        target_time_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
        jpl_group="sb",
        persistent_keep_marker=True,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=current_target,
        persistent_search_next_refresh_utc=datetime(
            2026, 4, 18, 13, 0, tzinfo=timezone.utc
        ),
        persistent_search_reference_time_utc=datetime(
            2026, 4, 18, 12, 0, tzinfo=timezone.utc
        ),
    )
    dummy.request_client_update = Mock()
    dummy._schedule_persistent_search_refresh = Mock()
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(
            window_module,
            "project_jpl_target_altaz_from_state_vector",
            lambda _target, **_kwargs: (49.1, 244.7),
        )

        with caplog.at_level("INFO"):
            SkyWindow._on_jpl_ready(
                dummy,
                {
                    "target": current_target,
                    "target_time_utc": datetime(
                        2026, 4, 18, 13, 0, tzinfo=timezone.utc
                    ),
                    "refreshed_at_utc": datetime(
                        2026, 4, 18, 13, 2, tzinfo=timezone.utc
                    ),
                    "horizons_epoch_utc": datetime(
                        2026, 4, 18, 13, 0, tzinfo=timezone.utc
                    ),
                    "horizons_position_km": (4.0, 5.0, 6.0),
                    "horizons_velocity_km_s": (0.4, 0.5, 0.6),
                    "rows": [["2026-Apr-18 13:00:00", "*", "m", "244.7", "49.1"]],
                    "reason": "timer",
                },
            )

    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.alt_deg == 49.1
    assert dummy.state.persistent_search_target.az_deg == 244.7
    assert dummy.state.persistent_search_target.horizons_position_km == (4.0, 5.0, 6.0)
    assert dummy.state.persistent_search_target.horizons_velocity_km_s == (
        0.4,
        0.5,
        0.6,
    )
    assert (
        "JPL persistent target refreshed: label=Voyager 1 kind=jpl_small_body group=sb"
        in caplog.text
    )
    assert "target_time_utc=2026-04-18T13:00:00+00:00" in caplog.text
    assert "alt=49.1 az=244.7 command=DES=-31;" in caplog.text


def test_refresh_projected_persistent_search_target_reprojects_state_vector(
    monkeypatch,
) -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._current_time_obj = lambda: astropy.time.Time("2026-04-18T13:00:00Z")
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        persistent_search_target=SearchJumpTarget(
            label="Voyager 1",
            kind="jpl_small_body",
            sort_key=(0.0, "voyager 1"),
            command="DES=-31;",
            alt_deg=48.6,
            az_deg=245.6,
            horizons_epoch_utc=datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc),
            horizons_position_km=(1.0, 2.0, 3.0),
            horizons_velocity_km_s=(0.1, 0.2, 0.3),
            persistent_keep_marker=True,
        ),
    )
    dummy.request_client_update = Mock()
    monkeypatch.setattr(
        "zstarview.gui.window_update_overlays.project_jpl_target_altaz_from_state_vector",
        lambda _target, **_kwargs: (12.5, 220.0),
    )

    window_updates_module.SkyWindowUpdatesMixin.refresh_projected_persistent_search_target(
        dummy
    )

    assert dummy.state.persistent_search_target is not None
    assert dummy.state.persistent_search_target.alt_deg == 12.5
    assert dummy.state.persistent_search_target.az_deg == 220.0
    assert dummy.request_client_update.called


def test_handle_client_resize_discards_cached_sky_disc_and_requests_refresh() -> None:
    dummy = _WindowStub(_startup_initial_data_loaded=True)
    sky_disc_image = QImage(4, 4, QImage.Format.Format_ARGB32_Premultiplied)
    dummy._frameless_frame = None
    dummy.menu_button = None
    dummy.size_grip = None
    dummy._compositor = _DummyCompositor()
    dummy._begin_viewport_interaction_mode = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.request_client_update = Mock()
    dummy.start_background_cloud_update = Mock()
    dummy._raise_overlay_widgets = Mock()
    dummy._discard_stale_disc_images = lambda: (
        window_module.SkyWindowCoreMixin._discard_stale_disc_images(  # type: ignore[attr-defined]
            dummy
        )
    )
    dummy.cloud_controller = None
    dummy._cloud_controller = None
    dummy.cloud_state = SimpleNamespace(
        image=np.zeros((2, 2, 4), dtype=np.uint8),
        missing_mask=np.ones((2, 2), dtype=np.uint8),
        cloud_amount_field=SimpleNamespace(),
        render_key=object(),
        request_id=1,
        missing_mask_key=2,
    )
    dummy.state = SkyWindowState(
        render_view_center=(20.0, 30.0),
        sky_disc_image=sky_disc_image,
    )
    dummy.width = lambda: 200
    dummy.height = lambda: 100
    dummy.client_width = lambda: 200
    dummy.client_height = lambda: 100

    event = QResizeEvent(QSize(220, 120), QSize(200, 100))

    SkyWindow._handle_client_resize(dummy, event)

    assert dummy._disc_generation == 1
    assert dummy.state.sky_disc_image is None
    assert dummy.cloud_state.image is None
    assert dummy.cloud_state.missing_mask is None
    assert dummy.cloud_state.cloud_amount_field is None
    assert dummy.cloud_state.render_key is None
    assert dummy.cloud_state.request_id is None
    assert dummy.cloud_state.missing_mask_key is None
    assert dummy._compositor.invalidated is True
    dummy._begin_viewport_interaction_mode.assert_called_once()
    dummy.request_sky_data_update.assert_called_once_with(
        reason="resize",
        allow_during_viewport_interaction=True,
    )
    dummy.request_client_update.assert_called_once()


def test_handle_client_resize_during_startup_defers_refresh() -> None:
    dummy = _WindowStub(_startup_initial_data_loaded=False)
    dummy._layout_startup_log_overlay = Mock()
    dummy._raise_overlay_widgets = Mock()
    dummy.request_sky_data_update = Mock()
    dummy.request_client_update = Mock()
    dummy._begin_viewport_interaction_mode = Mock()
    dummy._compositor = _DummyCompositor()
    dummy.width = lambda: 200
    dummy.height = lambda: 100
    dummy.client_width = lambda: 200
    dummy.client_height = lambda: 100

    event = QResizeEvent(QSize(220, 120), QSize(200, 100))

    SkyWindow._handle_client_resize(dummy, event)

    assert dummy._startup_resize_pending is True
    dummy._begin_viewport_interaction_mode.assert_not_called()
    dummy.request_sky_data_update.assert_not_called()
    dummy.request_client_update.assert_called_once()


def test_search_satellite_targets_resolves_known_artificial_satellites() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    targets = SkyWindow._search_satellite_targets(dummy, "ISS")

    assert len(targets) == 1
    assert targets[0].label == "ISS"
    assert targets[0].kind == "satellite"
    assert targets[0].alt_deg is None
    assert targets[0].az_deg is None


def test_search_jpl_targets_skips_solar_system_bodies(monkeypatch) -> None:
    lookup_calls: list[str] = []

    def fake_lookup(_query: str, *, group: str):
        lookup_calls.append(group)
        return {"count": 0, "result": []}

    monkeypatch.setattr(window_module, "fetch_horizons_lookup", fake_lookup)

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    targets = SkyWindow._search_jpl_targets(dummy, "Mars")

    assert targets == []
    assert lookup_calls == []


def test_search_jpl_targets_limits_candidates_to_500(monkeypatch) -> None:
    lookup_calls: list[str] = []

    def fake_lookup(_query: str, *, group: str):
        lookup_calls.append(group)
        if group == "mb":
            return {
                "count": 600,
                "result": [
                    {
                        "name": f"PANSTARRS-{idx}",
                        "pdes": str(idx),
                        "spkid": str(1000 + idx),
                        "type": "asteroid",
                    }
                    for idx in range(600)
                ],
            }
        return {"count": 0, "result": []}

    monkeypatch.setattr(window_module, "fetch_horizons_lookup", fake_lookup)

    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )
    dummy._target_time_utc = lambda: datetime(2026, 4, 18, 12, 0, tzinfo=timezone.utc)

    targets = SkyWindow._search_jpl_targets(dummy, "PANSTARRS")

    assert lookup_calls == ["sct", "mb", "sb"]
    assert len(targets) == 500
    assert targets[0].label == "PANSTARRS-0"
    assert targets[-1].label == "PANSTARRS-499"


def test_search_jpl_targets_skips_solar_system_bodies_directly() -> None:
    dummy = _WindowStub()
    dummy.viewer_data = ViewerData(
        location=(35.0, 139.0),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        view_center=(20.0, 30.0),
        observer_height_m=1.7,
    )

    assert SkyWindow._search_jpl_targets(dummy, "Sun") == []
    assert SkyWindow._search_jpl_targets(dummy, "Moon") == []
    assert SkyWindow._search_jpl_targets(dummy, "Mars") == []
