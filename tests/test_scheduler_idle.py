from __future__ import annotations

from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

from astropy.time import Time

from zstarview.gui.window import _calendar_second_delay_ms
from zstarview.gui.window_update_sky import SkyWindowSkyUpdatesMixin
from zstarview.gui.window_updates import SkyWindowUpdatesMixin


class _FakeBusyController:
    def __init__(self, busy: bool = False) -> None:
        self._busy = busy

    def has_in_flight_update(self) -> bool:
        return self._busy


class _SchedulerProbe(SkyWindowUpdatesMixin):
    def __init__(self) -> None:
        self.state = SimpleNamespace(
            jump_highlight_name=None,
            jump_highlight_altaz=None,
            jump_highlight_until_ms=0.0,
            sky_update_pending=False,
            pending_star_vmag_limit=None,
            viewport_interaction_mode=False,
            celestial_data=None,
            sky_next_refresh_utc=None,
            sky_disc_next_refresh_utc=None,
            cloud_next_refresh_utc=None,
            cloud_projection_next_refresh_utc=None,
            meteor_next_refresh_utc=None,
            persistent_search_next_refresh_utc=None,
            satellite_next_refresh_utc=None,
            satellite_projection_next_refresh_utc=None,
            aircraft_next_refresh_utc=None,
            aircraft_projection_next_refresh_utc=None,
            tropical_cyclone_projection_next_refresh_utc=None,
            dynamic_display_second=None,
            dynamic_display_time=None,
            dynamic_planets=None,
            dynamic_planet_bucket=None,
            prepared_dynamic_planets=None,
            prepared_dynamic_planet_bucket=None,
            twinkle_bucket=None,
            twinkle_targets=(),
            simplified_view_enabled=False,
            simplified_view_labels_enabled=True,
        )
        self._is_shutting_down = False
        self.sky_update_interval = 600
        self._sky_worker = _FakeBusyController(False)
        self._geo_satellite_enabled = False
        self._clouddisc = object()
        self.cloud_disc_alpha = 1.0
        self._cloud_controller = _FakeBusyController(False)
        self._geosatellite_controller = None
        self._satellite_controller = _FakeBusyController(False)
        self._aircraft_controller = _FakeBusyController(False)
        self._meteor_controller = _FakeBusyController(False)
        self._precipitation_controller = None
        self._road_night_lights_controller = None
        self._tropical_cyclone_controller = _FakeBusyController(False)
        self._jpl_small_body_controller = _FakeBusyController(False)
        self._terrain_horizon_controller = _FakeBusyController(False)
        self._water_overlay_controller = _FakeBusyController(False)
        self._urban_outline_controller = _FakeBusyController(False)
        self.presentation_id = "scenic"
        self.twinkle_enabled = True
        self.twinkle_count = 30
        self.meteor_opacity = 0.0
        self.tropical_cyclone_state = SimpleNamespace(
            projection_next_refresh_utc=None,
            next_check_utc=None,
            next_refresh_utc=None,
        )
        self.start_calls: list[tuple[str, dict[str, object]]] = []
        self.client_updates = 0
        self.current_time = Time(101.2, format="unix")

    def _current_time_obj(self) -> Time:
        return self.current_time

    def _request_dynamic_planet_update(self, target_time: Time | None = None) -> None:
        self.requested_planet_time = target_time

    def request_client_update(self) -> None:
        self.client_updates += 1

    def start_background_sky_data_update(self, **kwargs: object) -> bool:
        self.start_calls.append(("sky", dict(kwargs)))
        return True

    def start_background_cloud_update(self, **kwargs: object) -> None:
        self.start_calls.append(("cloud", dict(kwargs)))
        self._cloud_controller._busy = True

    def start_background_satellite_update(self, **kwargs: object) -> bool:
        self.start_calls.append(("satellite", dict(kwargs)))
        return True

    def start_background_aircraft_update(self, **kwargs: object) -> bool:
        self.start_calls.append(("aircraft", dict(kwargs)))
        return True

    def _satellite_layer_enabled(self) -> bool:
        return True

    def _aircraft_layer_enabled(self) -> bool:
        return True

    def reproject_satellite_overlay(self) -> None:
        self.start_calls.append(("satellite_projection", {}))
        self.state.satellite_projection_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=2
        )

    def reproject_aircraft_overlay(self) -> None:
        self.start_calls.append(("aircraft_projection", {}))
        self.state.aircraft_projection_next_refresh_utc = datetime.now(timezone.utc) + timedelta(
            seconds=2
        )

    def _tropical_cyclone_layer_enabled(self) -> bool:
        return True

    def reproject_tropical_cyclone_overlay(self) -> None:
        self.start_calls.append(("tropical_projection", {}))
        next_refresh = datetime.now(timezone.utc) + timedelta(seconds=2)
        self.state.tropical_cyclone_projection_next_refresh_utc = next_refresh
        self.tropical_cyclone_state.projection_next_refresh_utc = next_refresh

    def _start_persistent_search_refresh(self, *, reason: str = "timer") -> bool:
        self.start_calls.append(("jpl", {"reason": reason}))
        return True


def test_calendar_second_delay_tracks_next_boundary() -> None:
    assert _calendar_second_delay_ms(100.0) == 1000
    assert _calendar_second_delay_ms(100.001) == 999
    assert _calendar_second_delay_ms(100.999) == 1


def test_planet_worker_result_is_buffered_without_repaint() -> None:
    probe = SimpleNamespace(
        _is_shutting_down=False,
        _disc_generation=7,
        state=SimpleNamespace(
            prepared_dynamic_planets=None,
            prepared_dynamic_planet_bucket=None,
            dynamic_planet_requested_bucket=51,
        ),
        request_client_update=lambda: (_ for _ in ()).throw(
            AssertionError("worker completion must not repaint")
        ),
    )
    planets = [object()]

    SkyWindowSkyUpdatesMixin._on_planet_data_calculated(
        probe,
        {
            "planets": planets,
            "time_unix": 102.0,
            "render_generation": 7,
        },
    )

    assert probe.state.prepared_dynamic_planets == planets
    assert probe.state.prepared_dynamic_planet_bucket == 51
    assert probe.state.dynamic_planet_requested_bucket is None


def test_scheduler_tick_runs_one_task_per_idle_tick() -> None:
    probe = _SchedulerProbe()
    now = datetime.now(timezone.utc)
    probe.state.sky_next_refresh_utc = now - timedelta(seconds=1)
    probe.state.cloud_next_refresh_utc = now - timedelta(seconds=1)

    probe._on_scheduler_tick()
    assert [name for name, _ in probe.start_calls] == ["sky"]
    assert probe.state.sky_next_refresh_utc is not None
    assert probe.state.cloud_next_refresh_utc is not None

    probe._on_scheduler_tick()
    assert [name for name, _ in probe.start_calls] == ["sky", "cloud"]


def test_scheduler_tick_skips_when_busy() -> None:
    probe = _SchedulerProbe()
    now = datetime.now(timezone.utc)
    probe.state.sky_next_refresh_utc = now - timedelta(seconds=1)
    probe.state.cloud_next_refresh_utc = now - timedelta(seconds=1)
    probe._sky_worker = _FakeBusyController(True)

    probe._on_scheduler_tick()

    assert probe.start_calls == []


def test_scheduler_tick_does_not_repaint_dynamic_layers_on_odd_second() -> None:
    probe = _SchedulerProbe()
    probe._clouddisc = None
    probe._sky_worker = _FakeBusyController(True)
    probe.state.aircraft_projection_next_refresh_utc = datetime.now(timezone.utc) - timedelta(
        seconds=1
    )
    probe.state.satellite_projection_next_refresh_utc = datetime.now(timezone.utc) - timedelta(
        seconds=1
    )

    probe._on_scheduler_tick()

    assert probe.start_calls == []
    assert probe.client_updates == 0


def test_scheduler_tick_publishes_dynamic_layers_once_on_even_second() -> None:
    probe = _SchedulerProbe()
    probe.current_time = Time(102.8, format="unix")
    planets = [object()]
    probe.state.prepared_dynamic_planets = planets
    probe.state.prepared_dynamic_planet_bucket = 51

    probe._on_scheduler_tick()
    assert probe.client_updates == 1
    assert probe.state.dynamic_display_second == 102
    assert abs(float(probe.state.dynamic_display_time.unix) - 102.0) < 1e-9
    assert probe.state.dynamic_planets == planets

    probe._on_scheduler_tick()
    assert probe.client_updates == 1


def test_scheduler_tick_projects_tropical_cyclone_even_when_busy() -> None:
    probe = _SchedulerProbe()
    probe._sky_worker = _FakeBusyController(True)
    probe.state.tropical_cyclone_projection_next_refresh_utc = datetime.now(timezone.utc) - timedelta(
        seconds=1
    )

    probe._on_scheduler_tick()

    assert [name for name, _ in probe.start_calls] == ["tropical_projection"]
    assert probe.state.tropical_cyclone_projection_next_refresh_utc is not None
    assert (
        probe.state.tropical_cyclone_projection_next_refresh_utc
        > datetime.now(timezone.utc)
    )


def test_scheduler_tick_groups_data_refreshes_together() -> None:
    probe = _SchedulerProbe()
    now = datetime.now(timezone.utc)
    probe.state.cloud_next_refresh_utc = now - timedelta(seconds=1)
    probe.state.satellite_next_refresh_utc = now - timedelta(seconds=1)
    probe.state.aircraft_next_refresh_utc = now - timedelta(seconds=1)

    probe._on_scheduler_tick()
    assert [name for name, _ in probe.start_calls] == ["cloud"]

    probe._cloud_controller._busy = False
    probe._on_scheduler_tick()
    assert [name for name, _ in probe.start_calls] == ["cloud", "satellite"]

    probe._on_scheduler_tick()
    assert [name for name, _ in probe.start_calls] == ["cloud", "satellite", "aircraft"]
