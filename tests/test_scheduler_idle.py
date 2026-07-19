from __future__ import annotations

from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

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
            cloud_next_refresh_utc=None,
            cloud_projection_next_refresh_utc=None,
            persistent_search_next_refresh_utc=None,
            satellite_next_refresh_utc=None,
            satellite_projection_next_refresh_utc=None,
            aircraft_next_refresh_utc=None,
            aircraft_projection_next_refresh_utc=None,
            tropical_cyclone_projection_next_refresh_utc=None,
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
        self._tropical_cyclone_controller = _FakeBusyController(False)
        self._jpl_small_body_controller = _FakeBusyController(False)
        self._terrain_horizon_controller = _FakeBusyController(False)
        self._water_overlay_controller = _FakeBusyController(False)
        self._urban_outline_controller = _FakeBusyController(False)
        self.tropical_cyclone_state = SimpleNamespace(
            projection_next_refresh_utc=None,
            next_check_utc=None,
            next_refresh_utc=None,
        )
        self.start_calls: list[tuple[str, dict[str, object]]] = []
        self.client_updates = 0

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


def test_scheduler_tick_keeps_overlay_projection_when_busy() -> None:
    probe = _SchedulerProbe()
    probe._sky_worker = _FakeBusyController(True)
    probe.state.aircraft_projection_next_refresh_utc = datetime.now(timezone.utc) - timedelta(
        seconds=1
    )
    probe.state.satellite_projection_next_refresh_utc = datetime.now(timezone.utc) - timedelta(
        seconds=1
    )

    probe._on_scheduler_tick()

    assert [name for name, _ in probe.start_calls] == ["aircraft_projection"]


def test_scheduler_tick_puts_overlay_projection_after_cloud() -> None:
    probe = _SchedulerProbe()
    now = datetime.now(timezone.utc)
    probe.state.cloud_next_refresh_utc = now - timedelta(seconds=1)
    probe.state.aircraft_projection_next_refresh_utc = datetime.now(timezone.utc) - timedelta(
        seconds=1
    )
    probe.state.satellite_projection_next_refresh_utc = datetime.now(timezone.utc) - timedelta(
        seconds=1
    )

    probe._on_scheduler_tick()
    assert [name for name, _ in probe.start_calls] == ["cloud"]

    probe._cloud_controller._busy = False
    probe._on_scheduler_tick()
    assert [name for name, _ in probe.start_calls] == [
        "cloud",
        "aircraft_projection",
    ]

    probe._on_scheduler_tick()
    assert [name for name, _ in probe.start_calls] == [
        "cloud",
        "aircraft_projection",
        "satellite_projection",
    ]


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
