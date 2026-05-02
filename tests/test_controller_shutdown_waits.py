from __future__ import annotations

import threading
import time
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

import astropy.time
import numpy as np

from zstarview.gui.aircraft_controller import AircraftController
from zstarview.gui.jpl_small_body_controller import JplSmallBodyController
from zstarview.gui.satellite_controller import SatelliteController
from zstarview.gui.sky_worker import SkyDataWorker
from zstarview.gui.terrain_controller import TerrainHorizonController
from zstarview.gui.urban_outline_controller import UrbanOutlineController
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.types import ViewerData


def _assert_shutdown_waits(monkeypatch, controller, trigger_update, *, patch_attr: str = "_run_update") -> None:
    started = threading.Event()
    release = threading.Event()

    def fake_worker(**kwargs):
        started.set()
        release.wait(timeout=2.0)

    monkeypatch.setattr(controller, patch_attr, fake_worker)

    trigger_update()

    assert started.wait(timeout=1.0)

    finished = threading.Event()

    def run_shutdown() -> None:
        controller.shutdown(wait_timeout_s=1.0)
        finished.set()

    shutdown_thread = threading.Thread(target=run_shutdown)
    shutdown_thread.start()

    time.sleep(0.05)
    assert not finished.is_set()

    release.set()
    assert finished.wait(timeout=1.0)
    shutdown_thread.join(timeout=1.0)


def test_satellite_controller_shutdown_waits(monkeypatch) -> None:
    controller = SatelliteController()

    def trigger_update() -> None:
        controller.update(
            observer_lat=35.0,
            observer_lon=139.0,
            observer_height_m=1.7,
            time_obj=astropy.time.Time(datetime(2026, 4, 25, 0, 0, tzinfo=timezone.utc)),
            enabled_groups=("ISS",),
            reason="manual",
        )

    _assert_shutdown_waits(monkeypatch, controller, trigger_update)


def test_aircraft_controller_shutdown_waits(monkeypatch) -> None:
    controller = AircraftController()

    def trigger_update() -> None:
        controller.update(
            observer_lat=35.0,
            observer_lon=139.0,
            observer_height_m=1.7,
            time_obj=astropy.time.Time(datetime(2026, 4, 25, 0, 0, tzinfo=timezone.utc)),
            reason="manual",
        )

    _assert_shutdown_waits(monkeypatch, controller, trigger_update)


def test_jpl_controller_shutdown_waits(monkeypatch) -> None:
    controller = JplSmallBodyController()
    target = SimpleNamespace(label="demo", command="COMMAND", object_key="")

    def trigger_update() -> None:
        controller.update(
            observer_lat=35.0,
            observer_lon=139.0,
            observer_height_m=1.7,
            target=target,
            target_time_utc=datetime(2026, 4, 25, 0, 0, tzinfo=timezone.utc),
            reason="manual",
        )

    _assert_shutdown_waits(monkeypatch, controller, trigger_update)


def test_terrain_controller_shutdown_waits(monkeypatch, tmp_path) -> None:
    controller = TerrainHorizonController(cache_dir=tmp_path)

    def trigger_update() -> None:
        controller.update(lat=35.0, lon=139.0, observer_height_m=1.7, reason="manual")

    _assert_shutdown_waits(monkeypatch, controller, trigger_update)


def test_urban_outline_controller_shutdown_waits(monkeypatch, tmp_path) -> None:
    controller = UrbanOutlineController(derived_root_dir=tmp_path)
    viewer_data = ViewerData(location=(35.0, 139.0), timezone_name="Asia/Tokyo", city_name="Tokyo")

    monkeypatch.setattr(controller, "_required_derived_dirs", lambda _viewer_data: ())
    monkeypatch.setattr(controller, "_selected_skyscraper_tiles", lambda _viewer_data: ())

    def trigger_update() -> None:
        controller.update(viewer_data=viewer_data, reason="manual")

    _assert_shutdown_waits(monkeypatch, controller, trigger_update)


def test_sky_worker_shutdown_waits(monkeypatch) -> None:
    controller = SkyDataWorker()

    def trigger_update() -> None:
        controller.update(
            lat=35.0,
            lon=139.0,
            observer_height_m=1.7,
            view_center=(90.0, 180.0),
            star_catalog=np.empty(0, dtype=object),
            dso_catalog=None,
            star_vmag_limit=None,
            star_subset_indices=None,
            delta_t=timedelta(0),
            sky_disc_alpha=0.0,
            sky_disc_base_size=16,
            edge_fov_deg=90.0,
            content_fov_deg=90.0,
            theme=THEME_STYLES_BY_PRESET["night"],
            star_catalog_meta=None,
            render_width_px=16,
            render_height_px=16,
            render_generation=0,
        )

    _assert_shutdown_waits(monkeypatch, controller, trigger_update)
