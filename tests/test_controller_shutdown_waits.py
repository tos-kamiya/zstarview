from __future__ import annotations

import threading
import time
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
from urllib.error import URLError

import astropy.time
import numpy as np

from zstarview.clouddisc.types import DownloadCancelledError
from zstarview.gui.aircraft_controller import AircraftController
from zstarview.gui.jpl_small_body_controller import JplSmallBodyController
from zstarview.gui.satellite_controller import SatelliteController
from zstarview.gui.sky_worker import SkyDataWorker
from zstarview.gui.terrain_controller import TerrainHorizonController
from zstarview.gui.urban_outline_controller import UrbanOutlineController
from zstarview.gui.water_overlay_controller import WaterOverlayController
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


def test_aircraft_controller_timeout_is_logged_without_traceback(monkeypatch, caplog) -> None:
    controller = AircraftController()
    fetch_started = threading.Event()

    def fake_fetcher(_bbox):
        fetch_started.set()
        raise URLError(TimeoutError("timed out"))

    monkeypatch.setattr(controller, "_fetcher", fake_fetcher)

    with caplog.at_level("WARNING", logger="zstarview.gui.aircraft_controller"):
        controller.update(
            observer_lat=35.0,
            observer_lon=139.0,
            observer_height_m=1.7,
            time_obj=astropy.time.Time(datetime(2026, 4, 25, 0, 0, tzinfo=timezone.utc)),
            reason="manual",
        )

        assert fetch_started.wait(timeout=1.0)
        deadline = time.time() + 1.0
        while controller._active_workers and time.time() < deadline:  # noqa: SLF001
            time.sleep(0.01)

    assert "Aircraft update failed: <urlopen error timed out>" in caplog.text
    assert "Traceback (most recent call last)" not in caplog.text


def test_aircraft_controller_does_not_project_in_fetch_stage() -> None:
    projector_calls: list[object] = []

    controller = AircraftController(
        fetcher=lambda bbox: [],
        projector=lambda *args, **kwargs: projector_calls.append((args, kwargs)),
    )

    controller._run_update(
        observer_lat=35.0,
        observer_lon=139.0,
        observer_height_m=1.7,
        time_obj=astropy.time.Time(datetime(2026, 4, 25, 0, 0, tzinfo=timezone.utc)),
        reason="manual",
        request_id=0,
    )

    assert projector_calls == []


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


def test_terrain_controller_shutdown_cancels_download(monkeypatch, tmp_path) -> None:
    controller = TerrainHorizonController(cache_dir=tmp_path)
    started = threading.Event()

    def fake_fetch_copernicus_dem(*, abort_event=None, **_kwargs):
        assert abort_event is not None
        started.set()
        while not abort_event.is_set():
            time.sleep(0.01)
        raise DownloadCancelledError("Cancelled while downloading Copernicus DEM tiles")

    monkeypatch.setattr("zstarview.gui.terrain_controller.fetch_copernicus_dem", fake_fetch_copernicus_dem)

    controller.update(lat=35.0, lon=139.0, observer_height_m=1.7, reason="manual")
    assert started.wait(timeout=1.0)

    done = threading.Event()

    def run_shutdown() -> None:
        controller.shutdown(wait_timeout_s=1.0)
        done.set()

    shutdown_thread = threading.Thread(target=run_shutdown)
    shutdown_thread.start()

    assert done.wait(timeout=1.0)
    shutdown_thread.join(timeout=1.0)


def test_urban_outline_controller_shutdown_cancels_download(monkeypatch, tmp_path) -> None:
    controller = UrbanOutlineController(derived_root_dir=tmp_path)
    viewer_data = ViewerData(location=(35.0, 139.0), timezone_name="Asia/Tokyo", city_name="Tokyo")
    monkeypatch.setattr(controller, "_required_derived_dirs", lambda _viewer_data: (("building", tmp_path / "normal" / "bldg"),))
    monkeypatch.setattr(controller, "_selected_skyscraper_tiles", lambda _viewer_data: ())
    monkeypatch.setattr("zstarview.gui.urban_outline_controller.resolve_overture_release_for_cache_root", lambda **_kwargs: None)
    started = threading.Event()

    def fake_import_overture_buildings(*, abort_event=None, **_kwargs):
        assert abort_event is not None
        started.set()
        while not abort_event.is_set():
            time.sleep(0.01)
        raise DownloadCancelledError("Cancelled while running overturemaps download")

    monkeypatch.setattr("zstarview.gui.urban_outline_controller.import_overture_buildings", fake_import_overture_buildings)

    controller.update(viewer_data=viewer_data, reason="manual")
    assert started.wait(timeout=1.0)

    done = threading.Event()

    def run_shutdown() -> None:
        controller.shutdown(wait_timeout_s=1.0)
        done.set()

    shutdown_thread = threading.Thread(target=run_shutdown)
    shutdown_thread.start()

    assert done.wait(timeout=1.0)
    shutdown_thread.join(timeout=1.0)


def test_water_overlay_controller_shutdown_cancels_download(monkeypatch) -> None:
    controller = WaterOverlayController()

    def trigger_update() -> None:
        controller.update(
            viewer_data=ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Tokyo"),
            observer_ground_m=0.0,
            use_dem_ground=False,
            reason="manual",
        )

    _assert_shutdown_waits(monkeypatch, controller, trigger_update)


def test_sky_worker_shutdown_waits(monkeypatch) -> None:
    controller = SkyDataWorker()

    def trigger_update() -> None:
        controller.update(
            ephemeris=object(),
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
