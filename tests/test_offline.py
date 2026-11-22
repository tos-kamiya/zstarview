from __future__ import annotations

import argparse
import logging
import socket
import urllib.error
import urllib.request
from collections.abc import Iterator

import pytest

from zstarview import offline
from zstarview.search.resolver import resolve_search_targets


@pytest.fixture
def offline_env(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    offline.install_network_guard()
    try:
        yield
    finally:
        offline.uninstall_network_guard()


def test_is_offline_follows_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(offline.OFFLINE_ENV_VAR, raising=False)
    assert offline.is_offline() is False
    offline.ensure_online("x")
    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    assert offline.is_offline() is True
    with pytest.raises(offline.OfflineError):
        offline.ensure_online("x")


def test_offline_error_is_an_oserror() -> None:
    assert isinstance(offline.OfflineError("example.com"), OSError)


@pytest.mark.usefixtures("offline_env")
def test_guard_blocks_remote_name_resolution_and_urlopen() -> None:
    with pytest.raises(offline.OfflineError):
        socket.getaddrinfo("example.com", 443)
    with pytest.raises(OSError):
        urllib.request.urlopen("https://example.com/", timeout=1)


@pytest.mark.usefixtures("offline_env")
def test_guard_blocks_direct_ip_connect() -> None:
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        with pytest.raises(offline.OfflineError):
            sock.connect(("93.184.216.34", 80))
    finally:
        sock.close()


@pytest.mark.usefixtures("offline_env")
def test_guard_allows_loopback() -> None:
    assert socket.getaddrinfo("localhost", 80)
    assert socket.getaddrinfo("127.0.0.1", 80)
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind(("127.0.0.1", 0))
    server.listen(1)
    client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        client.connect(server.getsockname())
    finally:
        client.close()
        server.close()


def test_uninstall_restores_socket_functions() -> None:
    original = socket.getaddrinfo
    offline.install_network_guard()
    assert socket.getaddrinfo is not original
    offline.uninstall_network_guard()
    assert socket.getaddrinfo is original


def test_apply_offline_overrides_turns_off_short_lived_layers() -> None:
    args = argparse.Namespace(
        aircraft_opacity=0.4,
        satellite_opacity=0.7,
        precipitation_opacity=0.0,
        cloud_opacity=0.175,
        tropical_cyclone_opacity=0.5,
        meteor_trails_opacity=0.5,
        terrain_horizon_opacity=0.3,
        night_light_opacity=0.22,
        moon_style="image",
    )

    disabled = offline.apply_offline_overrides(args)

    for name in offline.SHORT_LIVED_OPACITY_ARGS:
        assert getattr(args, name) == 0.0
    assert args.terrain_horizon_opacity == 0.3
    assert args.night_light_opacity == 0.22
    assert args.moon_style == "sphere"
    assert "precipitation_opacity" not in disabled
    assert "aircraft_opacity" in disabled
    assert "moon_style=image" in disabled


def test_offline_cli_flag_is_available_in_both_parsers() -> None:
    from zstarview.cli import args as cli_args

    assert cli_args.parse_args(["--offline"]).offline is True
    assert cli_args.parse_args([]).offline is False
    export_args = cli_args.build_export_image_argument_parser().parse_args(
        ["--offline", "-o", "x.png"]
    )
    assert export_args.offline is True


def test_search_skips_online_callbacks_offline(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")

    def forbidden(_query: str) -> list[object]:
        raise AssertionError("online search must not run offline")

    resolution = resolve_search_targets(
        "ISS",
        (),
        satellite_search_callback=forbidden,  # type: ignore[arg-type]
        jpl_search_callback=forbidden,  # type: ignore[arg-type]
    )

    assert resolution.candidates == ()
    assert "offline" in resolution.status


def test_place_search_uses_cache_without_request_offline(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    from datetime import datetime, timezone

    from zstarview.location_resolver import nominatim, place_search
    from zstarview.location_resolver.place_cache import (
        build_place_cache_key,
        save_place_cache,
    )

    key = build_place_cache_key("Matsue Station", None, "en", limit=5)
    save_place_cache(
        key,
        (
            {
                "name": "Matsue Station",
                "display_name": "Matsue Station, Matsue, Japan",
                "lat": "35.4643",
                "lon": "133.0632",
                "category": "railway",
                "type": "station",
                "importance": 0.5,
                "osm_id": 1,
                "osm_type": "node",
            },
        ),
        original_query="Matsue Station",
        fetched_at_utc=datetime(2026, 1, 1, tzinfo=timezone.utc),
        cache_root=tmp_path,
    )
    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    monkeypatch.setattr(
        nominatim,
        "_fetch",
        lambda *_a, **_k: pytest.fail("Nominatim must not be called offline"),
    )

    candidates = place_search.search_place_candidates(
        "Matsue Station", cache_root=tmp_path
    )

    assert candidates
    assert candidates[0].cache_fetched_at_utc is not None
    assert not (tmp_path / place_search.RATE_LIMIT_FILE_NAME).exists()


def test_place_search_without_cache_fails_offline(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    from zstarview.location_resolver import place_search

    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")

    with pytest.raises(place_search.PlaceSearchNetworkError):
        place_search.search_place_candidates("Nowhere", cache_root=tmp_path)


def test_omitted_location_skips_ip_lookup_offline(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from zstarview.location_resolver import resolve as resolver

    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    monkeypatch.setattr(resolver, "_resolve_ground_elevation_m", lambda **_k: 0.0)
    monkeypatch.setattr(
        resolver,
        "fetch_location_by_ip",
        lambda: pytest.fail("ip-api must not be called offline"),
    )

    location = resolver.resolve_launch_location(
        None,
        load_last_city_func=lambda: "JP/Osaka",
        save_last_city_func=lambda _value: None,
    )

    assert "Osaka" in location.display_name


def test_overture_release_check_reuses_last_release_offline(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    from datetime import datetime, timezone

    from zstarview.data import import_overture_buildings as overture

    overture.write_overture_release_check_metadata(
        payload={
            "last_release_check_utc": "2020-01-01T00:00:00+00:00",
            "last_seen_overture_release": "2025-01-01.0",
        },
        cache_root_dir=tmp_path,
    )
    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    monkeypatch.setattr(
        overture,
        "fetch_latest_overture_release",
        lambda: pytest.fail("release catalog must not be fetched offline"),
    )

    release = overture.resolve_overture_release_for_cache_root(
        cache_root_dir=tmp_path, now_utc=datetime(2026, 1, 1, tzinfo=timezone.utc)
    )

    assert release == "2025-01-01.0"
    metadata = overture.read_overture_release_check_metadata(cache_root_dir=tmp_path)
    assert metadata is not None
    assert metadata["last_release_check_utc"] == "2020-01-01T00:00:00+00:00"


def test_overturemaps_subprocess_is_refused_offline(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from zstarview.data import import_overture_buildings as overture

    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")

    with pytest.raises(offline.OfflineError):
        overture._run_download_command(["overturemaps", "download"])


def test_hover_controllers_do_not_fetch_offline(monkeypatch: pytest.MonkeyPatch) -> None:
    from datetime import datetime, timezone

    from zstarview.gui.moon_hover_controller import MoonHoverController
    from zstarview.gui.solar_hover_controller import SolarHoverController

    class _Services:
        def submit(self, *_args: object, **_kwargs: object) -> object:
            raise AssertionError("hover images must not be fetched offline")

    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    now = datetime(2026, 1, 1, tzinfo=timezone.utc)

    assert SolarHoverController(_Services()).request(now) is None  # type: ignore[arg-type]
    assert MoonHoverController(_Services()).request(now) is None  # type: ignore[arg-type]


def test_urlerror_from_guard_is_oserror() -> None:
    assert issubclass(urllib.error.URLError, OSError)


def test_is_offline_failure_follows_exception_chain(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    try:
        try:
            raise urllib.error.URLError(offline.OfflineError("example.com"))
        except urllib.error.URLError as exc:
            raise RuntimeError("road data request failed") from exc
    except RuntimeError as wrapped:
        assert offline.is_offline_failure(wrapped) is True
        assert offline.failure_log_level(wrapped) == logging.INFO
    assert offline.is_offline_failure(RuntimeError("other")) is False
    assert offline.failure_log_level(RuntimeError("other")) == logging.WARNING
    monkeypatch.delenv(offline.OFFLINE_ENV_VAR)
    assert offline.is_offline_failure(offline.OfflineError("x")) is False


def test_road_lights_offline_use_stale_cache_without_fetch(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    from datetime import datetime, timedelta, timezone

    from zstarview import road_night_lights as road

    stale = road.RoadNightLightCacheSnapshot(
        ways=(road.RoadNightLightWay(1, "primary", ((139.7, 35.65), (139.71, 35.66))),),
        fetched_at_utc=datetime.now(timezone.utc) - timedelta(days=90),
    )
    key = road.road_night_lights_scope_key(
        observer_lat_deg=35.6580, observer_lon_deg=139.7016, radius_km=5.0
    )
    road.save_road_night_lights_cache(key, stale, cache_root=tmp_path)
    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    monkeypatch.setattr(
        road,
        "fetch_road_night_lights",
        lambda **_k: pytest.fail("road data must not be fetched offline"),
    )

    result, cache_hit = road.load_or_fetch_road_night_lights_with_source(
        observer_lat_deg=35.6580, observer_lon_deg=139.7016, cache_root=tmp_path
    )

    assert result == stale
    assert cache_hit is True


def test_road_lights_offline_without_cache_raise_cache_miss(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    from zstarview import road_night_lights as road

    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    monkeypatch.setattr(
        road,
        "fetch_road_night_lights",
        lambda **_k: pytest.fail("road data must not be fetched offline"),
    )

    with pytest.raises(offline.OfflineCacheMissError):
        road.load_or_fetch_road_night_lights_with_source(
            observer_lat_deg=10.0, observer_lon_deg=10.0, cache_root=tmp_path
        )


def test_export_layer_skipped_offline_is_not_a_failure(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    import threading

    from zstarview.cli import export_image_layers as layers

    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    done = threading.Event()
    done.set()
    thread = threading.Thread(target=lambda: None)
    thread.start()
    thread.join()
    failures: list[str] = []

    with caplog.at_level(logging.INFO):
        result = layers._await_background_task_result(
            label="urban",
            thread=thread,
            done=done,
            state={"error": offline.OfflineCacheMissError("urban outline data")},
            deadline=None,
            layer_failures=failures,
            allow_partial_data=False,
        )

    assert result is None
    assert failures == []
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_terrain_controller_reports_not_cached_offline(
    monkeypatch: pytest.MonkeyPatch, tmp_path, caplog: pytest.LogCaptureFixture
) -> None:
    from zstarview.gui import terrain_controller as terrain_controller_module

    def _raise_no_tiles(**_kwargs: object) -> None:
        raise RuntimeError(
            "No Copernicus DEM tiles were downloaded for the requested area."
        )

    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    monkeypatch.setattr(
        terrain_controller_module, "fetch_copernicus_dem", _raise_no_tiles
    )
    controller = terrain_controller_module.TerrainHorizonController(cache_dir=tmp_path)
    failed: list[object] = []
    controller.terrain_failed.connect(failed.append)

    with caplog.at_level(logging.INFO):
        controller._run_update(lat=-25.0, lon=134.0, observer_height_m=1.7, reason="initial")

    assert failed == [{"banner": "Terrain horizon: not cached"}]
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_export_water_footprints_offline_use_cache_only(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from types import SimpleNamespace

    from zstarview.cli import export_image_layers as layers

    monkeypatch.setenv(offline.OFFLINE_ENV_VAR, "1")
    monkeypatch.setattr(layers, "load_water_overlay_cache", lambda *_a, **_k: None)
    monkeypatch.setattr(
        layers,
        "host",
        lambda: pytest.fail("water polygons must not be fetched offline"),
    )

    footprints = layers._load_or_fetch_water_overlay_footprints(
        viewer_data=SimpleNamespace(lat_deg=10.0, lon_deg=10.0),  # type: ignore[arg-type]
        scan_radius_km=32.0,
        deadline=None,
    )

    assert footprints == ()
