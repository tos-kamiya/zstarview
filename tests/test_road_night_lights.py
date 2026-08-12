from __future__ import annotations

import threading
from concurrent.futures import Future
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

import numpy as np
import pytest

from zstarview.gui.road_night_lights_controller import RoadNightLightsController
from zstarview.render.terrain import _road_distance_attenuation
from zstarview.road_night_lights import (
    RoadNightLightCacheSnapshot,
    ROAD_NIGHT_LIGHT_LAMP_HEIGHT_M,
    ROAD_NIGHT_LIGHT_POINT_SPACING_M,
    ROAD_NIGHT_LIGHT_STROKE_HEIGHT_M,
    RoadNightLightWay,
    build_road_night_light_ground_sampler,
    build_road_night_lights_query,
    clip_road_night_light_way_to_annulus,
    is_road_night_light_lamp_enabled,
    load_or_fetch_road_night_lights,
    load_or_fetch_road_night_lights_with_source,
    load_road_night_lights_cache,
    parse_road_night_lights_payload,
    project_road_night_lights,
    road_night_lights_cache_is_recent,
    road_night_light_lamp_strength_factor,
    road_night_lights_cache_path,
    road_night_lights_scope_key,
    save_road_night_lights_cache,
    simplify_road_night_light_way_for_observer,
)
from zstarview.types import ViewerData


def test_road_lamp_spacing_reduces_point_density_by_half() -> None:
    assert ROAD_NIGHT_LIGHT_POINT_SPACING_M == 240.0


def test_build_ground_sampler_uses_vectorized_dem_grid(monkeypatch, tmp_path) -> None:
    calls: dict[str, object] = {}

    class FakeGrid:
        def sample_lonlat(self, lon_deg, lat_deg, *, method):
            calls["sample"] = (lon_deg.tolist(), lat_deg.tolist(), method)
            return np.asarray([123.0, 124.0], dtype=np.float64)

    class FakeDem:
        def __init__(self, paths, *, default_elevation_m):
            calls["dem_init"] = (paths, default_elevation_m)

        def build_grid(self, bbox):
            calls["bbox"] = bbox
            return FakeGrid()

        def close(self):
            calls["closed"] = True

    def fake_fetch(**kwargs):
        calls["fetch"] = kwargs
        return SimpleNamespace(paths=(tmp_path / "dem.tif",))

    monkeypatch.setattr("zstarview.road_night_lights.GeoTiffDem", FakeDem)
    monkeypatch.setattr("zstarview.road_night_lights.fetch_copernicus_dem", fake_fetch)

    sampler = build_road_night_light_ground_sampler(
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        cache_dir=tmp_path,
    )

    assert sampler([35.0, 35.1], [139.0, 139.1]).tolist() == [123.0, 124.0]
    assert calls["sample"] == ([139.0, 139.1], [35.0, 35.1], "bilinear")
    assert calls["closed"] is True
    assert calls["fetch"]["cache_dir"] == tmp_path


def test_build_ground_sampler_propagates_dem_failure(monkeypatch, tmp_path) -> None:
    def fail_fetch(**_kwargs):
        raise RuntimeError("DEM unavailable")

    monkeypatch.setattr("zstarview.road_night_lights.fetch_copernicus_dem", fail_fetch)

    with pytest.raises(RuntimeError, match="DEM unavailable"):
        build_road_night_light_ground_sampler(
            observer_lat_deg=35.0,
            observer_lon_deg=139.0,
            cache_dir=tmp_path,
        )


def test_project_road_lights_uses_distinct_dem_relative_heights(monkeypatch) -> None:
    projection_calls: list[dict[str, object]] = []

    def fake_project(**kwargs):
        projection_calls.append(kwargs)
        return tuple(
            SimpleNamespace(alt_deg=1.0, az_deg=2.0, distance_km=3.0)
            for _ in kwargs["target_latitude_deg"]
        )

    monkeypatch.setattr(
        "zstarview.road_night_lights.project_place_targets_to_altaz", fake_project
    )

    def ground_sampler(latitudes, longitudes):
        assert len(latitudes) == len(longitudes)
        return [100.0] * len(latitudes)

    polylines = project_road_night_lights(
        (RoadNightLightWay(1, "primary", ((139.0, 35.0), (139.002, 35.0))),),
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        observer_height_m=1.7,
        ground_elevation_m_sampler=ground_sampler,
    )

    assert len(polylines) == 1
    assert len(projection_calls) == 2
    stroke_call, lamp_call = projection_calls
    assert stroke_call["observer_height_m"] == 101.7
    assert lamp_call["observer_height_m"] == 101.7
    assert set(stroke_call["target_height_m"]) == {
        100.0 + ROAD_NIGHT_LIGHT_STROKE_HEIGHT_M
    }
    assert set(lamp_call["target_height_m"]) == {
        100.0 + ROAD_NIGHT_LIGHT_LAMP_HEIGHT_M
    }


def test_controller_builds_dem_sampler_before_projection(monkeypatch) -> None:
    way = RoadNightLightWay(1, "primary", ((139.0, 35.0), (139.01, 35.0)))
    snapshot = RoadNightLightCacheSnapshot((way,))
    sampler = object()
    projection_calls: list[dict[str, object]] = []

    monkeypatch.setattr(
        "zstarview.gui.road_night_lights_controller.load_or_fetch_road_night_lights_with_source",
        lambda **_kwargs: (snapshot, True),
    )
    monkeypatch.setattr(
        "zstarview.gui.road_night_lights_controller.simplify_road_night_light_way_for_observer",
        lambda source_way, **_kwargs: source_way,
    )
    monkeypatch.setattr(
        "zstarview.gui.road_night_lights_controller.clip_road_night_lights_to_annulus",
        lambda ways, **_kwargs: ways,
    )
    monkeypatch.setattr(
        "zstarview.gui.road_night_lights_controller.build_road_night_light_ground_sampler",
        lambda **_kwargs: sampler,
    )

    def fake_project(ways, **kwargs):
        projection_calls.append({"ways": ways, **kwargs})
        return ()

    monkeypatch.setattr(
        "zstarview.gui.road_night_lights_controller.project_road_night_lights",
        fake_project,
    )

    controller = RoadNightLightsController()
    payloads: list[dict[str, object]] = []
    controller.road_ready.connect(payloads.append)
    controller._run(
        ViewerData(location=(35.0, 139.0), timezone_name="UTC", city_name="Test")
    )

    assert projection_calls[0]["ground_elevation_m_sampler"] is sampler
    assert payloads == [{"polylines": [], "source": "Road: cache"}]


def test_controller_allows_same_location_retry_after_failure(monkeypatch) -> None:
    controller = RoadNightLightsController()
    viewer = ViewerData(
        location=(35.0, 139.0), timezone_name="UTC", city_name="Test"
    )
    failed_future: Future[None] = Future()

    def fail_run(_viewer_data: ViewerData) -> None:
        raise RuntimeError("temporary failure")

    monkeypatch.setattr(controller, "_run", fail_run)
    controller._run_in_thread(failed_future, viewer, (35.0, 139.0))

    with pytest.raises(RuntimeError, match="temporary failure"):
        failed_future.result()
    assert controller._key is None

    monkeypatch.setattr(threading.Thread, "start", lambda _thread: None)
    assert controller.update(viewer_data=viewer, reason="retry") is True


def test_build_query_fetches_all_supported_types_once() -> None:
    query = build_road_night_lights_query(
        observer_lat_deg=35.6580,
        observer_lon_deg=139.7016,
    )
    assert "around:10000.0,35.65800000,139.70160000" in query
    assert '"highway"~"^(motorway|trunk|primary|secondary|tertiary)$"' in query
    assert query.count('"highway"') == 1


def test_parse_payload_keeps_supported_ways_and_geometry() -> None:
    ways = parse_road_night_lights_payload(
        {
            "elements": [
                {
                    "type": "way",
                    "id": 1,
                    "tags": {"highway": "primary"},
                    "geometry": [{"lon": 0.0, "lat": 0.0}, {"lon": 0.1, "lat": 0.0}],
                },
                {
                    "type": "way",
                    "id": 2,
                    "tags": {"highway": "residential"},
                    "geometry": [{"lon": 0.0, "lat": 0.0}, {"lon": 0.1, "lat": 0.0}],
                },
            ]
        }
    )
    assert ways == (RoadNightLightWay(1, "primary", ((0.0, 0.0), (0.1, 0.0))),)


def test_cache_round_trip(tmp_path) -> None:
    key = road_night_lights_scope_key(
        observer_lat_deg=35.6580, observer_lon_deg=139.7016, radius_km=10.0
    )
    snapshot = RoadNightLightCacheSnapshot(
        ways=(RoadNightLightWay(1, "trunk", ((139.70, 35.65), (139.71, 35.66))),),
        fetched_at_utc=datetime(2026, 8, 9, tzinfo=timezone.utc),
    )
    save_road_night_lights_cache(key, snapshot, cache_root=tmp_path)
    assert road_night_lights_cache_path(key, cache_root=tmp_path).exists()
    assert load_road_night_lights_cache(key, cache_root=tmp_path) == snapshot


def test_clip_way_excludes_inner_half_km_and_keeps_outer_ten_km() -> None:
    way = RoadNightLightWay(1, "primary", ((0.0, 0.0), (0.12, 0.0)))
    fragments = clip_road_night_light_way_to_annulus(
        way,
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
    )
    assert len(fragments) == 1
    assert fragments[0].coordinates_lonlat[0][0] > 0.004
    assert fragments[0].coordinates_lonlat[-1][0] < 0.1


def test_simplify_way_removes_sub_grid_vertices() -> None:
    way = RoadNightLightWay(
        1,
        "primary",
        ((0.01, 0.0), (0.01000001, 0.0), (0.01000002, 0.0), (0.011, 0.0)),
    )
    simplified = simplify_road_night_light_way_for_observer(
        way,
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
    )
    assert len(simplified.coordinates_lonlat) < len(way.coordinates_lonlat)
    assert len(simplified.coordinates_lonlat) >= 2


def test_load_or_fetch_reuses_existing_cache(monkeypatch, tmp_path) -> None:
    key = road_night_lights_scope_key(
        observer_lat_deg=35.6580, observer_lon_deg=139.7016, radius_km=10.0
    )
    snapshot = RoadNightLightCacheSnapshot(
        ways=(RoadNightLightWay(1, "primary", ((139.7, 35.65), (139.71, 35.66))),),
        fetched_at_utc=datetime.now(timezone.utc),
    )
    save_road_night_lights_cache(key, snapshot, cache_root=tmp_path)
    monkeypatch.setattr(
        "zstarview.road_night_lights.fetch_road_night_lights",
        lambda **_kwargs: (_ for _ in ()).throw(AssertionError("cache miss")),
    )
    assert (
        load_or_fetch_road_night_lights(
            observer_lat_deg=35.6580,
            observer_lon_deg=139.7016,
            cache_root=tmp_path,
        )
        == snapshot
    )


def test_load_or_fetch_falls_back_to_five_km_after_timeout(monkeypatch, tmp_path) -> None:
    snapshot = RoadNightLightCacheSnapshot(
        ways=(RoadNightLightWay(1, "primary", ((139.7, 35.65), (139.71, 35.66))),)
    )
    calls: list[float] = []

    def fetch(**kwargs):
        calls.append(float(kwargs["radius_km"]))
        if len(calls) == 1:
            raise RuntimeError("HTTP 504")
        return snapshot

    monkeypatch.setattr("zstarview.road_night_lights.fetch_road_night_lights", fetch)

    result, cache_hit = load_or_fetch_road_night_lights_with_source(
        observer_lat_deg=35.6580,
        observer_lon_deg=139.7016,
        cache_root=tmp_path,
    )

    assert result == snapshot
    assert cache_hit is False
    assert calls == [10.0, 5.0]
    assert load_road_night_lights_cache(
        road_night_lights_scope_key(
            observer_lat_deg=35.6580,
            observer_lon_deg=139.7016,
            radius_km=5.0,
        ),
        cache_root=tmp_path,
    ) == snapshot


def test_load_or_fetch_returns_fresh_empty_snapshot(monkeypatch, tmp_path) -> None:
    snapshot = RoadNightLightCacheSnapshot(ways=())
    monkeypatch.setattr(
        "zstarview.road_night_lights.fetch_road_night_lights",
        lambda **_kwargs: snapshot,
    )

    result, cache_hit = load_or_fetch_road_night_lights_with_source(
        observer_lat_deg=24.643738,
        observer_lon_deg=110.612524,
        cache_root=tmp_path,
    )

    assert result == snapshot
    assert cache_hit is False


def test_road_cache_recent_and_expired_boundaries() -> None:
    now = datetime.now(timezone.utc)
    recent = RoadNightLightCacheSnapshot(
        ways=(), fetched_at_utc=now - timedelta(days=29, hours=23)
    )
    expired = RoadNightLightCacheSnapshot(
        ways=(), fetched_at_utc=now - timedelta(days=30, seconds=1)
    )

    assert road_night_lights_cache_is_recent(recent, now_utc=now) is True
    assert road_night_lights_cache_is_recent(expired, now_utc=now) is False
    assert (
        road_night_lights_cache_is_recent(
            RoadNightLightCacheSnapshot(ways=()), now_utc=now
        )
        is False
    )


def test_load_or_fetch_uses_stale_cache_when_refreshes_fail(monkeypatch, tmp_path) -> None:
    stale = RoadNightLightCacheSnapshot(
        ways=(RoadNightLightWay(1, "primary", ((139.7, 35.65), (139.71, 35.66))),),
        fetched_at_utc=datetime.now(timezone.utc) - timedelta(days=31),
    )
    key = road_night_lights_scope_key(
        observer_lat_deg=35.6580, observer_lon_deg=139.7016, radius_km=10.0
    )
    save_road_night_lights_cache(key, stale, cache_root=tmp_path)
    calls: list[float] = []

    def fetch(**kwargs):
        calls.append(float(kwargs["radius_km"]))
        raise RuntimeError("HTTP 504")

    monkeypatch.setattr("zstarview.road_night_lights.fetch_road_night_lights", fetch)

    result, cache_hit = load_or_fetch_road_night_lights_with_source(
        observer_lat_deg=35.6580,
        observer_lon_deg=139.7016,
        cache_root=tmp_path,
    )

    assert result == stale
    assert cache_hit is True
    assert calls == [10.0, 5.0]


def test_road_distance_attenuation_fades_smoothly() -> None:
    near = _road_distance_attenuation(0.5)
    middle = _road_distance_attenuation(5.25)
    far = _road_distance_attenuation(10.0)

    assert near == 1.0
    assert near > middle > far
    assert far == 0.2


def test_road_lamps_remain_visible_until_sun_altitude_reaches_minus_one() -> None:
    assert road_night_light_lamp_strength_factor(0.0) == 0.0
    assert road_night_light_lamp_strength_factor(-2.0) == 0.5
    assert road_night_light_lamp_strength_factor(-4.0) == 1.0
    assert is_road_night_light_lamp_enabled(-0.1) is True
    assert is_road_night_light_lamp_enabled(0.0) is False
