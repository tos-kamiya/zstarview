from __future__ import annotations

from datetime import datetime, timezone

from zstarview.render.terrain import _road_distance_attenuation
from zstarview.road_night_lights import (
    RoadNightLightCacheSnapshot,
    RoadNightLightWay,
    build_road_night_lights_query,
    clip_road_night_light_way_to_annulus,
    is_road_night_light_lamp_enabled,
    load_or_fetch_road_night_lights,
    load_or_fetch_road_night_lights_with_source,
    load_road_night_lights_cache,
    parse_road_night_lights_payload,
    road_night_light_lamp_strength_factor,
    road_night_lights_cache_path,
    road_night_lights_scope_key,
    save_road_night_lights_cache,
    simplify_road_night_light_way_for_observer,
)


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
        ways=(RoadNightLightWay(1, "primary", ((139.7, 35.65), (139.71, 35.66))),)
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
