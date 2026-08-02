from __future__ import annotations

from datetime import datetime, timezone
from types import SimpleNamespace

import zstarview.gui.water_overlay_controller as mod
from zstarview.gui.water_overlay_cache import WaterOverlayCacheSnapshot
from zstarview.gui.water_overlay_controller import WaterOverlayController
from zstarview.water_mask_interface import WaterSurfaceBandStats
from zstarview.water_overlay import (
    DEFAULT_WATER_AZIMUTH_STEP_DEG,
    WaterOverlayPoint,
    WaterPolygonFootprint,
)


def test_water_overlay_controller_uses_compact_failure_banner_and_log(
    monkeypatch,
    caplog,
) -> None:
    failures: list[dict[str, object]] = []
    controller = WaterOverlayController()
    controller.water_failed.connect(failures.append)

    monkeypatch.setattr(
        mod,
        "sample_water_surface_interface_points_with_stats",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("HTTP 504")),
    )
    monkeypatch.setattr(mod, "load_water_overlay_cache_for_location", lambda **kwargs: None)
    monkeypatch.setattr(mod, "load_water_overlay_cache", lambda *_args, **_kwargs: None)
    controller._active_key = (35.0, 139.0, 1.7, 0.0, False)

    with caplog.at_level("WARNING", logger="zstarview.gui.water_overlay_controller"):
        controller._run_update(
            lat_deg=35.0,
            lon_deg=139.0,
            observer_height_m=1.7,
            observer_ground_m=0.0,
            use_dem_ground=False,
            reason="manual",
            key=(35.0, 139.0, 1.7, 0.0, False),
            scope_key="scope",
            scan_radius_km=2.0,
            cached_scope=None,
        )

    assert failures == [{"banner": "Water: unavailable"}]
    assert "Water surface update failed: HTTP 504" in caplog.text
    assert "Traceback" not in caplog.text


def test_water_overlay_controller_uses_recent_disk_snapshot(monkeypatch) -> None:
    controller = WaterOverlayController()
    footprint = WaterPolygonFootprint(
        water_id="water-1",
        kind="water_polygon",
        outer_rings_lonlat=(((139.0, 35.0), (139.1, 35.0), (139.1, 35.1), (139.0, 35.0)),),
        inner_rings_lonlat=(),
        source="test",
        tags={},
    )
    snapshot = WaterOverlayCacheSnapshot(
        footprints=(footprint,),
        water_polygon_count=1,
        fetched_at_utc=datetime.now(timezone.utc),
    )

    monkeypatch.setattr(mod, "load_water_overlay_cache", lambda *_args, **_kwargs: snapshot)
    monkeypatch.setattr(
        mod,
        "fetch_overpass_json",
        lambda **_kwargs: (_ for _ in ()).throw(AssertionError("should not fetch overpass")),
    )

    scope_cache = controller._ensure_scope_cache(
        scope_key="scope",
        lat_deg=35.0,
        lon_deg=139.0,
        scan_radius_km=2.0,
        cached_scope=None,
        now_utc=datetime.now(timezone.utc),
    )

    assert scope_cache.footprints == snapshot.footprints
    assert scope_cache.fetched_at_utc == snapshot.fetched_at_utc


def test_water_overlay_controller_saves_fresh_disk_snapshot(monkeypatch) -> None:
    controller = WaterOverlayController()
    footprint = WaterPolygonFootprint(
        water_id="water-2",
        kind="water_polygon",
        outer_rings_lonlat=(((139.0, 35.0), (139.2, 35.0), (139.2, 35.2), (139.0, 35.0)),),
        inner_rings_lonlat=(),
        source="test",
        tags={},
    )
    saved: dict[str, object] = {}

    monkeypatch.setattr(mod, "load_water_overlay_cache", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        mod,
        "fetch_overpass_json",
        lambda **_kwargs: {"elements": []},
    )
    monkeypatch.setattr(
        mod,
        "extract_water_polygons",
        lambda *_args, **_kwargs: (footprint,),
    )
    monkeypatch.setattr(
        mod,
        "save_water_overlay_cache",
        lambda scope_key, snapshot, **_kwargs: saved.setdefault(
            "payload",
            (scope_key, snapshot),
        ),
    )

    scope_cache = controller._ensure_scope_cache(
        scope_key="scope",
        lat_deg=35.0,
        lon_deg=139.0,
        scan_radius_km=2.0,
        cached_scope=None,
        now_utc=datetime.now(timezone.utc),
    )

    from zstarview.water_overlay import simplify_water_footprints_for_observer

    expected = simplify_water_footprints_for_observer(
        (footprint,),
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
    )

    assert scope_cache.footprints == expected
    assert saved["payload"][0] == "scope"
    assert saved["payload"][1].footprints == expected


def test_water_overlay_controller_uses_other_radius_cache_after_fetch_failure(
    monkeypatch,
) -> None:
    controller = WaterOverlayController()
    footprint = WaterPolygonFootprint(
        water_id="cached-river",
        kind="natural_water",
        outer_rings_lonlat=(((139.0, 35.0), (139.1, 35.0), (139.1, 35.1), (139.0, 35.0)),),
        inner_rings_lonlat=(),
        source="test",
        tags={},
    )
    fallback = WaterOverlayCacheSnapshot(
        footprints=(footprint,),
        water_polygon_count=1,
        fetched_at_utc=datetime.now(timezone.utc),
    )
    monkeypatch.setattr(mod, "load_water_overlay_cache", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        mod,
        "load_water_overlay_cache_for_location",
        lambda **_kwargs: fallback,
    )
    monkeypatch.setattr(
        mod,
        "fetch_overpass_json",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("HTTP 504")),
    )

    scope_cache = controller._ensure_scope_cache(
        scope_key="earth_+35.0000_+139.0000_r5.00",
        lat_deg=35.0,
        lon_deg=139.0,
        scan_radius_km=5.0,
        cached_scope=None,
        now_utc=datetime.now(timezone.utc),
    )

    assert scope_cache.footprints == fallback.footprints


def test_water_overlay_controller_fast_mode_uses_sparsest_sampling(monkeypatch) -> None:
    controller = WaterOverlayController()
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        controller,
        "_spawn_worker",
        lambda *, target, kwargs, label: captured.setdefault("kwargs", kwargs),
    )

    controller.update(
        viewer_data=SimpleNamespace(lat_deg=35.0, lon_deg=139.0, observer_height_m=1.7),
        observer_ground_m=0.0,
        use_dem_ground=False,
        reason="viewport",
        fast_mode=True,
        surface_size_px=(4096, 4096),
    )

    assert float(captured["kwargs"]["azimuth_step_deg"]) == DEFAULT_WATER_AZIMUTH_STEP_DEG


def test_water_overlay_controller_uses_recent_cached_footprints_as_is(
    monkeypatch,
) -> None:
    controller = WaterOverlayController()
    footprint = WaterPolygonFootprint(
        water_id="river",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (0.0200, 0.0000),
                (0.0201, 0.0000),
                (0.0202, 0.0000),
                (0.0203, 0.0000),
                (0.0210, 0.0000),
                (0.0210, 0.0010),
                (0.0200, 0.0010),
                (0.0200, 0.0000),
            ),
        ),
        inner_rings_lonlat=(),
        source="way",
        tags={"natural": "water", "water": "river"},
    )
    snapshot = WaterOverlayCacheSnapshot(
        footprints=(footprint,),
        water_polygon_count=1,
        fetched_at_utc=datetime.now(timezone.utc),
    )
    saved: dict[str, object] = {}

    monkeypatch.setattr(mod, "load_water_overlay_cache", lambda *_args, **_kwargs: snapshot)
    monkeypatch.setattr(
        mod,
        "save_water_overlay_cache",
        lambda scope_key, snap, **_kwargs: saved.setdefault("payload", (scope_key, snap)),
    )

    scope_cache = controller._ensure_scope_cache(
        scope_key="scope",
        lat_deg=0.0,
        lon_deg=0.0,
        scan_radius_km=2.0,
        cached_scope=None,
        now_utc=datetime.now(timezone.utc),
    )

    assert scope_cache.footprints == snapshot.footprints
    assert saved == {}


def test_water_overlay_controller_does_not_refetch_empty_cached_footprints(
    monkeypatch,
) -> None:
    controller = WaterOverlayController()
    sample_stats = (
        WaterSurfaceBandStats("125m", 0, 0, 0, 0),
        WaterSurfaceBandStats("250m", 0, 0, 0, 0),
        WaterSurfaceBandStats("500m", 0, 0, 0, 0),
    )

    monkeypatch.setattr(
        mod,
        "sample_water_surface_interface_points_with_stats",
        lambda **kwargs: ((), sample_stats),
    )
    monkeypatch.setattr(
        mod,
        "fetch_water_overlay_footprints",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("should not fetch overpass")),
    )
    monkeypatch.setattr(mod, "sample_water_overlay_points", lambda *args, **kwargs: ())

    scope_cache = mod._WaterOverlayScopeCache(
        footprints=(),
        fetched_at_utc=datetime.now(timezone.utc),
    )

    active_dots, sea_mask_dots, sea_dots, inland_dots, dem_dots, band_stats, footprints = (
        controller._build_requested_variants(
            scope_cache,
            observer_lat_deg=0.0,
            observer_lon_deg=0.0,
            observer_height_m=0.0,
            observer_ground_m=0.0,
            use_dem_ground=False,
            scan_radius_km=2.0,
            target_ground_sampler=None,
            key=(0.0, 0.0, 0.0, 0.0, False, DEFAULT_WATER_AZIMUTH_STEP_DEG),
            scope_key="scope",
        )
    )

    assert active_dots == ()
    assert sea_mask_dots == ()
    assert sea_dots == ()
    assert inland_dots == ()
    assert dem_dots is None
    assert band_stats == sample_stats
    assert footprints == ()


def test_water_overlay_controller_uses_sea_mask_points_before_sampling(monkeypatch, caplog) -> None:
    controller = WaterOverlayController()
    observed: dict[str, object] = {}

    monkeypatch.setattr(
        mod,
        "sample_water_surface_interface_points_with_stats",
        lambda **kwargs: (
            observed.setdefault("kwargs", kwargs),
            (
                (WaterOverlayPoint("water-mask", 10.0, 20.0, 0.5, water_category="sea-250"),),
                (
                    WaterSurfaceBandStats("125m", 0, 0, 0, 0),
                    WaterSurfaceBandStats("250m", 1, 12, 1, 1),
                    WaterSurfaceBandStats("500m", 0, 0, 0, 0),
                ),
            ),
        )[1],
    )
    monkeypatch.setattr(mod, "fetch_water_overlay_footprints", lambda *args, **_kwargs: ())
    monkeypatch.setattr(mod, "sample_water_overlay_points", lambda *args, **_kwargs: ())

    scope_cache = mod._WaterOverlayScopeCache(
        footprints=(),
        fetched_at_utc=None,
    )

    with caplog.at_level("INFO", logger="zstarview.gui.water_overlay_controller"):
        controller._run_update(
            lat_deg=0.0,
            lon_deg=0.0,
            observer_height_m=0.0,
            observer_ground_m=0.0,
            use_dem_ground=False,
            reason="manual",
            key=(0.0, 0.0, 0.0, 0.0, False),
            scope_key="scope",
            scan_radius_km=2.0,
            cached_scope=scope_cache,
        )

    assert float(observed["kwargs"]["observer_height_m"]) == 0.0
    assert float(observed["kwargs"]["max_distance_km"]) == 2.0
    assert "Water band stats: 125m opened_tiles=0 raw=0 collapsed=0 visible=0" in caplog.text
    assert "Water band stats: 250m opened_tiles=1 raw=12 collapsed=1 visible=1" in caplog.text
    assert "Water band stats: 500m opened_tiles=0 raw=0 collapsed=0 visible=0" in caplog.text
    assert "Water mask dots: 1 visible, nearest sea dot 0.500 km, bands: 125m=0 250m=1 500m=0" in caplog.text


def test_build_requested_variants_combines_sea_and_inland_points(monkeypatch) -> None:
    controller = WaterOverlayController()
    footprint = object()
    scope_cache = mod._WaterOverlayScopeCache(
        footprints=(),
        fetched_at_utc=None,
    )

    monkeypatch.setattr(
        mod,
        "sample_water_surface_interface_points_with_stats",
        lambda **_kwargs: (
            (WaterOverlayPoint("sea", 1.0, 10.0, 0.5, water_category="sea-125"),),
            (WaterSurfaceBandStats("125m", 1, 1, 1, 1),),
        ),
    )
    monkeypatch.setattr(
        mod,
        "fetch_water_overlay_footprints",
        lambda **_kwargs: (footprint,),
    )
    monkeypatch.setattr(
        mod,
        "sample_water_overlay_points",
        lambda footprints, **_kwargs: (
            WaterOverlayPoint("inland", 2.0, 20.0, 1.5, water_category="river"),
            WaterOverlayPoint("inland", 3.0, 30.0, 1.8, water_category="lake"),
        )
        if footprints
        else (),
    )

    (
        active_dots,
        sea_mask_dots,
        sea_dots,
        inland_dots,
        dem_dots,
        band_stats,
        footprints,
    ) = controller._build_requested_variants(
        scope_cache,
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        observer_height_m=1.7,
        observer_ground_m=12.0,
        use_dem_ground=False,
        scan_radius_km=2.0,
        target_ground_sampler=None,
        key=(35.0, 139.0, 1.7, 12.0, False),
        scope_key="scope",
    )

    assert len(active_dots) == 3
    assert [point.water_category for point in sea_mask_dots] == ["sea-125"]
    assert [point.water_category for point in sea_dots] == ["sea-125", "river", "lake"]
    assert [point.water_category for point in inland_dots] == ["river", "lake"]
    assert dem_dots is None
    assert band_stats == (WaterSurfaceBandStats("125m", 1, 1, 1, 1),)
    assert footprints == (footprint,)


def test_run_update_emits_sea_then_combined_water_layers(monkeypatch) -> None:
    controller = WaterOverlayController()
    footprint = object()
    scope_cache = mod._WaterOverlayScopeCache(
        footprints=(footprint,),
        fetched_at_utc=None,
    )
    emitted: list[dict[str, object]] = []

    monkeypatch.setattr(
        mod,
        "sample_water_surface_interface_points_with_stats",
        lambda **_kwargs: (
            (
                WaterOverlayPoint("sea", 1.0, 10.0, 0.5, water_category="sea-125"),
            ),
            (WaterSurfaceBandStats("125m", 1, 1, 1, 1),),
        ),
    )
    monkeypatch.setattr(
        mod,
        "sample_water_overlay_points",
        lambda footprints, **_kwargs: (
            WaterOverlayPoint("inland", 2.0, 20.0, 1.5, water_category="river"),
            WaterOverlayPoint("inland", 3.0, 30.0, 1.8, water_category="lake"),
        )
        if footprints
        else (),
    )
    monkeypatch.setattr(
        controller,
        "_emit_variant",
        lambda dots, **payload: emitted.append(
            {
                "dots": tuple(dots),
                **payload,
            }
        ),
    )
    monkeypatch.setattr(controller, "_store_scope_cache", lambda *args, **_kwargs: None)
    controller._active_key = (35.0, 139.0, 1.7, 12.0, False)

    controller._run_update(
        lat_deg=35.0,
        lon_deg=139.0,
        observer_height_m=1.7,
        observer_ground_m=12.0,
        use_dem_ground=False,
        reason="manual",
        key=(35.0, 139.0, 1.7, 12.0, False),
        scope_key="scope",
        scan_radius_km=2.0,
        cached_scope=scope_cache,
    )

    assert len(emitted) == 2
    assert [point.water_category for point in emitted[0]["dots"]] == ["sea-125"]
    assert emitted[0]["mode"] == "sea"
    assert emitted[0]["inland_dots"] is None
    assert [point.water_category for point in emitted[1]["dots"]] == ["sea-125", "river", "lake"]
    assert emitted[1]["mode"] == "sea"
    assert [point.water_category for point in emitted[1]["sea_dots"]] == ["sea-125"]
    assert [point.water_category for point in emitted[1]["inland_dots"]] == ["river", "lake"]


def test_run_update_does_not_refetch_dem_for_gui_mode(monkeypatch) -> None:
    controller = WaterOverlayController()
    scope_cache = mod._WaterOverlayScopeCache(
        footprints=(),
        fetched_at_utc=None,
        sea_mask_dots=(WaterOverlayPoint("sea", 1.0, 10.0, 0.5, water_category="sea-125"),),
        sea_band_stats=(WaterSurfaceBandStats("125m", 1, 1, 1, 1),),
    )
    emitted: list[dict[str, object]] = []

    monkeypatch.setattr(
        mod,
        "fetch_water_overlay_footprints",
        lambda *args, **_kwargs: (WaterPolygonFootprint(
            water_id="poly",
            kind="water_polygon",
            outer_rings_lonlat=(((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)),),
            inner_rings_lonlat=(),
            source="test",
            tags={},
        ),),
    )
    monkeypatch.setattr(
        mod,
        "sample_water_overlay_points",
        lambda footprints, **_kwargs: (
            WaterOverlayPoint("inland", 2.0, 20.0, 1.5, water_category="river"),
        )
        if footprints
        else (),
    )
    monkeypatch.setattr(
        mod,
        "sample_water_surface_interface_points_with_stats",
        lambda **_kwargs: (
            (
                WaterOverlayPoint("sea", 1.0, 10.0, 0.5, water_category="sea-125"),
            ),
            (WaterSurfaceBandStats("125m", 1, 1, 1, 1),),
        ),
    )
    monkeypatch.setattr(
        controller,
        "_emit_variant",
        lambda dots, **payload: emitted.append({"dots": tuple(dots), **payload}),
    )
    monkeypatch.setattr(controller, "_store_scope_cache", lambda *args, **_kwargs: None)
    controller._active_key = (35.0, 139.0, 1.7, 12.0, False)

    controller._run_update(
        lat_deg=35.0,
        lon_deg=139.0,
        observer_height_m=1.7,
        observer_ground_m=12.0,
        use_dem_ground=True,
        reason="manual",
        key=(35.0, 139.0, 1.7, 12.0, False),
        scope_key="scope",
        scan_radius_km=2.0,
        cached_scope=scope_cache,
    )

    assert emitted
    assert emitted[-1]["mode"] == "dem"
