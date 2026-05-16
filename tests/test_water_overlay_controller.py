from __future__ import annotations

import zstarview.gui.water_overlay_controller as mod
from zstarview.gui.water_overlay_controller import WaterOverlayController
from zstarview.water_mask_interface import WaterSurfaceBandStats
from zstarview.water_overlay import WaterOverlayPoint


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
    controller._active_key = (35.0, 139.0, 1.7, 0.0, False)  # noqa: SLF001

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

    assert failures == [{"banner": "Water: HTTP 504"}]
    assert "Water surface update failed: HTTP 504" in caplog.text
    assert "Traceback" not in caplog.text


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

    scope_cache = mod._WaterOverlayScopeCache(  # noqa: SLF001
        footprints=(),
        fetched_at_utc=None,
    )

    with caplog.at_level("INFO", logger="zstarview.gui.water_overlay_controller"):
        controller._run_update(  # noqa: SLF001
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
    assert "Water band stats: 125m tiles=0 raw=0 collapsed=0 visible=0" in caplog.text
    assert "Water band stats: 250m tiles=1 raw=12 collapsed=1 visible=1" in caplog.text
    assert "Water band stats: 500m tiles=0 raw=0 collapsed=0 visible=0" in caplog.text
    assert "Water mask points: 1 visible, nearest sea point 0.500 km, bands: 125m=0 250m=1 500m=0" in caplog.text


def test_build_requested_variants_combines_sea_and_inland_points(monkeypatch) -> None:
    controller = WaterOverlayController()
    footprint = object()
    scope_cache = mod._WaterOverlayScopeCache(  # noqa: SLF001
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

    active_points, sea_points, dem_points, band_stats, footprints = controller._build_requested_variants(  # noqa: SLF001
        scope_cache,
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        observer_height_m=1.7,
        observer_ground_m=12.0,
        use_dem_ground=False,
        scan_radius_km=2.0,
        target_ground_sampler=None,
    )

    assert len(active_points) == 3
    assert [point.water_category for point in sea_points] == ["sea-125", "river", "lake"]
    assert dem_points is None
    assert band_stats == (WaterSurfaceBandStats("125m", 1, 1, 1, 1),)
    assert footprints == (footprint,)
