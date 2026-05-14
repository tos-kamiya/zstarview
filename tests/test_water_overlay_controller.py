from __future__ import annotations

from datetime import datetime, timezone

import zstarview.gui.water_overlay_controller as mod
from zstarview.gui.water_overlay_controller import WaterOverlayController
from zstarview.water_overlay import WaterPolygonFootprint


def test_water_overlay_controller_uses_compact_failure_banner_and_log(
    monkeypatch,
    caplog,
) -> None:
    failures: list[dict[str, object]] = []
    controller = WaterOverlayController()
    controller.water_failed.connect(failures.append)

    monkeypatch.setattr(controller, "_load_scope_snapshot", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        controller,
        "_select_cached_variant",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        "zstarview.gui.water_overlay_controller.fetch_overpass_json",
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


def test_water_overlay_controller_simplifies_footprints_before_sampling(monkeypatch) -> None:
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
    simplified = WaterPolygonFootprint(
        water_id="river",
        kind="natural_water",
        outer_rings_lonlat=(
            (
                (0.0200, 0.0000),
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
    observed: dict[str, object] = {}

    def _simplify_water_footprints_for_observer(footprints, **kwargs):
        observed["input"] = tuple(footprints)
        return (simplified,)

    def _sample_water_overlay_points(footprints, **kwargs):
        observed["sampled"] = tuple(footprints)
        return ()

    monkeypatch.setattr(
        mod,
        "simplify_water_footprints_for_observer",
        _simplify_water_footprints_for_observer,
    )
    monkeypatch.setattr(mod, "sample_water_overlay_points", _sample_water_overlay_points)
    scope_cache = mod._WaterOverlayScopeCache(  # noqa: SLF001
        footprints=(footprint,),
        fetched_at_utc=datetime.now(timezone.utc),
    )

    controller._build_requested_variants(  # noqa: SLF001
        scope_cache,
        observer_lat_deg=0.0,
        observer_lon_deg=0.0,
        observer_height_m=0.0,
        observer_ground_m=0.0,
        use_dem_ground=False,
        scan_radius_km=2.0,
    )

    assert observed["input"] == (footprint,)
    assert observed["sampled"] == (simplified,)
