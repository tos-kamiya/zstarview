from __future__ import annotations

import zstarview.gui.water_overlay_controller as mod
from zstarview.gui.water_overlay_controller import WaterOverlayController
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
        "sample_water_surface_interface_points",
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
        "sample_water_surface_interface_points",
        lambda **kwargs: (
            observed.setdefault("kwargs", kwargs),
            WaterOverlayPoint("water-mask", 10.0, 20.0, 0.5, water_category="sea"),
        )[1:],
    )

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
    assert "Water mask points: 1 visible, nearest sea point 0.500 km" in caplog.text
