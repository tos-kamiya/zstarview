from __future__ import annotations

from zstarview.gui.water_overlay_controller import WaterOverlayController


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
