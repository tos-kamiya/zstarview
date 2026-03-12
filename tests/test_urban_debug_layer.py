from __future__ import annotations

from pathlib import Path

from zstarview.types import ViewerData
from zstarview.urban_debug_layer import resolve_urban_debug_layer_for_viewer


def test_resolve_urban_debug_layer_for_viewer_builds_dynamic_layer(monkeypatch, tmp_path: Path) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Lat: 35.67, Lon: 139.76",
        observer_height_m=1.7,
    )

    monkeypatch.setattr(
        "zstarview.urban_debug_layer.select_derived_tile_envelopes",
        lambda *_args, **_kwargs: (type("Envelope", (), {"path": tmp_path / "tile.json"})(),),
    )
    monkeypatch.setattr(
        "zstarview.urban_debug_layer.parse_derived_tile_buildings",
        lambda *_args, **_kwargs: ("building",),
    )
    monkeypatch.setattr(
        "zstarview.urban_debug_layer.compute_debug_outlines",
        lambda *_args, **_kwargs: type(
            "Result",
            (),
            {
                "outlines": (
                    (
                        type("Point", (), {"altitude_deg": -1.0, "azimuth_deg": 10.0})(),
                        type("Point", (), {"altitude_deg": -2.0, "azimuth_deg": 12.0})(),
                    ),
                )
            },
        )(),
    )

    got = resolve_urban_debug_layer_for_viewer(
        viewer,
        tokyo23_derived_dir=tmp_path,
    )

    assert got == [[(-1.0, 10.0), (-2.0, 12.0)]]


def test_resolve_urban_debug_layer_for_viewer_returns_none_when_outside_coverage(tmp_path: Path) -> None:
    viewer = ViewerData(
        location=(34.69, 135.50),
        timezone_name="Asia/Tokyo",
        city_name="Osaka",
        observer_height_m=1.7,
    )

    got = resolve_urban_debug_layer_for_viewer(
        viewer,
        tokyo23_derived_dir=tmp_path,
    )

    assert got is None
