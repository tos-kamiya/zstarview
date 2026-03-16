from __future__ import annotations

from pathlib import Path

from zstarview.types import UrbanOutlinePolyline, ViewerData
from zstarview.urban_outline_layer import resolve_urban_outline_layer_for_viewer


def test_resolve_urban_outline_layer_for_viewer_builds_dynamic_layer(monkeypatch, tmp_path: Path) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Lat: 35.67, Lon: 139.76",
        observer_height_m=1.7,
    )
    derived_dir = tmp_path / "13100_tokyo23" / "bldg"
    derived_dir.mkdir(parents=True)

    monkeypatch.setattr(
        "zstarview.urban_outline_layer.select_derived_tile_envelopes",
        lambda *_args, **_kwargs: (type("Envelope", (), {"path": tmp_path / "tile.json"})(),),
    )
    parse_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.parse_derived_tile_buildings",
        lambda *args, **kwargs: parse_calls.append((args, kwargs)) or ("building",),
    )
    compute_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.compute_urban_outlines",
        lambda *args, **_kwargs: compute_calls.append(args) or type(
            "Result",
            (),
            {
                "outlines": (
                    type(
                        "Outline",
                        (),
                        {
                            "height_m": 45.0,
                            "points": (
                                type("Point", (), {"altitude_deg": -1.0, "azimuth_deg": 10.0})(),
                                type("Point", (), {"altitude_deg": -2.0, "azimuth_deg": 12.0})(),
                            ),
                        },
                    )(),
                )
            },
        )(),
    )

    got = resolve_urban_outline_layer_for_viewer(
        viewer,
        derived_root_dir=tmp_path,
    )

    assert got == [UrbanOutlinePolyline(points=[(-1.0, 10.0), (-2.0, 12.0)], height_m=45.0)]
    assert parse_calls[0][1] == {}
    assert compute_calls[0][0].viewpoint_height_m == 1.7
    assert compute_calls[0][0].observer_height_m == 1.7


def test_resolve_urban_outline_layer_for_viewer_returns_none_when_outside_coverage(tmp_path: Path) -> None:
    viewer = ViewerData(
        location=(34.69, 135.50),
        timezone_name="Asia/Tokyo",
        city_name="Osaka",
        observer_height_m=1.7,
    )
    derived_dir = tmp_path / "13100_tokyo23" / "bldg"
    derived_dir.mkdir(parents=True)

    got = resolve_urban_outline_layer_for_viewer(
        viewer,
        derived_root_dir=tmp_path,
    )

    assert got is None


def test_resolve_urban_outline_layer_for_viewer_prefers_explicit_derived_dir(
    monkeypatch,
    tmp_path: Path,
) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Lat: 35.67, Lon: 139.76",
        observer_height_m=1.7,
    )
    tokyo_dir = tmp_path / "13100_tokyo23" / "bldg"
    kyoto_dir = tmp_path / "26100_kyoto" / "bldg"
    tokyo_dir.mkdir(parents=True)
    kyoto_dir.mkdir(parents=True)
    seen_dirs = []

    def fake_select(derived_dir: Path, **_kwargs):
        seen_dirs.append(derived_dir)
        if derived_dir == tokyo_dir:
            return (type("Envelope", (), {"path": tmp_path / "tokyo.json"})(),)
        raise ValueError("outside coverage")

    monkeypatch.setattr("zstarview.urban_outline_layer.select_derived_tile_envelopes", fake_select)
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.parse_derived_tile_buildings",
        lambda *_args, **_kwargs: ("building",),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.compute_urban_outlines",
        lambda *_args, **_kwargs: type(
            "Result",
            (),
            {
                "outlines": (
                    type(
                        "Outline",
                        (),
                        {
                            "height_m": 45.0,
                            "points": (
                                type("Point", (), {"altitude_deg": -1.0, "azimuth_deg": 10.0})(),
                            ),
                        },
                    )(),
                )
            },
        )(),
    )

    got = resolve_urban_outline_layer_for_viewer(
        viewer,
        derived_root_dir=tmp_path,
        derived_dir=tokyo_dir,
    )

    assert got == [UrbanOutlinePolyline(points=[(-1.0, 10.0)], height_m=45.0)]
    assert seen_dirs == [tokyo_dir]
