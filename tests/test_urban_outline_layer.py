from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import numpy as np

from zstarview.data.skyscraper_tiles import SKYSCRAPER_OUTER_RADIUS_KM
from zstarview.data.urban_outline_common import BuildingFootprint
from zstarview.types import UrbanOutlinePolyline, ViewerData
from zstarview.urban_outline_layer import (
    _merge_building_footprints,
    _resolve_building_ground_elevations,
    _build_observer_centric_urban_outline_result,
    resolve_urban_outline_layer_for_viewer,
)


def test_resolve_urban_outline_layer_for_viewer_builds_dynamic_layer(monkeypatch, tmp_path: Path) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Lat: 35.67, Lon: 139.76",
        observer_height_m=1.7,
    )
    derived_dir = tmp_path / "13100_tokyo23" / "bldg"
    derived_dir.mkdir(parents=True)

    select_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.select_derived_tile_envelopes",
        lambda *args, **kwargs: select_calls.append((args, kwargs))
        or (type("Envelope", (), {"path": tmp_path / "tile.json"})(),),
    )
    parse_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.parse_derived_tile_buildings",
        lambda *args, **kwargs: parse_calls.append((args, kwargs))
        or (
            BuildingFootprint(
                building_id="building-1",
                height_m=45.0,
                rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.0, 35.0)),),
            ),
        ),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer._resolve_building_ground_elevations",
        lambda **kwargs: (12.0, kwargs["buildings"]),
    )
    compute_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.compute_urban_outlines",
        lambda *args, **kwargs: compute_calls.append((args, kwargs)) or type(
            "Result",
            (),
            {
                "outlines": (
                    type(
                        "Outline",
                        (),
                        {
                            "height_m": 45.0,
                            "distance_km": 0.85,
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

    assert got == [
        UrbanOutlinePolyline(
            points=[(-2.0, 12.0), (-1.0, 10.0)],
            height_m=45.0,
            distance_km=0.85,
        )
    ]
    assert select_calls[0][1]["radius_km"] == 2.5
    assert parse_calls[0][1] == {}
    assert compute_calls[0][0][0].viewpoint_height_m == 1.7
    assert compute_calls[0][0][0].observer_height_m == 1.7
    assert compute_calls[0][1]["observer_ground_elevation_m"] == 12.0


def test_resolve_urban_outline_layer_for_viewer_passes_far_range_distance_filters(
    monkeypatch,
    tmp_path: Path,
) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Lat: 35.67, Lon: 139.76",
        observer_height_m=1.7,
    )
    derived_dir = tmp_path / "skyscrapers" / "bldg"
    derived_dir.mkdir(parents=True)

    select_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.select_derived_tile_envelopes",
        lambda *args, **kwargs: select_calls.append((args, kwargs))
        or (type("Envelope", (), {"path": tmp_path / "tile.json"})(),),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.parse_derived_tile_buildings",
        lambda *_args, **_kwargs: (
            BuildingFootprint(
                building_id="skyscraper-1",
                height_m=180.0,
                rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.0, 35.0)),),
            ),
        ),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer._resolve_building_ground_elevations",
        lambda **kwargs: (0.0, kwargs["buildings"]),
    )
    compute_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.compute_urban_outlines",
        lambda *args, **kwargs: compute_calls.append((args, kwargs))
        or type("Result", (), {"outlines": ()})(),
    )

    got = resolve_urban_outline_layer_for_viewer(
        viewer,
        derived_root_dir=tmp_path,
        derived_dir=derived_dir,
        radius_km=SKYSCRAPER_OUTER_RADIUS_KM,
        min_distance_km=2.5,
    )

    assert got is None
    assert select_calls[0][1]["radius_km"] == SKYSCRAPER_OUTER_RADIUS_KM
    assert compute_calls[0][1]["radius_km"] == SKYSCRAPER_OUTER_RADIUS_KM
    assert compute_calls[0][1]["min_distance_km"] == 2.5


def test_resolve_urban_outline_layer_for_viewer_applies_runtime_min_height_filter(
    monkeypatch,
    tmp_path: Path,
) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Lat: 35.67, Lon: 139.76",
        observer_height_m=1.7,
    )
    derived_dir = tmp_path / "skyscrapers" / "bldg"
    derived_dir.mkdir(parents=True)

    monkeypatch.setattr(
        "zstarview.urban_outline_layer.select_derived_tile_envelopes",
        lambda *_args, **_kwargs: (type("Envelope", (), {"path": tmp_path / "tile.json"})(),),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.parse_derived_tile_buildings",
        lambda *_args, **_kwargs: (
            BuildingFootprint(
                building_id="lower",
                height_m=170.0,
                rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.0, 35.0)),),
            ),
            BuildingFootprint(
                building_id="higher",
                height_m=210.0,
                rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.0, 35.0)),),
            ),
        ),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer._resolve_building_ground_elevations",
        lambda **kwargs: (0.0, kwargs["buildings"]),
    )
    compute_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.compute_urban_outlines",
        lambda *args, **kwargs: compute_calls.append((args, kwargs))
        or type("Result", (), {"outlines": ()})(),
    )

    got = resolve_urban_outline_layer_for_viewer(
        viewer,
        derived_root_dir=tmp_path,
        derived_dir=derived_dir,
        radius_km=SKYSCRAPER_OUTER_RADIUS_KM,
        min_distance_km=2.5,
        min_height_m=200.0,
    )

    assert got is None
    assert tuple(building.building_id for building in compute_calls[0][0][1]) == ("higher",)


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
        lambda *_args, **_kwargs: (
            BuildingFootprint(
                building_id="building-1",
                height_m=45.0,
                rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.0, 35.0)),),
            ),
        ),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer._resolve_building_ground_elevations",
        lambda **kwargs: (0.0, kwargs["buildings"]),
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
                            "distance_km": 0.85,
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

    assert got == [
        UrbanOutlinePolyline(
            points=[(-1.0, 10.0)],
            height_m=45.0,
            distance_km=0.85,
        )
    ]
    assert seen_dirs == [tokyo_dir]


def test_resolve_urban_outline_layer_for_viewer_reprojects_cached_topocentric_result(
    monkeypatch,
    tmp_path: Path,
) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Lat: 35.67, Lon: 139.76",
        observer_height_m=1.7,
        view_center=(1.0, 20.0),
    )
    derived_dir = tmp_path / "13100_tokyo23" / "bldg"
    derived_dir.mkdir(parents=True)

    _build_observer_centric_urban_outline_result.cache_clear()

    compute_calls = []
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.select_derived_tile_envelopes",
        lambda *args, **kwargs: (type("Envelope", (), {"path": tmp_path / "tile.json"})(),),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.parse_derived_tile_buildings",
        lambda *_args, **_kwargs: (
            BuildingFootprint(
                building_id="building-1",
                height_m=45.0,
                rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.0, 35.0)),),
            ),
        ),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer._resolve_building_ground_elevations",
        lambda **kwargs: (0.0, kwargs["buildings"]),
    )
    def fake_linearize(run_points, *, view_center, edge_fov_deg):
        if float(view_center[0]) < 5.0:
            return [run_points[0], run_points[-1]]
        return list(run_points)

    monkeypatch.setattr(
        "zstarview.urban_outline_layer._maybe_linearize_run_points",
        fake_linearize,
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.compute_urban_outlines",
        lambda *args, **kwargs: compute_calls.append((args, kwargs))
        or type(
            "Result",
            (),
            {
                "outlines": (
                    type(
                        "Outline",
                        (),
                        {
                            "height_m": 45.0,
                            "distance_km": 0.85,
                            "points": (
                                type("Point", (), {"altitude_deg": 1.0, "azimuth_deg": 10.0})(),
                                type("Point", (), {"altitude_deg": 1.0, "azimuth_deg": 20.0})(),
                                type("Point", (), {"altitude_deg": 1.0, "azimuth_deg": 30.0})(),
                            ),
                        },
                    )(),
                )
            },
        )(),
    )

    first = resolve_urban_outline_layer_for_viewer(
        viewer,
        derived_root_dir=tmp_path,
        derived_dir=derived_dir,
    )
    second = resolve_urban_outline_layer_for_viewer(
        replace(viewer, view_center=(10.0, 20.0)),
        derived_root_dir=tmp_path,
        derived_dir=derived_dir,
    )

    assert len(compute_calls) == 1
    assert first is not None
    assert second is not None
    assert len(first[0].points) == 2
    assert len(second[0].points) == 3


def test_resolve_building_ground_elevations_falls_back_to_zero_when_dem_missing(
    monkeypatch,
) -> None:
    building = BuildingFootprint(
        building_id="building-1",
        height_m=45.0,
        rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.001, 35.001), (139.0, 35.0)),),
    )
    monkeypatch.setattr(
        "zstarview.urban_outline_layer.fetch_copernicus_dem",
        lambda **_kwargs: (_ for _ in ()).throw(
            RuntimeError("No Copernicus DEM tiles were downloaded for the requested area.")
        ),
    )

    observer_ground_m, buildings = _resolve_building_ground_elevations(
        lat_deg=35.67,
        lon_deg=139.76,
        buildings=(building,),
        radius_km=2.5,
        dem_cache_dir=Path("/tmp/unused"),
    )

    assert observer_ground_m == 0.0
    assert buildings == (building,)


def test_merge_building_footprints_prefers_parts_over_parent_buildings() -> None:
    base_building = BuildingFootprint(
        building_id="building-1",
        height_m=120.0,
        rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.0, 35.0)),),
    )
    building_part = BuildingFootprint(
        building_id="part-1",
        height_m=130.0,
        rings_lonlat=(((139.0, 35.0), (139.0005, 35.0), (139.0, 35.0)),),
        parent_building_id="building-1",
    )
    untouched_building = BuildingFootprint(
        building_id="building-2",
        height_m=80.0,
        rings_lonlat=(((139.0, 35.0), (139.002, 35.0), (139.0, 35.0)),),
    )

    got = _merge_building_footprints((base_building, building_part, untouched_building))

    assert got == (building_part, untouched_building)


def test_resolve_building_ground_elevations_preserves_min_height(monkeypatch) -> None:
    building = BuildingFootprint(
        building_id="floating-part",
        height_m=45.0,
        min_height_m=12.0,
        rings_lonlat=(((139.0, 35.0), (139.001, 35.0), (139.001, 35.001), (139.0, 35.0)),),
    )

    class _FakeGrid:
        def sample_lonlat(self, lon, lat, method="bilinear"):
            return np.full(lon.shape, 100.0, dtype=np.float64)

    class _FakeDem:
        def __init__(self, *_args, **_kwargs):
            pass

        def build_grid(self, _bbox):
            return _FakeGrid()

        def close(self):
            return None

    monkeypatch.setattr(
        "zstarview.urban_outline_layer.fetch_copernicus_dem",
        lambda **_kwargs: type("Download", (), {"paths": (Path("/tmp/fake.tif"),)})(),
    )
    monkeypatch.setattr("zstarview.urban_outline_layer.GeoTiffDem", _FakeDem)

    observer_ground_m, buildings = _resolve_building_ground_elevations(
        lat_deg=35.67,
        lon_deg=139.76,
        buildings=(building,),
        radius_km=2.5,
        dem_cache_dir=Path("/tmp/unused"),
    )

    assert observer_ground_m == 100.0
    assert buildings[0].ground_elevation_m == 100.0
    assert buildings[0].min_height_m == 12.0
