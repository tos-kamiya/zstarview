from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "urban_skyline_demo.py"
    spec = importlib.util.spec_from_file_location("urban_skyline_demo", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_list_japan_towers_contains_tokyo_skytree() -> None:
    mod = _load_module()

    towers = mod.list_japan_towers()
    names = {tower.name for tower in towers}

    assert "Tokyo Skytree" in names
    assert "CN Tower" not in names


def test_compute_urban_skyline_detects_eastern_building() -> None:
    mod = _load_module()
    tower = next(tower for tower in mod.list_japan_towers() if tower.name == "Tokyo Skytree")
    tower = tower.__class__(
        id=tower.id,
        qid=tower.qid,
        kind=tower.kind,
        name=tower.name,
        labels=tower.labels,
        names=tower.names,
        latitude_deg=35.0,
        longitude_deg=139.0,
        height_m=100.0,
        observer_height_m=100.0,
        meta=tower.meta,
    )
    buildings = (
        mod.BuildingFootprint(
            building_id="east-block",
            height_m=150.0,
            rings_lonlat=(
                (
                    (139.00110, 35.00010),
                    (139.00110, 34.99990),
                    (139.00130, 34.99990),
                    (139.00130, 35.00010),
                    (139.00110, 35.00010),
                ),
            ),
        ),
    )

    result = mod.compute_urban_skyline(
        tower,
        buildings,
        radius_km=1.0,
        azimuth_step_deg=1.0,
        edge_sample_step_m=2.0,
    )

    samples_by_az = {round(sample.azimuth_deg): sample.altitude_deg for sample in result.samples}
    east_alt = samples_by_az[90]

    assert result.buildings_considered == 1
    assert result.buildings_contributing == 1
    assert east_alt > 15.0
    assert east_alt < 35.0
    assert math.isclose(result.peak_azimuth_deg, 90.0, abs_tol=10.0)


def test_skyline_result_to_payload_is_keyed_for_app_import() -> None:
    mod = _load_module()
    tower = next(tower for tower in mod.list_japan_towers() if tower.name == "Tokyo Skytree")
    result = mod.SkylineResult(
        tower=tower,
        samples=(
            mod.SkylineSample(azimuth_deg=0.0, altitude_deg=1.25),
            mod.SkylineSample(azimuth_deg=0.5, altitude_deg=1.5),
        ),
        buildings_considered=2,
        buildings_contributing=1,
        peak_altitude_deg=1.5,
        peak_azimuth_deg=0.5,
    )

    payload = mod.skyline_result_to_payload(result)

    assert payload == {
        "name": "Tokyo Skytree",
        "profile": [
            {"az": 0.0, "alt": 1.25},
            {"az": 0.5, "alt": 1.5},
        ],
    }


def test_load_building_footprints_filters_low_buildings(tmp_path: Path) -> None:
    mod = _load_module()
    path = tmp_path / "buildings.geojson"
    path.write_text(
        """
        {
          "type": "FeatureCollection",
          "features": [
            {
              "type": "Feature",
              "properties": {"id": "low", "height": 20},
              "geometry": {
                "type": "Polygon",
                "coordinates": [[[139.0, 35.0], [139.0, 35.001], [139.001, 35.001], [139.001, 35.0], [139.0, 35.0]]]
              }
            },
            {
              "type": "Feature",
              "properties": {"id": "high", "height": 80},
              "geometry": {
                "type": "Polygon",
                "coordinates": [[[139.0, 35.0], [139.0, 35.001], [139.001, 35.001], [139.001, 35.0], [139.0, 35.0]]]
              }
            }
          ]
        }
        """,
        encoding="utf-8",
    )

    buildings = mod.load_building_footprints(
        path,
        height_fields=("height",),
        min_building_height_m=50.0,
    )

    assert [building.building_id for building in buildings] == ["high"]


def test_skyline_radius_results_to_payload_keeps_multiple_radii() -> None:
    mod = _load_module()
    tower = next(tower for tower in mod.list_japan_towers() if tower.name == "Tokyo Skytree")
    radius_results = (
        mod.SkylineRadiusResult(
            radius_km=0.1,
            result=mod.SkylineResult(
                tower=tower,
                samples=(mod.SkylineSample(0.0, 1.25),),
                buildings_considered=1,
                buildings_contributing=1,
                peak_altitude_deg=1.25,
                peak_azimuth_deg=0.0,
            ),
        ),
        mod.SkylineRadiusResult(
            radius_km=1.0,
            result=mod.SkylineResult(
                tower=tower,
                samples=(mod.SkylineSample(0.0, 1.5),),
                buildings_considered=2,
                buildings_contributing=1,
                peak_altitude_deg=1.5,
                peak_azimuth_deg=0.0,
            ),
        ),
    )

    payload = mod.skyline_radius_results_to_payload(radius_results)

    assert payload == {
        "name": "Tokyo Skytree",
        "profiles": [
            {"radius_km": 0.1, "profile": [{"az": 0.0, "alt": 1.25}]},
            {"radius_km": 1.0, "profile": [{"az": 0.0, "alt": 1.5}]},
        ],
    }


def test_update_altitude_bins_from_polyline_fills_intermediate_bins() -> None:
    mod = _load_module()
    import numpy as np

    altitudes = np.full(10, -90.0, dtype=np.float64)
    updated = mod.update_altitude_bins_from_polyline(
        altitudes,
        azimuth_deg=np.array([1.0, 4.0], dtype=np.float64),
        altitude_deg=np.array([-10.0, -4.0], dtype=np.float64),
        azimuth_step_deg=1.0,
    )

    assert updated is True
    assert altitudes[1] > -90.0
    assert altitudes[2] > -90.0
    assert altitudes[3] > -90.0
    assert altitudes[4] > -90.0


def test_update_altitude_bins_from_polyline_wraps_across_north() -> None:
    mod = _load_module()
    import numpy as np

    altitudes = np.full(360, -90.0, dtype=np.float64)
    mod.update_altitude_bins_from_polyline(
        altitudes,
        azimuth_deg=np.array([359.0, 1.0], dtype=np.float64),
        altitude_deg=np.array([-5.0, -5.0], dtype=np.float64),
        azimuth_step_deg=1.0,
    )

    assert altitudes[359] > -90.0
    assert altitudes[0] > -90.0
    assert altitudes[1] > -90.0
    assert altitudes[180] == -90.0


def test_update_altitude_bins_from_polyline_open_path_does_not_close_loop() -> None:
    mod = _load_module()
    import numpy as np

    altitudes = np.full(360, -90.0, dtype=np.float64)
    mod.update_altitude_bins_from_polyline(
        altitudes,
        azimuth_deg=np.array([10.0, 20.0], dtype=np.float64),
        altitude_deg=np.array([-5.0, -5.0], dtype=np.float64),
        azimuth_step_deg=1.0,
        closed=False,
    )

    assert altitudes[10] > -90.0
    assert altitudes[15] > -90.0
    assert altitudes[20] > -90.0
    assert altitudes[200] == -90.0


def test_iter_true_runs_splits_disconnected_mask_regions() -> None:
    mod = _load_module()
    import numpy as np

    mask = np.array([False, True, True, False, True, False], dtype=bool)

    got = list(mod.iter_true_runs(mask))

    assert got == [slice(1, 3), slice(4, 5)]


def test_update_altitude_bins_from_interval_fills_constant_altitude() -> None:
    mod = _load_module()
    import numpy as np

    altitudes = np.full(360, -90.0, dtype=np.float64)

    updated = mod.update_altitude_bins_from_interval(
        altitudes,
        start_azimuth_deg=10.0,
        end_azimuth_deg=20.0,
        altitude_deg=-3.0,
        azimuth_step_deg=1.0,
    )

    assert updated is True
    assert altitudes[10] == -3.0
    assert altitudes[15] == -3.0
    assert altitudes[20] == -3.0
    assert altitudes[200] == -90.0


def test_compute_band_ends_m_uses_next_band_start() -> None:
    mod = _load_module()
    import numpy as np

    starts = np.array([888.8888889, 1333.3333333, 2000.0], dtype=np.float64)

    got = mod.compute_band_ends_m(starts, fallback_band_width_m=90.0)

    assert np.allclose(got, np.array([978.8888889, 1423.3333333, 2090.0]))


def test_compute_cumulative_urban_skyline_uses_radius_band_not_cumulative() -> None:
    mod = _load_module()
    tower = next(tower for tower in mod.list_japan_towers() if tower.name == "Tokyo Skytree")
    tower = tower.__class__(
        id=tower.id,
        qid=tower.qid,
        kind=tower.kind,
        name=tower.name,
        labels=tower.labels,
        names=tower.names,
        latitude_deg=35.0,
        longitude_deg=139.0,
        height_m=100.0,
        observer_height_m=100.0,
        meta=tower.meta,
    )
    buildings = (
        mod.BuildingFootprint(
            building_id="near-east",
            height_m=150.0,
            rings_lonlat=(
                (
                    (139.00500, 35.00010),
                    (139.00500, 34.99990),
                    (139.00520, 34.99990),
                    (139.00520, 35.00010),
                    (139.00500, 35.00010),
                ),
            ),
        ),
        mod.BuildingFootprint(
            building_id="far-east",
            height_m=180.0,
            rings_lonlat=(
                (
                    (139.01000, 35.00010),
                    (139.01000, 34.99990),
                    (139.01020, 34.99990),
                    (139.01020, 35.00010),
                    (139.01000, 35.00010),
                ),
            ),
        ),
    )

    results = mod.compute_cumulative_urban_skyline(
        tower,
        buildings,
        radii_km=(0.45, 0.9),
        radius_band_width_m=120.0,
        azimuth_step_deg=1.0,
        edge_sample_step_m=2.0,
    )

    near = results[0].result
    far = results[1].result
    near_samples = {round(sample.azimuth_deg): sample.altitude_deg for sample in near.samples}
    far_samples = {round(sample.azimuth_deg): sample.altitude_deg for sample in far.samples}

    assert near.buildings_considered >= 1
    assert far.buildings_considered >= 1
    assert near_samples[90] != far_samples[90]
