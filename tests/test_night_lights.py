from __future__ import annotations

import json

import numpy as np
import pytest

from zstarview import night_lights
from zstarview import night_light_source


def test_terrain_visibility_threshold_curve_uses_distance_order() -> None:
    distances_m = np.asarray([500.0, 1000.0, 1500.0, 2500.0], dtype=np.float64)
    threshold = night_lights._terrain_visibility_threshold_curve(
        azimuth_deg=180.0,
        distances_m=distances_m,
        terrain_profile_altaz=[(1.0, 180.0), (2.0, 180.0)],
        terrain_profile_distances_m=[1000.0, 2000.0],
        terrain_secondary_ridges_altaz_layers=[
            [(3.0, 180.0)],
        ],
        terrain_secondary_ridges_distances_m_layers=[
            [1500.0],
        ],
    )

    assert threshold is not None
    assert np.allclose(
        threshold,
        np.asarray([-np.inf, 1.0, 3.0, 3.0], dtype=np.float64),
        equal_nan=False,
    )


def test_terrain_visibility_threshold_grid_stacks_distance_columns(monkeypatch) -> None:
    def fake_visibility_threshold_curve(*, azimuth_deg: float, distances_m: np.ndarray, **_kwargs) -> np.ndarray:
        return np.asarray(
            [float(azimuth_deg) + float(distances_m[0]) / 1000.0, float(azimuth_deg) + float(distances_m[1]) / 1000.0],
            dtype=np.float64,
        )

    monkeypatch.setattr(night_lights, "_terrain_visibility_threshold_curve", fake_visibility_threshold_curve)

    grid = night_lights._terrain_visibility_threshold_grid(
        az_grid=np.asarray([10.0, 20.0], dtype=np.float64),
        distances_m=np.asarray([1_000.0, 2_000.0], dtype=np.float64),
        terrain_profile_altaz=[(1.0, 10.0)],
        terrain_profile_distances_m=[1_000.0],
        terrain_secondary_ridges_altaz_layers=[[(2.0, 10.0)]],
        terrain_secondary_ridges_distances_m_layers=[[1_000.0]],
    )

    assert grid.shape == (2, 2)
    assert np.allclose(
        grid,
        np.asarray(
            [
                [11.0, 12.0],
                [21.0, 22.0],
            ],
            dtype=np.float64,
        ),
    )


def test_surface_point_apparent_altitudes_decrease_with_distance() -> None:
    distances_m = np.asarray([1000.0, 5000.0, 10000.0], dtype=np.float64)

    altitudes = night_lights._surface_point_apparent_altitudes(
        distances_m,
        observer_height_m=120.0,
        refraction_coefficient=0.13,
    )

    assert altitudes.shape == distances_m.shape
    assert np.all(np.isfinite(altitudes))
    assert altitudes[0] < altitudes[1] < altitudes[2] < 0.0


def test_apply_night_light_sample_floor_keeps_masked_samples_dark() -> None:
    samples = np.asarray([0.0, 1.0, 2.0, 0.0], dtype=np.float64)
    visibility_mask = np.asarray([True, True, False, True], dtype=bool)

    got = night_lights._apply_night_light_sample_floor(
        samples,
        visibility_mask,
        floor_value=0.2,
    )

    assert np.allclose(got, np.asarray([0.2, 1.2, 0.0, 0.2], dtype=np.float64))


def test_night_light_terrain_context_collects_inputs() -> None:
    terrain_sample_distances_m = np.asarray([1_000.0, 2_000.0], dtype=np.float64)
    terrain_sample_terrain_elevation_m = np.asarray([[100.0, 200.0]], dtype=np.float64)
    terrain_sample_azimuths_deg = [180.0]

    context = night_lights.NightLightTerrainContext.from_inputs(
        terrain_profile_altaz=[(1.0, 180.0)],
        terrain_profile_distances_m=[1_000.0],
        terrain_secondary_ridges_altaz_layers=[[(2.0, 190.0)]],
        terrain_secondary_ridges_distances_m_layers=[[2_000.0]],
        terrain_sample_azimuths_deg=terrain_sample_azimuths_deg,
        terrain_sample_distances_m=terrain_sample_distances_m,
        terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
    )

    assert context.terrain_profile_key == ((1.0, 180.0),)
    assert context.terrain_profile_distances_key == (1_000.0,)
    assert context.terrain_secondary_ridges_key == (((2.0, 190.0),),)
    assert context.terrain_secondary_ridges_distances_key == ((2_000.0,),)
    assert context.has_sample_grid
    assert context.terrain_sample_grid_key == (
        id(terrain_sample_azimuths_deg),
        id(terrain_sample_distances_m),
        id(terrain_sample_terrain_elevation_m),
    )
    assert context.terrain_sample_distances_m is terrain_sample_distances_m
    assert context.terrain_sample_terrain_elevation_m is terrain_sample_terrain_elevation_m
    assert context.terrain_sample_distances_key == (1_000.0, 2_000.0)
    assert context.terrain_sample_terrain_elevation_key == ((100.0, 200.0),)


def test_terrain_sample_edge_strength_rows_are_uniform_strengths() -> None:
    terrain_sample_distances_m = np.asarray([1_000.0, 3_000.0], dtype=np.float64)
    terrain_sample_terrain_elevation_m = np.asarray(
        [
            [100.0, 250.0],
            [50.0, 75.0],
        ],
        dtype=np.float64,
    )
    source_distances_m = np.asarray([1_000.0, 2_000.0, 3_000.0], dtype=np.float64)

    rows = night_lights._terrain_sample_edge_strength_rows(
        terrain_sample_distances_m=terrain_sample_distances_m,
        terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
        source_distances_m=source_distances_m,
    )

    assert rows is not None
    assert rows.shape == (2, 3)
    assert np.allclose(rows, np.ones((2, 3), dtype=np.float64))


def test_build_night_light_glow_fields_keeps_source_strengths_above_threshold(monkeypatch) -> None:
    monkeypatch.setattr(
        night_lights,
        "_accumulate_local_glow_strengths",
        lambda **_kwargs: np.asarray([0.8], dtype=np.float64),
    )
    monkeypatch.setattr(
        night_lights,
        "_accumulate_local_glow_field",
        lambda **_kwargs: np.full((night_lights._target_altitude_bins().size, 1), 0.2, dtype=np.float64),
    )

    fields = night_lights._build_night_light_glow_fields_from_samples(
        az_grid=np.asarray([180.0], dtype=np.float64),
        horizon_alt_values=np.asarray([-1.0], dtype=np.float64),
        distances_m=np.asarray([1_000.0], dtype=np.float64),
        source_matrix=np.asarray([[1.0]], dtype=np.float64),
        source_altitudes=np.asarray([[-10.0]], dtype=np.float64),
        terrain_context=night_lights.NightLightTerrainContext.from_inputs(
            terrain_profile_altaz=[(0.0, 180.0)],
            terrain_profile_distances_m=[1_000.0],
            terrain_secondary_ridges_altaz_layers=None,
            terrain_secondary_ridges_distances_m_layers=None,
        ),
        max_distance_km=3.0,
        terrain_visibility_threshold_grid=np.asarray([[0.0]], dtype=np.float64),
    )

    assert fields is not None
    strengths, field = fields
    assert strengths.shape == (1,)
    assert strengths[0] > 0.0
    assert field.shape[1] == 1


def test_night_light_distance_boost_grows_linearly() -> None:
    distances_m = np.asarray([0.0, 64_000.0, 128_000.0], dtype=np.float64)

    boost = night_lights._night_light_distance_boost(distances_m)

    assert np.allclose(boost, np.asarray([1.0, 1.5, 2.0], dtype=np.float64))


def test_ridge_glow_distance_gain_maps_far_edge_to_fifty() -> None:
    gain = night_lights._ridge_glow_distance_gain(max_distance_km=128.0)
    far_edge_boost = night_lights._night_light_distance_boost(
        np.asarray([128_000.0], dtype=np.float64),
        max_distance_km=128.0,
    )[0]

    assert np.isclose(gain, 127.5)
    assert np.isclose(far_edge_boost * gain, 255.0)


def test_night_light_distance_sigma_keeps_three_km_reference() -> None:
    base_sigma = night_lights.NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG

    sigma_3km = night_lights._night_light_distance_sigma_deg(3_000.0)
    sigma_6km = night_lights._night_light_distance_sigma_deg(6_000.0)

    assert np.isclose(sigma_3km, base_sigma)
    assert np.isclose(
        sigma_6km,
        base_sigma * (2.0 ** (-night_lights.NIGHT_LIGHTS_DISTANCE_SIGMA_GAMMA)),
    )
    assert sigma_6km < sigma_3km


def test_normalize_night_light_values_blends_linear_and_logarithmic_values() -> None:
    normalized = night_lights._normalize_night_light_values(
        np.asarray([0.5, 1.0], dtype=np.float64),
        1.0,
    )

    expected_log = np.log1p(0.5) / np.log1p(1.0)
    assert np.isclose(normalized[0], (0.5 + expected_log) / 2.0)
    assert np.isclose(normalized[1], 1.0)


def test_night_light_strength_factor_uses_minus_nine_to_minus_four_blend() -> None:
    assert np.isclose(night_lights.night_light_strength_factor(-9.0), 1.0)
    assert np.isclose(night_lights.night_light_strength_factor(-6.5), 0.5)
    assert np.isclose(night_lights.night_light_strength_factor(-4.0), 0.0)
    assert night_lights.is_night_light_enabled(-4.1)
    assert not night_lights.is_night_light_enabled(-4.0)


def test_release_manifest_schema_accepts_all_expected_tiles() -> None:
    manifest = {
        "dataset_version": night_light_source.NIGHT_LIGHTS_DATASET_VERSION,
        "tiles": {
            tile_name: {
                "path": f"{tile_name}.tif",
                "url": f"https://github.com/example/release/{tile_name}.tif",
                "sha256": "0" * 64,
                "width": 21600,
                "height": 21600,
                "resolution_degrees": [15.0 / 3600.0, 15.0 / 3600.0],
            }
            for tile_name in night_light_source.NIGHT_LIGHTS_TILE_NAMES
        },
    }

    got = night_light_source._validate_manifest(manifest)

    assert got["dataset_version"] == night_light_source.NIGHT_LIGHTS_DATASET_VERSION


def test_release_manifest_schema_rejects_non_https_tile_url() -> None:
    manifest = {
        "dataset_version": night_light_source.NIGHT_LIGHTS_DATASET_VERSION,
        "tiles": {
            tile_name: {
                "path": f"{tile_name}.tif",
                "url": f"https://github.com/example/release/{tile_name}.tif",
                "sha256": "0" * 64,
                "width": 21600,
                "height": 21600,
                "resolution_degrees": [15.0 / 3600.0, 15.0 / 3600.0],
            }
            for tile_name in night_light_source.NIGHT_LIGHTS_TILE_NAMES
        },
    }
    manifest["tiles"]["A1"]["url"] = "http://example.invalid/A1.tif"

    with pytest.raises(night_light_source.NightLightsManifestError):
        night_light_source._validate_manifest(manifest)


def test_release_tiles_download_only_requested_tiles(tmp_path, monkeypatch) -> None:
    manifest = {
        "dataset_version": night_light_source.NIGHT_LIGHTS_DATASET_VERSION,
        "tiles": {
            tile_name: {
                "path": f"{tile_name}.tif",
                "url": f"https://github.com/example/release/{tile_name}.tif",
                "sha256": "0" * 64,
                "width": 21600,
                "height": 21600,
                "resolution_degrees": [15.0 / 3600.0, 15.0 / 3600.0],
            }
            for tile_name in night_light_source.NIGHT_LIGHTS_TILE_NAMES
        },
    }
    dataset_root = tmp_path / night_light_source.NIGHT_LIGHTS_DATASET_VERSION
    dataset_root.mkdir()
    (dataset_root / "manifest.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )
    downloaded: list[str] = []

    def fake_download(url, destination, *, timeout_s):
        downloaded.append(destination.name)
        destination.write_bytes(b"tile")
        return "0" * 64

    monkeypatch.setattr(night_light_source, "_download_file", fake_download)
    monkeypatch.setattr(night_light_source, "_validate_geotiff", lambda *_args: None)

    paths = night_light_source._ensure_night_light_tiles(
        cache_root=tmp_path,
        tile_names={"A1", "D2"},
    )

    assert set(paths) == {"A1", "D2"}
    assert downloaded == ["A1.tif", "D2.tif"]


def test_gaussian_weight_lut_uses_half_degree_bins() -> None:
    sigma = night_lights.NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG
    weights = night_lights._lookup_gaussian_weights(
        np.asarray([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float64),
        sigma_deg=sigma,
        step_deg=0.5,
    )

    assert weights.shape == (5,)
    assert np.isclose(weights[0], 1.0)
    assert np.isclose(weights[1], weights[0])
    assert np.isclose(weights[2], night_lights._gaussian_weight_lut(sigma, 0.5)[1])
    assert weights[2] > weights[3] >= weights[4]


def test_gaussian_weights_clip_beyond_azimuth_limit() -> None:
    sigma = night_lights.NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG
    weights = night_lights._lookup_gaussian_weights(
        np.asarray([24.5, 25.0, 25.5, 30.0], dtype=np.float64),
        sigma_deg=sigma,
        step_deg=0.5,
        max_delta_deg=night_lights.NIGHT_LIGHTS_NEIGHBORHOOD_MAX_AZ_DELTA_DEG,
    )

    assert weights[0] > 0.0
    assert weights[1] > 0.0
    assert weights[2] == 0.0
    assert weights[3] == 0.0


def test_local_glow_supports_separate_altitude_and_azimuth_sigma() -> None:
    field = night_lights._accumulate_local_glow_field(
        source_azimuths_deg=np.asarray([0.0], dtype=np.float64),
        source_altitudes_deg=np.asarray([0.0], dtype=np.float64),
        source_strengths=np.asarray([1.0], dtype=np.float64),
        target_azimuths_deg=np.asarray([0.0, 4.0], dtype=np.float64),
        target_altitudes_deg=np.asarray([0.0, 4.0], dtype=np.float64),
        sigma_deg=4.0,
        azimuth_sigma_deg=2.0,
    )

    assert field[0, 1] < field[1, 0]


def test_target_altitude_bins_use_one_degree_step() -> None:
    bins = night_lights._target_altitude_bins()

    assert bins.size > 1
    assert np.isclose(bins[1] - bins[0], 1.0)
    assert np.isclose(bins[0], -90.0)


def test_compute_night_light_glow_profile_can_skip_night_light_tiles(monkeypatch) -> None:
    def fail_if_tiles_are_requested(**_kwargs):
        raise AssertionError("night light tiles should not be requested")

    monkeypatch.setattr(night_lights, "_ensure_night_light_tiles", fail_if_tiles_are_requested)

    profile = night_lights.compute_night_light_glow_profile(
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        sun_alt_deg=-5.0,
        terrain_profile_altaz=[(0.0, 180.0)],
        terrain_profile_distances_m=[1_000.0],
        terrain_secondary_ridges_altaz_layers=[[(1.0, 180.0)]],
        terrain_secondary_ridges_distances_m_layers=[[1_000.0]],
        terrain_sample_distances_m=np.asarray([1_000.0, 2_000.0], dtype=np.float64),
        terrain_sample_terrain_elevation_m=np.asarray([[100.0, 200.0]], dtype=np.float64),
        include_night_light_tiles=False,
    )

    assert profile is not None
    assert len(profile.samples) == 1
    assert np.any(np.asarray(profile.edge_alpha_grid, dtype=np.float64) > 0.0)
    assert np.allclose(np.asarray(profile.alpha_grid, dtype=np.float64), 0.0)
