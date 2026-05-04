from __future__ import annotations

import numpy as np
import pytest
from rasterio.transform import Affine

from zstarview.terrain.dem import WGS84_GEOD
from zstarview.terrain.horizon import (
    DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM,
    ObserverLocation,
    compute_horizon_layers,
    _prune_secondary_peak_indices_by_visibility,
    _prune_secondary_peak_indices_by_main_profile,
    _select_distance_band_peak_index,
    _should_break_secondary_ridge,
    build_distance_samples,
    compute_apparent_altitudes,
)
from zstarview.render.terrain import _distance_band_alpha, _distance_band_underlay_alpha, _distance_band_underlay_width, _distance_band_widths
from zstarview.render.terrain import draw_terrain_secondary_ridges
from zstarview.terrain.dem import DemGrid, sample_ground_elevation


class _IdentityTransformer:
    def transform(self, lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return lon_deg, lat_deg


class _FlatDemGrid:
    def sample_lonlat(
        self,
        lon_deg: np.ndarray,
        lat_deg: np.ndarray,
        *,
        method: str = "bilinear",
    ) -> np.ndarray:
        return np.zeros_like(lon_deg, dtype=np.float64)


def test_build_distance_samples_covers_requested_distance() -> None:
    samples = build_distance_samples(1.0, 300.0)

    assert len(samples) == 4
    assert samples[0] == 300.0
    assert samples[-1] == 1000.0


def test_compute_apparent_altitudes_is_zero_for_equal_height_flat_short_range() -> None:
    got = compute_apparent_altitudes(
        observer_elevation_m=100.0,
        target_elevation_m=np.array([100.0]),
        surface_distance_m=np.array([1.0]),
        earth_radius_m=6_371_008.8,
        refraction_coefficient=0.0,
    )

    assert abs(float(got[0])) < 1e-3


def test_should_break_secondary_ridge_on_large_distance_jump() -> None:
    assert _should_break_secondary_ridge(
        3_000.0,
        9_500.0,
        max_distance_jump_ratio=0.25,
    )


def test_should_keep_secondary_ridge_on_small_distance_change() -> None:
    assert not _should_break_secondary_ridge(
        3_000.0,
        3_700.0,
        max_distance_jump_ratio=0.25,
    )


def test_prune_secondary_peak_indices_by_visibility_drops_farther_lower_peak() -> None:
    altitudes = np.array([0.2, 0.7, 0.6, 0.9], dtype=np.float64)

    got = _prune_secondary_peak_indices_by_visibility(altitudes, [0, 1, 2, 3])

    assert got == [0, 1, 3]


def test_prune_secondary_peak_indices_by_visibility_keeps_farther_higher_peak() -> None:
    altitudes = np.array([0.2, 0.1, 0.4, 0.3], dtype=np.float64)

    got = _prune_secondary_peak_indices_by_visibility(altitudes, [0, 1, 2, 3])

    assert got == [0, 2]


def test_prune_secondary_peak_indices_by_main_profile_drops_lower_than_main() -> None:
    altitudes = np.array([0.2, 0.7, 0.6, 0.9], dtype=np.float64)
    distances = np.array([100.0, 200.0, 300.0, 400.0], dtype=np.float64)

    got = _prune_secondary_peak_indices_by_main_profile(
        altitudes,
        distances,
        [0, 1, 2, 3],
        main_peak_distance_m=250.0,
    )

    assert got == [0, 1]


def test_prune_secondary_peak_indices_by_main_profile_keeps_higher_than_main() -> None:
    altitudes = np.array([0.2, 0.7, 0.6, 0.9], dtype=np.float64)
    distances = np.array([100.0, 200.0, 300.0, 400.0], dtype=np.float64)

    got = _prune_secondary_peak_indices_by_main_profile(
        altitudes,
        distances,
        [0, 1, 2, 3],
        main_peak_distance_m=450.0,
    )

    assert got == [0, 1, 2, 3]


def test_select_distance_band_peak_index_prefers_highest_in_band() -> None:
    altitudes = np.array([0.1, 0.6, 0.3, 0.9, 0.4], dtype=np.float64)
    distances = np.array([5_000.0, 20_000.0, 35_000.0, 42_000.0, 70_000.0], dtype=np.float64)

    got = _select_distance_band_peak_index(
        altitudes,
        distances,
        band_min_m=10_000.0,
        band_max_m=50_000.0,
        include_upper_bound=False,
    )

    assert got == 3


def test_select_distance_band_peak_index_includes_last_upper_bound() -> None:
    altitudes = np.array([0.1, 0.6, 0.3, 0.9, 0.4], dtype=np.float64)
    distances = np.array([5_000.0, 20_000.0, 35_000.0, 42_000.0, 120_000.0], dtype=np.float64)

    got = _select_distance_band_peak_index(
        altitudes,
        distances,
        band_min_m=80_000.0,
        band_max_m=120_000.0,
        include_upper_bound=True,
    )

    assert got == 4


def test_distance_band_alpha_drops_with_distance() -> None:
    near = _distance_band_alpha(distance_km=0.5, band_count=9, opacity=0.38)
    mid = _distance_band_alpha(distance_km=4.0, band_count=9, opacity=0.38)
    far = _distance_band_alpha(distance_km=128.0, band_count=9, opacity=0.38)

    assert near > mid > far
    assert far < 0.08


def test_secondary_ridge_render_uses_four_km_alpha_for_all_bands(monkeypatch) -> None:
    seen_alphas: list[float] = []

    def fake_solid_pen(_color_rgb, alpha, _width):
        seen_alphas.append(float(alpha))
        return object()

    def fake_draw_glow(*_args, **_kwargs) -> None:
        return None

    class _FakePainter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

        def drawLine(self, *_args, **_kwargs) -> None:
            pass

    monkeypatch.setattr("zstarview.render.terrain._solid_pen", fake_solid_pen)
    monkeypatch.setattr("zstarview.render.terrain._draw_terrain_secondary_ridge_glow", fake_draw_glow)

    draw_terrain_secondary_ridges(
        _FakePainter(),  # type: ignore[arg-type]
        geometry=type("Geometry", (), {"center": (0, 0), "radius": 100})(),
        terrain_secondary_profile_layers=[
            [(1.0, 10.0), (2.0, 20.0)],
            [(1.0, 10.0), (2.0, 20.0)],
        ],
        terrain_secondary_profile_distances_m_layers=[
            [500.0, 500.0],
            [128000.0, 128000.0],
        ],
        view_center=(0.0, 0.0),
        opacity=0.38,
        line_width_scale=1.0,
        fast_mode=False,
        is_in_fov_func=lambda *_args, **_kwargs: True,
        altaz_to_normalized_xy_func=lambda alt_deg, az_deg, _view_center, *, edge_fov_deg=95.0: (
            alt_deg / 90.0,
            az_deg / 180.0,
        ),
        normalized_to_screen_xy_func=lambda nx, ny, _geometry: (nx, ny),
        split_by_gaps_func=lambda points: [points],
    )

    assert len(seen_alphas) == 4
    expected_band_alpha = _distance_band_alpha(distance_km=4.0, band_count=2, opacity=0.38)
    expected_underlay_alpha = _distance_band_underlay_alpha(distance_km=4.0, band_count=2, opacity=0.38)
    assert seen_alphas[0] == pytest.approx(expected_underlay_alpha)
    assert seen_alphas[1] == pytest.approx(expected_band_alpha)
    assert seen_alphas[2] == pytest.approx(expected_underlay_alpha)
    assert seen_alphas[3] == pytest.approx(expected_band_alpha)


def test_distance_band_underlay_blur_increases_while_alpha_drops() -> None:
    near_width = _distance_band_underlay_width(distance_km=0.5, band_count=9)
    mid_width = _distance_band_underlay_width(distance_km=4.0, band_count=9)
    far_width = _distance_band_underlay_width(distance_km=128.0, band_count=9)
    near_foreground = _distance_band_widths(distance_km=0.5, band_count=9)
    mid_foreground = _distance_band_widths(distance_km=4.0, band_count=9)
    far_foreground = _distance_band_widths(distance_km=128.0, band_count=9)
    near_alpha = _distance_band_underlay_alpha(distance_km=0.5, band_count=9, opacity=0.38)
    mid_alpha = _distance_band_underlay_alpha(distance_km=4.0, band_count=9, opacity=0.38)
    far_alpha = _distance_band_underlay_alpha(distance_km=128.0, band_count=9, opacity=0.38)

    assert near_width > near_foreground
    assert mid_width > mid_foreground
    assert far_width > far_foreground
    assert near_alpha > mid_alpha > far_alpha
    assert far_alpha < 0.01


def test_compute_horizon_layers_adds_nearest_secondary_band() -> None:
    layers = compute_horizon_layers(
        dem_grid=_FlatDemGrid(),
        geod=WGS84_GEOD,
        observer=ObserverLocation(
            latitude_deg=35.0,
            longitude_deg=139.0,
            observer_ground_m=0.0,
            observer_eye_m=0.0,
        ),
        azimuth_step_deg=180.0,
        distance_samples_m=build_distance_samples(128.0, 500.0),
        dem_resampling="bilinear",
        earth_radius_m=6_371_008.8,
        refraction_coefficient=0.0,
    )

    assert DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM[:2] == (0.5, 1.0)
    assert len(layers.secondary_layers) >= 2
    assert [point.distance_m for point in layers.secondary_layers[0]] == [500.0, 500.0]
    assert [point.distance_m for point in layers.secondary_layers[1]] == [1000.0, 1000.0]


def test_distance_band_widths_start_thickest_and_taper_down() -> None:
    widths = [
        _distance_band_widths(distance_km=distance_km, band_count=9)
        for distance_km in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0)
    ]

    assert widths[0] > widths[1] > widths[2]
    assert widths[0] == max(widths)


def test_distance_band_widths_follow_distance_even_with_more_bands() -> None:
    near = _distance_band_widths(distance_km=0.5, band_count=9)
    mid = _distance_band_widths(distance_km=1.0, band_count=9)
    farther = _distance_band_widths(distance_km=2.0, band_count=9)

    assert near > mid > farther


def test_sample_ground_elevation_uses_default_elevation_for_nodata() -> None:
    grid = DemGrid(
        array=np.array([[np.nan]], dtype=np.float64),
        transform=Affine.translation(0.0, 1.0) * Affine.scale(1.0, -1.0),
        crs="EPSG:4326",
        default_elevation_m=0.0,
        _to_grid=_IdentityTransformer(),
    )

    got = sample_ground_elevation(
        grid,
        latitude_deg=0.5,
        longitude_deg=0.5,
    )

    assert got == 0.0
