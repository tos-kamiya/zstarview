from __future__ import annotations

import numpy as np
import pytest
from PySide6.QtGui import QImage, QPainter
from PySide6.QtWidgets import QApplication

from zstarview import night_lights
from zstarview.render.night_lights import draw_night_light_glow
from zstarview.render import night_lights as night_lights_render
from zstarview.types import ScreenGeometry, ViewerData

app = QApplication.instance() or QApplication([])


def test_read_or_build_manifest_parses_tile_urls(tmp_path, monkeypatch) -> None:
    html = "\n".join(
        f'<a href="https://assets.science.nasa.gov/content/dam/science/esd/eo/images/'
        f'imagerecords/144000/144897/BlackMarble_2016_{tile}_geo_gray.tif">{tile}</a>'
        for tile in night_lights.NIGHT_LIGHTS_TILE_NAMES
    )
    monkeypatch.setattr(night_lights, "_read_url", lambda *_args, **_kwargs: html)

    manifest = night_lights._read_or_build_manifest(cache_root=tmp_path)

    assert manifest["dataset_version"] == night_lights.NIGHT_LIGHTS_DATASET_VERSION
    assert manifest["source_page_url"] == night_lights.NIGHT_LIGHTS_PAGE_URL
    assert set(manifest["tile_urls"]) == set(night_lights.NIGHT_LIGHTS_TILE_NAMES)
    assert (tmp_path / night_lights.NIGHT_LIGHTS_DATASET_VERSION / "manifest.json").exists()


def test_compute_night_light_glow_profile_returns_none_for_daytime() -> None:
    profile = night_lights.compute_night_light_glow_profile(
        observer_lat_deg=35.0,
        observer_lon_deg=139.0,
        sun_alt_deg=1.0,
    )
    assert profile is None


def test_night_light_strength_factor_is_continuous() -> None:
    assert night_lights.night_light_strength_factor(1.0) == 0.0
    assert night_lights.night_light_strength_factor(0.0) == 0.0
    assert 0.0 < night_lights.night_light_strength_factor(-3.0) < 1.0
    assert night_lights.night_light_strength_factor(-6.0) == 1.0
    assert night_lights.night_light_strength_factor(-10.0) == 1.0


def test_night_light_band_width_scale_decreases_with_distance() -> None:
    near = night_lights_render._night_light_band_width_scale(0.5)
    mid = night_lights_render._night_light_band_width_scale(1.0)
    far = night_lights_render._night_light_band_width_scale(128.0)

    assert near == 1.0
    assert near > mid > far
    assert far >= night_lights_render.NIGHT_LIGHTS_MIN_WIDTH_SCALE


def test_night_light_distance_scale_uses_inverse_sqrt() -> None:
    near = night_lights_render._night_light_distance_scale(0.5)
    mid = night_lights_render._night_light_distance_scale(2.0)
    far = night_lights_render._night_light_distance_scale(8.0)

    assert near == 1.0
    assert np.isclose(mid, 0.5)
    assert np.isclose(far, 0.25)


def test_layer_distance_km_treats_zero_as_first_secondary_ridge() -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=0.0, horizon_alt_deg=0.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
        band_profiles=(
            night_lights.NightLightDistanceBandProfile(
                min_distance_km=0.5,
                max_distance_km=3.0,
                samples=(
                    night_lights.NightLightGlowSample(
                        azimuth_deg=0.0,
                        horizon_alt_deg=0.0,
                        strength=1.0,
                    ),
                ),
            ),
            night_lights.NightLightDistanceBandProfile(
                min_distance_km=3.0,
                max_distance_km=9.0,
                samples=(
                    night_lights.NightLightGlowSample(
                        azimuth_deg=0.0,
                        horizon_alt_deg=0.0,
                        strength=1.0,
                    ),
                ),
            ),
        ),
    )

    assert np.isclose(night_lights_render._layer_distance_km(profile, 0), 3.0)
    assert np.isclose(night_lights_render._layer_distance_km(profile, 1), 9.0)
    with pytest.raises(AssertionError):
        night_lights_render._layer_distance_km(profile, 2)


def test_night_light_distance_attenuation_uses_inverse_square() -> None:
    attenuation = night_lights._night_light_distance_attenuation(
        np.asarray([1000.0, 2000.0, 4000.0], dtype=np.float64)
    )

    assert np.allclose(attenuation, np.asarray([1.0, 0.25, 0.0625], dtype=np.float64))


def test_band_lower_edge_altitudes_use_previous_layer_internal_division() -> None:
    current = np.asarray([2.0, 4.0], dtype=np.float64)
    previous = np.asarray([0.0, 2.0], dtype=np.float64)

    assert np.allclose(
        night_lights_render._band_lower_edge_altitudes(current, None),
        np.asarray([0.2, 0.4], dtype=np.float64),
    )
    assert np.allclose(
        night_lights_render._band_lower_edge_altitudes(current, previous),
        np.asarray([0.2, 2.2], dtype=np.float64),
    )


def test_compute_night_light_glow_profile_uses_mocked_sampling(tmp_path, monkeypatch) -> None:
    tile_paths = {}
    for tile in night_lights.NIGHT_LIGHTS_TILE_NAMES:
        path = tmp_path / f"{tile}.tif"
        path.write_text("dummy", encoding="utf-8")
        tile_paths[tile] = path

    monkeypatch.setattr(night_lights, "_ensure_night_light_tiles", lambda **_kwargs: tile_paths)
    monkeypatch.setattr(
        night_lights,
        "_build_azimuth_grid",
        lambda *_args, **_kwargs: (
            np.asarray([0.0, 90.0, 180.0], dtype=np.float64),
            np.zeros(3, dtype=np.float64),
        ),
    )
    monkeypatch.setattr(
        night_lights,
        "_sample_ray_brightness_curve",
        lambda **kwargs: np.cumsum(
            np.full(
                np.asarray(kwargs["distances_m"], dtype=np.float64).shape,
                float(kwargs["azimuth_deg"]) / 180.0,
                dtype=np.float64,
            )
        ),
    )

    night_lights._compute_night_light_base_profile.cache_clear()
    try:
        profile = night_lights.compute_night_light_glow_profile(
            observer_lat_deg=35.0,
            observer_lon_deg=139.0,
            sun_alt_deg=-5.0,
        )
    finally:
        night_lights._compute_night_light_base_profile.cache_clear()

    assert profile is not None
    assert len(profile.samples) == 3
    assert len(profile.band_profiles) > 0
    strengths = [sample.strength for sample in profile.samples]
    assert strengths[0] <= strengths[1] <= strengths[2]
    band_strengths = [band.samples[0].strength for band in profile.band_profiles]
    assert band_strengths[0] <= band_strengths[-1]


def test_compute_night_light_glow_profile_has_band_profiles(tmp_path, monkeypatch) -> None:
    tile_paths = {}
    for tile in night_lights.NIGHT_LIGHTS_TILE_NAMES:
        path = tmp_path / f"{tile}.tif"
        path.write_text("dummy", encoding="utf-8")
        tile_paths[tile] = path

    monkeypatch.setattr(night_lights, "_ensure_night_light_tiles", lambda **_kwargs: tile_paths)
    monkeypatch.setattr(
        night_lights,
        "_build_azimuth_grid",
        lambda *_args, **_kwargs: (
            np.asarray([0.0, 90.0, 180.0], dtype=np.float64),
            np.zeros(3, dtype=np.float64),
        ),
    )

    def _sample_ray_brightness_curve(**kwargs) -> np.ndarray:
        distances = np.asarray(kwargs["distances_m"], dtype=np.float64)
        return np.cumsum(np.ones_like(distances, dtype=np.float64))

    monkeypatch.setattr(night_lights, "_sample_ray_brightness_curve", _sample_ray_brightness_curve)
    night_lights._compute_night_light_base_profile.cache_clear()
    try:
        profile = night_lights.compute_night_light_glow_profile(
            observer_lat_deg=35.0,
            observer_lon_deg=139.0,
            sun_alt_deg=-5.0,
        )
    finally:
        night_lights._compute_night_light_base_profile.cache_clear()

    assert profile is not None
    assert len(profile.band_profiles) >= 2
    first_band = profile.band_profiles[0].samples[0].strength
    last_band = profile.band_profiles[-1].samples[0].strength
    assert first_band <= last_band


def test_compute_night_light_glow_profile_reuses_location_cache(tmp_path, monkeypatch) -> None:
    tile_paths = {}
    for tile in night_lights.NIGHT_LIGHTS_TILE_NAMES:
        path = tmp_path / f"{tile}.tif"
        path.write_text("dummy", encoding="utf-8")
        tile_paths[tile] = path

    monkeypatch.setattr(night_lights, "_ensure_night_light_tiles", lambda **_kwargs: tile_paths)
    monkeypatch.setattr(
        night_lights,
        "_build_azimuth_grid",
        lambda *_args, **_kwargs: (
            np.asarray([0.0, 90.0, 180.0], dtype=np.float64),
            np.zeros(3, dtype=np.float64),
        ),
    )
    call_count = {"count": 0}

    def _sample_ray_brightness(**kwargs) -> np.ndarray:
        call_count["count"] += 1
        distances = np.asarray(kwargs["distances_m"], dtype=np.float64)
        return np.cumsum(np.full(distances.shape, float(kwargs["azimuth_deg"]) / 180.0, dtype=np.float64))

    monkeypatch.setattr(night_lights, "_sample_ray_brightness_curve", _sample_ray_brightness)
    night_lights._compute_night_light_base_profile.cache_clear()
    try:
        profile1 = night_lights.compute_night_light_glow_profile(
            observer_lat_deg=35.0,
            observer_lon_deg=139.0,
            sun_alt_deg=-5.0,
        )
        profile2 = night_lights.compute_night_light_glow_profile(
            observer_lat_deg=35.0,
            observer_lon_deg=139.0,
            sun_alt_deg=-10.0,
        )
    finally:
        night_lights._compute_night_light_base_profile.cache_clear()

    assert profile1 is not None
    assert profile2 is not None
    assert call_count["count"] == 3


def test_draw_night_light_glow_smoke() -> None:
    image = QImage(200, 100, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        profile = night_lights.NightLightGlowProfile(
            samples=(
                night_lights.NightLightGlowSample(azimuth_deg=170.0, horizon_alt_deg=0.0, strength=0.4),
                night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
                night_lights.NightLightGlowSample(azimuth_deg=190.0, horizon_alt_deg=0.0, strength=0.4),
            ),
            sun_alt_deg=-5.0,
            band_profiles=(
                night_lights.NightLightDistanceBandProfile(
                    min_distance_km=0.5,
                    max_distance_km=3.0,
                    samples=(
                        night_lights.NightLightGlowSample(
                            azimuth_deg=170.0,
                            horizon_alt_deg=0.0,
                            strength=0.4,
                        ),
                        night_lights.NightLightGlowSample(
                            azimuth_deg=180.0,
                            horizon_alt_deg=0.0,
                            strength=1.0,
                        ),
                        night_lights.NightLightGlowSample(
                            azimuth_deg=190.0,
                            horizon_alt_deg=0.0,
                            strength=0.4,
                        ),
                    ),
                ),
            ),
        )
        viewer_data = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(0.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        )
        draw_night_light_glow(
            painter,
            geometry=ScreenGeometry(center=(100, 50), radius=45),
            profile=profile,
            terrain_profile_altaz=None,
            terrain_secondary_ridges_altaz_layers=[
                [
                    (0.0, 170.0),
                    (0.0, 180.0),
                    (0.0, 190.0),
                ]
            ],
            viewer_data=viewer_data,
            sun_alt_deg=-5.0,
        )
    finally:
        painter.end()

    assert any(
        image.pixelColor(x, y).alpha() > 0
        for x in range(image.width())
        for y in range(image.height())
    )
    assert any(
        (
            (color := image.pixelColor(x, y)).alpha() > 0
            and max(color.red(), color.green(), color.blue()) - min(color.red(), color.green(), color.blue()) <= 20
        )
        for x in range(image.width())
        for y in range(image.height())
    )


def test_draw_night_light_glow_draws_sky_glow() -> None:
    image = QImage(200, 100, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        profile = night_lights.NightLightGlowProfile(
            samples=(
                night_lights.NightLightGlowSample(azimuth_deg=170.0, horizon_alt_deg=0.0, strength=0.4),
                night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
                night_lights.NightLightGlowSample(azimuth_deg=190.0, horizon_alt_deg=0.0, strength=0.4),
            ),
            sun_alt_deg=-5.0,
        )
        viewer_data = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(0.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        )
        draw_night_light_glow(
            painter,
            geometry=ScreenGeometry(center=(100, 50), radius=45),
            profile=profile,
            terrain_profile_altaz=[
                (0.0, 170.0),
                (0.0, 180.0),
                (0.0, 190.0),
            ],
            terrain_secondary_ridges_altaz_layers=None,
            viewer_data=viewer_data,
            sun_alt_deg=-5.0,
        )
    finally:
        painter.end()

    assert any(
        image.pixelColor(x, y).alpha() > 0
        for x in range(90, 111)
        for y in range(45, 56)
    )


def test_draw_night_light_glow_fades_toward_zenith() -> None:
    image = QImage(200, 100, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        profile = night_lights.NightLightGlowProfile(
            samples=(
                night_lights.NightLightGlowSample(azimuth_deg=170.0, horizon_alt_deg=0.0, strength=0.4),
                night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
                night_lights.NightLightGlowSample(azimuth_deg=190.0, horizon_alt_deg=0.0, strength=0.4),
            ),
            sun_alt_deg=-5.0,
            band_half_width_deg=10.0,
        )
        viewer_data = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(0.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        )
        draw_night_light_glow(
            painter,
            geometry=ScreenGeometry(center=(100, 50), radius=45),
            profile=profile,
            terrain_profile_altaz=[
                (0.0, 170.0),
                (0.0, 180.0),
                (0.0, 190.0),
            ],
            terrain_secondary_ridges_altaz_layers=None,
            viewer_data=viewer_data,
            sun_alt_deg=-5.0,
        )
    finally:
        painter.end()

    x = 100
    ys = [y for y in range(image.height()) if image.pixelColor(x, y).alpha() > 0]
    assert ys
    lower = image.pixelColor(x, max(ys))
    upper = image.pixelColor(x, min(ys))
    assert lower.alpha() > upper.alpha()


def test_sky_glow_uses_three_expanding_steps() -> None:
    boundaries = night_lights_render._sky_glow_step_boundaries()
    widths = np.diff(boundaries)
    assert len(widths) == 3
    assert np.isclose(boundaries[0], 0.0)
    assert np.isclose(boundaries[-1], 1.0)
    assert np.allclose(widths, np.asarray([10.0, 20.0, 40.0]) / 70.0)

    alphas = night_lights_render._sky_glow_step_alpha_scales()
    assert len(alphas) == 3
    assert np.isclose(alphas[0], 1.0)
    assert np.all(alphas[1:] < alphas[:-1])
    assert np.isclose(night_lights_render.NIGHT_LIGHTS_STREET_LIGHT_GLOW_ALPHA_BASE, 1.0)
    assert np.isclose(night_lights_render.NIGHT_LIGHTS_SKY_GLOW_ALPHA_SCALE, 0.8)
    assert np.isclose(night_lights_render.NIGHT_LIGHTS_SKY_GLOW_ALPHA_FLOOR, 0.01)
    assert night_lights_render.NIGHT_LIGHTS_SKY_GLOW_RGB == night_lights_render.NIGHT_LIGHTS_GLOW_RGB
    assert night_lights_render.NIGHT_LIGHTS_SKY_GLOW_WIDTH_WEIGHTS == (10.0, 20.0, 40.0)
    assert (
        night_lights_render.NIGHT_LIGHTS_SKY_GLOW_ALPHA_SCALE
        != night_lights_render.NIGHT_LIGHTS_STREET_LIGHT_GLOW_ALPHA_BASE
    )


def test_sky_glow_directional_altitudes_limit_downward_steps() -> None:
    raw_altitudes = [10.0, 9.8, 9.6, 9.4, 9.2]
    boundary8 = night_lights_render._sky_glow_directional_altitudes(raw_altitudes, 8)
    boundary32 = night_lights_render._sky_glow_directional_altitudes(raw_altitudes, 32)

    assert len(boundary8) == len(raw_altitudes)
    assert len(boundary32) == len(raw_altitudes)
    assert np.allclose(boundary8, raw_altitudes)
    assert np.allclose(boundary32, [10.0, 9.95, 9.9, 9.85, 9.8])
    assert all((prev - cur) <= 0.2 + 1.0e-9 for prev, cur in zip(boundary8, boundary8[1:]))
    assert all((prev - cur) <= 0.05 + 1.0e-9 for prev, cur in zip(boundary32, boundary32[1:]))


def test_sky_glow_directional_altitudes_preserve_rises() -> None:
    raw_altitudes = [10.0, 9.8, 10.2, 9.7]
    boundary8 = night_lights_render._sky_glow_directional_altitudes(raw_altitudes, 8)

    assert len(boundary8) == len(raw_altitudes)
    assert boundary8[2] >= raw_altitudes[2]
    assert boundary8[2] >= boundary8[1]
    assert tuple(night_lights_render.NIGHT_LIGHTS_SKY_GLOW_WINDOW_SIZES) == (8, 16, 32)


def test_draw_night_light_glow_respects_opacity() -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(
                azimuth_deg=175.0,
                horizon_alt_deg=0.0,
                strength=0.7,
            ),
            night_lights.NightLightGlowSample(
                azimuth_deg=185.0,
                horizon_alt_deg=0.0,
                strength=1.0,
            ),
        ),
        sun_alt_deg=-5.0,
        band_profiles=(
            night_lights.NightLightDistanceBandProfile(
                min_distance_km=0.5,
                max_distance_km=3.0,
                samples=(
                    night_lights.NightLightGlowSample(
                        azimuth_deg=175.0,
                        horizon_alt_deg=0.0,
                        strength=0.7,
                    ),
                    night_lights.NightLightGlowSample(
                        azimuth_deg=185.0,
                        horizon_alt_deg=0.0,
                        strength=1.0,
                    ),
                ),
            ),
        )
    )

    full = QImage(120, 80, QImage.Format.Format_ARGB32_Premultiplied)
    full.fill(0)
    p_full = QPainter(full)
    try:
        viewer_data = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(0.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        )
        draw_night_light_glow(
            p_full,
            geometry=ScreenGeometry(center=(60, 40), radius=36),
            profile=profile,
            terrain_profile_altaz=None,
            terrain_secondary_ridges_altaz_layers=[[(0.0, 175.0), (0.0, 185.0)]],
            viewer_data=viewer_data,
            opacity=1.0,
            sun_alt_deg=-5.0,
        )
    finally:
        p_full.end()

    dim = QImage(120, 80, QImage.Format.Format_ARGB32_Premultiplied)
    dim.fill(0)
    p_dim = QPainter(dim)
    try:
        viewer_data = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(0.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        )
        draw_night_light_glow(
            p_dim,
            geometry=ScreenGeometry(center=(60, 40), radius=36),
            profile=profile,
            terrain_profile_altaz=None,
            terrain_secondary_ridges_altaz_layers=[[(0.0, 175.0), (0.0, 185.0)]],
            viewer_data=viewer_data,
            opacity=0.25,
            sun_alt_deg=-5.0,
        )
    finally:
        p_dim.end()

    independent = QImage(120, 80, QImage.Format.Format_ARGB32_Premultiplied)
    independent.fill(0)
    p_independent = QPainter(independent)
    try:
        viewer_data = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(0.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        )
        draw_night_light_glow(
            p_independent,
            geometry=ScreenGeometry(center=(60, 40), radius=36),
            profile=profile,
            terrain_profile_altaz=[
                (0.0, 175.0),
                (0.0, 180.0),
                (0.0, 185.0),
            ],
            terrain_secondary_ridges_altaz_layers=None,
            viewer_data=viewer_data,
            opacity=0.0,
            sky_glow_opacity=0.25,
            sun_alt_deg=-5.0,
        )
    finally:
        p_independent.end()

    zero = QImage(120, 80, QImage.Format.Format_ARGB32_Premultiplied)
    zero.fill(0)
    p_zero = QPainter(zero)
    try:
        viewer_data = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(0.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
        )
        draw_night_light_glow(
            p_zero,
            geometry=ScreenGeometry(center=(60, 40), radius=36),
            profile=profile,
            terrain_profile_altaz=None,
            terrain_secondary_ridges_altaz_layers=[[(0.0, 175.0), (0.0, 185.0)]],
            viewer_data=viewer_data,
            opacity=0.0,
            sky_glow_opacity=0.0,
            sun_alt_deg=-5.0,
        )
    finally:
        p_zero.end()

    assert any(full.pixelColor(x, y).alpha() > 0 for x in range(full.width()) for y in range(full.height()))
    assert any(dim.pixelColor(x, y).alpha() > 0 for x in range(dim.width()) for y in range(dim.height()))
    assert any(independent.pixelColor(x, y).alpha() > 0 for x in range(independent.width()) for y in range(independent.height()))
    assert not any(zero.pixelColor(x, y).alpha() > 0 for x in range(zero.width()) for y in range(zero.height()))


def test_draw_night_light_glow_does_not_draw_without_secondary_ridges() -> None:
    image = QImage(200, 100, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        profile = night_lights.NightLightGlowProfile(
            samples=(
                night_lights.NightLightGlowSample(azimuth_deg=170.0, horizon_alt_deg=0.0, strength=0.4),
                night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
                night_lights.NightLightGlowSample(azimuth_deg=190.0, horizon_alt_deg=0.0, strength=0.4),
            ),
            sun_alt_deg=-5.0,
        )
        draw_night_light_glow(
            painter,
            geometry=ScreenGeometry(center=(100, 50), radius=45),
            profile=profile,
            terrain_profile_altaz=None,
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(0.0, 180.0),
                edge_fov_deg=95.0,
                content_fov_deg=110.0,
            ),
            sun_alt_deg=-5.0,
        )
    finally:
        painter.end()

    assert not any(
        image.pixelColor(x, y).alpha() > 0
        for x in range(image.width())
        for y in range(image.height())
    )
