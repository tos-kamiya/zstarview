from __future__ import annotations

import numpy as np
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
            terrain_profile_altaz=[(0.0, 175.0), (0.0, 185.0)],
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
            terrain_profile_altaz=[(0.0, 175.0), (0.0, 185.0)],
            viewer_data=viewer_data,
            opacity=0.25,
            sun_alt_deg=-5.0,
        )
    finally:
        p_dim.end()

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
            terrain_profile_altaz=[(0.0, 175.0), (0.0, 185.0)],
            viewer_data=viewer_data,
            opacity=0.0,
            sun_alt_deg=-5.0,
        )
    finally:
        p_zero.end()

    assert any(full.pixelColor(x, y).alpha() > 0 for x in range(full.width()) for y in range(full.height()))
    assert any(dim.pixelColor(x, y).alpha() > 0 for x in range(dim.width()) for y in range(dim.height()))
    assert not any(zero.pixelColor(x, y).alpha() > 0 for x in range(zero.width()) for y in range(zero.height()))


def test_draw_night_light_glow_uses_flat_horizon_when_missing_terrain() -> None:
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

    assert any(
        image.pixelColor(x, y).alpha() > 0
        for x in range(image.width())
        for y in range(image.height())
    )
