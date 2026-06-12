from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage, QPainter
from PySide6.QtWidgets import QApplication

from zstarview import night_lights
from zstarview.render.ridge_glow import (
    RIDGE_GLOW_SKY_ALPHA_FLOOR,
    RIDGE_GLOW_SKY_ALPHA_SCALE,
    RIDGE_GLOW_SKY_LAYER_SPECS,
    RIDGE_GLOW_SKY_RGB,
    _ridge_glow_directional_altitudes,
    draw_ridge_glow_normal,
)
from zstarview.types import ScreenGeometry, ViewerData

app = QApplication.instance() or QApplication([])


def test_ridge_glow_uses_five_expanding_steps() -> None:
    assert np.isclose(RIDGE_GLOW_SKY_ALPHA_SCALE, 0.85)
    assert np.isclose(RIDGE_GLOW_SKY_ALPHA_FLOOR, 0.02)
    assert RIDGE_GLOW_SKY_RGB == (255, 170, 48)
    assert len(RIDGE_GLOW_SKY_LAYER_SPECS) == 5
    assert [spec.upper_alt_offset_deg for spec in RIDGE_GLOW_SKY_LAYER_SPECS] == [0.3, 0.7, 1.6, 3.5, 7.4]
    assert [spec.alpha_scale for spec in RIDGE_GLOW_SKY_LAYER_SPECS] == [1.0, 0.5, 0.25, 0.125, 0.12]
    assert [spec.focus_scale for spec in RIDGE_GLOW_SKY_LAYER_SPECS] == [1.8, 1.6, 1.4, 1.2, 1.1]
    assert [spec.window_size for spec in RIDGE_GLOW_SKY_LAYER_SPECS] == [2, 4, 8, 16, 32]
    alpha_scales = np.asarray([spec.alpha_scale for spec in RIDGE_GLOW_SKY_LAYER_SPECS], dtype=np.float64)
    assert len(alpha_scales) == 5
    assert np.allclose(alpha_scales, np.asarray([1.0, 0.5, 0.25, 0.125, 0.12]), atol=0.01)
    assert np.all(alpha_scales[1:] < alpha_scales[:-1])
    width_offsets = np.asarray([spec.upper_alt_offset_deg for spec in RIDGE_GLOW_SKY_LAYER_SPECS], dtype=np.float64)
    assert len(width_offsets) == 5
    assert np.allclose(width_offsets, np.asarray([0.3, 0.7, 1.6, 3.5, 7.4]), atol=0.01)
    focus_scales = np.asarray([spec.focus_scale for spec in RIDGE_GLOW_SKY_LAYER_SPECS], dtype=np.float64)
    assert len(focus_scales) == 5
    assert np.allclose(focus_scales, np.asarray([1.8, 1.6, 1.4, 1.2, 1.1]), atol=0.01)
    assert np.all(focus_scales[:-1] >= focus_scales[1:])


def test_ridge_glow_directional_altitudes_use_center_weighted_average() -> None:
    raw_altitudes = [1.0, 3.0, 2.0, 5.0, 4.0, 2.0, 1.0]
    boundary3 = _ridge_glow_directional_altitudes(raw_altitudes, 3)
    boundary5 = _ridge_glow_directional_altitudes(raw_altitudes, 5)

    assert len(boundary3) == len(raw_altitudes)
    assert len(boundary5) == len(raw_altitudes)
    assert np.allclose(
        boundary3,
        [1.6666666667, 2.25, 3.0, 4.0, 3.75, 2.25, 1.3333333333],
    )
    assert np.allclose(
        boundary5,
        [1.8333333333, 2.5, 3.0, 3.5555555556, 3.2222222222, 2.625, 1.8333333333],
    )
    assert boundary3[2] > boundary3[1]
    assert boundary5[3] > boundary5[2]
    assert tuple(spec.window_size for spec in RIDGE_GLOW_SKY_LAYER_SPECS) == (2, 4, 8, 16, 32)


def test_ridge_glow_directional_altitudes_preserve_rises() -> None:
    raw_altitudes = [10.0, 9.8, 10.2, 9.7]
    boundary3 = _ridge_glow_directional_altitudes(raw_altitudes, 3)

    assert len(boundary3) == len(raw_altitudes)
    assert boundary3[2] >= boundary3[1]
    assert boundary3[2] >= boundary3[3]
    assert tuple(spec.window_size for spec in RIDGE_GLOW_SKY_LAYER_SPECS) == (2, 4, 8, 16, 32)


def test_draw_ridge_glow_smoke() -> None:
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
        draw_ridge_glow_normal(
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


def test_draw_ridge_glow_fades_toward_zenith() -> None:
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
        draw_ridge_glow_normal(
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


def test_draw_ridge_glow_respects_opacity() -> None:
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
        ),
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
        draw_ridge_glow_normal(
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
        draw_ridge_glow_normal(
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
        draw_ridge_glow_normal(
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
            ridge_glow_opacity=0.25,
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
        draw_ridge_glow_normal(
            p_zero,
            geometry=ScreenGeometry(center=(60, 40), radius=36),
            profile=profile,
            terrain_profile_altaz=None,
            terrain_secondary_ridges_altaz_layers=[[(0.0, 175.0), (0.0, 185.0)]],
            viewer_data=viewer_data,
            opacity=0.0,
            ridge_glow_opacity=0.0,
            sun_alt_deg=-5.0,
        )
    finally:
        p_zero.end()

    assert any(full.pixelColor(x, y).alpha() > 0 for x in range(full.width()) for y in range(full.height()))
    assert any(dim.pixelColor(x, y).alpha() > 0 for x in range(dim.width()) for y in range(dim.height()))
    assert any(independent.pixelColor(x, y).alpha() > 0 for x in range(independent.width()) for y in range(independent.height()))
    assert not any(zero.pixelColor(x, y).alpha() > 0 for x in range(zero.width()) for y in range(zero.height()))


def test_draw_ridge_glow_does_not_draw_without_secondary_ridges() -> None:
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
        draw_ridge_glow_normal(
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
