from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication

import zstarview.gui.composite as render_composite
from zstarview import night_lights
from zstarview.types import ScreenGeometry

app = QApplication.instance() or QApplication([])


def test_glow_mask_to_qimage_uses_bright_base_color() -> None:
    mask = render_composite.GlowMask(
        alpha=np.asarray([[0.0, 0.5], [1.0, 0.25]], dtype=np.float32),
        scale=0.25,
    )

    image = render_composite._glow_mask_to_qimage(mask, (200, 100, 50))

    assert image.format() == QImage.Format.Format_ARGB32_Premultiplied
    center = image.pixelColor(0, 1)
    assert center.alpha() > 0
    assert center.red() >= center.green() >= center.blue()
    assert center.red() == 255
    assert center.green() == 128
    assert center.blue() == 64


def test_glow_mask_to_qimage_applies_stable_noise() -> None:
    mask = render_composite.GlowMask(
        alpha=np.full((6, 6), 0.5, dtype=np.float32),
        scale=0.25,
    )

    image1 = render_composite._glow_mask_to_qimage(mask, (200, 100, 50))
    image2 = render_composite._glow_mask_to_qimage(mask, (200, 100, 50))

    alpha_grid1 = np.asarray(
        [[image1.pixelColor(x, y).alpha() for x in range(image1.width())] for y in range(image1.height())],
        dtype=np.uint8,
    )
    alpha_grid2 = np.asarray(
        [[image2.pixelColor(x, y).alpha() for x in range(image2.width())] for y in range(image2.height())],
        dtype=np.uint8,
    )

    assert np.array_equal(alpha_grid1, alpha_grid2)
    assert len(np.unique(alpha_grid1)) > 1


def test_night_light_ray_alpha_field_decays_above_horizon() -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=-20.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )

    alpha = render_composite._night_light_ray_alpha_field(
        profile=profile,
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        view_center=(0.0, 180.0),
        terrain_profile_altaz=None,
        opacity=0.5,
        sun_alt_deg=-5.0,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
    )

    assert alpha.shape == (80, 80)
    assert alpha.dtype == np.float32
    assert np.any(alpha > 0.0)
    assert np.all((alpha >= 0.0) & (alpha <= 1.0))
    center = alpha[40, 40]
    upper = alpha[20, 40]
    assert upper < center


def test_night_light_ray_alpha_field_is_soft_below_horizon() -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )

    alpha = render_composite._night_light_ray_alpha_field(
        profile=profile,
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        view_center=(0.0, 180.0),
        terrain_profile_altaz=[(10.0, 180.0)],
        opacity=0.5,
        sun_alt_deg=-5.0,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
    )

    assert alpha[40, 40] > 0.0


def test_ridge_glow_height_is_lower_than_night_light_height() -> None:
    assert render_composite.GLOW_MASK_RIDGE_GLOW_HEIGHT_DEG < render_composite.GLOW_MASK_NIGHT_LIGHT_HEIGHT_DEG


def test_cumulative_max_ridge_altitude_uses_previous_layers() -> None:
    azimuths = np.asarray([170.0, 180.0, 190.0], dtype=np.float32)
    layers = [
        [(1.0, 170.0), (2.0, 180.0), (1.0, 190.0)],
        [(3.0, 170.0), (1.0, 180.0), (4.0, 190.0)],
        [(2.0, 170.0), (5.0, 180.0), (3.0, 190.0)],
    ]

    cumulative = render_composite._cumulative_max_ridge_altitude(layers, azimuths)

    assert np.allclose(cumulative, np.asarray([3.0, 5.0, 4.0], dtype=np.float32))


def test_build_glow_mask_uses_ray_alpha_field(monkeypatch) -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=-20.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )

    mask = render_composite._build_glow_mask(
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        view_center=(0.0, 180.0),
        terrain_profile_altaz=None,
        terrain_secondary_ridges_altaz_layers=None,
        night_light_glow_profile=profile,
        night_light_opacity=0.5,
        ridge_glow_opacity=0.0,
        night_light_sun_alt_deg=-5.0,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
    )

    assert mask is not None
    assert mask.alpha.shape == (20, 20)
    assert np.any(mask.alpha > 0.0)


def test_build_glow_mask_routes_secondary_layers_into_night_light_mask(monkeypatch) -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=-20.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )
    observed: dict[str, object] = {}

    def fake_night_light_ray_alpha_field(**kwargs) -> np.ndarray:
        observed["secondary_layers"] = kwargs.get("terrain_secondary_ridges_altaz_layers")
        return np.full((20, 20), 0.5, dtype=np.float32)

    monkeypatch.setattr(render_composite, "_night_light_ray_alpha_field", fake_night_light_ray_alpha_field)

    mask = render_composite._build_glow_mask(
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        view_center=(0.0, 180.0),
        terrain_profile_altaz=[(0.0, 180.0)],
        terrain_secondary_ridges_altaz_layers=[
            [(1.0, 180.0)],
            [(2.0, 180.0)],
        ],
        night_light_glow_profile=profile,
        night_light_opacity=0.5,
        ridge_glow_opacity=0.0,
        night_light_sun_alt_deg=-5.0,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
    )

    assert mask is not None
    assert observed["secondary_layers"] == [[(1.0, 180.0)], [(2.0, 180.0)]]


def test_build_glow_mask_skips_fast_mode(monkeypatch) -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=0.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )

    mask = render_composite._build_glow_mask(
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        view_center=(0.0, 180.0),
        terrain_profile_altaz=[(0.0, 180.0)],
        terrain_secondary_ridges_altaz_layers=[[(0.0, 180.0)]],
        night_light_glow_profile=profile,
        night_light_opacity=0.5,
        ridge_glow_opacity=0.5,
        night_light_sun_alt_deg=-5.0,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        fast_mode=True,
    )

    assert mask is None
