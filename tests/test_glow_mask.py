from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage, QPainter
from PySide6.QtWidgets import QApplication

import zstarview.gui.composite as render_composite
from zstarview import night_lights
from zstarview.types import ScreenGeometry, ViewerData

app = QApplication.instance() or QApplication([])


def _make_viewer_data(
    *,
    view_center: tuple[float, float] = (0.0, 180.0),
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = 90.0,
) -> ViewerData:
    return ViewerData(
        location=(35.0, 135.0),
        timezone_name="UTC",
        city_name="",
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
        observer_height_m=1.7,
    )


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


def test_crop_night_light_alpha_grid_altitude_bins_trims_inactive_rows() -> None:
    altitude_bins = np.asarray([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0], dtype=np.float64)
    alpha_grid = np.asarray(
        [
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.75],
            [0.0, 0.0],
            [0.0, 0.0],
        ],
        dtype=np.float64,
    )

    got_altitude_bins, got_alpha_grid = render_composite._crop_night_light_alpha_grid_altitude_bins(
        altitude_bins,
        alpha_grid,
    )

    assert np.allclose(got_altitude_bins, np.asarray([0.0, 1.0, 2.0], dtype=np.float64))
    assert got_alpha_grid.shape == (3, 2)
    assert np.allclose(got_alpha_grid[1], np.asarray([0.0, 0.75], dtype=np.float64))


def test_night_light_ray_alpha_field_decays_above_horizon() -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=-20.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )

    alpha = render_composite._night_light_ray_alpha_field(
        profile=profile,
        viewer_data=_make_viewer_data(),
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        opacity=0.5,
        sun_alt_deg=-5.0,
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
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=10.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )

    alpha = render_composite._night_light_ray_alpha_field(
        profile=profile,
        viewer_data=_make_viewer_data(),
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        opacity=0.5,
        sun_alt_deg=-5.0,
    )

    assert alpha[40, 40] > 0.0


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
        viewer_data=_make_viewer_data(),
        night_light_glow_profile=profile,
        night_light_opacity=0.5,
        night_light_sun_alt_deg=-5.0,
    )

    assert mask is not None
    assert mask.alpha.shape == (20, 20)
    assert np.any(mask.alpha > 0.0)


def test_build_glow_mask_uses_profile_only_for_night_light_mask(monkeypatch) -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=-20.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )
    observed: dict[str, object] = {}

    def fake_night_light_ray_alpha_field(**kwargs) -> np.ndarray:
        observed["kwargs"] = kwargs
        return np.full((20, 20), 0.5, dtype=np.float32)

    monkeypatch.setattr(render_composite, "_night_light_ray_alpha_field", fake_night_light_ray_alpha_field)

    mask = render_composite._build_glow_mask(
        width=80,
        height=80,
        geometry=ScreenGeometry(center=(40, 40), radius=36),
        viewer_data=_make_viewer_data(),
        night_light_glow_profile=profile,
        night_light_opacity=0.5,
        night_light_sun_alt_deg=-5.0,
    )

    assert mask is not None
    assert isinstance(observed["kwargs"]["viewer_data"], ViewerData)
    assert tuple(observed["kwargs"]["viewer_data"].view_center) == (0.0, 180.0)
    assert observed["kwargs"]["alpha_grid"] == ()


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
        viewer_data=_make_viewer_data(),
        night_light_glow_profile=profile,
        night_light_opacity=0.5,
        night_light_sun_alt_deg=-5.0,
        fast_mode=True,
    )

    assert mask is None


def test_compositor_reuses_cached_glow_mask_across_unrelated_frame_changes(monkeypatch) -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=-10.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )
    calls: list[int] = []

    def fake_build_glow_mask(**kwargs):
        calls.append(int(kwargs["width"]))
        return render_composite.GlowMask(
            alpha=np.full((8, 8), 0.5, dtype=np.float32),
            scale=0.25,
        )

    monkeypatch.setattr(render_composite, "_build_glow_mask", fake_build_glow_mask)

    compositor = render_composite.SkyCompositorCache()
    geom = ScreenGeometry(center=(16, 16), radius=16)
    sky1 = np.zeros((32, 32, 4), dtype=np.uint8)
    sky1[..., :3] = 80
    sky1[..., 3] = 255
    sky2 = np.zeros((32, 32, 4), dtype=np.uint8)
    sky2[..., :3] = 120
    sky2[..., 3] = 255

    canvas1 = QImage(32, 32, QImage.Format_ARGB32_Premultiplied)
    canvas1.fill(0)
    painter1 = QPainter(canvas1)
    compositor.draw(
        painter1,
        geom,
        render_composite.np_rgba_to_qimage(sky1),
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        terrain_profile_altaz=[(0.0, 180.0)],
        night_light_glow_profile=profile,
        night_light_opacity=0.2,
        night_light_sun_alt_deg=-5.0,
        content_fov_deg=90.0,
    )
    painter1.end()

    canvas2 = QImage(32, 32, QImage.Format_ARGB32_Premultiplied)
    canvas2.fill(0)
    painter2 = QPainter(canvas2)
    compositor.draw(
        painter2,
        geom,
        render_composite.np_rgba_to_qimage(sky2),
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        terrain_profile_altaz=[(0.0, 180.0)],
        night_light_glow_profile=profile,
        night_light_opacity=0.2,
        night_light_sun_alt_deg=-5.0,
        content_fov_deg=90.0,
    )
    painter2.end()

    assert len(calls) == 1


def test_compositor_builds_edge_glow_mask_separately(monkeypatch) -> None:
    calls: list[str] = []
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=-10.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )

    def fake_build_glow_mask(**kwargs):
        calls.append("night")
        return render_composite.GlowMask(
            alpha=np.full((8, 8), 0.5, dtype=np.float32),
            scale=0.25,
        )

    def fake_build_edge_glow_mask(**kwargs):
        calls.append("edge")
        return render_composite.GlowMask(
            alpha=np.full((8, 8), 0.25, dtype=np.float32),
            scale=0.25,
        )

    monkeypatch.setattr(render_composite, "_build_glow_mask", fake_build_glow_mask)
    monkeypatch.setattr(render_composite, "_build_edge_glow_mask", fake_build_edge_glow_mask)

    compositor = render_composite.SkyCompositorCache()
    geom = ScreenGeometry(center=(16, 16), radius=16)
    sky = np.zeros((32, 32, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    canvas = QImage(32, 32, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        render_composite.np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        terrain_profile_altaz=[(0.0, 180.0)],
        night_light_glow_profile=profile,
        night_light_opacity=0.2,
        night_light_sun_alt_deg=-5.0,
        content_fov_deg=90.0,
    )
    painter.end()

    assert calls == ["night", "edge"]


def test_compositor_passes_ridge_glow_opacity_to_edge_mask(monkeypatch) -> None:
    profile = night_lights.NightLightGlowProfile(
        samples=(
            night_lights.NightLightGlowSample(azimuth_deg=180.0, horizon_alt_deg=-10.0, strength=1.0),
        ),
        sun_alt_deg=-5.0,
    )
    captured: dict[str, float] = {}

    def fake_edge_ray_alpha_field(**kwargs):
        captured["opacity"] = float(kwargs["opacity"])
        return np.full((8, 8), 0.25, dtype=np.float32)

    monkeypatch.setattr(render_composite, "_night_light_edge_ray_alpha_field", fake_edge_ray_alpha_field)

    compositor = render_composite.SkyCompositorCache()
    geom = ScreenGeometry(center=(16, 16), radius=16)
    sky = np.zeros((32, 32, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    canvas = QImage(32, 32, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        render_composite.np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        terrain_profile_altaz=[(0.0, 180.0)],
        night_light_glow_profile=profile,
        night_light_opacity=0.2,
        ridge_glow_opacity=0.8,
        night_light_sun_alt_deg=-5.0,
        content_fov_deg=90.0,
    )
    painter.end()

    assert captured["opacity"] == 0.8
