from __future__ import annotations

import numpy as np
import pytest
from PySide6.QtCore import Qt
from PySide6.QtGui import QImage, QPainter

import zstarview.gui.composite as render_composite
from zstarview import night_lights as night_lights_module
from zstarview.gui.composite import (
    NEVER_RISES_GUIDE_WIDTH_SCALE,
    SkyCompositorCache,
    mask_cloud_alpha_by_missing,
    overlay_missing_tint,
)
from zstarview.render.earth_guide import earth_guide_line_alpha
from zstarview.render.guides import (
    REFERENCE_LINE_FG_WIDTH,
    REFERENCE_LINE_MID_ALPHA,
    REFERENCE_LINE_MID_WIDTH,
    REFERENCE_LINE_OUTER_ALPHA,
    REFERENCE_LINE_OUTER_WIDTH,
)
from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba
from zstarview.types import ScreenGeometry


def test_overlay_missing_with_hatch_tints_only_missing_region() -> None:
    base = np.zeros((16, 16, 4), dtype=np.uint8)
    base[..., :3] = 30
    base[..., 3] = 255

    missing = np.zeros((16, 16), dtype=np.uint8)
    missing[4:12, 4:12] = 255

    out = qimage_to_np_rgba(
        overlay_missing_tint(
            np_rgba_to_qimage(base),
            missing,
            tint_rgba=(255, 220, 80, 90),
        )
    )

    assert int(out[1, 1, 0]) == 30
    assert int(out[8, 8, 0]) > 30


def test_compositor_cache_key_includes_missing_mask() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 20
    sky[..., 3] = 255
    cloud = np.zeros((64, 64, 4), dtype=np.uint8)
    cloud[..., :3] = 255
    cloud[..., 3] = 80

    missing_none = np.zeros((64, 64), dtype=np.uint8)
    missing_half = np.zeros((64, 64), dtype=np.uint8)
    missing_half[:, :32] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache(ground_tint_opacity=1.0)

    canvas1 = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas1.fill(0)
    p1 = QPainter(canvas1)
    compositor.draw(
        p1,
        geom,
        np_rgba_to_qimage(sky),
        np_rgba_to_qimage(cloud),
        cloud_alpha=0.4,
        view_center=(0.0, 0.0),
        cloud_amount_field=None,
        missing_mask=missing_none,
        content_fov_deg=90.0,
    )
    p1.end()

    canvas2 = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas2.fill(0)
    p2 = QPainter(canvas2)
    compositor.draw(
        p2,
        geom,
        np_rgba_to_qimage(sky),
        np_rgba_to_qimage(cloud),
        cloud_alpha=0.4,
        view_center=(0.0, 0.0),
        cloud_amount_field=None,
        missing_mask=missing_half,
        content_fov_deg=90.0,
    )
    p2.end()

    arr1 = qimage_to_np_rgba(canvas1)
    arr2 = qimage_to_np_rgba(canvas2)
    assert np.any(arr1 != arr2)


def test_mask_cloud_alpha_by_missing_cuts_cloud_pixels() -> None:
    cloud = np.zeros((12, 12, 4), dtype=np.uint8)
    cloud[..., :3] = 255
    cloud[..., 3] = 180
    missing = np.zeros((12, 12), dtype=np.uint8)
    missing[3:9, 4:10] = 255

    out = qimage_to_np_rgba(mask_cloud_alpha_by_missing(np_rgba_to_qimage(cloud), missing))
    assert int(out[1, 1, 3]) == 180
    assert int(out[5, 6, 3]) == 0


def test_compositor_terrain_profile_does_not_add_ground_fill() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache(ground_tint_opacity=1.0)
    terrain_profile = [(30.0, float(az)) for az in range(360)]

    canvas_flat = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas_flat.fill(0)
    p_flat = QPainter(canvas_flat)
    compositor.draw(
        p_flat,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 0.0),
        content_fov_deg=90.0,
    )
    p_flat.end()

    canvas_terrain = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas_terrain.fill(0)
    p_terrain = QPainter(canvas_terrain)
    compositor.draw(
        p_terrain,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 0.0),
        terrain_profile_altaz=terrain_profile,
        content_fov_deg=90.0,
    )
    p_terrain.end()

    arr_flat = qimage_to_np_rgba(canvas_flat)
    arr_terrain = qimage_to_np_rgba(canvas_terrain)

    assert np.array_equal(arr_flat[24, 32, :3], np.array([100, 100, 100], dtype=np.uint8))
    assert np.array_equal(arr_terrain[24, 32, :3], np.array([100, 100, 100], dtype=np.uint8))


def test_compositor_fast_mode_matches_normal_mode_without_ground_fill() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache(ground_tint_opacity=1.0)
    terrain_profile = [(30.0, float(az)) for az in range(360)]

    canvas_normal = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas_normal.fill(0)
    painter_normal = QPainter(canvas_normal)
    compositor.draw(
        painter_normal,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        observer_lat_deg=35.0,
        terrain_profile_altaz=terrain_profile,
        content_fov_deg=90.0,
        fast_mode=False,
    )
    painter_normal.end()

    canvas_fast = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas_fast.fill(0)
    painter_fast = QPainter(canvas_fast)
    compositor.draw(
        painter_fast,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        observer_lat_deg=35.0,
        terrain_profile_altaz=terrain_profile,
        content_fov_deg=90.0,
        fast_mode=True,
    )
    painter_fast.end()

    arr_normal = qimage_to_np_rgba(canvas_normal)
    arr_fast = qimage_to_np_rgba(canvas_fast)
    assert np.array_equal(arr_normal, arr_fast)


def test_compositor_fast_mode_skips_night_light_overlay(monkeypatch) -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache(ground_tint_opacity=1.0)
    night_profile = night_lights_module.NightLightGlowProfile(
        samples=(
            night_lights_module.NightLightGlowSample(
                azimuth_deg=180.0,
                horizon_alt_deg=0.0,
                strength=1.0,
            ),
        ),
        sun_alt_deg=-5.0,
    )

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        observer_lat_deg=35.0,
        terrain_profile_altaz=[(0.0, float(az)) for az in range(360)],
        night_light_glow_profile=night_profile,
        content_fov_deg=90.0,
        fast_mode=True,
    )
    painter.end()


def test_compositor_ground_reset_replaces_lower_disc_with_background() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache(ground_tint_opacity=1.0)
    terrain_profile = [(45.0, float(az)) for az in range(360)]

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 0.0),
        terrain_profile_altaz=terrain_profile,
        ground_reset_rgba=(12, 34, 56, 255),
        content_fov_deg=90.0,
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    assert np.array_equal(arr[32, 32, :3], np.array([12, 34, 56], dtype=np.uint8))
    assert int(arr[32, 32, 3]) == 255


def test_compositor_ground_tint_opacity_no_longer_changes_output() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache(ground_tint_opacity=0.0)
    terrain_profile = [(30.0, float(az)) for az in range(360)]

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 0.0),
        terrain_profile_altaz=terrain_profile,
        content_fov_deg=90.0,
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    assert np.array_equal(arr[24, 32, :3], np.array([100, 100, 100], dtype=np.uint8))


def test_compositor_observer_latitude_draws_never_rises_outline() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache(ground_tint_opacity=1.0)

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        observer_lat_deg=35.0,
        content_fov_deg=90.0,
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    rgb = arr[..., :3]
    outline_mask = np.any(rgb != np.array([100, 100, 100], dtype=np.uint8)[None, None, :], axis=2)
    assert np.any(outline_mask)
    changed_rgb = rgb[outline_mask]
    assert np.all(changed_rgb[:, 0] >= changed_rgb[:, 1])
    assert np.all(changed_rgb[:, 1] >= changed_rgb[:, 2])
    assert np.all(arr[..., 3][outline_mask] == 255)


def test_compositor_guidelines_toggle_hides_never_rises_outline() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache(ground_tint_opacity=1.0)

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 180.0),
        observer_lat_deg=35.0,
        show_guidelines=False,
        earth_guide_opacity=0.0,
        content_fov_deg=90.0,
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    assert np.array_equal(arr[..., :3], sky[..., :3])


def test_never_rises_outline_uses_double_width_scale(monkeypatch) -> None:
    base = QImage(16, 16, QImage.Format_ARGB32_Premultiplied)
    base.fill(0)
    geom = ScreenGeometry(center=(8, 8), radius=8)

    class _DummyColor:
        def __init__(self, alpha: float) -> None:
            self._alpha = alpha

        def alphaF(self) -> float:
            return self._alpha

    class _DummyPen:
        def __init__(self, width: float, alpha: float) -> None:
            self._width = width
            self._color = _DummyColor(alpha)

        def widthF(self) -> float:
            return self._width

        def color(self) -> _DummyColor:
            return self._color

    class _DummyPainter:
        def __init__(self) -> None:
            self.pen_widths: list[float] = []
            self.pen_styles: list[object] = []
            self.pen_alphas: list[float] = []

        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def setRenderHint(self, *_args, **_kwargs) -> None:  # noqa: N802 - Qt naming
            pass

        def setPen(self, pen, *_args, **_kwargs) -> None:
            self.pen_widths.append(float(pen.widthF()))
            self.pen_styles.append(pen.style())
            self.pen_alphas.append(float(pen.color().alphaF()))

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

        def end(self) -> None:
            pass

    class _DummyQPainter:
        class RenderHint:
            Antialiasing = object()

        def __new__(cls, *_args, **_kwargs):
            return dummy_painter

    dummy_painter = _DummyPainter()
    monkeypatch.setattr(render_composite, "_never_rises_circle_altaz", lambda _lat_deg: [(0.0, 0.0), (0.0, 10.0)])
    monkeypatch.setattr(render_composite, "split_by_gaps", lambda projected: [projected])
    monkeypatch.setattr(render_composite, "_clip_polyline_to_radius", lambda fragment, _ratio: [fragment])
    monkeypatch.setattr(render_composite, "normalized_to_screen_xy", lambda x, y, _geometry: (x, y))
    monkeypatch.setattr(render_composite, "QPainter", _DummyQPainter)

    render_composite._overlay_never_rises_outline(
        base,
        geometry=geom,
        view_center=(0.0, 180.0),
        observer_lat_deg=35.0,
        never_rises_opacity=0.2,
        content_fov_deg=90.0,
    )

    assert NEVER_RISES_GUIDE_WIDTH_SCALE == pytest.approx(4.5)
    assert dummy_painter.pen_widths == [
        pytest.approx(REFERENCE_LINE_OUTER_WIDTH * NEVER_RISES_GUIDE_WIDTH_SCALE),
        pytest.approx(REFERENCE_LINE_MID_WIDTH * NEVER_RISES_GUIDE_WIDTH_SCALE),
        pytest.approx(REFERENCE_LINE_FG_WIDTH * NEVER_RISES_GUIDE_WIDTH_SCALE),
    ]
    assert dummy_painter.pen_alphas == [
        pytest.approx((REFERENCE_LINE_OUTER_ALPHA / 255.0) * 0.5),
        pytest.approx((REFERENCE_LINE_MID_ALPHA / 255.0) * 0.5),
        pytest.approx(round(255.0 * earth_guide_line_alpha(0.2) * 0.5) / 255.0),
    ]
    assert dummy_painter.pen_styles == [
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
        Qt.PenStyle.SolidLine,
    ]
