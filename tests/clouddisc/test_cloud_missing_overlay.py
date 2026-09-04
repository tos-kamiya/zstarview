from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pytest
from PySide6.QtCore import Qt
from PySide6.QtGui import QImage, QPainter

import zstarview.gui.composite as render_composite
from zstarview import night_lights as night_lights_module
from zstarview.clouddisc.altaz_grid import CloudAltAzGrid
from zstarview.clouddisc.types import SourceKey
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
from zstarview.types import ScreenGeometry, ViewerData, ViewProjection


def _viewer_data(
    view_center: tuple[float, float] = (0.0, 0.0),
    latitude: float = 0.0,
) -> ViewerData:
    return ViewerData(
        location=(latitude, 0.0),
        timezone_name="UTC",
        city_name="",
        view_center=view_center,
        edge_fov_deg=90.0,
        content_fov_deg=90.0,
        observer_height_m=0.0,
    )


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

    missing_none = np.zeros((64, 64), dtype=np.uint8)
    missing_half = np.zeros((64, 64), dtype=np.uint8)
    missing_half[:, :32] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()

    canvas1 = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas1.fill(0)
    p1 = QPainter(canvas1)
    compositor.draw(
        p1,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.4,
        viewer_data=_viewer_data(),
        missing_mask=missing_none,
    )
    p1.end()

    canvas2 = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas2.fill(0)
    p2 = QPainter(canvas2)
    compositor.draw(
        p2,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.4,
        viewer_data=_viewer_data(),
        missing_mask=missing_half,
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


def test_compositor_clips_sky_layers_below_terrain_horizon() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()
    terrain_profile = [(30.0, float(az)) for az in range(360)]

    canvas_flat = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas_flat.fill(0)
    p_flat = QPainter(canvas_flat)
    compositor.draw(
        p_flat,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        viewer_data=_viewer_data(),
    )
    p_flat.end()

    canvas_terrain = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas_terrain.fill(0)
    p_terrain = QPainter(canvas_terrain)
    compositor.draw(
        p_terrain,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        viewer_data=_viewer_data(),
        terrain_profile_altaz=terrain_profile,
    )
    p_terrain.end()

    arr_flat = qimage_to_np_rgba(canvas_flat)
    arr_terrain = qimage_to_np_rgba(canvas_terrain)

    assert np.array_equal(arr_flat[24, 32, :3], np.array([100, 100, 100], dtype=np.uint8))
    assert np.array_equal(arr_terrain[24, 32, :3], np.array([0, 0, 0], dtype=np.uint8))


def test_terrain_clip_inverse_projects_only_ridge_guard(monkeypatch) -> None:
    width, height = 512, 288
    sky = np.full((height, width, 4), 160, dtype=np.uint8)
    sky[..., 3] = 255
    geometry = ScreenGeometry(center=(width // 2, height // 2), radius=height // 2)
    terrain_profile = [(12.0, float(az)) for az in range(360)]
    selected_pixel_counts: list[int] = []
    original = render_composite._inverse_project_pixel_coordinates

    def record_selected_pixels(x, y, **kwargs):
        selected_pixel_counts.append(int(np.asarray(x).size))
        return original(x, y, **kwargs)

    monkeypatch.setattr(
        render_composite,
        "_inverse_project_pixel_coordinates",
        record_selected_pixels,
    )

    clipped = render_composite._clip_below_terrain_horizon(
        np_rgba_to_qimage(sky),
        geometry=geometry,
        projection=ViewProjection(
            view_center=(0.0, 180.0),
            edge_fov_deg=90.0,
            content_fov_deg=90.75,
        ),
        terrain_profile_altaz=terrain_profile,
    )

    assert selected_pixel_counts
    assert selected_pixel_counts[0] < width * height // 2
    clipped_rgba = qimage_to_np_rgba(clipped)
    assert int(clipped_rgba[height - 1, width // 2, 3]) == 0
    assert int(clipped_rgba[height // 4, width // 2, 3]) == 255


def test_compositor_fast_mode_matches_normal_mode_without_ground_fill() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()
    terrain_profile = [(30.0, float(az)) for az in range(360)]

    canvas_normal = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas_normal.fill(0)
    painter_normal = QPainter(canvas_normal)
    compositor.draw(
        painter_normal,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        viewer_data=_viewer_data((0.0, 180.0), latitude=35.0),
        terrain_profile_altaz=terrain_profile,
        earth_guide_opacity=0.0,
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
        cloud_alpha=0.0,
        viewer_data=_viewer_data((0.0, 180.0), latitude=35.0),
        terrain_profile_altaz=terrain_profile,
        earth_guide_opacity=0.0,
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
    compositor = SkyCompositorCache()
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
        cloud_alpha=0.0,
        viewer_data=_viewer_data((0.0, 180.0), latitude=35.0),
        terrain_profile_altaz=[(0.0, float(az)) for az in range(360)],
        night_light_glow_profile=night_profile,
        fast_mode=True,
    )
    painter.end()


def test_compositor_renders_cloud_grid_without_cloud_image() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 80
    sky[..., 3] = 255

    grid = CloudAltAzGrid(
        amount=np.ones((90, 720), dtype=np.float32),
        missing_mask=np.zeros((90, 720), dtype=np.uint8),
        alt_min_deg=0.0,
        alt_max_deg=90.0,
        az_min_deg=0.0,
        az_max_deg=360.0,
        observer_lat=35.0,
        observer_lon=135.0,
        satellite="Geo-sat",
        product="infrared",
        time_utc=datetime(2026, 6, 22, tzinfo=timezone.utc),
        shells_km=(),
        source_key=SourceKey(
            satellite="Geo-sat",
            provider="infrared",
            timeslot_utc=datetime(2026, 6, 22, tzinfo=timezone.utc),
        ),
        coverage_ratio=1.0,
        grid_resolution_deg=0.5,
    )

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.4,
        viewer_data=_viewer_data((45.0, 180.0)),
        cloud_altaz_grid=grid,
    )
    painter.end()

    out = qimage_to_np_rgba(canvas)
    assert np.any(out != sky)


def test_compositor_does_not_apply_full_window_ground_reset() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()
    terrain_profile = [(45.0, float(az)) for az in range(360)]

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        viewer_data=_viewer_data(),
        terrain_profile_altaz=terrain_profile,
        ground_reset_rgba=(12, 34, 56, 255),
        earth_guide_opacity=0.0,
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    assert np.array_equal(arr[32, 32, :3], np.array([0, 0, 0], dtype=np.uint8))
    assert int(arr[32, 32, 3]) == 0


def test_compositor_has_no_ground_tint_option() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()
    terrain_profile = [(30.0, float(az)) for az in range(360)]

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        viewer_data=_viewer_data(),
        terrain_profile_altaz=terrain_profile,
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    assert np.array_equal(arr[24, 32, :3], np.array([0, 0, 0], dtype=np.uint8))


def test_compositor_observer_latitude_draws_never_rises_outline() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        viewer_data=_viewer_data((0.0, 180.0), latitude=35.0),
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    rgb = arr[..., :3]
    changed_mask = np.any(
        rgb != np.array([100, 100, 100], dtype=np.uint8)[None, None, :],
        axis=2,
    )
    outline_mask = changed_mask & (arr[..., 3] > 0)
    assert np.any(outline_mask)
    changed_rgb = rgb[outline_mask]
    assert np.all(changed_rgb[:, 0] >= changed_rgb[:, 1])
    assert np.all(changed_rgb[:, 1] >= changed_rgb[:, 2])
    assert np.all(arr[..., 3][outline_mask] > 0)


def test_compositor_guidelines_toggle_hides_never_rises_outline() -> None:
    sky = np.zeros((64, 64, 4), dtype=np.uint8)
    sky[..., :3] = 100
    sky[..., 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()

    canvas = QImage(64, 64, QImage.Format_ARGB32_Premultiplied)
    canvas.fill(0)
    painter = QPainter(canvas)
    compositor.draw(
        painter,
        geom,
        np_rgba_to_qimage(sky),
        cloud_alpha=0.0,
        viewer_data=_viewer_data((0.0, 180.0), latitude=35.0),
        show_guidelines=False,
        earth_guide_opacity=0.0,
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    yy, xx = np.ogrid[:64, :64]
    distance = np.hypot(xx - 32.0, yy - 32.0)
    inner_disc = distance < 31.0
    outside_disc = distance > 35.0
    assert np.all(arr[..., 3][inner_disc] == 255)
    assert np.all(arr[..., :3][inner_disc] == 100)
    assert np.all(arr[..., 3][outside_disc] == 0)


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

        def setRenderHint(self, *_args, **_kwargs) -> None:
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
