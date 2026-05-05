from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage, QPainter

from zstarview.gui.composite import SkyCompositorCache, mask_cloud_alpha_by_missing, overlay_missing_tint
from zstarview.types import ScreenGeometry
from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


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
    assert np.array_equal(arr_normal[40, 32, :3], np.array([100, 100, 100], dtype=np.uint8))
    assert np.array_equal(arr_fast[40, 32, :3], np.array([100, 100, 100], dtype=np.uint8))


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
    outline_rgb = np.array([240, 173, 122], dtype=np.int16)
    rgb = arr[..., :3].astype(np.int16)
    outline_mask = np.all(np.abs(rgb - outline_rgb[None, None, :]) <= 3, axis=2)
    assert np.any(outline_mask)
    assert np.all(arr[..., 3][outline_mask] == 51)
    assert np.array_equal(arr[40, 32, :3], np.array([100, 100, 100], dtype=np.uint8))
