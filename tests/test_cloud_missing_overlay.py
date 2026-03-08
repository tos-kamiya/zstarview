from __future__ import annotations

import numpy as np
from PySide6.QtGui import QImage, QPainter

from zstarview.ui.composite import SkyCompositorCache, mask_cloud_alpha_by_missing, overlay_missing_tint
from zstarview.types import ScreenGeometry
from zstarview.utils.qt import np_rgba_to_qimage, qimage_to_np_rgba


def test_overlay_missing_with_hatch_tints_only_missing_region() -> None:
    base = np.zeros((16, 16, 4), dtype=np.uint8)
    base[..., :3] = 30
    base[..., 3] = 255

    missing = np.zeros((16, 16, 4), dtype=np.uint8)
    missing[4:12, 4:12, 3] = 255

    out = qimage_to_np_rgba(
        overlay_missing_tint(
            np_rgba_to_qimage(base),
            np_rgba_to_qimage(missing),
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

    missing_none = np.zeros((64, 64, 4), dtype=np.uint8)
    missing_half = np.zeros((64, 64, 4), dtype=np.uint8)
    missing_half[:, :32, 3] = 255

    geom = ScreenGeometry(center=(32, 32), radius=32)
    compositor = SkyCompositorCache()

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
        stripe_density=None,
        missing_mask=np_rgba_to_qimage(missing_none),
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
        stripe_density=None,
        missing_mask=np_rgba_to_qimage(missing_half),
    )
    p2.end()

    arr1 = qimage_to_np_rgba(canvas1)
    arr2 = qimage_to_np_rgba(canvas2)
    assert np.any(arr1 != arr2)


def test_mask_cloud_alpha_by_missing_cuts_cloud_pixels() -> None:
    cloud = np.zeros((12, 12, 4), dtype=np.uint8)
    cloud[..., :3] = 255
    cloud[..., 3] = 180
    missing = np.zeros((12, 12, 4), dtype=np.uint8)
    missing[3:9, 4:10, 3] = 255

    out = qimage_to_np_rgba(mask_cloud_alpha_by_missing(np_rgba_to_qimage(cloud), np_rgba_to_qimage(missing)))
    assert int(out[1, 1, 3]) == 180
    assert int(out[5, 6, 3]) == 0


def test_compositor_terrain_profile_tints_ground_below_terrain_horizon() -> None:
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
        None,
        cloud_alpha=0.0,
        view_center=(0.0, 0.0),
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
    )
    p_terrain.end()

    arr_flat = qimage_to_np_rgba(canvas_flat)
    arr_terrain = qimage_to_np_rgba(canvas_terrain)

    assert np.array_equal(arr_flat[24, 32, :3], np.array([100, 100, 100], dtype=np.uint8))
    assert np.array_equal(arr_terrain[24, 32, :3], np.array([93, 76, 33], dtype=np.uint8))


def test_compositor_ground_tint_opacity_zero_makes_ground_fill_black() -> None:
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
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    assert np.array_equal(arr[24, 32, :3], np.array([0, 0, 0], dtype=np.uint8))


def test_compositor_reapplies_never_rises_tint_after_ground_fill() -> None:
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
    )
    painter.end()

    arr = qimage_to_np_rgba(canvas)
    assert np.array_equal(arr[40, 32, :3], np.array([112, 79, 36], dtype=np.uint8))
