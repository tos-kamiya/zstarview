import numpy as np
from PySide6.QtCore import QRect

from zstarview.paths import HatchConfig
from zstarview.gui.composite import (
    build_cloud_amount_field,
    compose_cloud_over_sky,
    render_variable_width_cloud_stripes,
)
from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


def test_cloud_amount_field_tracks_alpha_strength() -> None:
    arr = np.zeros((64, 64, 4), dtype=np.uint8)
    arr[:, :32, 3] = 48
    arr[:, 32:, 3] = 220

    field = build_cloud_amount_field(np_rgba_to_qimage(arr), bins=64)

    left_mean = float(field.amount[:, :32].mean())
    right_mean = float(field.amount[:, 32:].mean())
    assert right_mean > left_mean + 0.12
    assert right_mean > left_mean * 1.8


def test_variable_width_cloud_stripes_keep_fixed_alpha() -> None:
    arr = np.zeros((64, 64, 4), dtype=np.uint8)
    arr[..., 3] = 180
    field = build_cloud_amount_field(np_rgba_to_qimage(arr), bins=64)
    cfg = HatchConfig(20, 19, 8, 173)

    out = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 128, 128, cfg))
    positive = out[..., 3] > 0
    assert np.any(positive)
    assert np.array_equal(np.unique(out[..., 3][positive]), np.array([173], dtype=np.uint8))


def test_compose_cloud_addition_is_weighted_by_cloud_alpha() -> None:
    sky = np.zeros((8, 8, 4), dtype=np.uint8)
    sky[..., :3] = 20
    sky[..., 3] = 255

    cloud = np.zeros((8, 8, 4), dtype=np.uint8)
    cloud[..., :3] = 255
    cloud[..., 3] = 0
    cloud[4, 4, 3] = 255

    out = qimage_to_np_rgba(
        compose_cloud_over_sky(
            sky_img=np_rgba_to_qimage(sky),
            cloud_img_rgba=np_rgba_to_qimage(cloud),
            dest_rect=QRect(0, 0, 8, 8),
            cloud_opacity=1.0,
            gray_mix=0.0,
        )
    )
    # Alpha=0 cloud pixels should not brighten sky.
    assert int(out[1, 1, 0]) == 20
    # Alpha=255 cloud pixel should be brighter due to additive cloud term.
    assert int(out[4, 4, 0]) > 20


def test_compose_cloud_without_sky_uses_opaque_black_disc_base() -> None:
    sky = np.zeros((8, 8, 4), dtype=np.uint8)

    cloud = np.zeros((8, 8, 4), dtype=np.uint8)
    cloud[..., :3] = 255
    cloud[..., 3] = 0
    cloud[4, 4, 3] = 255

    out = qimage_to_np_rgba(
        compose_cloud_over_sky(
            sky_img=np_rgba_to_qimage(sky),
            cloud_img_rgba=np_rgba_to_qimage(cloud),
            dest_rect=QRect(0, 0, 8, 8),
            cloud_opacity=1.0,
            gray_mix=0.0,
        )
    )

    assert int(out[1, 1, 3]) == 255
    assert int(out[4, 4, 3]) == 255
    assert int(out[1, 1, 0]) == 0
    assert int(out[4, 4, 0]) > 0


def test_variable_width_cloud_stripes_make_dense_regions_wider() -> None:
    base = np.zeros((128, 128, 4), dtype=np.uint8)
    base[:, :64, 3] = 48
    base[:, 64:, 3] = 220

    field = build_cloud_amount_field(np_rgba_to_qimage(base), bins=96)
    cfg = HatchConfig(20, 19, 8, 255)

    out = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 256, 256, cfg))
    left_ratio = float(np.mean(out[:, :128, 3] > 0))
    right_ratio = float(np.mean(out[:, 128:, 3] > 0))
    assert right_ratio > left_ratio + 0.08


def test_variable_width_cloud_stripes_stay_sparse_on_larger_canvas() -> None:
    base = np.zeros((128, 128, 4), dtype=np.uint8)
    base[..., :3] = 255
    base[..., 3] = 180

    field = build_cloud_amount_field(np_rgba_to_qimage(base), bins=96)
    cfg = HatchConfig(20, 19, 8, 255)

    small = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 200, 200, cfg))
    large = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 800, 800, cfg))

    # Larger canvas should not increase stripe occupancy; it should stay sparse or become sparser.
    small_ratio = float(np.mean(small[..., 3] > 0))
    large_ratio = float(np.mean(large[..., 3] > 0))
    assert large_ratio <= small_ratio + 0.05
