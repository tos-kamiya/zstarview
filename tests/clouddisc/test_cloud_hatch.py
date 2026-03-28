import numpy as np
from PySide6.QtCore import QRect

from zstarview.paths import HatchConfig
from zstarview.gui.composite import (
    build_stripe_density_field,
    cloud_with_hatched_alpha,
    compose_cloud_over_sky,
    render_hatched_cloud_from_density,
)
from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


def test_cloud_hatch_reduces_rgb_and_alpha_together() -> None:
    arr = np.zeros((64, 64, 4), dtype=np.uint8)
    arr[..., :3] = 200
    arr[..., 3] = 255
    out = qimage_to_np_rgba(
        cloud_with_hatched_alpha(np_rgba_to_qimage(arr), HatchConfig(20, 19, 8, 255))
    )
    # Some pixels should be fully cut by hatch.
    assert int(out[..., 3].min()) == 0
    # Where alpha is fully cut, RGB should also be cut.
    cut = out[..., 3] == 0
    assert np.any(cut)
    assert int(out[..., :3][cut].max()) == 0


def test_cloud_hatch_flattens_alpha_within_each_cloud_band() -> None:
    h, w = 64, 64
    arr = np.zeros((h, w, 4), dtype=np.uint8)
    arr[..., :3] = 255
    # Deliberately include variation so flattening effect is observable.
    arr[..., 3] = np.tile(np.linspace(0, 255, w, dtype=np.uint8), (h, 1))

    cfg = HatchConfig(20, 19, 8, 255)
    out = qimage_to_np_rgba(cloud_with_hatched_alpha(np_rgba_to_qimage(arr), cfg))

    xs = np.arange(w, dtype=np.int32)[None, :]
    ys = np.arange(h, dtype=np.int32)[:, None]
    period = max(2, int(round(np.hypot(cfg.tile_w_px, cfg.tile_h_px) * 0.5)))
    band = max(1, int(round(cfg.line_px)))
    u = xs - ys
    u_mod = np.mod(u, period)
    dist = np.minimum(u_mod, period - u_mod)
    keep_mask = dist > (band / 2.0)
    stripe_id = np.floor_divide(u, period)

    # In cloud bands (keep_mask), alpha should be constant per stripe.
    for sid in np.unique(stripe_id[keep_mask]):
        m = keep_mask & (stripe_id == sid)
        if int(np.count_nonzero(m)) <= 1:
            continue
        uniq_alpha = np.unique(out[..., 3][m])
        assert len(uniq_alpha) == 1


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


def test_density_hatch_keeps_stripes_sparse_on_larger_canvas() -> None:
    base = np.zeros((128, 128, 4), dtype=np.uint8)
    base[..., :3] = 255
    base[..., 3] = 180

    density = build_stripe_density_field(np_rgba_to_qimage(base), bins=96)
    cfg = HatchConfig(20, 19, 8, 255)

    small = qimage_to_np_rgba(render_hatched_cloud_from_density(density, 200, 200, cfg))
    large = qimage_to_np_rgba(render_hatched_cloud_from_density(density, 800, 800, cfg))

    # Larger canvas should not increase stripe occupancy; it should stay sparse or become sparser.
    small_ratio = float(np.mean(small[..., 3] > 0))
    large_ratio = float(np.mean(large[..., 3] > 0))
    assert large_ratio <= small_ratio + 0.05
