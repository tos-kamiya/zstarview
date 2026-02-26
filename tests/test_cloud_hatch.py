import numpy as np
from PySide6.QtCore import QRect

from zstarview.paths import HatchConfig
from zstarview.ui.composite import cloud_with_hatched_alpha, compose_cloud_over_sky
from zstarview.ui.composite import make_hatch_tile_qimage
from zstarview.utils.qt import np_rgba_to_qimage, qimage_to_np_rgba


def test_make_hatch_tile_contains_stripes() -> None:
    tile = qimage_to_np_rgba(make_hatch_tile_qimage(20, 19, 8, 255))
    alpha = tile[..., 3]
    assert int(alpha.min()) == 0
    assert int(alpha.max()) == 255
    assert 0.10 < float((alpha == 255).mean()) < 0.90


def test_make_hatch_tile_respects_strength() -> None:
    tile = qimage_to_np_rgba(make_hatch_tile_qimage(20, 19, 8, 90))
    alpha = tile[..., 3]
    assert int(alpha.max()) == 90


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
