import numpy as np
from PySide6.QtCore import QRect

from zstarview.paths import HatchConfig
from zstarview.gui.composite import (
    CloudAmountField,
    _cloud_stripe_fade_factor,
    _render_alpha_scaled_cloud_stripes_rgba,
    _stripe_render_grids,
    build_cloud_amount_field,
    compose_cloud_over_sky,
    render_variable_width_cloud_stripes,
)
from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


def _alpha_runs(row: np.ndarray) -> list[tuple[int, int]]:
    positive = np.flatnonzero(row > 0)
    if positive.size == 0:
        return []
    runs: list[tuple[int, int]] = []
    start = int(positive[0])
    prev = int(positive[0])
    for idx in positive[1:]:
        current = int(idx)
        if current != prev + 1:
            runs.append((start, prev))
            start = current
        prev = current
    runs.append((start, prev))
    return runs


def test_cloud_amount_field_tracks_alpha_strength() -> None:
    arr = np.zeros((64, 64, 4), dtype=np.uint8)
    arr[:, :32, 3] = 48
    arr[:, 32:, 3] = 220

    field = build_cloud_amount_field(np_rgba_to_qimage(arr), bins=64)

    left_mean = float(field.amount[:, :32].mean())
    right_mean = float(field.amount[:, 32:].mean())
    assert right_mean > left_mean + 0.12
    assert right_mean > left_mean * 1.8


def test_variable_width_cloud_stripes_keep_peak_alpha_at_base_line() -> None:
    arr = np.zeros((64, 64, 4), dtype=np.uint8)
    arr[..., 3] = 180
    field = build_cloud_amount_field(np_rgba_to_qimage(arr), bins=64)
    cfg = HatchConfig(20, 19, 8, 173)

    out = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 128, 128, cfg, content_fov_deg=90.0))
    positive = out[..., 3] > 0
    assert np.any(positive)
    assert int(out[..., 3][positive].max()) == 173
    assert int(out[..., 3][positive].min()) < 173


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
            content_fov_deg=90.0,
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
            content_fov_deg=90.0,
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

    out = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 256, 256, cfg, content_fov_deg=90.0))
    left_ratio = float(np.mean(out[:, :128, 3] > 0))
    right_ratio = float(np.mean(out[:, 128:, 3] > 0))
    assert right_ratio > left_ratio + 0.08


def test_variable_width_cloud_stripes_stay_sparse_on_larger_canvas() -> None:
    base = np.zeros((128, 128, 4), dtype=np.uint8)
    base[..., :3] = 255
    base[..., 3] = 180

    field = build_cloud_amount_field(np_rgba_to_qimage(base), bins=96)
    cfg = HatchConfig(20, 19, 8, 255)

    small = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 200, 200, cfg, content_fov_deg=90.0))
    large = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 800, 800, cfg, content_fov_deg=90.0))

    # Larger canvas should not increase stripe occupancy; it should stay sparse or become sparser.
    small_ratio = float(np.mean(small[..., 3] > 0))
    large_ratio = float(np.mean(large[..., 3] > 0))
    assert large_ratio <= small_ratio + 0.05


def test_variable_width_cloud_stripes_fade_away_from_base_line() -> None:
    cfg = HatchConfig(20, 19, 8, 255)
    field = CloudAmountField(
        amount=np.full((96, 96), 1.0, dtype=np.float32),
        u_min=-2.0,
        u_max=2.0,
        v_min=-2.0,
        v_max=2.0,
        nonzero_lo=0.0,
        nonzero_hi=1.0,
        source_cache_key=7,
    )

    out = qimage_to_np_rgba(
        render_variable_width_cloud_stripes(
            field,
            192,
            192,
            cfg,
            width_factor=0.5,
            content_fov_deg=90.0,
        )
    )
    row_index = 96
    runs = _alpha_runs(out[row_index, :, 3])
    assert runs

    start, end = runs[len(runs) // 2]
    start_alpha = int(out[row_index, start, 3])
    end_alpha = int(out[row_index, end, 3])
    assert start_alpha > end_alpha + 40
    assert end_alpha <= 170


def test_cloud_stripe_fade_factor_is_ease_out() -> None:
    fade_span = 10.0
    fade_start = float(_cloud_stripe_fade_factor(np.array([0.5]), fade_span)[0])
    fade_mid = float(_cloud_stripe_fade_factor(np.array([0.5 + fade_span * 0.5]), fade_span)[0])
    fade_end = float(_cloud_stripe_fade_factor(np.array([0.5 + fade_span]), fade_span)[0])

    assert fade_start == 1.0
    assert fade_mid > 0.75
    assert fade_end == 0.5


def test_stripe_render_grids_use_baseline_projection_for_sampling() -> None:
    phase, line_mask, inside_disc, sample_idx = _stripe_render_grids(
        16,
        16,
        10,
        8.5,
        8.0,
        8.0,
        8.0,
        16,
        16,
        90.0,
    )

    assert bool(line_mask[2, 6])
    assert bool(line_mask[1, 7])
    assert bool(inside_disc[2, 6])
    assert bool(inside_disc[1, 7])
    assert int(sample_idx[2, 6]) == int(sample_idx[1, 7])
    assert float(phase[2, 6]) != float(phase[1, 7])


def test_variable_width_cloud_stripes_anchor_lower_left_edge() -> None:
    cfg = HatchConfig(20, 19, 8, 255)
    low_field = CloudAmountField(
        amount=np.full((96, 96), 0.2, dtype=np.float32),
        u_min=-2.0,
        u_max=2.0,
        v_min=-2.0,
        v_max=2.0,
        nonzero_lo=0.0,
        nonzero_hi=1.0,
        source_cache_key=1,
    )
    high_field = CloudAmountField(
        amount=np.full((96, 96), 1.0, dtype=np.float32),
        u_min=-2.0,
        u_max=2.0,
        v_min=-2.0,
        v_max=2.0,
        nonzero_lo=0.0,
        nonzero_hi=1.0,
        source_cache_key=2,
    )
    low_out = qimage_to_np_rgba(render_variable_width_cloud_stripes(low_field, 192, 192, cfg, content_fov_deg=90.0))
    high_out = qimage_to_np_rgba(render_variable_width_cloud_stripes(high_field, 192, 192, cfg, content_fov_deg=90.0))

    row_index = 96
    low_runs = _alpha_runs(low_out[row_index, :, 3])
    high_runs = _alpha_runs(high_out[row_index, :, 3])
    assert low_runs
    assert high_runs

    low_start, low_end = low_runs[1]
    matching_high = next((run for run in high_runs if run[0] == low_start), None)
    assert matching_high is not None
    high_start, high_end = matching_high
    assert high_start == low_start
    assert high_end > low_end


def test_variable_width_cloud_stripes_drop_to_zero_for_tiny_cloud_amount() -> None:
    cfg = HatchConfig(20, 19, 8, 255)
    tiny_field = CloudAmountField(
        amount=np.full((96, 96), 0.01, dtype=np.float32),
        u_min=-2.0,
        u_max=2.0,
        v_min=-2.0,
        v_max=2.0,
        nonzero_lo=0.0,
        nonzero_hi=1.0,
        source_cache_key=3,
    )

    out = qimage_to_np_rgba(render_variable_width_cloud_stripes(tiny_field, 192, 192, cfg, content_fov_deg=90.0))
    assert not np.any(out[..., 3] > 0)


def test_alpha_scaled_cloud_stripes_encode_cloud_amount_in_alpha() -> None:
    cfg = HatchConfig(20, 19, 8, 255)
    field = CloudAmountField(
        amount=np.full((96, 96), 0.21, dtype=np.float32),
        u_min=-2.0,
        u_max=2.0,
        v_min=-2.0,
        v_max=2.0,
        nonzero_lo=0.2,
        nonzero_hi=0.22,
        source_cache_key=4,
    )
    field.amount[:, 48:] = 0.22

    out = _render_alpha_scaled_cloud_stripes_rgba(
        field,
        192,
        192,
        cfg,
        target_stripes=12,
        width_factor=0.2,
        content_fov_deg=90.0,
    )
    positive = out[..., 3] > 0
    assert np.any(positive)
    left_mean = float(out[:, :96, 3][out[:, :96, 3] > 0].mean())
    right_mean = float(out[:, 96:, 3][out[:, 96:, 3] > 0].mean())
    assert right_mean > left_mean + 40.0


def test_variable_width_cloud_stripes_use_content_fov_for_sampling_extent() -> None:
    cfg = HatchConfig(20, 19, 8, 255)
    y, x = np.ogrid[:96, :96]
    center = 47.5
    field = CloudAmountField(
        amount=np.clip(1.0 - (np.hypot(x - center, y - center) / (96 / 2)), 0.0, 1.0).astype(np.float32),
        u_min=-2.0,
        u_max=2.0,
        v_min=-2.0,
        v_max=2.0,
        nonzero_lo=0.0,
        nonzero_hi=1.0,
        source_cache_key=5,
    )

    narrow = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 192, 192, cfg, content_fov_deg=90.0))
    wide = qimage_to_np_rgba(render_variable_width_cloud_stripes(field, 192, 192, cfg, content_fov_deg=120.0))

    assert int(np.count_nonzero(wide[8, :, 3])) > int(np.count_nonzero(narrow[8, :, 3])) + 20


def test_variable_width_cloud_stripes_use_fractional_alpha_for_partial_line() -> None:
    cfg = HatchConfig(20, 19, 8, 200)
    field = CloudAmountField(
        amount=np.full((96, 96), 0.42, dtype=np.float32),
        u_min=-2.0,
        u_max=2.0,
        v_min=-2.0,
        v_max=2.0,
        nonzero_lo=0.0,
        nonzero_hi=1.0,
        source_cache_key=6,
    )

    out = qimage_to_np_rgba(
        render_variable_width_cloud_stripes(
            field,
            192,
            192,
            cfg,
            target_stripes=12,
            width_factor=0.5,
            content_fov_deg=90.0,
        )
    )
    positive = out[..., 3][out[..., 3] > 0]
    assert np.any(positive)
    assert int(positive.max()) == 200
    assert np.any(positive < 200)
