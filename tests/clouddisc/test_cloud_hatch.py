import datetime as dt

import numpy as np
from PySide6.QtCore import QRect

from zstarview.clouddisc.altaz_grid import CloudAltAzGrid
from zstarview.clouddisc.types import SourceKey
from zstarview.gui.composite import (
    CloudAmountField,
    _cloud_render_content_fov_deg,
    _cloud_stripe_fade_factor,
    _render_alpha_scaled_cloud_stripes_rgba_from_altaz_grid,
    _scale_qimage_preserving_aspect,
    _render_variable_width_cloud_stripes_rgba_from_altaz_grid,
    _scaled_cloud_target_stripes,
    _stripe_render_grids,
    build_cloud_amount_field,
    compose_cloud_over_sky,
    render_variable_width_cloud_stripes,
)
from zstarview.paths import HatchConfig
from zstarview.render.pipeline import compute_star_render_surface_size
from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


def _make_altaz_grid(amount_value: float = 1.0) -> CloudAltAzGrid:
    amount = np.full((90, 720), amount_value, dtype=np.float32)
    missing = np.zeros((90, 720), dtype=np.uint8)
    source_key = SourceKey(
        satellite="G19",
        provider="GOES",
        timeslot_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        sat_priority=("AUTO",),
    )
    return CloudAltAzGrid(
        amount=amount,
        missing_mask=missing,
        alt_min_deg=0.0,
        alt_max_deg=90.0,
        az_min_deg=0.0,
        az_max_deg=360.0,
        observer_lat=35.0,
        observer_lon=135.0,
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        shells_km=(6374.0, 6376.0, 6378.0),
        source_key=source_key,
        coverage_ratio=1.0,
    )


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


def test_compose_cloud_low_opacity_preserves_sky_color() -> None:
    sky = np.zeros((8, 8, 4), dtype=np.uint8)
    sky[..., 0] = 100
    sky[..., 1] = 50
    sky[..., 2] = 0
    sky[..., 3] = 255

    cloud = np.zeros((8, 8, 4), dtype=np.uint8)
    cloud[..., 3] = 0
    cloud[4, 4, 3] = 255

    out = qimage_to_np_rgba(
        compose_cloud_over_sky(
            sky_img=np_rgba_to_qimage(sky),
            cloud_img_rgba=np_rgba_to_qimage(cloud),
            dest_rect=QRect(0, 0, 8, 8),
            cloud_opacity=0.1,
            gray_mix=1.0,
            content_fov_deg=90.0,
        )
    )

    assert int(out[4, 4, 0]) > 80
    assert int(out[4, 4, 1]) > 40
    assert int(out[4, 4, 2]) < 20


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


def test_scale_qimage_preserving_aspect_centers_without_stretching() -> None:
    src = np.zeros((2, 4, 4), dtype=np.uint8)
    src[..., 3] = 255

    out = qimage_to_np_rgba(_scale_qimage_preserving_aspect(np_rgba_to_qimage(src), 8, 8))

    assert int(out[:, :, 3].max()) == 255
    assert np.all(out[:2, :, 3] == 0)
    assert np.all(out[6:, :, 3] == 0)
    assert np.all(out[2:6, :, 3] == 255)


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


def test_scaled_cloud_target_stripes_tracks_buffer_size() -> None:
    assert _scaled_cloud_target_stripes(50, 600, 600) == 50
    assert _scaled_cloud_target_stripes(50, 848, 848) == 71
    assert _scaled_cloud_target_stripes(50, 300, 300) == 25


def test_variable_width_cloud_stripes_use_star_surface_density() -> None:
    base = np.zeros((128, 128, 4), dtype=np.uint8)
    base[..., :3] = 255
    base[..., 3] = 180

    field = build_cloud_amount_field(np_rgba_to_qimage(base), bins=96)
    cfg = HatchConfig(20, 19, 8, 255)

    large_w, large_h = 800, 800
    low_w, low_h = compute_star_render_surface_size(large_w, large_h, 800, 600)
    assert (low_w, low_h) == (693, 693)

    small = qimage_to_np_rgba(
        render_variable_width_cloud_stripes(
            field,
            200,
            200,
            cfg,
            content_fov_deg=90.0,
            density_reference_size=(200, 200),
        )
    )
    large = qimage_to_np_rgba(
        render_variable_width_cloud_stripes(
            field,
            large_w,
            large_h,
            cfg,
            content_fov_deg=90.0,
            density_reference_size=(low_w, low_h),
        )
    )

    small_runs = _alpha_runs(small[100, :, 3])
    large_runs = _alpha_runs(large[400, :, 3])
    assert len(large_runs) > len(small_runs)


def test_variable_width_cloud_stripes_fade_away_from_center_line() -> None:
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
    center = out.shape[1] // 2
    assert start < center < end
    assert abs((center - start) - (end - center)) <= 1
    center_alpha = int(out[row_index, center, 3])
    left_alpha = int(out[row_index, start, 3])
    right_alpha = int(out[row_index, end, 3])
    assert center_alpha > left_alpha + 40
    assert center_alpha > right_alpha + 40


def test_cloud_stripe_fade_factor_is_ease_out() -> None:
    fade_span = 10.0
    fade_start = float(_cloud_stripe_fade_factor(np.array([0.5]), fade_span)[0])
    fade_mid = float(_cloud_stripe_fade_factor(np.array([0.5 + fade_span * 0.5]), fade_span)[0])
    fade_end = float(_cloud_stripe_fade_factor(np.array([0.5 + fade_span]), fade_span)[0])

    assert fade_start == 1.0
    assert fade_mid > 0.75
    assert fade_end == 0.5


def test_cloud_render_content_fov_adds_sky_margin() -> None:
    assert _cloud_render_content_fov_deg(110.0) == 122.0
    assert _cloud_render_content_fov_deg(170.0) == 180.0


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


def test_variable_width_cloud_stripes_anchor_center_line() -> None:
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

    center = low_out.shape[1] // 2
    low_match = next((run for run in low_runs if run[0] <= center <= run[1]), None)
    high_match = next((run for run in high_runs if run[0] <= center <= run[1]), None)
    assert low_match is not None
    assert high_match is not None
    low_start, low_end = low_match
    high_start, high_end = high_match
    assert low_start < center < low_end
    assert high_start < center < high_end
    assert (high_end - high_start) > (low_end - low_start)


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


def test_alpha_scaled_cloud_stripes_encode_cloud_amount_in_alpha_from_altaz_grid() -> None:
    cfg = HatchConfig(20, 19, 8, 255)
    amount = np.full((90, 720), 0.2, dtype=np.float32)
    amount[:, 360:] = 0.8
    grid = CloudAltAzGrid(
        amount=amount,
        missing_mask=np.zeros((90, 720), dtype=np.uint8),
        alt_min_deg=0.0,
        alt_max_deg=90.0,
        az_min_deg=0.0,
        az_max_deg=360.0,
        observer_lat=35.0,
        observer_lon=135.0,
        satellite="G19",
        product="CMIPF-C13",
        time_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
        shells_km=(6374.0, 6376.0, 6378.0),
        source_key=SourceKey(
            satellite="G19",
            provider="GOES",
            timeslot_utc=dt.datetime(2026, 6, 22, 12, 0, 0, tzinfo=dt.timezone.utc),
            sat_priority=("AUTO",),
        ),
        coverage_ratio=1.0,
    )

    out = _render_alpha_scaled_cloud_stripes_rgba_from_altaz_grid(
        grid,
        192,
        192,
        cfg,
        view_center=(45.0, 180.0),
        target_stripes=12,
        width_factor=0.2,
        content_fov_deg=90.0,
    )
    positive = out[..., 3] > 0
    assert np.any(positive)
    left_alpha = out[:, :96, 3]
    right_alpha = out[:, 96:, 3]
    left_mean = float(left_alpha[left_alpha > 0].mean())
    right_mean = float(right_alpha[right_alpha > 0].mean())
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


def test_variable_width_cloud_stripes_can_render_directly_from_altaz_grid() -> None:
    cfg = HatchConfig(20, 19, 8, 255)
    grid = _make_altaz_grid(1.0)

    out = _render_variable_width_cloud_stripes_rgba_from_altaz_grid(
        grid,
        192,
        192,
        cfg,
        view_center=(45.0, 180.0),
        content_fov_deg=90.0,
    )
    assert out.shape == (192, 192, 4)
    assert np.any(out[..., 3] > 0)


def test_alpha_scaled_cloud_stripes_can_render_directly_from_altaz_grid() -> None:
    cfg = HatchConfig(20, 19, 8, 255)
    grid = _make_altaz_grid(1.0)

    out = _render_alpha_scaled_cloud_stripes_rgba_from_altaz_grid(
        grid,
        192,
        192,
        cfg,
        view_center=(45.0, 180.0),
        content_fov_deg=90.0,
    )
    assert out.shape == (192, 192, 4)
    assert np.any(out[..., 3] > 0)


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
