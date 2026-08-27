"""Cloud field, rasterization, and sky/cloud compositing helpers."""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import cast

import numpy as np
from PySide6.QtCore import QPointF, QRect, Qt
from PySide6.QtGui import (
    QBrush,
    QColor,
    QImage,
    QPainter,
    QPainterPath,
)

from ..clouddisc.altaz_grid import CloudAltAzGrid
from ..paths import (
    HatchConfig,
    ThemeStyle,
)
from ..render.ground_mask import inverse_project_disc as _shared_inverse_project_disc
from ..render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba
from ..render.sky_disc import sky_color_near_solar_horizon
from ..types import ScreenGeometry, ViewProjection

CLOUD_DAY_RGB = (255, 255, 255)
CLOUD_NIGHT_RGB = (184, 196, 224)
CLOUD_TINT_START_SUN_ALT_DEG = -6.0
CLOUD_TINT_FULL_SUN_ALT_DEG = -12.0
CLOUD_NIGHT_BOOST = 0.3
CLOUD_SUNSET_WARM_END_SUN_ALT_DEG = -8.0
CLOUD_SUNSET_SHELL_MIX = (0.22, 0.15, 0.08)
CLOUD_SUNSET_WARM_RB_OFFSET = 10.0
CLOUD_SUNSET_WARM_RB_SCALE = 80.0
CLOUD_SUNSET_WARM_RG_OFFSET = 5.0
CLOUD_SUNSET_WARM_RG_SCALE = 60.0
HALFTONE_MIN_GRID_DELTA_PX = 22.0
HALFTONE_LEVEL_DIAMETER_BASE_SCALE = 0.675


def _inverse_project_disc(
    width: int,
    height: int,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
    content_fov_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return _shared_inverse_project_disc(
        width,
        height,
        geometry,
        view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )


@dataclass(frozen=True)
class CloudAmountField:
    """Compact cloud-amount field in normalized 45-degree (u, v) space."""

    amount: np.ndarray  # float32 in [0,1], shape=(bins_u, bins_v)
    u_min: float
    u_max: float
    v_min: float
    v_max: float
    nonzero_lo: float
    nonzero_hi: float
    source_cache_key: int


def _sample_altaz_grid_amount(
    grid: CloudAltAzGrid,
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
) -> np.ndarray:
    """Sample a `CloudAltAzGrid.amount` field at arbitrary alt/az positions."""
    amount = np.asarray(grid.amount, dtype=np.float32)
    if amount.ndim != 2 or amount.size == 0:
        return np.zeros_like(np.asarray(alt_deg, dtype=np.float32), dtype=np.float32)

    alt = np.asarray(alt_deg, dtype=np.float64)
    az = np.asarray(az_deg, dtype=np.float64)
    if alt.shape != az.shape:
        raise ValueError("alt_deg and az_deg must have the same shape")

    alt_bins, az_bins = amount.shape
    alt_span = max(1.0e-6, float(grid.alt_max_deg) - float(grid.alt_min_deg))
    az_span = max(1.0e-6, float(grid.az_max_deg) - float(grid.az_min_deg))

    valid = np.isfinite(alt) & np.isfinite(az)
    valid &= alt >= float(grid.alt_min_deg)
    valid &= alt <= float(grid.alt_max_deg)
    if not np.any(valid):
        return np.zeros_like(alt, dtype=np.float32)

    alt_pos = (alt - float(grid.alt_min_deg)) / alt_span * float(alt_bins) - 0.5
    az_pos = ((az - float(grid.az_min_deg)) % az_span) / az_span * float(az_bins) - 0.5

    alt0 = np.floor(alt_pos).astype(np.int64, copy=False)
    az0 = np.floor(az_pos).astype(np.int64, copy=False)
    alt1 = np.clip(alt0 + 1, 0, alt_bins - 1)
    az1 = (az0 + 1) % az_bins
    alt0 = np.clip(alt0, 0, alt_bins - 1)
    az0 = az0 % az_bins

    wa = np.clip(alt_pos - np.floor(alt_pos), 0.0, 1.0).astype(np.float32, copy=False)
    wb = np.clip(az_pos - np.floor(az_pos), 0.0, 1.0).astype(np.float32, copy=False)

    sampled = (
        (1.0 - wa) * (1.0 - wb) * amount[alt0, az0]
        + wa * (1.0 - wb) * amount[alt1, az0]
        + (1.0 - wa) * wb * amount[alt0, az1]
        + wa * wb * amount[alt1, az1]
    ).astype(np.float32, copy=False)
    sampled[~valid] = 0.0
    return sampled


def _sample_altaz_grid_to_screen_map(
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    geometry: ScreenGeometry | None,
    projection: ViewProjection,
) -> np.ndarray:
    """Project a `CloudAltAzGrid` into a per-pixel sampled amount map."""
    w = max(1, int(width))
    h = max(1, int(height))
    sampled = np.zeros((h, w), dtype=np.float32)
    alt_deg, az_deg, inside = _inverse_project_disc(
        w,
        h,
        geometry
        if geometry is not None
        else ScreenGeometry(center=((w - 1) // 2, (h - 1) // 2), radius=max(1, min(w, h) // 2)),
        tuple(float(value) for value in projection.view_center),
        edge_fov_deg=float(projection.edge_fov_deg),
        content_fov_deg=float(projection.content_fov_deg),
    )
    if alt_deg.size == 0 or not np.any(inside):
        return sampled

    inside_idx = np.flatnonzero(inside)
    sampled_inside = _sample_altaz_grid_amount(grid, alt_deg, az_deg)
    sampled.reshape(-1)[inside_idx] = sampled_inside
    return sampled


def _scale_qimage_preserving_aspect(
    image: QImage,
    width: int,
    height: int,
) -> QImage:
    """Scale an image without stretching it to a new aspect ratio."""
    target_w = max(1, int(width))
    target_h = max(1, int(height))
    if image.width() == target_w and image.height() == target_h:
        return image
    scaled = image.scaled(
        target_w,
        target_h,
        Qt.KeepAspectRatio,
        Qt.SmoothTransformation,
    )
    if scaled.width() == target_w and scaled.height() == target_h:
        return scaled
    canvas = QImage(target_w, target_h, QImage.Format.Format_ARGB32_Premultiplied)
    canvas.fill(Qt.transparent)
    painter = QPainter(canvas)
    try:
        x = (target_w - scaled.width()) // 2
        y = (target_h - scaled.height()) // 2
        painter.drawImage(x, y, scaled)
    finally:
        painter.end()
    return canvas


CLOUD_DAY_RGB = (255, 255, 255)
CLOUD_NIGHT_RGB = (184, 196, 224)
CLOUD_TINT_START_SUN_ALT_DEG = -6.0
CLOUD_TINT_FULL_SUN_ALT_DEG = -12.0
CLOUD_NIGHT_BOOST = 0.3


def _smoothstep_scalar(edge0: float, edge1: float, x: float) -> float:
    t = float(np.clip((float(x) - float(edge0)) / (float(edge1) - float(edge0)), 0.0, 1.0))
    return t * t * (3.0 - 2.0 * t)


def _cloud_tint_rgb_for_sun_alt(sun_alt_deg: float | None) -> tuple[int, int, int]:
    """Return the display tint for white cloud marks at the current sun altitude."""
    if sun_alt_deg is None:
        return CLOUD_DAY_RGB
    night_mix = _smoothstep_scalar(CLOUD_TINT_START_SUN_ALT_DEG, CLOUD_TINT_FULL_SUN_ALT_DEG, float(sun_alt_deg))
    return cast(
        tuple[int, int, int],
        tuple(
            int(round(float(day) + (float(night) - float(day)) * night_mix))
            for day, night in zip(CLOUD_DAY_RGB, CLOUD_NIGHT_RGB)
        ),
    )


def _cloud_tint_rgb_for_theme(
    theme: ThemeStyle | None,
    sun_alt_deg: float | None,
) -> tuple[int, int, int]:
    if theme is not None and theme.window_background.base_rgb == (255, 255, 255):
        return (132, 146, 164)
    return _cloud_tint_rgb_for_sun_alt(sun_alt_deg)


@lru_cache(maxsize=8)
def _solar_horizon_rgb(
    sun_alt_deg: float,
    sun_az_deg: float,
    aerosol_optical_depth: float | None,
) -> tuple[int, int, int]:
    """Cache the modelled 0..10 degree colour used by shell tinting."""
    return sky_color_near_solar_horizon(
        (float(sun_alt_deg), float(sun_az_deg)),
        alpha=1.0,
        exposure=1.0,
        aerosol_optical_depth=aerosol_optical_depth,
    )[:3]


def _sunset_cloud_tint_rgb(
    sun_altaz: tuple[float, float] | None,
    *,
    base_rgb: tuple[int, int, int],
    shell_index: int,
    aerosol_optical_depth: float | None = None,
) -> tuple[int, int, int]:
    """Return a shell tint using the modelled colour near the solar horizon."""
    if sun_altaz is None or not (0 <= shell_index < len(CLOUD_SUNSET_SHELL_MIX)):
        return tuple(int(np.clip(value, 0, 255)) for value in base_rgb)

    sun_alt_deg = float(sun_altaz[0])
    if sun_alt_deg < CLOUD_SUNSET_WARM_END_SUN_ALT_DEG:
        return tuple(int(np.clip(value, 0, 255)) for value in base_rgb)

    horizon_rgb = np.asarray(
        _solar_horizon_rgb(sun_alt_deg, float(sun_altaz[1]), aerosol_optical_depth),
        dtype=np.float32,
    )
    red_blue_warmth = float(
        np.clip(
            (horizon_rgb[0] - horizon_rgb[2] - CLOUD_SUNSET_WARM_RB_OFFSET)
            / CLOUD_SUNSET_WARM_RB_SCALE,
            0.0,
            1.0,
        )
    )
    red_green_warmth = float(
        np.clip(
            (horizon_rgb[0] - horizon_rgb[1] - CLOUD_SUNSET_WARM_RG_OFFSET)
            / CLOUD_SUNSET_WARM_RG_SCALE,
            0.0,
            1.0,
        )
    )
    warm_factor = min(red_blue_warmth, red_green_warmth)
    if warm_factor <= 0.0:
        return tuple(int(np.clip(value, 0, 255)) for value in base_rgb)

    base = np.asarray(base_rgb, dtype=np.float32)
    mix = warm_factor * float(CLOUD_SUNSET_SHELL_MIX[shell_index])
    tinted = base + (horizon_rgb - base) * mix
    return tuple(int(np.clip(round(float(value)), 0, 255)) for value in tinted)


def _smooth_cloud_amount_grid(values: np.ndarray) -> np.ndarray:
    """Apply a small edge-preserving blur to keep stripe widths from flickering."""
    padded = np.pad(values.astype(np.float32, copy=False), ((1, 1), (1, 1)), mode="edge")
    smoothed = (
        padded[:-2, :-2]
        + 2.0 * padded[:-2, 1:-1]
        + padded[:-2, 2:]
        + 2.0 * padded[1:-1, :-2]
        + 4.0 * padded[1:-1, 1:-1]
        + 2.0 * padded[1:-1, 2:]
        + padded[2:, :-2]
        + 2.0 * padded[2:, 1:-1]
        + padded[2:, 2:]
    ) / 16.0
    return np.clip(smoothed, 0.0, 1.0).astype(np.float32, copy=False)



@lru_cache(maxsize=4)
def _cloud_amount_bin_index_grids(
    height: int,
    width: int,
    bins_u: int,
    bins_v: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Cache normalized cloud bin indices for a given source image shape."""
    h = max(1, int(height))
    w = max(1, int(width))
    cy, cx = (h - 1) * 0.5, (w - 1) * 0.5
    r = max(1.0, min(cx, cy))
    y, x = np.ogrid[:h, :w]
    xn = (x - cx) / r
    yn = (y - cy) / r
    u_idx = np.clip((xn - yn + 2.0) * (bins_u / 4.0), 0.0, bins_u - 1).astype(np.int32)
    v_idx = np.clip((xn + yn + 2.0) * (bins_v / 4.0), 0.0, bins_v - 1).astype(np.int32)
    return (u_idx, v_idx)


@lru_cache(maxsize=4)
def _stripe_render_grids(
    width: int,
    height: int,
    period: int,
    max_band: float,
    cx: float,
    cy: float,
    rr: float,
    bins_u: int,
    bins_v: int,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = 90.0,
    centered: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cache stripe geometry and baseline-projected field sampling grids."""
    w = max(1, int(width))
    h = max(1, int(height))
    xs = np.arange(w, dtype=np.int32)[None, :]
    ys = np.arange(h, dtype=np.int32)[:, None]
    u_pix = xs - ys
    period_i = max(1, int(period))
    if centered:
        center_offset = period_i // 2
        u_mod = np.mod(u_pix + center_offset, period_i)
        # Anchor the center line on a pixel coordinate. This keeps equal
        # pixel distances on both sides of the stroke at the same level.
        phase = np.abs(u_mod.astype(np.float32, copy=False) - float(center_offset))
        u_base = np.floor_divide(u_pix + center_offset, period_i) * period_i + center_offset
    else:
        u_mod = np.mod(u_pix, period_i)
        phase = u_mod.astype(np.float32, copy=False) + 0.5
        u_base = np.floor_divide(u_pix, period_i) * period_i
    line_mask = phase <= float(max_band)

    max_r = max(0.0, float(content_fov_deg) / max(1.0e-6, float(edge_fov_deg)))
    y, x = np.ogrid[:h, :w]
    inside_disc = ((x - cx) ** 2 + (y - cy) ** 2) <= ((rr * max_r) + 0.25) ** 2
    sample_radius = max(1.0, rr * max_r)
    v_pix = xs + ys
    x_base = (v_pix + u_base).astype(np.float32, copy=False) * 0.5
    y_base = (v_pix - u_base).astype(np.float32, copy=False) * 0.5
    xn = (x_base - cx) / sample_radius
    yn = (y_base - cy) / sample_radius
    u_idx = np.clip((xn - yn + 2.0) * (bins_u / 4.0), 0.0, bins_u - 1).astype(np.int32)
    v_idx = np.clip((xn + yn + 2.0) * (bins_v / 4.0), 0.0, bins_v - 1).astype(np.int32)
    return (phase, line_mask, inside_disc, u_idx * bins_v + v_idx)


def _scaled_cloud_target_stripes(
    target_stripes: int,
    reference_width: int,
    reference_height: int,
) -> int:
    """Return the requested stripe count as an absolute value.

    Previously this scaled the stripe count with the reference render
    surface size so that density stayed constant across window sizes.
    The caller now treats `target_stripes` as an absolute stripe count,
    so this helper simply clamps and returns it unchanged.  The
    `reference_width`/`reference_height` parameters are kept for API
    compatibility.
    """
    return max(1, int(target_stripes))


def _cloud_stripe_fade_factor(phase: np.ndarray, fade_span: float) -> np.ndarray:
    """Return a gentle fade curve for variable-width cloud stripes."""
    progress = np.clip((phase - 0.5) / max(1.0, float(fade_span)), 0.0, 1.0)
    return 1.0 - 0.5 * np.square(progress)


def _quantize_cloud_width_levels(normalized: np.ndarray) -> np.ndarray:
    """Quantize normalized width-mode cloud amount to five levels."""
    clipped = np.clip(np.asarray(normalized, dtype=np.float32), 0.0, 1.0)
    return np.floor(clipped * 4.0 + 0.5) / 4.0


def _cloud_render_content_fov_deg(content_fov_deg: float) -> float:
    """Return a slightly expanded FOV for cloud rendering before final clipping."""
    return min(180.0, max(0.0, float(content_fov_deg) + 12.0))


_HALFTONE_GRID_REFERENCE_DIAMETER = 600.0


def _halftone_grid_delta(
    output_diameter: float,
    target_stripes: int,
    *,
    min_grid_delta_px: float = HALFTONE_MIN_GRID_DELTA_PX,
) -> float:
    """Return the halftone grid spacing in pixels.

    The grid spacing follows the same base-plus-upscale rule as the star
    layer: the spacing scales proportionally with the output diameter and is
    clamped to a minimum so compact windows do not collapse the halftone into
    clutter.
    """
    diameter = max(1.0, float(output_diameter))
    return max(float(min_grid_delta_px), diameter / max(1, int(target_stripes)))


def _halftone_level_diameters(delta: float, width_factor: float) -> tuple[float, ...]:
    """Return halftone dot diameters for the 8 quantization levels."""
    wf = max(0.01, float(width_factor))
    diam_scale = max(1.0, float(delta)) / 30.0 * wf * 0.5 * HALFTONE_LEVEL_DIAMETER_BASE_SCALE
    return (
        0.0,
        4.0 * diam_scale,
        8.0 * diam_scale,
        12.0 * diam_scale,
        16.0 * diam_scale,
        20.0 * diam_scale,
        24.0 * diam_scale,
        28.0 * diam_scale,
    )


def build_cloud_amount_field_from_rgba(
    cloud: np.ndarray,
    *,
    bins: int = 320,
    source_cache_key: int = 0,
) -> CloudAmountField:
    """Build a compact cloud-amount field from an RGBA image in normalized (u, v) space."""
    h, w = cloud.shape[:2]

    alpha01 = cloud[..., 3].astype(np.float32) / 255.0
    inside = alpha01 > 0.0
    if not np.any(inside):
        amount = np.zeros((bins, bins), dtype=np.float32)
        return CloudAmountField(
            amount=amount,
            u_min=-2.0,
            u_max=2.0,
            v_min=-2.0,
            v_max=2.0,
            nonzero_lo=0.0,
            nonzero_hi=1.0,
            source_cache_key=int(source_cache_key),
        )

    u_min, u_max = -2.0, 2.0
    v_min, v_max = -2.0, 2.0
    bins_u = max(32, int(bins))
    bins_v = bins_u

    u_idx, v_idx = _cloud_amount_bin_index_grids(h, w, bins_u, bins_v)

    ids = (u_idx[inside] * bins_v + v_idx[inside]).astype(np.int64, copy=False)
    vals = alpha01[inside].astype(np.float64, copy=False)
    size = bins_u * bins_v
    sums = np.bincount(ids, weights=vals, minlength=size)
    counts = np.bincount(ids, minlength=size)
    means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0).astype(np.float32, copy=False)
    amount = cast(np.ndarray, np.asarray(means, dtype=np.float32).reshape((bins_u, bins_v)))
    amount = _smooth_cloud_amount_grid(amount)
    positive = amount[amount > 0.0]
    if positive.size > 0:
        nonzero_lo = float(np.percentile(positive, 12.0))
        nonzero_hi = float(np.percentile(positive, 92.0))
        if nonzero_hi <= nonzero_lo + 1e-6:
            nonzero_lo = float(positive.min())
            nonzero_hi = float(positive.max())
    else:
        nonzero_lo = 0.0
        nonzero_hi = 1.0

    return CloudAmountField(
        amount=amount,
        u_min=u_min,
        u_max=u_max,
        v_min=v_min,
        v_max=v_max,
        nonzero_lo=nonzero_lo,
        nonzero_hi=nonzero_hi,
        source_cache_key=int(source_cache_key),
    )


def build_cloud_amount_field(source_img: QImage, *, bins: int = 320) -> CloudAmountField:
    """Build a compact cloud-amount field from a cloud image in normalized (u, v) space."""
    cloud = qimage_to_np_rgba(
        source_img
        if source_img.format() == QImage.Format_RGBA8888
        else source_img.convertToFormat(QImage.Format_RGBA8888)
    )
    return build_cloud_amount_field_from_rgba(
        cloud,
        bins=bins,
        source_cache_key=int(source_img.cacheKey()),
    )


def _render_variable_width_cloud_stripes_rgba(
    cloud_amount: CloudAmountField,
    width: int,
    height: int,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry | None = None,
    *,
    target_stripes: int = 50,
    width_factor: float = 0.85,
    density_reference_size: tuple[int, int] | None = None,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = 90.0,
) -> np.ndarray:
    """Render fixed-opacity cloud stripes whose width increases with cloud amount."""
    w = max(1, int(width))
    h = max(1, int(height))
    ref_w, ref_h = (
        (w, h)
        if density_reference_size is None
        else (max(1, int(density_reference_size[0])), max(1, int(density_reference_size[1])))
    )

    diameter_px = float(min(w, h))
    stripes = _scaled_cloud_target_stripes(target_stripes, ref_w, ref_h)
    wf = float(np.clip(width_factor, 0.1, 0.95))
    base_period = int(np.clip(round(diameter_px / stripes), 14, 64))
    period = base_period
    max_band = max(1.0, float(base_period) * wf * 0.5)

    if geometry is None:
        cx = (w - 1) * 0.5
        cy = (h - 1) * 0.5
        rr = max(1.0, min(cx, cy))
    else:
        cx = float(geometry.center[0])
        cy = float(geometry.center[1])
        rr = max(1.0, float(geometry.radius))

    out = np.zeros((h, w, 4), dtype=np.uint8)

    bins_u, bins_v = cloud_amount.amount.shape
    phase, line_mask, inside_disc, sample_idx = _stripe_render_grids(
        w,
        h,
        period,
        max_band,
        cx,
        cy,
        rr,
        bins_u,
        bins_v,
        edge_fov_deg,
        _cloud_render_content_fov_deg(content_fov_deg),
        centered=True,
    )
    sampled = np.clip(cloud_amount.amount.reshape(-1)[sample_idx], 0.0, 1.0)
    if cloud_amount.nonzero_hi > cloud_amount.nonzero_lo + 1e-6:
        normalized = (sampled - cloud_amount.nonzero_lo) / (cloud_amount.nonzero_hi - cloud_amount.nonzero_lo)
    else:
        normalized = sampled
    normalized = _quantize_cloud_width_levels(normalized)
    present = sampled > 0.03
    line_index = np.floor(phase).astype(np.int32, copy=False)
    local_levels = np.where(present, normalized * max_band, 0.0)
    whole_levels = np.floor(local_levels).astype(np.int32, copy=False)
    frac_levels = np.clip(local_levels - whole_levels, 0.0, 1.0)
    full_mask = inside_disc & line_mask & (line_index >= 0) & (line_index < whole_levels)
    partial_mask = inside_disc & line_mask & (line_index >= 0) & (line_index == whole_levels) & (frac_levels > 1e-6)
    if not np.any(full_mask) and not np.any(partial_mask):
        return out

    alpha_u8 = int(np.clip(round(float(hatch_cfg.strength)), 1, 255))
    fade_span = max(1.0, float(max_band) - 0.5)
    fade = _cloud_stripe_fade_factor(phase, fade_span)
    if np.any(full_mask):
        out[..., :3][full_mask] = 255
        out[..., 3][full_mask] = np.clip(np.round(alpha_u8 * fade[full_mask]), 0, alpha_u8).astype(np.uint8)
    if np.any(partial_mask):
        out[..., :3][partial_mask] = 255
        partial_alpha = np.clip(np.round(frac_levels * alpha_u8), 0, alpha_u8).astype(np.uint8)
        out[..., 3][partial_mask] = np.clip(
            np.round(partial_alpha[partial_mask].astype(np.float32) * fade[partial_mask]),
            0,
            alpha_u8,
        ).astype(np.uint8)
    return out


def _render_variable_width_cloud_stripes_rgba_from_amount_map(
    sampled_amount: np.ndarray,
    width: int,
    height: int,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry | None = None,
    *,
    target_stripes: int = 50,
    width_factor: float = 0.85,
    density_reference_size: tuple[int, int] | None = None,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = 90.0,
    stripe_rgb: tuple[int, int, int] = (255, 255, 255),
) -> np.ndarray:
    """Render variable-width cloud stripes from a per-pixel sampled amount map."""
    w = max(1, int(width))
    h = max(1, int(height))
    sampled = np.clip(np.asarray(sampled_amount, dtype=np.float32), 0.0, 1.0)
    if sampled.shape != (h, w):
        raise ValueError("sampled_amount shape must match the requested output size")

    ref_w, ref_h = (
        (w, h)
        if density_reference_size is None
        else (max(1, int(density_reference_size[0])), max(1, int(density_reference_size[1])))
    )

    diameter_px = float(min(w, h))
    stripes = _scaled_cloud_target_stripes(target_stripes, ref_w, ref_h)
    wf = float(np.clip(width_factor, 0.1, 0.95))
    base_period = int(np.clip(round(diameter_px / stripes), 14, 64))
    max_band = max(1.0, float(base_period) * wf * 0.5)

    if geometry is None:
        cx = (w - 1) * 0.5
        cy = (h - 1) * 0.5
        rr = max(1.0, min(cx, cy))
    else:
        cx = float(geometry.center[0])
        cy = float(geometry.center[1])
        rr = max(1.0, float(geometry.radius))

    out = np.zeros((h, w, 4), dtype=np.uint8)
    phase, line_mask, inside_disc, _ = _stripe_render_grids(
        w,
        h,
        base_period,
        max_band,
        cx,
        cy,
        rr,
        1,
        1,
        edge_fov_deg,
        _cloud_render_content_fov_deg(content_fov_deg),
        centered=True,
    )
    present = sampled > 0.03
    line_index = np.floor(phase).astype(np.int32, copy=False)
    if np.any(present):
        nonzero = sampled[present]
        nonzero_lo = float(np.percentile(nonzero, 12.0)) if nonzero.size > 0 else 0.0
        nonzero_hi = float(np.percentile(nonzero, 92.0)) if nonzero.size > 0 else 1.0
        if nonzero_hi <= nonzero_lo + 1e-6 and nonzero.size > 0:
            nonzero_lo = float(nonzero.min())
            nonzero_hi = float(nonzero.max())
        if nonzero_hi > nonzero_lo + 1e-6:
            normalized = (sampled - nonzero_lo) / (nonzero_hi - nonzero_lo)
        else:
            normalized = sampled
    else:
        normalized = sampled
    # Keep width-mode cloud density readable as five stable levels rather
    # than allowing every source-value fluctuation to change the band width.
    normalized = _quantize_cloud_width_levels(normalized)
    local_levels = np.where(present, normalized * max_band, 0.0)
    whole_levels = np.floor(local_levels).astype(np.int32, copy=False)
    frac_levels = np.clip(local_levels - whole_levels, 0.0, 1.0)
    full_mask = inside_disc & line_mask & (line_index >= 0) & (line_index < whole_levels)
    partial_mask = inside_disc & line_mask & (line_index >= 0) & (line_index == whole_levels) & (frac_levels > 1e-6)
    if not np.any(full_mask) and not np.any(partial_mask):
        return out

    alpha_u8 = int(np.clip(round(float(hatch_cfg.strength)), 1, 255))
    fade_span = max(1.0, float(max_band) - 0.5)
    fade = _cloud_stripe_fade_factor(phase, fade_span)
    if np.any(full_mask):
        out[..., :3][full_mask] = np.asarray(stripe_rgb, dtype=np.uint8)
        out[..., 3][full_mask] = np.clip(np.round(alpha_u8 * fade[full_mask]), 0, alpha_u8).astype(np.uint8)
    if np.any(partial_mask):
        out[..., :3][partial_mask] = np.asarray(stripe_rgb, dtype=np.uint8)
        partial_alpha = np.clip(np.round(frac_levels * alpha_u8), 0, alpha_u8).astype(np.uint8)
        out[..., 3][partial_mask] = np.clip(
            np.round(partial_alpha[partial_mask].astype(np.float32) * fade[partial_mask]),
            0,
            alpha_u8,
        ).astype(np.uint8)
    return out


def _render_alpha_scaled_cloud_stripes_rgba_from_amount_map(
    sampled_amount: np.ndarray,
    width: int,
    height: int,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry | None = None,
    *,
    target_stripes: int = 50,
    width_factor: float = 0.2,
    density_reference_size: tuple[int, int] | None = None,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = 90.0,
) -> np.ndarray:
    """Render alpha-scaled cloud stripes from a per-pixel sampled amount map."""
    w = max(1, int(width))
    h = max(1, int(height))
    sampled = np.clip(np.asarray(sampled_amount, dtype=np.float32), 0.0, 1.0)
    if sampled.shape != (h, w):
        raise ValueError("sampled_amount shape must match the requested output size")

    ref_w, ref_h = (
        (w, h)
        if density_reference_size is None
        else (max(1, int(density_reference_size[0])), max(1, int(density_reference_size[1])))
    )

    diameter_px = float(min(w, h))
    stripes = _scaled_cloud_target_stripes(target_stripes, ref_w, ref_h)
    wf = max(0.01, float(width_factor))
    period = int(np.clip(round(diameter_px / stripes), 14, 64))
    max_band = max(1.0, float(period) * wf)

    if geometry is None:
        cx = (w - 1) * 0.5
        cy = (h - 1) * 0.5
        rr = max(1.0, min(cx, cy))
    else:
        cx = float(geometry.center[0])
        cy = float(geometry.center[1])
        rr = max(1.0, float(geometry.radius))

    out = np.zeros((h, w, 4), dtype=np.uint8)
    phase, line_mask, inside_disc, _ = _stripe_render_grids(
        w,
        h,
        period,
        max_band,
        cx,
        cy,
        rr,
        1,
        1,
        edge_fov_deg,
        _cloud_render_content_fov_deg(content_fov_deg),
    )
    draw_mask = inside_disc & line_mask & (phase <= max_band)
    if not np.any(draw_mask):
        return out

    nonzero = sampled[sampled > 0.0]
    if nonzero.size > 0:
        nonzero_lo = float(np.percentile(nonzero, 12.0))
        nonzero_hi = float(np.percentile(nonzero, 92.0))
        if nonzero_hi <= nonzero_lo + 1e-6:
            nonzero_lo = float(nonzero.min())
            nonzero_hi = float(nonzero.max())
        if nonzero_hi > nonzero_lo + 1e-6:
            normalized = (sampled - nonzero_lo) / (nonzero_hi - nonzero_lo)
        else:
            normalized = sampled
    else:
        normalized = sampled
    normalized = np.clip(normalized, 0.0, 1.0)

    alpha_scale = float(np.clip(hatch_cfg.strength, 0, 255)) / 255.0
    alpha = np.zeros((h, w), dtype=np.float32)
    alpha[draw_mask] = normalized[draw_mask] * 255.0 * alpha_scale
    positive = alpha > 0.5
    if not np.any(positive):
        return out

    out[..., :3][positive] = 255
    out[..., 3] = np.clip(np.round(alpha), 0, 255).astype(np.uint8)
    return out


def _render_variable_width_cloud_stripes_rgba_from_altaz_grid(
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry | None = None,
    *,
    projection: ViewProjection,
    target_stripes: int = 50,
    width_factor: float = 0.85,
    density_reference_size: tuple[int, int] | None = None,
    stripe_rgb: tuple[int, int, int] = (255, 255, 255),
) -> np.ndarray:
    """Render variable-width cloud stripes directly from a `CloudAltAzGrid`."""
    sampled_amount = _sample_altaz_grid_to_screen_map(
        grid,
        width,
        height,
        geometry,
        projection,
    )
    return _render_variable_width_cloud_stripes_rgba_from_amount_map(
        sampled_amount,
        width,
        height,
        hatch_cfg,
        geometry=geometry,
        target_stripes=target_stripes,
        width_factor=width_factor,
        density_reference_size=density_reference_size,
        edge_fov_deg=float(projection.edge_fov_deg),
        content_fov_deg=float(projection.content_fov_deg),
        stripe_rgb=stripe_rgb,
    )


def _render_alpha_scaled_cloud_stripes_rgba_from_altaz_grid(
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry | None = None,
    *,
    projection: ViewProjection,
    target_stripes: int = 50,
    width_factor: float = 0.2,
    density_reference_size: tuple[int, int] | None = None,
) -> np.ndarray:
    """Render alpha-scaled cloud stripes directly from a `CloudAltAzGrid`."""
    sampled_amount = _sample_altaz_grid_to_screen_map(
        grid,
        width,
        height,
        geometry,
        projection,
    )
    return _render_alpha_scaled_cloud_stripes_rgba_from_amount_map(
        sampled_amount,
        width,
        height,
        hatch_cfg,
        geometry=geometry,
        target_stripes=target_stripes,
        width_factor=width_factor,
        density_reference_size=density_reference_size,
        edge_fov_deg=float(projection.edge_fov_deg),
        content_fov_deg=float(projection.content_fov_deg),
    )


def _inverse_project_points(
    xs: np.ndarray,
    ys: np.ndarray,
    cx: float,
    cy: float,
    radius: float,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project screen points to (alt, az).

    Returns (alts_deg, azs_deg, inside_mask) for points within the content FOV.
    Points outside the disc receive NaN.
    """
    nx = (np.asarray(xs, dtype=np.float64) - cx) / radius
    ny = (np.asarray(ys, dtype=np.float64) - cy) / radius
    rr2 = nx * nx + ny * ny
    max_r = max(0.0, float(content_fov_deg) / max(1.0e-6, float(edge_fov_deg)))
    inside = rr2 <= (max_r * max_r)

    alts = np.full_like(xs, np.nan, dtype=np.float64)
    azs = np.full_like(ys, np.nan, dtype=np.float64)

    if not np.any(inside):
        return alts.astype(np.float32), azs.astype(np.float32), inside

    r = np.sqrt(rr2[inside])
    theta = np.radians(r * max(1.0e-6, float(edge_fov_deg)))
    psi = np.arctan2(nx[inside], -ny[inside])

    alt_c, az_c = view_center
    eps = 1e-3
    phi1 = np.float64(math.radians(np.clip(float(alt_c), -90.0 + eps, 90.0 - eps)))
    lam1 = np.float64(math.radians(float(az_c)))

    sin_phi1 = np.sin(phi1)
    cos_phi1 = np.cos(phi1)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    sin_phi2 = sin_phi1 * cos_theta + cos_phi1 * sin_theta * np.cos(psi)
    sin_phi2 = np.clip(sin_phi2, -1.0, 1.0)
    phi2 = np.arcsin(sin_phi2)

    y_val = np.sin(psi) * sin_theta * cos_phi1
    x_val = cos_theta - sin_phi1 * sin_phi2
    lam2 = lam1 + np.arctan2(y_val, x_val)

    alts[inside] = np.degrees(phi2)
    azs[inside] = np.degrees(lam2) % 360.0
    return alts.astype(np.float32), azs.astype(np.float32), inside


def _render_halftone_cloud_rgba_from_altaz_grid(
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry | None = None,
    *,
    projection: ViewProjection,
    target_stripes: int = 100,
    width_factor: float = 1.0,
    density_reference_size: tuple[int, int] | None = None,
    grid_phase: tuple[float, float] = (0.0, 0.0),
    grid_scale: float = 1.0,
) -> np.ndarray:
    """Render quantized halftone cloud circles/chains from a `CloudAltAzGrid`.

    Cloud amount is linearly quantized into 8 levels (0/4/8/12/16/20/24/28 px, scaled by delta/30/√2) on a
    screen-fixed square 2D grid rotated 45 degrees (u = x - y, v = x + y).
    Each grid cell with level > 0 is drawn as a circle whose diameter
    equals the quantized level.  Each non-empty grid cell is rendered
    as an individual circle rather than connecting adjacent cells,
    producing a field of round dots (halftone).
    """
    w = max(1, int(width))
    h = max(1, int(height))

    if geometry is None:
        cx = (w - 1) * 0.5
        cy = (h - 1) * 0.5
        rr = max(1.0, min(cx, cy))
    else:
        cx = float(geometry.center[0])
        cy = float(geometry.center[1])
        rr = max(1.0, float(geometry.radius))

    # Grid spacing — square cells (delta_u == delta_v)
    # Use the actual output diameter for grid spacing, but enforce a minimum
    # spacing so compact windows do not collapse the halftone into clutter.
    output_diameter = float(min(w, h))
    delta = _halftone_grid_delta(
        output_diameter,
        target_stripes,
        min_grid_delta_px=HALFTONE_MIN_GRID_DELTA_PX * max(0.01, float(grid_scale)),
    )
    delta_u = delta
    delta_v = delta
    phase_u, phase_v = (float(grid_phase[0]), float(grid_phase[1]))

    # Circle diameters per quantized level scale with grid spacing.
    level_diameters = _halftone_level_diameters(delta, width_factor)

    # Disc geometry
    edge_fov = float(projection.edge_fov_deg)
    content_fov = float(projection.content_fov_deg)
    max_r = max(0.0, content_fov / max(1.0e-6, edge_fov))
    disc_radius = rr * max_r
    # Disc center in (u,v) = (x-y, x+y) space
    u_disc = cx - cy
    v_disc = cx + cy

    # Margin for circle radius extending beyond cell center.
    margin = max(level_diameters) * 0.5 + 2.0
    # A phase-shifted grid can leave its nearest cell center farther from the
    # edge than the circle radius. Generate a small extra ring of centers and
    # let the final painter clip it to the actual empty-sky disc.
    sample_margin = max(margin, delta * 0.75)
    sampling_disc_radius = disc_radius + sample_margin
    sampling_disc_radius_sq = sampling_disc_radius * sampling_disc_radius
    disc_radius_uv = sampling_disc_radius * math.sqrt(2.0)

    # Grid index ranges clipped to viewport extent.
    # Expand the accepted sampling area slightly so circles near the screen edge
    # are still generated even when their centers fall just outside the FOV disc.
    u_min = max(-(h - 1), u_disc - disc_radius_uv - sample_margin)
    u_max = min(w - 1, u_disc + disc_radius_uv + sample_margin)
    v_min = max(0.0, v_disc - disc_radius_uv - sample_margin)
    v_max = min(w + h - 2, v_disc + disc_radius_uv + sample_margin)

    i_min = int(math.floor(u_min / delta_u))
    i_max = int(math.ceil(u_max / delta_u))
    j_min = int(math.floor(v_min / delta_v))
    j_max = int(math.ceil(v_max / delta_v))

    # Collect cell centers inside the disc
    cell_xs: list[float] = []
    cell_ys: list[float] = []
    cell_meta: list[tuple[int, int]] = []  # (i, j)

    for i in range(i_min, i_max + 1):
        u_line = (i + phase_u) * delta_u
        for j in range(j_min, j_max + 1):
            v_center = (j + 0.5 + phase_v) * delta_v
            x = (u_line + v_center) * 0.5
            y = (v_center - u_line) * 0.5

            dx = x - cx
            dy = y - cy
            if dx * dx + dy * dy > sampling_disc_radius_sq:
                continue

            cell_xs.append(x)
            cell_ys.append(y)
            cell_meta.append((i, j))

    if not cell_xs:
        out = np.zeros((h, w, 4), dtype=np.uint8)
        return out

    # Batch inverse-project
    view_center = (float(projection.view_center[0]), float(projection.view_center[1]))
    cloud_content_fov = _cloud_render_content_fov_deg(content_fov)
    alts, azs, inside = _inverse_project_points(
        np.array(cell_xs, dtype=np.float64),
        np.array(cell_ys, dtype=np.float64),
        cx,
        cy,
        rr,
        view_center,
        edge_fov,
        cloud_content_fov,
    )

    # Batch sample cloud amount
    amounts = _sample_altaz_grid_amount(grid, alts, azs)

    # Percentile normalization (matching width-style approach).
    # Collect non-zero amounts from cells inside the disc, compute
    # 12%ile / 92%ile, then normalize so thresholds adapt to the
    # actual data range instead of using absolute values.
    _raw = np.asarray(amounts, dtype=np.float64)
    _inside_mask = inside & (_raw > 0.03)
    if np.any(_inside_mask):
        _nonzero = _raw[_inside_mask]
        _lo = float(np.percentile(_nonzero, 12.0))
        _hi = float(np.percentile(_nonzero, 92.0))
        if _hi <= _lo + 1e-6:
            _lo = float(_nonzero.min())
            _hi = float(_nonzero.max())
        if _hi > _lo + 1e-6:
            _span = _hi - _lo
            _normalized = np.clip((_raw - _lo) / _span, 0.0, 1.0)
        else:
            _normalized = np.clip(_raw, 0.0, 1.0)
    else:
        _normalized = np.clip(_raw, 0.0, 1.0)

    # Quantization thresholds: nonlinear, powers of 1.4 (1.4 ** -i).
    # Suppresses low-amount noise while spreading high amounts across
    # more levels.
    # level 0: [0.000, 0.095), 1: [0.095, 0.133), 2: [0.133, 0.186),
    #       3: [0.186, 0.260), 4: [0.260, 0.364), 5: [0.364, 0.510),
    #       6: [0.510, 0.714), 7: [0.714, 1.000]
    t0 = 1.4 ** -7
    t1 = 1.4 ** -6
    t2 = 1.4 ** -5
    t3 = 1.4 ** -4
    t4 = 1.4 ** -3
    t5 = 1.4 ** -2
    t6 = 1.4 ** -1

    # Build cell-level map: (i, j) -> level (0-7)
    cell_level: dict[tuple[int, int], int] = {}
    for k in range(len(cell_xs)):
        if not inside[k]:
            continue
        amount = float(_normalized[k])
        if amount < t0:
            level = 0
        elif amount < t1:
            level = 1
        elif amount < t2:
            level = 2
        elif amount < t3:
            level = 3
        elif amount < t4:
            level = 4
        elif amount < t5:
            level = 5
        elif amount < t6:
            level = 6
        else:
            level = 7
        if level > 0:
            cell_level[cell_meta[k]] = level

    if not cell_level:
        out = np.zeros((h, w, 4), dtype=np.uint8)
        return out

    # Draw each non-empty grid cell as an individual circle.
    circles: list[tuple[float, float, float]] = []  # (x, y, diameter)

    for k in range(len(cell_xs)):
        key = cell_meta[k]
        if key not in cell_level:
            continue
        i, j = key
        level = cell_level[key]
        diam = level_diameters[level]
        x = ((i + phase_u) * delta_u + (j + 0.5 + phase_v) * delta_v) * 0.5
        y = ((j + 0.5 + phase_v) * delta_v - (i + phase_u) * delta_u) * 0.5
        circles.append((x, y, diam))

    # Render with QPainter
    image = QImage(w, h, QImage.Format_ARGB32_Premultiplied)
    image.fill(Qt.GlobalColor.transparent)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

    clip_path = QPainterPath()
    clip_path.addEllipse(QPointF(cx, cy), disc_radius, disc_radius)
    painter.setClipPath(clip_path)

    alpha = int(np.clip(hatch_cfg.strength, 0, 255))
    base_color = QColor(255, 255, 255, alpha)

    # Draw each cloud cell as a uniform circle.  Cloud amount is represented
    # by the quantized diameter, not by a color or alpha gradient.
    for x, y, diam in circles:
        radius = max(0.5, diam * 0.5)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QBrush(base_color))
        painter.drawEllipse(QPointF(x, y), radius, radius)

    painter.end()
    return qimage_to_np_rgba(image)

def render_variable_width_cloud_stripes(
    cloud_amount: CloudAmountField,
    width: int,
    height: int,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry | None = None,
    *,
    target_stripes: int = 50,
    width_factor: float = 0.85,
    density_reference_size: tuple[int, int] | None = None,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float,
) -> QImage:
    out = _render_variable_width_cloud_stripes_rgba(
        cloud_amount,
        width,
        height,
        hatch_cfg,
        geometry=geometry,
        target_stripes=target_stripes,
        width_factor=width_factor,
        density_reference_size=density_reference_size,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )
    return np_rgba_to_qimage(out)




def compose_cloud_over_sky(
    sky_img: QImage,
    cloud_img_rgba: QImage | np.ndarray,
    dest_rect: QRect,
    geometry: ScreenGeometry | None = None,
    *,
    cloud_opacity: float = 1.0,
    gray_mix: float = 1.0,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float = 90.0,
    sun_alt_deg: float | None = None,
    cloud_tint_rgb: tuple[int, int, int] | None = None,
    transparent_sky_rgb: tuple[int, int, int] | None = None,
) -> QImage:
    """Composite cloud over sky with optional gray desaturation behind clouds.

    - sky is blended with a grayscale version where cloud alpha is present
    - cloud color is additively applied with ``cloud_opacity``
    - final image is clipped to a circle
    """
    w, h = dest_rect.width(), dest_rect.height()

    if sky_img.width() != w or sky_img.height() != h:
        sky_img = _scale_qimage_preserving_aspect(sky_img, w, h)

    sky_np = qimage_to_np_rgba(sky_img)
    if isinstance(cloud_img_rgba, QImage):
        if cloud_img_rgba.width() != w or cloud_img_rgba.height() != h:
            cloud_img_rgba = _scale_qimage_preserving_aspect(cloud_img_rgba, w, h)
        cloud_np = qimage_to_np_rgba(cloud_img_rgba)
    else:
        cloud_np = cloud_img_rgba

    if geometry is None:
        cx = (w - 1) * 0.5
        cy = (h - 1) * 0.5
        rr = max(1.0, min(cx, cy))
    else:
        cx = float(geometry.center[0]) - float(dest_rect.x())
        cy = float(geometry.center[1]) - float(dest_rect.y())
        rr = max(1.0, float(geometry.radius))
    max_r = max(0.0, float(content_fov_deg) / max(1.0e-6, float(edge_fov_deg)))
    y, x = np.ogrid[:h, :w]
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    disc_mask = r2 <= ((rr * max_r) + 0.25) ** 2

    if not np.any(sky_np[..., 3]):
        if transparent_sky_rgb is not None:
            sky_np[..., :3][disc_mask] = np.asarray(
                transparent_sky_rgb,
                dtype=np.uint8,
            )
        sky_np[..., 3][disc_mask] = 255

    cop = float(np.clip(cloud_opacity, 0.0, 1.0))
    sky_rgb_u16 = sky_np[..., :3].astype(np.uint16)
    sky_alpha_u16 = sky_np[..., 3].astype(np.uint16)
    r = sky_rgb_u16[..., 0]
    g = sky_rgb_u16[..., 1]
    b = sky_rgb_u16[..., 2]
    gray_u8 = ((77 * r + 150 * g + 29 * b) >> 8).astype(np.uint8)

    # Keep the grayscale underlay tied to the visible cloud strength so low-opacity
    # clouds do not still flatten the sky into a fully gray disc.
    # Use a stronger, independently clipped coefficient for desaturation so
    # low cloud opacity still gives the cloud pattern a visible color shift.
    gray_opacity = float(np.clip(cop * 4.0, 0.0, 1.0))
    a = (
        cloud_np[..., 3].astype(np.float32)
        / 255.0
        * float(np.clip(gray_mix, 0.0, 1.0))
        * gray_opacity
    )
    a8 = (np.clip(a, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint16)
    inv_a8 = 255 - a8

    gray_u16 = gray_u8.astype(np.uint16)
    base_u16 = (inv_a8[:, :, None] * sky_rgb_u16 + a8[:, :, None] * gray_u16[:, :, None]) // 255

    out_alpha_u16 = sky_alpha_u16.copy()
    if cop > 0.0:
        cop_u16 = int(round(cop * 255))
        cloud_rgb_u32 = cloud_np[..., :3].astype(np.uint32)
        tint_rgb = np.asarray(
            cloud_tint_rgb
            if cloud_tint_rgb is not None
            else _cloud_tint_rgb_for_sun_alt(sun_alt_deg),
            dtype=np.uint32,
        )
        cloud_luma = np.mean(cloud_np[..., :3].astype(np.float32), axis=2)
        # A zero-RGB cloud pixel is a transparent/placeholder cloud sample,
        # not the deliberately shaded blue-grey center circle.
        cloud_darkness_u32 = np.where(
            cloud_luma > 1.0,
            np.clip(255.0 - cloud_luma, 0.0, 255.0),
            0.0,
        ).astype(np.uint32)
        cloud_rgb_u32 = (cloud_rgb_u32 * tint_rgb[None, None, :]) // np.uint32(255)
        cloud_a_u32 = cloud_np[..., 3].astype(np.uint32)[:, :, None]
        darken_u32 = (
            base_u16.astype(np.uint32)
            * cloud_darkness_u32[:, :, None]
            * cloud_a_u32
            * np.uint32(cop_u16)
            // np.uint32(255 * 255 * 255)
        )
        base_u32 = np.maximum(
            0,
            base_u16.astype(np.int64) - darken_u32.astype(np.int64),
        ).astype(np.uint32)
        # A shaded cell must not be brightened back to white by the additive
        # cloud pass. Reduce the white-light contribution in proportion to the
        # same darkness that was subtracted from the sky.
        cloud_light_u32 = 255 - cloud_darkness_u32[:, :, None]
        add_u32 = (
            cloud_rgb_u32
            * cloud_a_u32
            * np.uint32(cop_u16)
            * cloud_light_u32
            // np.uint32(255 * 255 * 255)
        )
        out_u16 = base_u32 + add_u32
        np.minimum(out_u16, 255, out=out_u16)
        cloud_alpha_u16 = ((cloud_np[..., 3].astype(np.uint32) * np.uint32(cop_u16)) // np.uint32(255)).astype(np.uint16)
        out_alpha_u16 = np.maximum(out_alpha_u16, cloud_alpha_u16)
    else:
        out_u16 = base_u16

    out = np.empty((h, w, 4), dtype=np.uint8)
    out[..., :3] = out_u16.astype(np.uint8)
    out[..., 3] = out_alpha_u16.astype(np.uint8)

    out[..., 3][~disc_mask] = 0
    out[..., :3][~disc_mask] = 0

    return np_rgba_to_qimage(out)


def _mask_cloud_alpha_by_missing_rgba(
    cloud: np.ndarray,
    missing_mask_alpha: np.ndarray,
) -> np.ndarray:
    """Set cloud alpha to zero where missing-mask alpha is present."""
    h, w = cloud.shape[:2]
    if missing_mask_alpha.shape != (h, w):
        y_idx = np.rint(np.linspace(0, missing_mask_alpha.shape[0] - 1, h)).astype(np.int32)
        x_idx = np.rint(np.linspace(0, missing_mask_alpha.shape[1] - 1, w)).astype(np.int32)
        missing_mask_alpha = missing_mask_alpha[y_idx][:, x_idx]
    cut = missing_mask_alpha > 0
    if np.any(cut):
        cloud[..., 3][cut] = 0
        cloud[..., :3][cut] = 0
    return cloud


def mask_cloud_alpha_by_missing(
    source_img: QImage,
    missing_mask_alpha: np.ndarray,
) -> QImage:
    cloud = qimage_to_np_rgba(
        source_img
        if source_img.format() == QImage.Format_RGBA8888
        else source_img.convertToFormat(QImage.Format_RGBA8888)
    )
    cloud = _mask_cloud_alpha_by_missing_rgba(cloud, missing_mask_alpha)
    return np_rgba_to_qimage(cloud)
