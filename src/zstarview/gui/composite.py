# -*- coding: utf-8 -*-
"""
Sky and cloud compositing utilities and cache.

This module provides:
- Pure helpers to convert rendered cloud RGBA into a compact cloud-amount field.
- A function to composite sky and cloud layers without relying on any global state.
- A small cache class that handles scaling, stripe rendering, compositing, and reuse.
"""
from __future__ import annotations

import math
import colorsys
from dataclasses import dataclass
from functools import lru_cache
from typing import Optional, Tuple, cast

import numpy as np
from PySide6.QtCore import QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..night_lights import NightLightGlowProfile
from ..night_lights import night_light_strength_factor
from ..paths import (
    CLOUD_HATCH_DEFAULT,
    CLOUD_MISSING_TINT_RGBA,
    PALETTE_NEVER_RISES_GUIDE_RGB,
    HatchConfig,
    ThemeStyle,
)
from ..render.background import (
    dimalt_ring_brightness_score_from_rgba,
    dimalt_ring_pen_color_from_color,
    draw_altitude_ring_overlay,
)
from ..render.earth_guide import (
    draw_earth_guide,
    earth_guide_line_alpha,
)
from ..render.geometry import normalized_to_screen_xy
from ..render.guides import (
    REFERENCE_LINE_FG_WIDTH,
    REFERENCE_LINE_MID_ALPHA,
    REFERENCE_LINE_MID_WIDTH,
    REFERENCE_LINE_OUTER_ALPHA,
    REFERENCE_LINE_OUTER_WIDTH,
    _clip_polyline_to_radius,
    split_by_gaps,
)
from ..render.night_lights import NIGHT_LIGHTS_GLOW_RGB
from ..render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba
from ..types import ScreenGeometry, ViewerData

NEVER_RISES_GUIDE_WIDTH_SCALE = 4.5
NEVER_RISES_GUIDE_ALPHA_SCALE = 0.5
ALT_RING_DIMALT_SAMPLE_AZ_STEP_DEG = 30.0


def _dimalt_ring_color_for_sky_image(
    sky_img: QImage,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    *,
    alt_deg: float,
    edge_fov_deg: float,
) -> QColor | None:
    """Estimate a dimalt stroke color for one altitude ring."""
    if sky_img.isNull() or geometry.radius < 1:
        return None

    width = sky_img.width()
    height = sky_img.height()
    if width <= 0 or height <= 0:
        return None

    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    radius = float(geometry.radius)
    samples: list[tuple[float, QColor]] = []
    az_deg = 0.0
    while az_deg <= 360.0 + 1.0e-6:
        nx, ny = altaz_to_normalized_xy(
            float(alt_deg),
            float(az_deg),
            view_center,
            edge_fov_deg=edge_fov_deg,
        )
        x = int(round(cx + (nx * radius)))
        y = int(round(cy + (ny * radius)))
        if 0 <= x < width and 0 <= y < height:
            color = sky_img.pixelColor(x, y)
            if color.alpha() > 0:
                samples.append(
                    (
                        dimalt_ring_brightness_score_from_rgba(
                            color.red(),
                            color.green(),
                            color.blue(),
                            color.alpha(),
                        ),
                        color,
                    )
                )
        az_deg += ALT_RING_DIMALT_SAMPLE_AZ_STEP_DEG

    if not samples:
        return None

    samples.sort(key=lambda item: item[0])
    mid = len(samples) // 2
    if len(samples) % 2 == 1:
        _, color = samples[mid]
    else:
        _, color = samples[mid - 1]
    return dimalt_ring_pen_color_from_color(color)


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


@dataclass(frozen=True)
class GlowMask:
    """Low-resolution alpha-only glow buffer in normalized screen space."""

    alpha: np.ndarray
    scale: float


GLOW_MASK_SCALE = 0.25
GLOW_MASK_TINT_RGB = NIGHT_LIGHTS_GLOW_RGB
GLOW_MASK_NOISE_VARIATION = 0.16
GLOW_MASK_NIGHT_LIGHT_HEIGHT_DEG = 36.0
GLOW_MASK_NIGHT_LIGHT_DECAY_RATE = 2.4
GLOW_MASK_NIGHT_LIGHT_HORIZON_SIGMA_DEG = 12.0


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


def _lift_rgb_value_to_max(rgb: tuple[int, int, int]) -> tuple[int, int, int]:
    """Raise a tint color to the brightest HSV value while preserving hue and saturation."""
    r, g, b = (max(0, min(255, int(component))) / 255.0 for component in rgb)
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    if v <= 0.0:
        return (0, 0, 0)
    lifted = colorsys.hsv_to_rgb(h, s, 1.0)
    return tuple(int(round(max(0.0, min(1.0, channel)) * 255.0)) for channel in lifted)


def _circular_interp_profile_samples(
    samples: tuple[object, ...] | list[object],
    azimuths_deg: np.ndarray,
    *,
    value_attr: str,
) -> np.ndarray:
    """Interpolate azimuth-sorted profile samples across the 0/360 seam."""
    if not samples:
        return np.zeros_like(np.asarray(azimuths_deg, dtype=np.float64), dtype=np.float64)

    ordered = sorted(
        samples,
        key=lambda sample: float(getattr(sample, "azimuth_deg")) % 360.0,
    )
    sample_az = np.asarray([float(getattr(sample, "azimuth_deg")) % 360.0 for sample in ordered], dtype=np.float64)
    sample_vals = np.asarray([float(getattr(sample, value_attr)) for sample in ordered], dtype=np.float64)
    if sample_az.size == 1:
        return np.full_like(np.asarray(azimuths_deg, dtype=np.float64), float(sample_vals[0]), dtype=np.float64)

    sample_az_ext = np.concatenate([sample_az[-1:] - 360.0, sample_az, sample_az[:1] + 360.0])
    sample_vals_ext = np.concatenate([sample_vals[-1:], sample_vals, sample_vals[:1]])
    return np.interp(np.asarray(azimuths_deg, dtype=np.float64) % 360.0, sample_az_ext, sample_vals_ext)


def _night_light_ray_alpha_field(
    *,
    profile: NightLightGlowProfile,
    width: int,
    height: int,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
    opacity: float,
    sun_alt_deg: float | None,
    edge_fov_deg: float,
    content_fov_deg: float,
) -> np.ndarray:
    """Build a ray-sampled glow alpha field from the night-light profile."""
    alpha = np.zeros((0, 0), dtype=np.float32)
    if not profile.samples:
        return alpha

    width = max(1, int(width))
    height = max(1, int(height))
    alt_deg, az_deg, inside = _inverse_project_disc(
        width,
        height,
        geometry,
        view_center,
        edge_fov_deg=float(edge_fov_deg),
        content_fov_deg=float(content_fov_deg),
    )
    if alt_deg.size == 0 or not np.any(inside):
        return np.zeros((height, width), dtype=np.float32)

    alpha = np.zeros((height, width), dtype=np.float32)
    inside_az = np.asarray(az_deg, dtype=np.float32)
    inside_alt = np.asarray(alt_deg, dtype=np.float32)
    if terrain_secondary_ridges_altaz_layers:
        horizon_alt = _cumulative_max_ridge_altitude(
            terrain_secondary_ridges_altaz_layers,
            inside_az,
        )
    else:
        horizon_source = terrain_profile_altaz if terrain_profile_altaz else [
            (float(sample.horizon_alt_deg), float(sample.azimuth_deg))
            for sample in profile.samples
        ]
        horizon_alt = _interpolate_terrain_horizon_altitude(inside_az, horizon_source)
    brightness = _circular_interp_profile_samples(profile.samples, inside_az, value_attr="strength")
    sun_factor = 1.0 if sun_alt_deg is None else float(night_light_strength_factor(sun_alt_deg))
    layer_opacity = max(0.0, min(1.0, float(opacity)))
    if layer_opacity <= 0.0 or sun_factor <= 0.0:
        return alpha

    above_horizon = inside_alt - horizon_alt
    horizon_sigma = max(1.0e-6, float(GLOW_MASK_NIGHT_LIGHT_HORIZON_SIGMA_DEG))
    horizon_factor = np.exp(-np.abs(above_horizon) / horizon_sigma)
    normalized_height = np.clip(np.maximum(above_horizon, 0.0) / float(GLOW_MASK_NIGHT_LIGHT_HEIGHT_DEG), 0.0, 1.0)
    vertical_falloff = np.exp(-float(GLOW_MASK_NIGHT_LIGHT_DECAY_RATE) * normalized_height)
    glow_alpha = np.clip(
        layer_opacity
        * sun_factor
        * np.clip(brightness, 0.0, 1.0)
        * horizon_factor
        * vertical_falloff,
        0.0,
        1.0,
    ).astype(np.float32, copy=False)
    inside_idx = np.flatnonzero(inside)
    alpha_flat = alpha.reshape(-1)
    alpha_flat[inside_idx] = glow_alpha
    return alpha


def _stable_glow_noise_grid(height: int, width: int) -> np.ndarray:
    """Return a deterministic, screen-fixed noise field in [0, 1]."""
    h = max(0, int(height))
    w = max(0, int(width))
    if h == 0 or w == 0:
        return np.empty((h, w), dtype=np.float32)

    y = np.arange(h, dtype=np.uint32)[:, None]
    x = np.arange(w, dtype=np.uint32)[None, :]
    value = x * np.uint32(374761393) + y * np.uint32(668265263) + np.uint32(362437)
    value ^= value >> np.uint32(13)
    value *= np.uint32(1274126177)
    value ^= value >> np.uint32(16)
    noise = value.astype(np.float32) * (1.0 / 4294967295.0)
    return _smooth_cloud_amount_grid(noise)


def _glow_mask_to_qimage(glow_mask: GlowMask, tint_rgb: tuple[int, int, int]) -> QImage:
    """Convert a low-resolution glow mask into a premultiplied RGBA image."""
    alpha = np.asarray(glow_mask.alpha, dtype=np.float32)
    if alpha.ndim != 2:
        raise ValueError("GlowMask.alpha must be a 2D array")
    alpha = np.clip(alpha, 0.0, 1.0)
    if alpha.size == 0:
        return QImage()

    noise = _stable_glow_noise_grid(alpha.shape[0], alpha.shape[1])
    if noise.shape == alpha.shape:
        alpha = np.clip(
            alpha * (1.0 - (0.5 * GLOW_MASK_NOISE_VARIATION) + (GLOW_MASK_NOISE_VARIATION * noise)),
            0.0,
            1.0,
        )

    base_rgb = _lift_rgb_value_to_max(tint_rgb)
    rgb = np.empty((alpha.shape[0], alpha.shape[1], 3), dtype=np.uint8)
    rgb[:, :, 0] = base_rgb[0]
    rgb[:, :, 1] = base_rgb[1]
    rgb[:, :, 2] = base_rgb[2]
    rgba = np.zeros((alpha.shape[0], alpha.shape[1], 4), dtype=np.uint8)
    rgba[:, :, :3] = rgb
    rgba[:, :, 3] = np.clip(np.round(alpha * 255.0), 0, 255).astype(np.uint8)
    return np_rgba_to_qimage(rgba).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)


def _build_glow_mask(
    *,
    width: int,
    height: int,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    terrain_profile_altaz: list[tuple[float, float]] | None,
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None,
    night_light_glow_profile: NightLightGlowProfile | None,
    night_light_opacity: float,
    ridge_glow_opacity: float,
    night_light_sun_alt_deg: float | None,
    edge_fov_deg: float,
    content_fov_deg: float,
    fast_mode: bool = False,
    scale: float = GLOW_MASK_SCALE,
) -> GlowMask | None:
    if (
        night_light_glow_profile is None
        or not night_light_glow_profile.samples
        or (float(night_light_opacity) <= 0.0 and float(ridge_glow_opacity) <= 0.0)
        or fast_mode
    ):
        return None
    if width <= 0 or height <= 0:
        return None

    mask_scale = max(0.05, float(scale))
    low_w = max(1, int(round(float(width) * mask_scale)))
    low_h = max(1, int(round(float(height) * mask_scale)))
    low_geometry = ScreenGeometry(
        center=(
            max(0, int(round(float(geometry.center[0]) * mask_scale))),
            max(0, int(round(float(geometry.center[1]) * mask_scale))),
        ),
        radius=max(1, int(round(float(geometry.radius) * mask_scale))),
    )
    alpha = _night_light_ray_alpha_field(
        profile=night_light_glow_profile,
        width=low_w,
        height=low_h,
        geometry=low_geometry,
        view_center=view_center,
        terrain_profile_altaz=terrain_profile_altaz,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        opacity=float(night_light_opacity),
        sun_alt_deg=night_light_sun_alt_deg,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )

    if not np.any(alpha > 0.0):
        return None
    return GlowMask(alpha=alpha, scale=mask_scale)


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
        phase = np.abs(u_mod.astype(np.float32, copy=False) + 0.5 - float(center_offset))
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
    """Scale the stripe count with the reference render surface size."""
    base_diameter = 600.0
    reference_diameter = max(
        1.0,
        float(min(max(1, int(reference_width)), max(1, int(reference_height)))),
    )
    scaled = float(max(1, int(target_stripes))) * reference_diameter / base_diameter
    return max(1, int(round(scaled)))


def _cloud_stripe_fade_factor(phase: np.ndarray, fade_span: float) -> np.ndarray:
    """Return a gentle fade curve for variable-width cloud stripes."""
    progress = np.clip((phase - 0.5) / max(1.0, float(fade_span)), 0.0, 1.0)
    return 1.0 - 0.5 * np.square(progress)


def _cloud_render_content_fov_deg(content_fov_deg: float) -> float:
    """Return a slightly expanded FOV for cloud rendering before final clipping."""
    return min(180.0, max(0.0, float(content_fov_deg) + 12.0))


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


def build_cloud_amount_field(cloud_img: QImage, *, bins: int = 320) -> CloudAmountField:
    """Build a compact cloud-amount field from a cloud image in normalized (u, v) space."""
    cloud = qimage_to_np_rgba(
        cloud_img if cloud_img.format() == QImage.Format_RGBA8888 else cloud_img.convertToFormat(QImage.Format_RGBA8888)
    )
    return build_cloud_amount_field_from_rgba(cloud, bins=bins, source_cache_key=int(cloud_img.cacheKey()))


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
    normalized = np.clip(normalized, 0.0, 1.0)
    present = sampled > 0.03
    line_index = np.floor(phase - 0.5).astype(np.int32, copy=False)
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


def _render_alpha_scaled_cloud_stripes_rgba(
    cloud_amount: CloudAmountField,
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
    """Render fixed-width cloud stripes whose alpha follows cloud amount."""
    w = max(1, int(width))
    h = max(1, int(height))
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
    )
    draw_mask = inside_disc & line_mask & (phase <= max_band)
    if not np.any(draw_mask):
        return out

    sampled = np.clip(cloud_amount.amount.reshape(-1)[sample_idx], 0.0, 1.0)
    if cloud_amount.nonzero_hi > cloud_amount.nonzero_lo + 1e-6:
        normalized = (sampled - cloud_amount.nonzero_lo) / (
            cloud_amount.nonzero_hi - cloud_amount.nonzero_lo
        )
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
) -> QImage:
    """Composite cloud over sky with optional gray desaturation behind clouds.

    - sky is blended with a grayscale version where cloud alpha is present
    - cloud color is additively applied with ``cloud_opacity``
    - final image is clipped to a circle
    """
    w, h = dest_rect.width(), dest_rect.height()

    if sky_img.width() != w or sky_img.height() != h:
        sky_img = sky_img.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

    sky_np = qimage_to_np_rgba(sky_img)
    if isinstance(cloud_img_rgba, QImage):
        if cloud_img_rgba.width() != w or cloud_img_rgba.height() != h:
            cloud_img_rgba = cloud_img_rgba.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)
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
    a = (
        cloud_np[..., 3].astype(np.float32)
        / 255.0
        * float(np.clip(gray_mix, 0.0, 1.0))
        * cop
    )
    a8 = (np.clip(a, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint16)
    inv_a8 = 255 - a8

    gray_u16 = gray_u8.astype(np.uint16)
    base_u16 = (inv_a8[:, :, None] * sky_rgb_u16 + a8[:, :, None] * gray_u16[:, :, None]) // 255

    out_alpha_u16 = sky_alpha_u16.copy()
    if cop > 0.0:
        cop_u16 = int(round(cop * 255))
        cloud_rgb_u32 = cloud_np[..., :3].astype(np.uint32)
        cloud_a_u32 = cloud_np[..., 3].astype(np.uint32)[:, :, None]
        add_u32 = (cloud_rgb_u32 * cloud_a_u32 * np.uint32(cop_u16)) // np.uint32(255 * 255)
        out_u16 = base_u16 + add_u32.astype(np.uint16)
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


def apply_altitude_ring_highlights(
    sky_img: QImage,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    *,
    theme: ThemeStyle | None = None,
    edge_fov_deg: float = 90.0,
    altaz_rings_mode: str = "dimalt",
) -> QImage:
    """Add a subtle Alt-ring highlight inside a sky-disc image."""
    if sky_img.isNull():
        return sky_img

    out = sky_img.copy()
    radius = float(geometry.radius)
    if radius < 1.0:
        return out
    if theme is None:
        return out

    painter = QPainter(out)
    try:
        if altaz_rings_mode == "dimalt":
            def dimalt_ring_color_for_alt_deg(alt_deg: float) -> QColor | None:
                return _dimalt_ring_color_for_sky_image(
                    out,
                    geometry,
                    view_center,
                    alt_deg=alt_deg,
                    edge_fov_deg=edge_fov_deg,
                )
            draw_altitude_ring_overlay(
                painter,
                QRectF(0.0, 0.0, out.width(), out.height()),
                geometry,
                view_center=view_center,
                theme=theme,
                edge_fov_deg=edge_fov_deg,
                ring_color_for_alt_deg=dimalt_ring_color_for_alt_deg,
            )
    finally:
        painter.end()
    return out


def overlay_missing_tint(
    base_img: QImage,
    missing_mask_alpha: np.ndarray,
    *,
    tint_rgba: Tuple[int, int, int, int] = CLOUD_MISSING_TINT_RGBA,
) -> QImage:
    """Overlay missing-data regions with a faint yellow solid tint."""
    w, h = base_img.width(), base_img.height()

    out = qimage_to_np_rgba(base_img if base_img.format() == QImage.Format_RGBA8888 else base_img.convertToFormat(QImage.Format_RGBA8888))
    if missing_mask_alpha.shape != (h, w):
        y_idx = np.rint(np.linspace(0, missing_mask_alpha.shape[0] - 1, h)).astype(np.int32)
        x_idx = np.rint(np.linspace(0, missing_mask_alpha.shape[1] - 1, w)).astype(np.int32)
        missing_mask_alpha = missing_mask_alpha[y_idx][:, x_idx]

    # Treat missing-mask alpha as coverage indicator in [0,1].
    m = missing_mask_alpha.astype(np.float32) / 255.0
    valid = m > 0.0
    if not np.any(valid):
        return np_rgba_to_qimage(out)

    tint_r, tint_g, tint_b, tint_a = (int(np.clip(v, 0, 255)) for v in tint_rgba)
    alpha01 = (float(tint_a) / 255.0) * m

    # Blend tint into RGB only where missing-mask is present; preserve source alpha.
    out_rgb = out[..., :3].astype(np.float32)
    tint_rgb = np.array([tint_r, tint_g, tint_b], dtype=np.float32)
    a = np.clip(alpha01[..., None], 0.0, 1.0)
    out_rgb[valid] = out_rgb[valid] * (1.0 - a[valid]) + tint_rgb * a[valid]
    out[..., :3] = np.clip(np.round(out_rgb), 0, 255).astype(np.uint8)
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


def _is_far_size_mismatch(source_w: int, source_h: int, target_w: int, target_h: int) -> bool:
    source_w = max(1, int(source_w))
    source_h = max(1, int(source_h))
    target_w = max(1, int(target_w))
    target_h = max(1, int(target_h))
    ratio = max(
        float(source_w) / float(target_w),
        float(target_w) / float(source_w),
        float(source_h) / float(target_h),
        float(target_h) / float(source_h),
    )
    return ratio >= 1.5


def mask_cloud_alpha_by_missing(
    cloud_img: QImage,
    missing_mask_alpha: np.ndarray,
) -> QImage:
    cloud = qimage_to_np_rgba(
        cloud_img if cloud_img.format() == QImage.Format_RGBA8888 else cloud_img.convertToFormat(QImage.Format_RGBA8888)
    )
    cloud = _mask_cloud_alpha_by_missing_rgba(cloud, missing_mask_alpha)
    return np_rgba_to_qimage(cloud)


def _inverse_project_disc(
    width: int,
    height: int,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    *,
    edge_fov_deg: float,
    content_fov_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project square composited pixels up to the requested content FOV."""
    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    radius = max(1.0, float(geometry.radius))
    ys = (np.arange(height, dtype=np.float32) - cy) / radius
    xs = (np.arange(width, dtype=np.float32) - cx) / radius
    nx, ny = np.meshgrid(xs, ys)

    rr2 = nx * nx + ny * ny
    max_r = max(0.0, float(content_fov_deg) / max(1.0e-6, float(edge_fov_deg)))
    inside = rr2 <= (max_r * max_r)
    if not np.any(inside):
        return np.array([], dtype=np.float32), np.array([], dtype=np.float32), inside

    r = np.sqrt(rr2[inside]).astype(np.float32)
    theta = np.radians(r * max(1.0e-6, float(edge_fov_deg)))
    psi = np.arctan2(nx[inside], -ny[inside])

    alt_c, az_c = view_center
    eps = 1e-3
    phi1 = np.float32(math.radians(np.clip(alt_c, -90.0 + eps, 90.0 - eps)))
    lam1 = np.float32(math.radians(az_c))

    sin_phi1 = np.sin(phi1)
    cos_phi1 = np.cos(phi1)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    sin_phi2 = sin_phi1 * cos_theta + cos_phi1 * sin_theta * np.cos(psi)
    sin_phi2 = np.clip(sin_phi2, -1.0, 1.0)
    phi2 = np.arcsin(sin_phi2)

    y = np.sin(psi) * sin_theta * cos_phi1
    x = cos_theta - sin_phi1 * sin_phi2
    lam2 = lam1 + np.arctan2(y, x)

    alt = np.degrees(phi2).astype(np.float32)
    az = (np.degrees(lam2) + 360.0) % 360.0
    return alt, az.astype(np.float32), inside


def _interpolate_terrain_horizon_altitude(
    azimuth_deg: np.ndarray,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> np.ndarray:
    """Interpolate terrain-horizon altitude for azimuth samples."""
    if not terrain_profile_altaz:
        return np.zeros_like(azimuth_deg, dtype=np.float32)

    profile = np.asarray(terrain_profile_altaz, dtype=np.float32)
    if profile.ndim != 2 or profile.shape[1] != 2 or profile.shape[0] == 0:
        return np.zeros_like(azimuth_deg, dtype=np.float32)

    altitudes = profile[:, 0]
    azimuths = np.mod(profile[:, 1], 360.0)
    order = np.argsort(azimuths)
    azimuths = azimuths[order]
    altitudes = altitudes[order]
    azimuths, unique_idx = np.unique(azimuths, return_index=True)
    altitudes = altitudes[unique_idx]
    if azimuths.size == 0:
        return np.zeros_like(azimuth_deg, dtype=np.float32)

    azimuths_ext = np.concatenate((azimuths[-1:] - 360.0, azimuths, azimuths[:1] + 360.0))
    altitudes_ext = np.concatenate((altitudes[-1:], altitudes, altitudes[:1]))
    return np.interp(azimuth_deg, azimuths_ext, altitudes_ext).astype(np.float32)


def _cumulative_max_ridge_altitude(
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None,
    azimuth_deg: np.ndarray,
) -> np.ndarray:
    """Return the highest secondary-ridge altitude encountered so far per azimuth."""
    if not terrain_secondary_ridges_altaz_layers:
        return np.zeros_like(azimuth_deg, dtype=np.float32)

    cumulative = np.full_like(np.asarray(azimuth_deg, dtype=np.float32), -np.inf, dtype=np.float32)
    saw_any = False
    for layer in terrain_secondary_ridges_altaz_layers:
        if not layer:
            continue
        layer_alt = _interpolate_terrain_horizon_altitude(
            np.asarray(azimuth_deg, dtype=np.float32),
            [(float(alt), float(az) % 360.0) for alt, az in layer],
        )
        cumulative = np.maximum(cumulative, layer_alt)
        saw_any = True
    if not saw_any:
        return np.zeros_like(azimuth_deg, dtype=np.float32)
    return np.where(np.isfinite(cumulative), cumulative, 0.0).astype(np.float32)


def _never_rises_mask(
    alt_deg: np.ndarray,
) -> np.ndarray:
    """Disabled placeholder."""
    return np.zeros_like(alt_deg, dtype=bool)


def _neu_unit_to_altaz(vec: np.ndarray) -> Tuple[float, float]:
    north = float(vec[0])
    east = float(vec[1])
    up = float(np.clip(float(vec[2]), -1.0, 1.0))
    alt_deg = math.degrees(math.asin(up))
    az_deg = math.degrees(math.atan2(east, north)) % 360.0
    return alt_deg, az_deg


def _never_rises_circle_altaz(
    observer_lat_deg: float | None,
) -> list[tuple[float, float]]:
    if observer_lat_deg is None:
        return []
    lat = float(np.clip(observer_lat_deg, -90.0, 90.0))
    if abs(lat) < 1.0e-6:
        return []

    dec_deg = lat - 90.0 if lat > 0.0 else lat + 90.0
    lat_rad = math.radians(lat)
    dec_rad = math.radians(dec_deg)
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_dec = math.sin(dec_rad)
    cos_dec = math.cos(dec_rad)

    circle: list[tuple[float, float]] = []
    step_deg = 4.0
    for theta_deg in range(0, 360 + int(step_deg), int(step_deg)):
        ha_rad = math.radians(float(theta_deg))
        sin_alt = (sin_lat * sin_dec) + (cos_lat * cos_dec * math.cos(ha_rad))
        sin_alt = max(-1.0, min(1.0, sin_alt))
        alt_deg = math.degrees(math.asin(sin_alt))
        y = -cos_dec * math.sin(ha_rad)
        x = (sin_dec * cos_lat) - (cos_dec * sin_lat * math.cos(ha_rad))
        az_deg = math.degrees(math.atan2(y, x)) % 360.0
        circle.append((float(alt_deg), float(az_deg)))
    return circle


def _apply_ground_reset(
    base_img: QImage,
    *,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
    ground_reset_rgba: tuple[int, int, int, int] | None = None,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float,
) -> QImage:
    """Reset the disc below the geometric or terrain horizon to a neutral background."""
    if ground_reset_rgba is None:
        return base_img
    out = qimage_to_np_rgba(
        base_img if base_img.format() == QImage.Format_RGBA8888 else base_img.convertToFormat(QImage.Format_RGBA8888)
    )
    alt, az, inside = _inverse_project_disc(
        out.shape[1],
        out.shape[0],
        geometry,
        view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
    )
    if alt.size == 0:
        return np_rgba_to_qimage(out)

    horizon_alt = _interpolate_terrain_horizon_altitude(az, terrain_profile_altaz)
    ground_mask = alt < horizon_alt
    if not np.any(ground_mask):
        return np_rgba_to_qimage(out)
    rgb = out[..., :3][inside].astype(np.float32) / 255.0
    alpha = out[..., 3][inside].astype(np.float32) / 255.0
    reset_rgb = np.array(ground_reset_rgba[:3], dtype=np.float32) / 255.0
    reset_alpha = float(np.clip(float(ground_reset_rgba[3]) / 255.0, 0.0, 1.0))
    rgb[ground_mask] = reset_rgb[None, :]
    alpha[ground_mask] = reset_alpha
    out[..., :3][inside] = np.clip(np.round(rgb * 255.0), 0, 255).astype(np.uint8)
    out[..., 3][inside] = np.clip(np.round(alpha * 255.0), 0, 255).astype(np.uint8)
    return np_rgba_to_qimage(out)


def _overlay_never_rises_outline(
    base_img: QImage,
    *,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    observer_lat_deg: float | None = None,
    never_rises_opacity: float = 0.2,
    edge_fov_deg: float = 90.0,
    content_fov_deg: float,
) -> QImage:
    """Draw a thin never-rises outline using the historic accent color."""
    if observer_lat_deg is None:
        return base_img
    alpha_u8 = int(np.clip(round(255.0 * float(never_rises_opacity)), 0, 255))
    if alpha_u8 <= 0:
        return base_img
    out = qimage_to_np_rgba(
        base_img if base_img.format() == QImage.Format_RGBA8888 else base_img.convertToFormat(QImage.Format_RGBA8888)
    )
    if observer_lat_deg is None:
        return np_rgba_to_qimage(out)

    circle_altaz = _never_rises_circle_altaz(observer_lat_deg)
    if len(circle_altaz) < 2:
        return np_rgba_to_qimage(out)

    projected: list[tuple[float, float]] = []
    for alt_deg, az_deg in circle_altaz:
        nx, ny = altaz_to_normalized_xy(
            float(alt_deg),
            float(az_deg),
            view_center,
            edge_fov_deg=edge_fov_deg,
        )
        projected.append((float(nx), float(ny)))

    paint_img = np_rgba_to_qimage(out).convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(paint_img)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    try:
        outline_color = np.array(PALETTE_NEVER_RISES_GUIDE_RGB, dtype=np.uint8)
        pen_color = QColor(int(outline_color[0]), int(outline_color[1]), int(outline_color[2]))
        outer = QColor(pen_color)
        outer.setAlpha(int(np.clip(round(REFERENCE_LINE_OUTER_ALPHA * NEVER_RISES_GUIDE_ALPHA_SCALE), 0, 255)))
        mid = QColor(pen_color)
        mid.setAlpha(int(np.clip(round(REFERENCE_LINE_MID_ALPHA * NEVER_RISES_GUIDE_ALPHA_SCALE), 0, 255)))
        fg = QColor(pen_color)
        fg.setAlpha(
            int(
                np.clip(
                    round(255.0 * earth_guide_line_alpha(never_rises_opacity) * NEVER_RISES_GUIDE_ALPHA_SCALE),
                    0,
                    255,
                )
            )
        )

        outer_pen = QPen(outer, REFERENCE_LINE_OUTER_WIDTH * NEVER_RISES_GUIDE_WIDTH_SCALE, Qt.PenStyle.SolidLine)
        outer_pen.setCosmetic(True)
        outer_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        outer_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

        mid_pen = QPen(mid, REFERENCE_LINE_MID_WIDTH * NEVER_RISES_GUIDE_WIDTH_SCALE, Qt.PenStyle.SolidLine)
        mid_pen.setCosmetic(True)
        mid_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        mid_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

        fg_pen = QPen(fg, REFERENCE_LINE_FG_WIDTH * NEVER_RISES_GUIDE_WIDTH_SCALE, Qt.PenStyle.SolidLine)
        fg_pen.setCosmetic(True)
        fg_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        fg_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

        for fragment in split_by_gaps(projected):
            clipped_frags = _clip_polyline_to_radius(
                fragment,
                content_fov_deg / max(1.0e-6, float(edge_fov_deg)),
            )
            for clipped_frag in clipped_frags:
                if len(clipped_frag) < 2:
                    continue
                poly = QPolygonF(
                    [
                        QPointF(*normalized_to_screen_xy(x, y, geometry))
                        for x, y in clipped_frag
                    ]
                )
                painter.setPen(outer_pen)
                painter.drawPolyline(poly)
                painter.setPen(mid_pen)
                painter.drawPolyline(poly)
                painter.setPen(fg_pen)
                painter.drawPolyline(poly)
    finally:
        painter.end()
    return paint_img


def _overlay_earth_guide(
    base_img: QImage,
    *,
    geometry: ScreenGeometry,
    viewer_data: ViewerData | None,
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
    earth_guide_opacity: float = 0.028,
    visibility_boost: float = 1.0,
    fast_mode: bool = False,
) -> QImage:
    if viewer_data is None:
        return base_img
    out = (
        base_img
        if base_img.format() == QImage.Format.Format_ARGB32_Premultiplied
        else base_img.convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
    )
    painter = QPainter(out)
    painter.setRenderHint(QPainter.Antialiasing, True)
    try:
        draw_earth_guide(
            painter,
            geometry=geometry,
            viewer_data=viewer_data,
            terrain_profile_altaz=terrain_profile_altaz,
            earth_guide_opacity=float(earth_guide_opacity),
            visibility_boost=float(visibility_boost),
            fast_mode=fast_mode,
        )
    finally:
        painter.end()
    return out


class SkyCompositorCache:
    """Manage compositing and reuse the last composited image via a cache key."""

    def __init__(
        self,
        *,
        hatch_cfg: HatchConfig = CLOUD_HATCH_DEFAULT,
        gray_mix: float = 1.0,
        cloud_target_stripes: int = 50,
        cloud_stripe_width_factor: float = 0.85,
        cloud_stripe_mode: str = "width",
        missing_tint_rgba: Tuple[int, int, int, int] = CLOUD_MISSING_TINT_RGBA,
        ground_tint_opacity: float = 0.04,
    ) -> None:
        self._hatch_cfg = hatch_cfg
        self._gray_mix = gray_mix
        self._cloud_target_stripes = max(1, int(cloud_target_stripes))
        self._cloud_stripe_width_factor = max(0.01, float(cloud_stripe_width_factor))
        self._cloud_stripe_mode = "alpha" if str(cloud_stripe_mode) == "alpha" else "width"
        self._missing_tint_rgba: Tuple[int, int, int, int] = (
            int(np.clip(missing_tint_rgba[0], 0, 255)),
            int(np.clip(missing_tint_rgba[1], 0, 255)),
            int(np.clip(missing_tint_rgba[2], 0, 255)),
            int(np.clip(missing_tint_rgba[3], 0, 255)),
        )
        self._ground_tint_opacity = float(np.clip(ground_tint_opacity, 0.0, 1.0))
        self._composited_img: Optional[QImage] = None
        self._composite_key: Optional[Tuple] = None

    def invalidate(self) -> None:
        self._composite_key = None
        self._composited_img = None

    def draw(
        self,
        painter: QPainter,
        geometry: ScreenGeometry,
        sky_img: Optional[QImage],
        cloud_img: Optional[np.ndarray | QImage],
        *,
        cloud_alpha: float,
        density_reference_size: tuple[int, int] | None = None,
        view_center: Tuple[float, float] = (0.0, 0.0),
        observer_lat_deg: float | None = None,
        observer_lon_deg: float | None = None,
        observer_height_m: float = 0.0,
        cloud_amount_field: Optional[CloudAmountField] = None,
        missing_mask: Optional[np.ndarray] = None,
        show_guidelines: bool = True,
        terrain_profile_altaz: list[tuple[float, float]] | None = None,
        terrain_profile_distances_m: list[float] | None = None,
        terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
        terrain_secondary_ridges_distances_m_layers: list[list[float]] | None = None,
        night_light_glow_profile: NightLightGlowProfile | None = None,
        terrain_horizon_opacity: float = 0.003,
        earth_guide_opacity: float = 0.028,
        earth_guide_visibility_boost: float = 1.0,
        night_light_opacity: float = 0.07,
        ridge_glow_opacity: float = 0.3,
        night_light_sun_alt_deg: float | None = None,
        never_rises_opacity: float = 0.2,
        ground_reset_rgba: tuple[int, int, int, int] | None = None,
        theme: ThemeStyle | None = None,
        edge_fov_deg: float = 90.0,
        content_fov_deg: float,
        fast_mode: bool = False,
        sky_disc_altaz_rings: str = "dimalt",
    ) -> None:
        """Composite the sky/cloud layers (with cache) and draw into painter."""
        viewport = painter.viewport()
        x = int(viewport.x())
        y = int(viewport.y())
        w = int(viewport.width())
        h = int(viewport.height())

        sky_ck = int(sky_img.cacheKey()) if sky_img else 0
        cloud_ck = (
            int(cloud_img.cacheKey())
            if isinstance(cloud_img, QImage)
            else id(cloud_img)
            if cloud_img is not None
            else 0
        )
        amount_ck = (
            int(cloud_amount_field.source_cache_key)
            if cloud_amount_field is not None
            else cloud_ck
        )
        missing_ck = id(missing_mask) if missing_mask is not None else 0
        terrain_key = (
            tuple((round(float(alt), 3), round(float(az) % 360.0, 3)) for alt, az in terrain_profile_altaz)
            if terrain_profile_altaz
            else ()
        )
        terrain_distance_key = (
            tuple(round(float(distance_m), 3) for distance_m in terrain_profile_distances_m)
            if terrain_profile_distances_m is not None
            else ()
        )
        terrain_secondary_key = (
            tuple(
                tuple((round(float(alt), 3), round(float(az) % 360.0, 3)) for alt, az in layer)
                for layer in terrain_secondary_ridges_altaz_layers
            )
            if terrain_secondary_ridges_altaz_layers
            else ()
        )
        terrain_secondary_distance_key = (
            tuple(
                tuple(round(float(distance_m), 3) for distance_m in layer)
                for layer in terrain_secondary_ridges_distances_m_layers
            )
            if terrain_secondary_ridges_distances_m_layers
            else ()
        )
        night_light_key = (
            (
                tuple(
                    (
                        round(float(sample.azimuth_deg) % 360.0, 3),
                        round(float(sample.horizon_alt_deg), 3),
                        round(float(sample.strength), 4),
                    )
                    for sample in getattr(night_light_glow_profile, "samples", ())
                ),
                tuple(
                    (
                        round(float(band_profile.min_distance_km), 3),
                        round(float(band_profile.max_distance_km), 3),
                        tuple(
                            (
                                round(float(sample.azimuth_deg) % 360.0, 3),
                                round(float(sample.horizon_alt_deg), 3),
                                round(float(sample.strength), 4),
                            )
                            for sample in getattr(band_profile, "samples", ())
                        ),
                    )
                    for band_profile in getattr(night_light_glow_profile, "band_profiles", ())
                ),
            )
            if night_light_glow_profile is not None
            else ()
        )
        hatch_key = (
            self._hatch_cfg.tile_w_px,
            self._hatch_cfg.tile_h_px,
            self._hatch_cfg.line_px,
            self._hatch_cfg.strength,
        )
        comp_key = (
            "comp",
            sky_ck,
            cloud_ck,
            amount_ck,
            missing_ck,
            terrain_key,
            terrain_distance_key,
            terrain_secondary_key,
            terrain_secondary_distance_key,
            w,
            h,
            tuple(geometry.center),
            int(geometry.radius),
            float(cloud_alpha),
            float(view_center[0]),
            float(view_center[1]),
            float(content_fov_deg),
            None if observer_lat_deg is None else float(observer_lat_deg),
            None if observer_lon_deg is None else float(observer_lon_deg),
            float(observer_height_m),
            bool(show_guidelines),
            float(terrain_horizon_opacity),
            float(earth_guide_opacity),
            float(night_light_opacity),
            float(ridge_glow_opacity),
            None if night_light_sun_alt_deg is None else round(float(night_light_sun_alt_deg), 3),
            float(never_rises_opacity),
            bool(fast_mode),
            hatch_key,
            night_light_key,
            self._missing_tint_rgba,
            self._gray_mix,
            None
            if ground_reset_rgba is None
            else tuple(int(np.clip(c, 0, 255)) for c in ground_reset_rgba),
            None
            if density_reference_size is None
            else (
                max(1, int(density_reference_size[0])),
                max(1, int(density_reference_size[1])),
            ),
            self._cloud_target_stripes,
            self._cloud_stripe_width_factor,
            self._cloud_stripe_mode,
            float(GLOW_MASK_SCALE),
            tuple(int(c) for c in GLOW_MASK_TINT_RGB),
            str(sky_disc_altaz_rings),
            None
            if theme is None
            else (
                tuple(int(c) for c in theme.window_background.base_rgb),
                tuple(int(c) for c in theme.window_background.delta_rgb),
                int(theme.window_background.outer_alpha),
                int(theme.window_background.edge_alpha),
                bool(theme.window_background.flat_background),
                bool(theme.window_background.draw_outer_border),
            ),
        )

        if self._composite_key != comp_key or self._composited_img is None:
            def _scaled(qimg: Optional[QImage], *, reject_far_mismatch: bool = False) -> Optional[QImage]:
                if qimg is None:
                    return None
                if qimg.width() == w and qimg.height() == h:
                    return qimg
                if reject_far_mismatch and _is_far_size_mismatch(qimg.width(), qimg.height(), w, h):
                    return None
                return qimg.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

            def _black_disc_image() -> QImage:
                img = QImage(w, h, QImage.Format_ARGB32_Premultiplied)
                img.fill(Qt.transparent)
                arr = qimage_to_np_rgba(img)
                cy, cx = (h - 1) * 0.5, (w - 1) * 0.5
                cx = float(geometry.center[0]) - float(x)
                cy = float(geometry.center[1]) - float(y)
                rr = max(1.0, float(geometry.radius))
                max_r = max(0.0, float(content_fov_deg) / max(1.0e-6, float(edge_fov_deg)))
                yy, xx = np.ogrid[:h, :w]
                disc_mask = ((xx - cx) ** 2 + (yy - cy) ** 2) <= ((rr * max_r) + 0.25) ** 2
                arr[..., 3][disc_mask] = 255
                return np_rgba_to_qimage(arr)

            sky_s = _scaled(sky_img, reject_far_mismatch=True)
            if sky_s is None:
                sky_s = _black_disc_image()
            cloud_s = cloud_img
            if isinstance(cloud_s, QImage):
                cloud_s = qimage_to_np_rgba(_scaled(cloud_s))
            missing_s = missing_mask

            if cloud_s is not None and cloud_alpha > 0.0:
                if cloud_amount_field is None:
                    cloud_amount_field = build_cloud_amount_field_from_rgba(
                        np.array(cloud_s, copy=False),
                        source_cache_key=cloud_ck,
                    )
                if self._cloud_stripe_mode == "alpha":
                    cloud_s = _render_alpha_scaled_cloud_stripes_rgba(
                        cloud_amount_field,
                        w,
                        h,
                        self._hatch_cfg,
                        geometry=geometry,
                        target_stripes=self._cloud_target_stripes,
                        width_factor=self._cloud_stripe_width_factor,
                        density_reference_size=density_reference_size,
                        edge_fov_deg=edge_fov_deg,
                        content_fov_deg=content_fov_deg,
                    )
                else:
                    cloud_s = _render_variable_width_cloud_stripes_rgba(
                        cloud_amount_field,
                        w,
                        h,
                        self._hatch_cfg,
                        geometry=geometry,
                        target_stripes=self._cloud_target_stripes,
                        width_factor=self._cloud_stripe_width_factor,
                        density_reference_size=density_reference_size,
                        edge_fov_deg=edge_fov_deg,
                        content_fov_deg=content_fov_deg,
                    )
                if missing_s is not None:
                    cloud_s = _mask_cloud_alpha_by_missing_rgba(cloud_s, missing_s)
            if sky_disc_altaz_rings == "dimalt" and sky_s is not None:
                sky_s = apply_altitude_ring_highlights(
                    sky_s,
                    geometry,
                    view_center,
                    theme=theme,
                    edge_fov_deg=edge_fov_deg,
                    altaz_rings_mode="dimalt",
                )
                if cloud_s is None or cloud_alpha <= 0.0:
                    composited = sky_s
                else:
                    composited = compose_cloud_over_sky(
                        sky_img=sky_s,
                        cloud_img_rgba=cloud_s,
                        dest_rect=QRect(0, 0, w, h),
                        geometry=geometry,
                        cloud_opacity=cloud_alpha,
                        gray_mix=self._gray_mix,
                        edge_fov_deg=edge_fov_deg,
                        content_fov_deg=content_fov_deg,
                    )
            if cloud_s is None or cloud_alpha <= 0.0:
                composited = sky_s
            else:
                composited = compose_cloud_over_sky(
                    sky_img=sky_s,
                    cloud_img_rgba=cloud_s,
                    dest_rect=QRect(0, 0, w, h),
                    geometry=geometry,
                    cloud_opacity=cloud_alpha,
                    gray_mix=self._gray_mix,
                    edge_fov_deg=edge_fov_deg,
                    content_fov_deg=content_fov_deg,
                )
            composited = _apply_ground_reset(
                composited,
                geometry=geometry,
                view_center=view_center,
                terrain_profile_altaz=terrain_profile_altaz,
                ground_reset_rgba=ground_reset_rgba,
                edge_fov_deg=edge_fov_deg,
                content_fov_deg=content_fov_deg,
            )
            earth_viewer_data = (
                None
                if observer_lat_deg is None or observer_lon_deg is None
                else ViewerData(
                    location=(float(observer_lat_deg), float(observer_lon_deg)),
                    timezone_name="UTC",
                    city_name="",
                    view_center=view_center,
                    edge_fov_deg=float(edge_fov_deg),
                    content_fov_deg=float(content_fov_deg),
                    observer_height_m=float(observer_height_m),
                )
            )
            composited = _overlay_earth_guide(
                composited,
                geometry=geometry,
                viewer_data=earth_viewer_data,
                terrain_profile_altaz=terrain_profile_altaz,
                earth_guide_opacity=earth_guide_opacity,
                visibility_boost=earth_guide_visibility_boost,
                fast_mode=fast_mode,
            )
            glow_mask = _build_glow_mask(
                width=w,
                height=h,
                geometry=geometry,
                view_center=view_center,
                terrain_profile_altaz=terrain_profile_altaz,
                terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
                night_light_glow_profile=night_light_glow_profile,
                night_light_opacity=float(night_light_opacity),
                ridge_glow_opacity=float(ridge_glow_opacity),
                night_light_sun_alt_deg=night_light_sun_alt_deg,
                edge_fov_deg=edge_fov_deg,
                content_fov_deg=content_fov_deg,
                fast_mode=fast_mode,
            )
            if glow_mask is not None:
                glow_image = _glow_mask_to_qimage(glow_mask, GLOW_MASK_TINT_RGB)
                if not glow_image.isNull():
                    glow_painter = QPainter(composited)
                    try:
                        glow_painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
                        glow_painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
                        glow_painter.drawImage(QRect(0, 0, w, h), glow_image)
                    finally:
                        glow_painter.end()
            if show_guidelines:
                composited = _overlay_never_rises_outline(
                    composited,
                    geometry=geometry,
                    view_center=view_center,
                    observer_lat_deg=observer_lat_deg,
                    never_rises_opacity=never_rises_opacity,
                    edge_fov_deg=edge_fov_deg,
                    content_fov_deg=content_fov_deg,
                )
            if missing_s is not None:
                composited = overlay_missing_tint(
                    composited,
                    missing_s,
                    tint_rgba=self._missing_tint_rgba,
                )

            self._composited_img = composited
            self._composite_key = comp_key

        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode_SourceOver)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        painter.drawImage(QRect(x, y, w, h), self._composited_img)
        painter.restore()
