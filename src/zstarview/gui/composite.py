# -*- coding: utf-8 -*-
"""
Sky and cloud compositing utilities and cache.

This module provides:
- Pure helpers to convert rendered cloud RGBA into a compact cloud-amount field.
- A function to composite sky and cloud layers without relying on any global state.
- A small cache class that handles scaling, stripe rendering, compositing, and reuse.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from typing import Optional, Tuple, cast

import numpy as np
from PySide6.QtCore import QRect, Qt
from PySide6.QtGui import QImage, QPainter

from ..paths import CLOUD_HATCH_DEFAULT, CLOUD_MISSING_TINT_RGBA, HatchConfig
from ..render.sky_disc import GROUND_TINT_RGB, NEVER_RISES_TINT_RGB, NEVER_RISES_TINT_STRENGTH
from ..types import ScreenGeometry
from ..render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


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
    content_fov_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cache stripe geometry and field sampling grids for a given viewport."""
    w = max(1, int(width))
    h = max(1, int(height))
    xs = np.arange(w, dtype=np.int32)[None, :]
    ys = np.arange(h, dtype=np.int32)[:, None]
    u_pix = xs - ys
    u_mod = np.mod(u_pix, int(period))
    phase = u_mod.astype(np.float32, copy=False) + 0.5
    line_mask = phase <= float(max_band)

    max_r = max(0.0, float(content_fov_deg) / 90.0)
    y, x = np.ogrid[:h, :w]
    inside_disc = ((x - cx) ** 2 + (y - cy) ** 2) <= ((rr * max_r) + 0.25) ** 2
    sample_radius = max(1.0, rr * max_r)
    xn = (x - cx) / sample_radius
    yn = (y - cy) / sample_radius
    u_idx = np.clip((xn - yn + 2.0) * (bins_u / 4.0), 0.0, bins_u - 1).astype(np.int32)
    v_idx = np.clip((xn + yn + 2.0) * (bins_v / 4.0), 0.0, bins_v - 1).astype(np.int32)
    return (phase, line_mask, inside_disc, u_idx * bins_v + v_idx)


def build_cloud_amount_field_from_rgba(
    cloud: np.ndarray,
    *,
    bins: int = 192,
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


def build_cloud_amount_field(cloud_img: QImage, *, bins: int = 192) -> CloudAmountField:
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
    content_fov_deg: float = 90.0,
) -> np.ndarray:
    """Render fixed-opacity cloud stripes whose width increases with cloud amount."""
    w = max(1, int(width))
    h = max(1, int(height))

    diameter_px = float(min(w, h))
    stripes = max(1, int(target_stripes))
    wf = float(np.clip(width_factor, 0.1, 0.95))
    period = int(np.clip(round(diameter_px / stripes), 14, 64))
    max_band = max(2.0, float(period) * wf)

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
        content_fov_deg,
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
    if np.any(full_mask):
        out[..., :3][full_mask] = 255
        out[..., 3][full_mask] = alpha_u8
    if np.any(partial_mask):
        out[..., :3][partial_mask] = 255
        partial_alpha = np.clip(np.round(frac_levels * alpha_u8), 0, alpha_u8).astype(np.uint8)
        out[..., 3][partial_mask] = partial_alpha[partial_mask]
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
    content_fov_deg: float = 90.0,
) -> np.ndarray:
    """Render fixed-width cloud stripes whose alpha follows cloud amount."""
    w = max(1, int(width))
    h = max(1, int(height))

    diameter_px = float(min(w, h))
    stripes = max(1, int(target_stripes))
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
        content_fov_deg,
    )
    draw_mask = inside_disc & line_mask & (phase <= max_band)
    if not np.any(draw_mask):
        return out

    sampled = np.clip(cloud_amount.amount.reshape(-1)[sample_idx], 0.0, 1.0)
    alpha_scale = float(np.clip(hatch_cfg.strength, 0, 255)) / 255.0
    alpha = np.zeros((h, w), dtype=np.float32)
    alpha[draw_mask] = sampled[draw_mask] * 255.0 * alpha_scale
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
    content_fov_deg: float = 90.0,
) -> QImage:
    out = _render_variable_width_cloud_stripes_rgba(
        cloud_amount,
        width,
        height,
        hatch_cfg,
        geometry=geometry,
        target_stripes=target_stripes,
        width_factor=width_factor,
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
    max_r = max(0.0, float(content_fov_deg) / 90.0)
    y, x = np.ogrid[:h, :w]
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    disc_mask = r2 <= ((rr * max_r) + 0.25) ** 2

    if not np.any(sky_np[..., 3]):
        sky_np[..., 3][disc_mask] = 255

    sky_rgb_u16 = sky_np[..., :3].astype(np.uint16)
    sky_alpha_u16 = sky_np[..., 3].astype(np.uint16)
    r = sky_rgb_u16[..., 0]
    g = sky_rgb_u16[..., 1]
    b = sky_rgb_u16[..., 2]
    gray_u8 = ((77 * r + 150 * g + 29 * b) >> 8).astype(np.uint8)

    a = (cloud_np[..., 3].astype(np.float32) / 255.0) * float(np.clip(gray_mix, 0.0, 1.0))
    a8 = (np.clip(a, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint16)
    inv_a8 = 255 - a8

    gray_u16 = gray_u8.astype(np.uint16)
    base_u16 = (inv_a8[:, :, None] * sky_rgb_u16 + a8[:, :, None] * gray_u16[:, :, None]) // 255

    cop = float(np.clip(cloud_opacity, 0.0, 1.0))
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
    content_fov_deg: float = 90.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project square composited pixels up to the requested content FOV."""
    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    radius = max(1.0, float(geometry.radius))
    ys = (np.arange(height, dtype=np.float32) - cy) / radius
    xs = (np.arange(width, dtype=np.float32) - cx) / radius
    nx, ny = np.meshgrid(xs, ys)

    rr2 = nx * nx + ny * ny
    max_r = max(0.0, float(content_fov_deg) / 90.0)
    inside = rr2 <= (max_r * max_r)
    if not np.any(inside):
        return np.array([], dtype=np.float32), np.array([], dtype=np.float32), inside

    r = np.sqrt(rr2[inside]).astype(np.float32)
    theta = r * (np.pi / 2.0)
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


def _never_rises_mask(
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    observer_lat_deg: float | None,
) -> np.ndarray:
    """Return a boolean mask for declinations that never rise at the observer latitude."""
    if observer_lat_deg is None or alt_deg.size == 0:
        return np.zeros_like(alt_deg, dtype=bool)

    lat = float(np.clip(observer_lat_deg, -90.0, 90.0))
    if lat == 0.0:
        return np.zeros_like(alt_deg, dtype=bool)

    lat_rad = math.radians(lat)
    alt_rad = np.radians(alt_deg)
    az_rad = np.radians(az_deg)
    sin_dec = np.sin(alt_rad) * math.sin(lat_rad) + np.cos(alt_rad) * math.cos(lat_rad) * np.cos(az_rad)
    sin_dec = np.clip(sin_dec, -1.0, 1.0)
    dec = np.degrees(np.arcsin(sin_dec))
    if lat > 0.0:
        return dec <= (lat - 90.0)
    return dec >= (lat + 90.0)


def _apply_ground_tint(
    base_img: QImage,
    *,
    geometry: ScreenGeometry,
    view_center: Tuple[float, float],
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
    ground_tint_opacity: float = 1.0,
    observer_lat_deg: float | None = None,
    content_fov_deg: float = 90.0,
) -> QImage:
    """Tint the composited disc below the geometric or terrain horizon."""
    out = qimage_to_np_rgba(
        base_img if base_img.format() == QImage.Format_RGBA8888 else base_img.convertToFormat(QImage.Format_RGBA8888)
    )
    alt, az, inside = _inverse_project_disc(
        out.shape[1],
        out.shape[0],
        geometry,
        view_center,
        content_fov_deg=content_fov_deg,
    )
    if alt.size == 0:
        return np_rgba_to_qimage(out)

    horizon_alt = _interpolate_terrain_horizon_altitude(az, terrain_profile_altaz)
    ground_mask = alt < horizon_alt
    if not np.any(ground_mask):
        return np_rgba_to_qimage(out)
    rgb = out[..., :3][inside].astype(np.float32) / 255.0
    opacity = np.float32(np.clip(ground_tint_opacity, 0.0, 1.0))
    rgb[ground_mask] = GROUND_TINT_RGB[None, :] * opacity
    never_rises = _never_rises_mask(alt, az, observer_lat_deg)
    never_rises_ground = ground_mask & never_rises
    if np.any(never_rises_ground):
        rgb[never_rises_ground] = np.clip(
            rgb[never_rises_ground] + NEVER_RISES_TINT_RGB[None, :] * np.float32(NEVER_RISES_TINT_STRENGTH),
            0.0,
            1.0,
        )
    out[..., :3][inside] = np.clip(np.round(rgb * 255.0), 0, 255).astype(np.uint8)
    return np_rgba_to_qimage(out)


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
        ground_tint_opacity: float = 1.0,
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
        view_center: Tuple[float, float] = (0.0, 0.0),
        observer_lat_deg: float | None = None,
        cloud_amount_field: Optional[CloudAmountField] = None,
        missing_mask: Optional[np.ndarray] = None,
        terrain_profile_altaz: list[tuple[float, float]] | None = None,
        content_fov_deg: float = 90.0,
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
            w,
            h,
            tuple(geometry.center),
            int(geometry.radius),
            float(cloud_alpha),
            float(view_center[0]),
            float(view_center[1]),
            float(content_fov_deg),
            None if observer_lat_deg is None else float(observer_lat_deg),
            hatch_key,
            self._missing_tint_rgba,
            self._ground_tint_opacity,
            self._gray_mix,
            self._cloud_target_stripes,
            self._cloud_stripe_width_factor,
            self._cloud_stripe_mode,
        )

        if self._composite_key != comp_key or self._composited_img is None:
            def _scaled(qimg: Optional[QImage]) -> Optional[QImage]:
                if qimg is None:
                    return None
                if qimg.width() == w and qimg.height() == h:
                    return qimg
                return qimg.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

            def _black_disc_image() -> QImage:
                img = QImage(w, h, QImage.Format_ARGB32_Premultiplied)
                img.fill(Qt.transparent)
                arr = qimage_to_np_rgba(img)
                cy, cx = (h - 1) * 0.5, (w - 1) * 0.5
                cx = float(geometry.center[0]) - float(x)
                cy = float(geometry.center[1]) - float(y)
                rr = max(1.0, float(geometry.radius))
                max_r = max(0.0, float(content_fov_deg) / 90.0)
                yy, xx = np.ogrid[:h, :w]
                disc_mask = ((xx - cx) ** 2 + (yy - cy) ** 2) <= ((rr * max_r) + 0.25) ** 2
                arr[..., 3][disc_mask] = 255
                return np_rgba_to_qimage(arr)

            sky_s = _scaled(sky_img)
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
                        content_fov_deg=content_fov_deg,
                    )
                if missing_s is not None:
                    cloud_s = _mask_cloud_alpha_by_missing_rgba(cloud_s, missing_s)

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
                    content_fov_deg=content_fov_deg,
                )
            composited = _apply_ground_tint(
                composited,
                geometry=geometry,
                view_center=view_center,
                terrain_profile_altaz=terrain_profile_altaz,
                ground_tint_opacity=self._ground_tint_opacity,
                observer_lat_deg=observer_lat_deg,
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
