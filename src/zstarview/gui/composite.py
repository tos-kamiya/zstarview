# -*- coding: utf-8 -*-
"""
Sky and cloud compositing utilities and cache.

This module provides:
- Pure helpers to convert rendered cloud RGBA into a compact cloud-amount field.
- A function to composite sky and cloud layers without relying on any global state.
- A small cache class that handles scaling, stripe rendering, compositing, and reuse.
"""
from __future__ import annotations

import colorsys
import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Optional, Tuple, cast

import numpy as np
from PySide6.QtCore import QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPainterPath, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..clouddisc.altaz_grid import CloudAltAzGrid
from ..night_lights import (
    NIGHT_LIGHTS_GLOW_RGB,
    NightLightGlowProfile,
    night_light_strength_factor,
)
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
from ..render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba
from ..types import ScreenGeometry, ViewerData, ViewProjection

NEVER_RISES_GUIDE_WIDTH_SCALE = 4.5
NEVER_RISES_GUIDE_ALPHA_SCALE = 0.5
ALT_RING_DIMALT_SAMPLE_AZ_STEP_DEG = 30.0
HALFTONE_MIN_GRID_DELTA_PX = 22.0
HALFTONE_LEVEL_DIAMETER_BASE_SCALE = 0.9


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


GLOW_MASK_SCALE = 0.25
GLOW_MASK_TINT_RGB = NIGHT_LIGHTS_GLOW_RGB
GLOW_MASK_NOISE_VARIATION = 0.16
GLOW_MASK_NIGHT_LIGHT_HEIGHT_DEG = 30.0
GLOW_MASK_NIGHT_LIGHT_DECAY_RATE = 2.4
GLOW_MASK_NIGHT_LIGHT_ALTITUDE_CROP_ALPHA_THRESHOLD = 1.0e-4
GLOW_MASK_NIGHT_LIGHT_ALTITUDE_CROP_PAD_ROWS = 1


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
    return cast(tuple[int, int, int], tuple(int(round(max(0.0, min(1.0, channel)) * 255.0)) for channel in lifted))


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


def _crop_night_light_alpha_grid_altitude_bins(
    altitude_bins: np.ndarray,
    alpha_grid: np.ndarray,
    *,
    alpha_threshold: float = GLOW_MASK_NIGHT_LIGHT_ALTITUDE_CROP_ALPHA_THRESHOLD,
    pad_rows: int = GLOW_MASK_NIGHT_LIGHT_ALTITUDE_CROP_PAD_ROWS,
) -> tuple[np.ndarray, np.ndarray]:
    """Crop inactive altitude rows while keeping a small pad around active ones."""
    altitude_bins_arr = np.asarray(altitude_bins, dtype=np.float64).reshape(-1)
    alpha_grid_arr = np.asarray(alpha_grid, dtype=np.float64)
    if altitude_bins_arr.ndim != 1 or alpha_grid_arr.ndim != 2:
        return altitude_bins_arr, alpha_grid_arr
    if altitude_bins_arr.size == 0 or alpha_grid_arr.shape[0] != altitude_bins_arr.size:
        return altitude_bins_arr, alpha_grid_arr
    alpha_grid_arr = np.clip(alpha_grid_arr, 0.0, None)
    row_strengths = np.max(alpha_grid_arr, axis=1)
    active_rows = np.flatnonzero(row_strengths > max(0.0, float(alpha_threshold)))
    if active_rows.size == 0:
        return altitude_bins_arr, alpha_grid_arr
    pad = max(0, int(pad_rows))
    start_index = max(0, int(active_rows[0]) - pad)
    end_index = min(alpha_grid_arr.shape[0] - 1, int(active_rows[-1]) + pad)
    cropped_altitude_bins = altitude_bins_arr[start_index : end_index + 1]
    cropped_alpha_grid = alpha_grid_arr[start_index : end_index + 1, :]
    return cropped_altitude_bins, cropped_alpha_grid


def _interp_night_light_alpha_grid(
    profile: NightLightGlowProfile,
    azimuths_deg: np.ndarray,
    altitudes_deg: np.ndarray,
    *,
    alpha_grid: tuple[tuple[float, ...], ...] | None = None,
) -> np.ndarray | None:
    """Interpolate a precomputed night-light alpha field in alt/az space."""
    altitude_bins = np.asarray(getattr(profile, "altitude_bins_deg", ()), dtype=np.float64)
    alpha_grid = np.asarray(
        getattr(profile, "alpha_grid", ()) if alpha_grid is None else alpha_grid,
        dtype=np.float64,
    )
    if (
        altitude_bins.ndim != 1
        or alpha_grid.ndim != 2
        or altitude_bins.size == 0
        or alpha_grid.shape[0] != altitude_bins.size
        or alpha_grid.shape[1] != len(getattr(profile, "samples", ()))
    ):
        return None
    altitude_bins, alpha_grid = _crop_night_light_alpha_grid_altitude_bins(altitude_bins, alpha_grid)
    ordered = sorted(
        getattr(profile, "samples", ()),
        key=lambda sample: float(sample.azimuth_deg) % 360.0,
    )
    if not ordered:
        return None
    sample_az = np.asarray([float(sample.azimuth_deg) % 360.0 for sample in ordered], dtype=np.float64)
    if sample_az.size == 1:
        az_values = alpha_grid[:, 0]
    else:
        order = np.argsort(sample_az)
        sample_az = sample_az[order]
        alpha_grid = alpha_grid[:, order]
        sample_az_ext = np.concatenate([sample_az[-1:] - 360.0, sample_az, sample_az[:1] + 360.0])
        grid_ext = np.concatenate([alpha_grid[:, -1:], alpha_grid, alpha_grid[:, :1]], axis=1)
        az = np.asarray(azimuths_deg, dtype=np.float64).reshape(-1) % 360.0
        az_interp = np.empty((altitude_bins.size, az.size), dtype=np.float64)
        for alt_index in range(altitude_bins.size):
            az_interp[alt_index, :] = np.interp(az, sample_az_ext, grid_ext[alt_index, :])
        alt = np.asarray(altitudes_deg, dtype=np.float64).reshape(-1)
        return np.asarray(
            [
                np.interp(
                    float(target_alt),
                    altitude_bins,
                    az_interp[:, index],
                    left=0.0,
                    right=0.0,
                )
                for index, target_alt in enumerate(alt.tolist())
            ],
            dtype=np.float64,
        )

    alt = np.asarray(altitudes_deg, dtype=np.float64).reshape(-1)
    return np.interp(alt, altitude_bins, az_values, left=0.0, right=0.0).astype(np.float64, copy=False)


def _night_light_ray_alpha_field(
    *,
    profile: NightLightGlowProfile,
    viewer_data: ViewerData,
    width: int,
    height: int,
    geometry: ScreenGeometry,
    opacity: float,
    sun_alt_deg: float | None,
    alpha_grid: tuple[tuple[float, ...], ...] | None = None,
    apply_vertical_decay_when_grid_brightness: bool = False,
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
        tuple(float(value) for value in viewer_data.view_center),
        edge_fov_deg=float(viewer_data.edge_fov_deg),
        content_fov_deg=float(viewer_data.content_fov_deg),
    )
    if alt_deg.size == 0 or not np.any(inside):
        return np.zeros((height, width), dtype=np.float32)

    alpha = np.zeros((height, width), dtype=np.float32)
    inside_az = np.asarray(az_deg, dtype=np.float32)
    inside_alt = np.asarray(alt_deg, dtype=np.float32)
    main_horizon_source = [
        (float(sample.horizon_alt_deg), float(sample.azimuth_deg))
        for sample in profile.samples
    ]
    main_horizon_alt = _interpolate_terrain_horizon_altitude(inside_az, main_horizon_source)
    night_horizon_alt = main_horizon_alt
    grid_brightness = _interp_night_light_alpha_grid(
        profile,
        inside_az,
        inside_alt,
        alpha_grid=alpha_grid,
    )
    brightness = (
        grid_brightness
        if grid_brightness is not None
        else _circular_interp_profile_samples(profile.samples, inside_az, value_attr="strength")
    )
    sun_factor = 1.0 if sun_alt_deg is None else float(night_light_strength_factor(sun_alt_deg))
    night_layer_opacity = max(0.0, min(1.0, float(opacity)))
    if night_layer_opacity <= 0.0 or sun_factor <= 0.0:
        return alpha

    night_above_horizon = inside_alt - night_horizon_alt
    main_height = max(1.0e-6, float(GLOW_MASK_NIGHT_LIGHT_HEIGHT_DEG))
    night_horizon_factor = np.ones_like(night_above_horizon, dtype=np.float32)
    main_height_ratio = np.clip(np.maximum(night_above_horizon, 0.0) / main_height, 0.0, 1.0)
    main_vertical_falloff = (
        np.exp(-float(GLOW_MASK_NIGHT_LIGHT_DECAY_RATE) * main_height_ratio)
        if grid_brightness is None or apply_vertical_decay_when_grid_brightness
        else np.ones_like(main_height_ratio, dtype=np.float32)
    )
    glow_alpha = np.clip(
        sun_factor
        * np.clip(brightness, 0.0, 1.0)
        * (night_layer_opacity * night_horizon_factor * main_vertical_falloff),
        0.0,
        1.0,
    ).astype(np.float32, copy=False)
    inside_idx = np.flatnonzero(inside)
    alpha_flat = alpha.reshape(-1)
    alpha_flat[inside_idx] = glow_alpha
    return alpha


def _night_light_edge_ray_alpha_field(
    *,
    profile: NightLightGlowProfile,
    viewer_data: ViewerData,
    width: int,
    height: int,
    geometry: ScreenGeometry,
    opacity: float,
    sun_alt_deg: float | None,
) -> np.ndarray:
    edge_grid = getattr(profile, "edge_alpha_grid", ())
    if not edge_grid:
        return np.zeros((max(1, int(height)), max(1, int(width))), dtype=np.float32)
    return _night_light_ray_alpha_field(
        profile=profile,
        viewer_data=viewer_data,
        width=width,
        height=height,
        geometry=geometry,
        opacity=opacity,
        sun_alt_deg=sun_alt_deg,
        alpha_grid=edge_grid,
        apply_vertical_decay_when_grid_brightness=True,
    )


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
    viewer_data: ViewerData,
    night_light_glow_profile: NightLightGlowProfile | None,
    night_light_opacity: float,
    night_light_sun_alt_deg: float | None,
    ridge_glow_opacity: float = 0.03,
    fast_mode: bool = False,
    scale: float = GLOW_MASK_SCALE,
) -> GlowMask | None:
    return _build_glow_mask_for_grid(
        width=width,
        height=height,
        geometry=geometry,
        viewer_data=viewer_data,
        night_light_glow_profile=night_light_glow_profile,
        night_light_opacity=night_light_opacity,
        night_light_sun_alt_deg=night_light_sun_alt_deg,
        ridge_glow_opacity=ridge_glow_opacity,
        fast_mode=fast_mode,
        scale=scale,
        alpha_grid_attr="alpha_grid",
    )


def _build_edge_glow_mask(
    *,
    width: int,
    height: int,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    night_light_glow_profile: NightLightGlowProfile | None,
    ridge_glow_opacity: float,
    night_light_sun_alt_deg: float | None,
    fast_mode: bool = False,
    scale: float = GLOW_MASK_SCALE,
) -> GlowMask | None:
    return _build_glow_mask_for_grid(
        width=width,
        height=height,
        geometry=geometry,
        viewer_data=viewer_data,
        night_light_glow_profile=night_light_glow_profile,
        night_light_opacity=ridge_glow_opacity,
        night_light_sun_alt_deg=night_light_sun_alt_deg,
        ridge_glow_opacity=ridge_glow_opacity,
        fast_mode=fast_mode,
        scale=scale,
        alpha_grid_attr="edge_alpha_grid",
    )


def _build_glow_mask_for_grid(
    *,
    width: int,
    height: int,
    geometry: ScreenGeometry,
    viewer_data: ViewerData,
    night_light_glow_profile: NightLightGlowProfile | None,
    night_light_opacity: float,
    night_light_sun_alt_deg: float | None,
    ridge_glow_opacity: float = 0.03,
    fast_mode: bool,
    scale: float,
    alpha_grid_attr: str,
) -> GlowMask | None:
    effective_opacity = (
        float(ridge_glow_opacity)
        if alpha_grid_attr == "edge_alpha_grid"
        else float(night_light_opacity)
    )
    if (
        night_light_glow_profile is None
        or not night_light_glow_profile.samples
        or effective_opacity <= 0.0
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
    alpha_grid = getattr(night_light_glow_profile, alpha_grid_attr, ())
    if alpha_grid_attr == "edge_alpha_grid":
        alpha = _night_light_edge_ray_alpha_field(
            profile=night_light_glow_profile,
            viewer_data=viewer_data,
            width=low_w,
            height=low_h,
            geometry=low_geometry,
            opacity=effective_opacity,
            sun_alt_deg=night_light_sun_alt_deg,
        )
    else:
        alpha = _night_light_ray_alpha_field(
            profile=night_light_glow_profile,
            viewer_data=viewer_data,
            width=low_w,
            height=low_h,
            geometry=low_geometry,
            opacity=float(night_light_opacity),
            sun_alt_deg=night_light_sun_alt_deg,
            alpha_grid=alpha_grid,
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


def _cloud_render_content_fov_deg(content_fov_deg: float) -> float:
    """Return a slightly expanded FOV for cloud rendering before final clipping."""
    return min(180.0, max(0.0, float(content_fov_deg) + 12.0))


_HALFTONE_GRID_REFERENCE_DIAMETER = 600.0


def _halftone_grid_delta(output_diameter: float, target_stripes: int) -> float:
    """Return the halftone grid spacing in pixels.

    The grid spacing follows the same base-plus-upscale rule as the star
    layer: the spacing scales proportionally with the output diameter and is
    clamped to a minimum so compact windows do not collapse the halftone into
    clutter.
    """
    min_grid_delta_px = HALFTONE_MIN_GRID_DELTA_PX
    diameter = max(1.0, float(output_diameter))
    return max(min_grid_delta_px, diameter / max(1, int(target_stripes)))


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
    line_index = np.floor(phase - 0.5).astype(np.int32, copy=False)
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
    normalized = np.clip(normalized, 0.0, 1.0)
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
    view_center: Tuple[float, float],
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
    delta = _halftone_grid_delta(output_diameter, target_stripes)
    delta_u = delta
    delta_v = delta

    # Circle diameters per quantized level scale with grid spacing.
    level_diameters = _halftone_level_diameters(delta, width_factor)

    # Disc geometry
    edge_fov = float(projection.edge_fov_deg)
    content_fov = float(projection.content_fov_deg)
    max_r = max(0.0, content_fov / max(1.0e-6, edge_fov))
    disc_radius = rr * max_r
    disc_radius_sq = disc_radius * disc_radius

    # Disc center in (u,v) = (x-y, x+y) space
    u_disc = cx - cy
    v_disc = cx + cy
    disc_radius_uv = disc_radius * math.sqrt(2.0)

    # Margin for circle radius extending beyond cell center.
    margin = max(level_diameters) * 0.5 + 2.0
    sample_margin = margin

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
        u_line = i * delta_u
        for j in range(j_min, j_max + 1):
            v_center = (j + 0.5) * delta_v
            x = (u_line + v_center) * 0.5
            y = (v_center - u_line) * 0.5

            dx = x - cx
            dy = y - cy
            if dx * dx + dy * dy > disc_radius_sq:
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
        x = (i * delta_u + (j + 0.5) * delta_v) * 0.5
        y = ((j + 0.5) * delta_v - i * delta_u) * 0.5
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

    # Draw circles — drawPoint with RoundCap pen produces a filled circle
    # whose diameter equals the pen width.
    for x, y, diam in circles:
        pen = QPen(base_color, diam, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin)
        painter.setPen(pen)
        painter.drawPoint(QPointF(x, y))

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
        mode = str(cloud_stripe_mode).strip().lower()
        self._cloud_stripe_mode = mode if mode in ("alpha", "halftone") else "width"
        self._missing_tint_rgba: Tuple[int, int, int, int] = (
            int(np.clip(missing_tint_rgba[0], 0, 255)),
            int(np.clip(missing_tint_rgba[1], 0, 255)),
            int(np.clip(missing_tint_rgba[2], 0, 255)),
            int(np.clip(missing_tint_rgba[3], 0, 255)),
        )
        self._ground_tint_opacity = float(np.clip(ground_tint_opacity, 0.0, 1.0))
        self._composited_img: Optional[QImage] = None
        self._composite_key: Optional[Tuple] = None
        self._glow_mask_cache_stamp: Optional[Tuple] = None
        self._glow_mask_cache: GlowMask | None = None
        self._edge_glow_mask_cache_stamp: Optional[Tuple] = None
        self._edge_glow_mask_cache: GlowMask | None = None

    def invalidate(self) -> None:
        self._composite_key = None
        self._composited_img = None
        self._glow_mask_cache_stamp = None
        self._glow_mask_cache = None
        self._edge_glow_mask_cache_stamp = None
        self._edge_glow_mask_cache = None

    @staticmethod
    def _night_light_glow_key(
        night_light_glow_profile: NightLightGlowProfile | None,
        *,
        alpha_grid_attr: str = "alpha_grid",
    ) -> tuple[
        tuple[int, int, tuple[int, int]],
        tuple[tuple[float, float, float], ...],
    ]:
        alpha_grid = np.asarray(getattr(night_light_glow_profile, alpha_grid_attr, ()), dtype=np.float32)
        altitude_bins = tuple(getattr(night_light_glow_profile, "altitude_bins_deg", ()))
        alpha_shape = (
            int(alpha_grid.shape[0]) if alpha_grid.ndim >= 1 else 0,
            int(alpha_grid.shape[1]) if alpha_grid.ndim >= 2 else 0,
        )
        return (
            (
                len(altitude_bins),
                id(getattr(night_light_glow_profile, alpha_grid_attr, ())),
                alpha_shape,
            ),
            tuple(
                (
                    round(float(sample.azimuth_deg) % 360.0, 3),
                    round(float(sample.horizon_alt_deg), 3),
                    round(float(sample.strength), 4),
                )
                for sample in getattr(night_light_glow_profile, "samples", ())
            ),
        )

    def _glow_mask_cache_key(
        self,
        *,
        width: int,
        height: int,
        geometry: ScreenGeometry,
        viewer_data: ViewerData,
        night_light_glow_profile: NightLightGlowProfile | None,
        night_light_opacity: float,
        night_light_sun_alt_deg: float | None,
        fast_mode: bool,
        alpha_grid_attr: str,
        glow_kind: str,
    ) -> tuple[object, ...]:
        return (
            glow_kind,
            int(width),
            int(height),
            tuple(geometry.center),
            int(geometry.radius),
            tuple(float(value) for value in viewer_data.view_center),
            float(viewer_data.content_fov_deg),
            float(viewer_data.edge_fov_deg),
            self._night_light_glow_key(
                night_light_glow_profile,
                alpha_grid_attr=alpha_grid_attr,
            ),
            float(night_light_opacity),
            None if night_light_sun_alt_deg is None else round(float(night_light_sun_alt_deg), 3),
            bool(fast_mode),
            float(GLOW_MASK_SCALE),
        )

    def _resolve_glow_mask(
        self,
        *,
        glow_key: tuple[object, ...],
        width: int,
        height: int,
        geometry: ScreenGeometry,
        viewer_data: ViewerData,
        night_light_glow_profile: NightLightGlowProfile | None,
        night_light_opacity: float,
        night_light_sun_alt_deg: float | None,
        fast_mode: bool,
    ) -> GlowMask | None:
        if self._glow_mask_cache_stamp != glow_key:
            self._glow_mask_cache = _build_glow_mask(
                width=width,
                height=height,
                geometry=geometry,
                viewer_data=viewer_data,
                night_light_glow_profile=night_light_glow_profile,
                night_light_opacity=night_light_opacity,
                night_light_sun_alt_deg=night_light_sun_alt_deg,
                fast_mode=fast_mode,
            )
            self._glow_mask_cache_stamp = glow_key
        return self._glow_mask_cache

    def _resolve_edge_glow_mask(
        self,
        *,
        glow_key: tuple[object, ...],
        width: int,
        height: int,
        geometry: ScreenGeometry,
        viewer_data: ViewerData,
        night_light_glow_profile: NightLightGlowProfile | None,
        ridge_glow_opacity: float,
        night_light_sun_alt_deg: float | None,
        fast_mode: bool,
    ) -> GlowMask | None:
        if self._edge_glow_mask_cache_stamp != glow_key:
            self._edge_glow_mask_cache = _build_edge_glow_mask(
                width=width,
                height=height,
                geometry=geometry,
                viewer_data=viewer_data,
                night_light_glow_profile=night_light_glow_profile,
                ridge_glow_opacity=ridge_glow_opacity,
                night_light_sun_alt_deg=night_light_sun_alt_deg,
                fast_mode=fast_mode,
            )
            self._edge_glow_mask_cache_stamp = glow_key
        return self._edge_glow_mask_cache

    def draw(
        self,
        painter: QPainter,
        geometry: ScreenGeometry,
        sky_img: Optional[QImage],
        *,
        cloud_alpha: float,
        density_reference_size: tuple[int, int] | None = None,
        view_center: Tuple[float, float] = (0.0, 0.0),
        observer_lat_deg: float | None = None,
        observer_lon_deg: float | None = None,
        observer_height_m: float = 0.0,
        cloud_altaz_grid: CloudAltAzGrid | None = None,
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
        night_light_opacity: float = 0.04,
        ridge_glow_opacity: float = 0.03,
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
        glow_viewer_data = ViewerData(
            location=(
                float(observer_lat_deg) if observer_lat_deg is not None else 0.0,
                float(observer_lon_deg) if observer_lon_deg is not None else 0.0,
            ),
            timezone_name="UTC",
            city_name="",
            view_center=tuple(float(value) for value in view_center),
            edge_fov_deg=float(edge_fov_deg),
            content_fov_deg=float(content_fov_deg),
            observer_height_m=float(observer_height_m),
        )
        cloud_projection = ViewProjection(
            view_center=tuple(float(value) for value in view_center),
            edge_fov_deg=float(edge_fov_deg),
            content_fov_deg=float(content_fov_deg),
        )

        sky_ck = int(sky_img.cacheKey()) if sky_img else 0
        altaz_ck = (
            (
                int(cloud_altaz_grid.amount.shape[0]),
                int(cloud_altaz_grid.amount.shape[1]),
                round(float(cloud_altaz_grid.observer_lat), 3),
                round(float(cloud_altaz_grid.observer_lon), 3),
                str(getattr(cloud_altaz_grid.source_key, "satellite", "")),
                str(getattr(cloud_altaz_grid.source_key, "provider", "")),
                str(getattr(cloud_altaz_grid.source_key, "timeslot_utc", "")),
                round(float(cloud_altaz_grid.coverage_ratio), 6),
                round(float(cloud_altaz_grid.grid_resolution_deg), 6),
            )
            if cloud_altaz_grid is not None
            else ()
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
                (
                    len(tuple(getattr(night_light_glow_profile, "altitude_bins_deg", ()))),
                    id(getattr(night_light_glow_profile, "alpha_grid", ())),
                ),
                tuple(
                    (
                        round(float(sample.azimuth_deg) % 360.0, 3),
                        round(float(sample.horizon_alt_deg), 3),
                        round(float(sample.strength), 4),
                    )
                    for sample in getattr(night_light_glow_profile, "samples", ())
                ),
            )
            if night_light_glow_profile is not None
            else ()
        )
        glow_key = self._glow_mask_cache_key(
            width=w,
            height=h,
            geometry=geometry,
            viewer_data=glow_viewer_data,
            night_light_glow_profile=night_light_glow_profile,
            night_light_opacity=night_light_opacity,
            night_light_sun_alt_deg=night_light_sun_alt_deg,
            fast_mode=fast_mode,
            alpha_grid_attr="alpha_grid",
            glow_kind="glow",
        )
        edge_glow_key = self._glow_mask_cache_key(
            width=w,
            height=h,
            geometry=geometry,
            viewer_data=glow_viewer_data,
            night_light_glow_profile=night_light_glow_profile,
            night_light_opacity=ridge_glow_opacity,
            night_light_sun_alt_deg=night_light_sun_alt_deg,
            fast_mode=fast_mode,
            alpha_grid_attr="edge_alpha_grid",
            glow_kind="edge_glow",
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
            glow_key,
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
            edge_glow_key,
            altaz_ck,
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
            def _scaled(qimg: Optional[QImage]) -> Optional[QImage]:
                if qimg is None:
                    return None
                if qimg.width() == w and qimg.height() == h:
                    return qimg
                return _scale_qimage_preserving_aspect(qimg, w, h)

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

            sky_s = _scaled(sky_img)
            if sky_s is None:
                sky_s = _black_disc_image()
            missing_s = missing_mask
            cloud_s: np.ndarray | None = None

            if cloud_alpha > 0.0:
                if cloud_altaz_grid is not None:
                    if self._cloud_stripe_mode == "alpha":
                        cloud_s = _render_alpha_scaled_cloud_stripes_rgba_from_altaz_grid(
                            cloud_altaz_grid,
                            w,
                            h,
                            self._hatch_cfg,
                            geometry=geometry,
                            projection=cloud_projection,
                            target_stripes=self._cloud_target_stripes,
                            width_factor=self._cloud_stripe_width_factor,
                            density_reference_size=density_reference_size,
                        )
                    elif self._cloud_stripe_mode == "halftone":
                        cloud_s = _render_halftone_cloud_rgba_from_altaz_grid(
                            cloud_altaz_grid,
                            w,
                            h,
                            self._hatch_cfg,
                            geometry=geometry,
                            projection=cloud_projection,
                            target_stripes=self._cloud_target_stripes,
                            width_factor=self._cloud_stripe_width_factor,
                            density_reference_size=density_reference_size,
                        )
                    else:
                        cloud_s = _render_variable_width_cloud_stripes_rgba_from_altaz_grid(
                            cloud_altaz_grid,
                            w,
                            h,
                            self._hatch_cfg,
                            geometry=geometry,
                            projection=cloud_projection,
                            target_stripes=self._cloud_target_stripes,
                            width_factor=self._cloud_stripe_width_factor,
                            density_reference_size=density_reference_size,
                        )
                if missing_s is not None and cloud_s is not None:
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
            glow_mask = self._resolve_glow_mask(
                glow_key=glow_key,
                width=w,
                height=h,
                geometry=geometry,
                viewer_data=glow_viewer_data,
                night_light_glow_profile=night_light_glow_profile,
                night_light_opacity=night_light_opacity,
                night_light_sun_alt_deg=night_light_sun_alt_deg,
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
            edge_glow_mask = self._resolve_edge_glow_mask(
                glow_key=edge_glow_key,
                width=w,
                height=h,
                geometry=geometry,
                viewer_data=glow_viewer_data,
                night_light_glow_profile=night_light_glow_profile,
                ridge_glow_opacity=float(ridge_glow_opacity),
                night_light_sun_alt_deg=night_light_sun_alt_deg,
                fast_mode=fast_mode,
            )
            if edge_glow_mask is not None:
                edge_glow_image = _glow_mask_to_qimage(edge_glow_mask, GLOW_MASK_TINT_RGB)
                if not edge_glow_image.isNull():
                    edge_glow_painter = QPainter(composited)
                    try:
                        edge_glow_painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
                        edge_glow_painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Plus)
                        edge_glow_painter.drawImage(QRect(0, 0, w, h), edge_glow_image)
                    finally:
                        edge_glow_painter.end()
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
