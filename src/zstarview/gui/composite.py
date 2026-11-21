"""
Sky compositing utilities and cache.

This module composites the sky disc with terrain, glow, and guide layers,
reuses the last composited image via a cache key, and draws the precomputed
voxel cloud raster over the sky.
"""
from __future__ import annotations

import colorsys
import math
from collections.abc import Sequence
from dataclasses import dataclass, replace
from functools import lru_cache
from typing import cast

import numpy as np
from PySide6.QtCore import QPointF, QRect, QRectF, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPainterPath, QPen, QPolygonF

from ..astro import altaz_to_normalized_xy
from ..night_lights import (
    NIGHT_LIGHTS_GLOW_RGB,
    NightLightGlowProfile,
    NightLightGlowSample,
    night_light_strength_factor,
)
from ..night_lights_constants import RIDGE_GLOW_RGB
from ..paths import (
    EARTH_GUIDE_DEFAULT_OPACITY,
    NIGHT_LIGHT_DEFAULT_OPACITY,
    PALETTE_NEVER_RISES_GUIDE_RGB,
    RIDGE_GLOW_DEFAULT_OPACITY,
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
from ..render.ground_mask import (
    interpolate_terrain_horizon_altitude as _shared_interpolate_terrain_horizon_altitude,
)
from ..render.ground_mask import (
    inverse_project_disc as _shared_inverse_project_disc,
)
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
from ..render.sky_disc import SKY_DISC_OVERSCAN_DEG, SKY_DISC_RENDER_SCALE
from ..render.stars import StarRenderCache
from ..types import ScreenGeometry, ViewerData, ViewProjection
@dataclass(frozen=True)
class GlowMask:
    """Low-resolution alpha-only glow buffer in normalized screen space."""

    alpha: np.ndarray
    scale: float

NEVER_RISES_GUIDE_WIDTH_SCALE = 4.5
NEVER_RISES_GUIDE_ALPHA_SCALE = 0.5
ALT_RING_DIMALT_SAMPLE_AZ_STEP_DEG = 30.0


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
        Qt.AspectRatioMode.KeepAspectRatio,
        Qt.TransformationMode.SmoothTransformation,
    )
    if scaled.width() == target_w and scaled.height() == target_h:
        return scaled
    canvas = QImage(target_w, target_h, QImage.Format.Format_ARGB32_Premultiplied)
    canvas.fill(Qt.GlobalColor.transparent)
    painter = QPainter(canvas)
    try:
        x = (target_w - scaled.width()) // 2
        y = (target_h - scaled.height()) // 2
        painter.drawImage(x, y, scaled)
    finally:
        painter.end()
    return canvas


def _smooth_unit_grid(values: np.ndarray) -> np.ndarray:
    """Apply a small 3x3 binomial blur and clip the result to [0, 1]."""
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

def _dimalt_ring_color_for_sky_image(
    sky_img: QImage,
    geometry: ScreenGeometry,
    projection: ViewProjection,
    *,
    alt_deg: float,
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
            projection.view_center,
            edge_fov_deg=projection.edge_fov_deg,
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


def _clip_sky_image_to_disc(
    image: QImage,
    *,
    geometry: ScreenGeometry,
    viewport_offset: tuple[int, int],
    edge_fov_deg: float,
    content_fov_deg: float,
) -> QImage:
    """Apply an antialiased final-resolution clip to the sky disc."""
    viewport_x, viewport_y = viewport_offset
    cx = float(geometry.center[0]) - float(viewport_x)
    cy = float(geometry.center[1]) - float(viewport_y)
    radius = float(geometry.radius) * max(
        0.0,
        float(content_fov_deg) / max(1.0e-6, float(edge_fov_deg)),
    )
    clipped = QImage(image.size(), QImage.Format.Format_ARGB32_Premultiplied)
    clipped.fill(Qt.GlobalColor.transparent)
    painter = QPainter(clipped)
    try:
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        clip_path = QPainterPath()
        clip_path.addEllipse(QPointF(cx, cy), radius, radius)
        painter.setClipPath(clip_path)
        painter.drawImage(0, 0, image)
    finally:
        painter.end()
    return clipped


GLOW_MASK_SCALE = 0.25
GLOW_MASK_TINT_RGB = NIGHT_LIGHTS_GLOW_RGB
EDGE_GLOW_MASK_TINT_RGB = RIDGE_GLOW_RGB
GLOW_MASK_NOISE_VARIATION = 0.16
GLOW_MASK_NIGHT_LIGHT_HEIGHT_DEG = 30.0
GLOW_MASK_NIGHT_LIGHT_DECAY_RATE = 2.4
GLOW_MASK_NIGHT_LIGHT_ALTITUDE_CROP_ALPHA_THRESHOLD = 1.0e-4
GLOW_MASK_NIGHT_LIGHT_ALTITUDE_CROP_PAD_ROWS = 1





def _scale_rgb(rgb: tuple[int, int, int], factor: float) -> tuple[int, int, int]:
    """Scale each RGB channel by `factor`, rounding to integers."""
    return (
        int(round(rgb[0] * factor)),
        int(round(rgb[1] * factor)),
        int(round(rgb[2] * factor)),
    )


def _lift_rgb_value_to_max(rgb: tuple[int, int, int]) -> tuple[int, int, int]:
    """Raise a tint color to the brightest HSV value while preserving hue and saturation."""
    r, g, b = (max(0, min(255, int(component))) / 255.0 for component in rgb)
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    if v <= 0.0:
        return (0, 0, 0)
    lifted = colorsys.hsv_to_rgb(h, s, 1.0)
    return cast(tuple[int, int, int], tuple(int(round(max(0.0, min(1.0, channel)) * 255.0)) for channel in lifted))


def _circular_interp_profile_samples(
    samples: Sequence[NightLightGlowSample],
    azimuths_deg: np.ndarray,
    *,
    value_attr: str,
) -> np.ndarray:
    """Interpolate azimuth-sorted profile samples across the 0/360 seam."""
    if not samples:
        return np.zeros_like(np.asarray(azimuths_deg, dtype=np.float64), dtype=np.float64)

    ordered = sorted(
        samples,
        key=lambda sample: float(sample.azimuth_deg) % 360.0,
    )
    sample_az = np.asarray([float(sample.azimuth_deg) % 360.0 for sample in ordered], dtype=np.float64)
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
    alpha_grid: np.ndarray | tuple[tuple[float, ...], ...] | None = None,
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
        ViewProjection(
            view_center=viewer_data.view_center,
            edge_fov_deg=viewer_data.edge_fov_deg,
            content_fov_deg=viewer_data.content_fov_deg,
        ),
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
    return _smooth_unit_grid(noise)


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
    ridge_glow_opacity: float = RIDGE_GLOW_DEFAULT_OPACITY,
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
    ridge_glow_opacity: float = RIDGE_GLOW_DEFAULT_OPACITY,
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


def apply_altitude_ring_highlights(
    sky_img: QImage,
    geometry: ScreenGeometry,
    projection: ViewProjection,
    *,
    theme: ThemeStyle | None = None,
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
                    projection,
                    alt_deg=alt_deg,
                )
            draw_altitude_ring_overlay(
                painter,
                geometry,
                view_center=projection.view_center,
                edge_fov_deg=projection.edge_fov_deg,
                ring_color_for_alt_deg=dimalt_ring_color_for_alt_deg,
            )
    finally:
        painter.end()
    return out


def _inverse_project_disc(
    width: int,
    height: int,
    geometry: ScreenGeometry,
    projection: ViewProjection,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project square composited pixels up to the requested content FOV."""
    return _shared_inverse_project_disc(
        width,
        height,
        geometry,
        projection.view_center,
        edge_fov_deg=projection.edge_fov_deg,
        content_fov_deg=projection.content_fov_deg,
    )

def _interpolate_terrain_horizon_altitude(
    azimuth_deg: np.ndarray,
    terrain_profile_altaz: list[tuple[float, float]] | None,
) -> np.ndarray:
    """Interpolate terrain-horizon altitude for azimuth samples."""
    return _shared_interpolate_terrain_horizon_altitude(
        azimuth_deg,
        terrain_profile_altaz,
    )


def _clip_below_terrain_horizon(
    base_img: QImage,
    *,
    geometry: ScreenGeometry,
    projection: ViewProjection,
    terrain_profile_altaz: list[tuple[float, float]] | None,
    softness_px: float = 1.5,
) -> QImage:
    """Clip sky layers below the terrain horizon with a narrow full-resolution pass."""
    if not terrain_profile_altaz:
        return base_img
    out = qimage_to_np_rgba(
        base_img
        if base_img.format() == QImage.Format.Format_RGBA8888
        else base_img.convertToFormat(QImage.Format.Format_RGBA8888)
    )
    height, width = out.shape[:2]
    coverage = _terrain_horizon_coverage(
        width,
        height,
        geometry.center,
        geometry.radius,
        (float(projection.view_center[0]), float(projection.view_center[1])),
        float(projection.edge_fov_deg),
        float(projection.content_fov_deg),
        tuple((float(alt), float(az)) for alt, az in terrain_profile_altaz),
        float(softness_px),
    )
    affected = coverage < 1.0
    if not np.any(affected):
        return base_img
    pixels = out[affected].astype(np.float32)
    pixels *= coverage[affected, None]
    out[affected] = np.clip(
        np.round(pixels), 0, 255
    ).astype(np.uint8)
    return np_rgba_to_qimage(out)


@lru_cache(maxsize=2)
def _terrain_horizon_coverage(
    width: int,
    height: int,
    geometry_center: tuple[int, int],
    geometry_radius: int,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    terrain_profile_altaz: tuple[tuple[float, float], ...],
    softness_px: float,
) -> np.ndarray:
    """Cache the static per-pixel coverage derived from the terrain horizon."""
    geometry = ScreenGeometry(geometry_center, geometry_radius)
    projection = ViewProjection(view_center, edge_fov_deg, content_fov_deg)
    scale = float(np.clip(SKY_DISC_RENDER_SCALE, 0.01, 1.0))
    low_width = max(1, int(math.ceil(width * scale)))
    low_height = max(1, int(math.ceil(height * scale)))
    scale_x = low_width / float(width)
    scale_y = low_height / float(height)
    low_geometry = ScreenGeometry(
        center=(
            int(round(float(geometry.center[0]) * scale_x)),
            int(round(float(geometry.center[1]) * scale_y)),
        ),
        radius=max(1, int(round(float(geometry.radius) * min(scale_x, scale_y)))),
    )
    low_pixel_guard_deg = (
        2.0
        * float(edge_fov_deg)
        / max(1.0, float(low_geometry.radius))
    )
    low_alt, low_az, low_inside = _inverse_project_disc(
        low_width,
        low_height,
        low_geometry,
        replace(projection, content_fov_deg=content_fov_deg + low_pixel_guard_deg),
    )
    if low_alt.size == 0:
        return np.ones((height, width), dtype=np.float64)
    low_horizon_alt = _interpolate_terrain_horizon_altitude(
        low_az, list(terrain_profile_altaz)
    )
    # Keep two low-resolution pixels on either side of the ridge.  Pixels well
    # below that guard can be discarded from the coarse surface; only the
    # guard is inverse-projected again at display resolution.
    guard_alt_deg = low_pixel_guard_deg
    coarse_ground = np.zeros((low_height, low_width), dtype=bool)
    coarse_guard = np.zeros((low_height, low_width), dtype=bool)
    coarse_ground[low_inside] = low_alt < (low_horizon_alt - guard_alt_deg)
    coarse_guard[low_inside] = (
        (low_alt >= (low_horizon_alt - guard_alt_deg))
        & (low_alt <= (low_horizon_alt + guard_alt_deg))
    )
    # Nearest-neighbour expansion maps one low-resolution cell to a block of
    # display pixels. Expand the guard before that mapping so a ridge crossing
    # a cell corner cannot leave a small unclassified ground wedge behind.
    expanded_guard = coarse_guard.copy()
    for _ in range(2):
        padded = np.pad(expanded_guard, 1, mode="constant", constant_values=False)
        expanded_guard = np.logical_or.reduce(
            tuple(
                padded[dy : dy + low_height, dx : dx + low_width]
                for dy in range(3)
                for dx in range(3)
            )
        )
    coarse_guard = expanded_guard
    coarse_ground &= ~coarse_guard
    low_y = np.minimum(
        low_height - 1,
        np.floor(np.arange(height, dtype=np.float32) * scale_y).astype(np.intp),
    )
    low_x = np.minimum(
        low_width - 1,
        np.floor(np.arange(width, dtype=np.float32) * scale_x).astype(np.intp),
    )
    ground = coarse_ground[low_y][:, low_x]
    guard = coarse_guard[low_y][:, low_x]
    result = np.ones((height, width), dtype=np.float64)
    result[ground] = 0.0

    guard_y, guard_x = np.nonzero(guard)
    if guard_x.size == 0:
        result.setflags(write=False)
        return result
    alt, az, inside = _inverse_project_pixel_coordinates(
        guard_x,
        guard_y,
        geometry=geometry,
        projection=projection,
    )
    horizon_alt = _interpolate_terrain_horizon_altitude(
        az, list(terrain_profile_altaz)
    )
    horizon_softness = max(
        1.0e-4,
        float(softness_px) * float(edge_fov_deg) / max(1.0, float(geometry.radius)),
    )
    coverage = np.clip(
        (alt - horizon_alt + horizon_softness) / (2.0 * horizon_softness),
        0.0,
        1.0,
    )
    coverage = coverage * coverage * (3.0 - 2.0 * coverage)
    coverage[~inside] = 1.0
    result[guard_y, guard_x] = coverage
    result.setflags(write=False)
    return result


def _inverse_project_pixel_coordinates(
    x: np.ndarray,
    y: np.ndarray,
    *,
    geometry: ScreenGeometry,
    projection: ViewProjection,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project selected display pixels without allocating a full grid."""
    radius = max(1.0, float(geometry.radius))
    nx = (np.asarray(x, dtype=np.float32) - float(geometry.center[0])) / radius
    ny = (np.asarray(y, dtype=np.float32) - float(geometry.center[1])) / radius
    rr2 = nx * nx + ny * ny
    edge_fov = max(1.0e-6, float(projection.edge_fov_deg))
    max_r = max(0.0, float(projection.content_fov_deg) / edge_fov)
    inside = rr2 <= max_r * max_r

    r = np.sqrt(rr2).astype(np.float32)
    theta = np.radians(r * edge_fov)
    psi = np.arctan2(nx, -ny)
    alt_c, az_c = projection.view_center
    eps = 1.0e-3
    phi1 = np.float32(math.radians(np.clip(alt_c, -90.0 + eps, 90.0 - eps)))
    lam1 = np.float32(math.radians(az_c))
    sin_phi1 = np.sin(phi1)
    cos_phi1 = np.cos(phi1)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    sin_phi2 = np.clip(
        sin_phi1 * cos_theta + cos_phi1 * sin_theta * np.cos(psi),
        -1.0,
        1.0,
    )
    phi2 = np.arcsin(sin_phi2)
    projected_y = np.sin(psi) * sin_theta * cos_phi1
    projected_x = cos_theta - sin_phi1 * sin_phi2
    lam2 = lam1 + np.arctan2(projected_y, projected_x)
    alt = np.degrees(phi2).astype(np.float32)
    az = ((np.degrees(lam2) + 360.0) % 360.0).astype(np.float32)
    return alt, az, inside

def _neu_unit_to_altaz(vec: np.ndarray) -> tuple[float, float]:
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
    projection: ViewProjection,
    terrain_profile_altaz: list[tuple[float, float]] | None = None,
    ground_reset_rgba: tuple[int, int, int, int] | None = None,
) -> QImage:
    """Reset the disc below the geometric or terrain horizon to a neutral background."""
    if ground_reset_rgba is None:
        return base_img
    out = qimage_to_np_rgba(
        base_img if base_img.format() == QImage.Format.Format_RGBA8888 else base_img.convertToFormat(QImage.Format.Format_RGBA8888)
    )
    alt, az, inside = _inverse_project_disc(
        out.shape[1],
        out.shape[0],
        geometry,
        projection,
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
    projection: ViewProjection,
    observer_lat_deg: float | None = None,
    never_rises_opacity: float = 0.2,
) -> QImage:
    """Draw a thin never-rises outline using the historic accent color."""
    if observer_lat_deg is None:
        return base_img
    alpha_u8 = int(np.clip(round(255.0 * float(never_rises_opacity)), 0, 255))
    if alpha_u8 <= 0:
        return base_img
    out = qimage_to_np_rgba(
        base_img if base_img.format() == QImage.Format.Format_RGBA8888 else base_img.convertToFormat(QImage.Format.Format_RGBA8888)
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
            projection.view_center,
            edge_fov_deg=projection.edge_fov_deg,
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
                float(projection.content_fov_deg)
                / max(1.0e-6, float(projection.edge_fov_deg)),
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
    earth_guide_opacity: float = EARTH_GUIDE_DEFAULT_OPACITY,
    visibility_boost: float = 1.0,
    fast_mode: bool = False,
    theme: ThemeStyle | None = None,
    sun_altaz: tuple[float, float] | None = None,
) -> QImage:
    if viewer_data is None or float(earth_guide_opacity) <= 0.0:
        return base_img
    out = (
        base_img
        if base_img.format() == QImage.Format.Format_ARGB32_Premultiplied
        else base_img.convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)
    )
    painter = QPainter(out)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    try:
        draw_earth_guide(
            painter,
            geometry=geometry,
            viewer_data=viewer_data,
            sun_altaz=sun_altaz,
            terrain_profile_altaz=terrain_profile_altaz,
            earth_guide_opacity=float(earth_guide_opacity),
            visibility_boost=float(visibility_boost),
            fast_mode=fast_mode,
            layer_style=None if theme is None else theme.overlays.earth_guide,
            target_size=(out.width(), out.height()),
        )
    finally:
        painter.end()
    return out


def _additive_rgb_overlay(base_img: QImage, overlay: np.ndarray) -> QImage:
    """Add an already-opacity-scaled RGB overlay to an image."""
    base = qimage_to_np_rgba(
        base_img
        if base_img.format() == QImage.Format.Format_RGBA8888
        else base_img.convertToFormat(QImage.Format.Format_RGBA8888)
    )
    if overlay.shape != base.shape:
        return base_img
    rgb = base[..., :3].astype(np.uint16) + overlay[..., :3].astype(np.uint16)
    base[..., :3] = np.minimum(rgb, 255).astype(np.uint8)
    return np_rgba_to_qimage(base)


class SkyCompositorCache:
    """Manage compositing and reuse the last composited image via a cache key."""

    def __init__(self) -> None:
        self._composited_img: QImage | None = None
        self._composite_key: tuple | None = None
        self._glow_mask_cache_stamp: tuple | None = None
        self._glow_mask_cache: GlowMask | None = None
        self._edge_glow_mask_cache_stamp: tuple | None = None
        self._edge_glow_mask_cache: GlowMask | None = None
        self.star_render_cache = StarRenderCache()

    def invalidate(self) -> None:
        self._composite_key = None
        self._composited_img = None
        self._glow_mask_cache_stamp = None
        self._glow_mask_cache = None
        self._edge_glow_mask_cache_stamp = None
        self._edge_glow_mask_cache = None
        self.star_render_cache.invalidate()
    def draw_cloud_voxel_overlay(
        self,
        painter: QPainter,
        *,
        geometry: ScreenGeometry,
        projection: ViewProjection,
        cloud_rgba: np.ndarray | None,
        cloud_alpha: float,
    ) -> None:
        """Draw a precomputed voxel raster over the sky and below all stars."""
        if cloud_rgba is None or cloud_alpha <= 0.0:
            return
        image = np_rgba_to_qimage(cloud_rgba)
        if image.isNull():
            return
        # Worker pixels span content_fov, while geometry.radius spans edge_fov.
        # Expand both the bounds and scale so the angular coordinates agree.
        content_radius = float(geometry.radius) * max(
            float(projection.edge_fov_deg), float(projection.content_fov_deg)
        ) / max(1.0e-6, float(projection.edge_fov_deg))
        painter.save()
        painter.setOpacity(float(np.clip(cloud_alpha, 0.0, 1.0)))
        painter.drawImage(
            QRectF(
                float(geometry.center[0]) - content_radius,
                float(geometry.center[1]) - content_radius,
                content_radius * 2 + 1,
                content_radius * 2 + 1,
            ),
            image,
        )
        painter.restore()

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
        sky_img: QImage | None,
        *,
        viewer_data: ViewerData,
        show_guidelines: bool = True,
        terrain_profile_altaz: list[tuple[float, float]] | None = None,
        terrain_profile_distances_m: list[float] | None = None,
        terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
        terrain_secondary_ridges_distances_m_layers: list[list[float]] | None = None,
        night_light_glow_profile: NightLightGlowProfile | None = None,
        terrain_horizon_opacity: float = 0.003,
        earth_guide_opacity: float = EARTH_GUIDE_DEFAULT_OPACITY,
        earth_guide_visibility_boost: float = 1.0,
        night_light_opacity: float = NIGHT_LIGHT_DEFAULT_OPACITY,
        ridge_glow_opacity: float = RIDGE_GLOW_DEFAULT_OPACITY,
        sun_altaz: tuple[float, float] | None = None,
        molecular_cloud_overlay: np.ndarray | None = None,
        never_rises_opacity: float = 0.2,
        ground_reset_rgba: tuple[int, int, int, int] | None = None,
        theme: ThemeStyle | None = None,
        fast_mode: bool = False,
        draw_sky_disc: bool = True,
        sky_disc_altaz_rings: str = "dimalt",
    ) -> None:
        """Composite the sky layers (with cache) and draw into painter."""
        view_center = viewer_data.view_center
        observer_lat_deg = viewer_data.location[0]
        observer_lon_deg = viewer_data.location[1]
        observer_height_m = viewer_data.observer_height_m
        edge_fov_deg = float(viewer_data.edge_fov_deg)
        content_fov_deg = float(viewer_data.content_fov_deg)
        viewport = painter.viewport()
        x = int(viewport.x())
        y = int(viewport.y())
        w = int(viewport.width())
        h = int(viewport.height())
        glow_viewer_data = viewer_data
        # Keep the solar position canonical: callers provide the complete
        # position, while consumers that only need altitude derive it here.
        night_light_sun_alt_deg = None if sun_altaz is None else float(sun_altaz[0])
        sky_projection = ViewProjection(
            view_center=viewer_data.view_center,
            edge_fov_deg=float(edge_fov_deg),
            content_fov_deg=float(content_fov_deg),
        )

        sky_ck = int(sky_img.cacheKey()) if sky_img else 0
        molecular_cloud_ck = id(molecular_cloud_overlay) if molecular_cloud_overlay is not None else 0
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
        comp_key = (
            "comp",
            sky_ck,
            molecular_cloud_ck,
            terrain_key,
            terrain_distance_key,
            terrain_secondary_key,
            terrain_secondary_distance_key,
            w,
            h,
            tuple(geometry.center),
            int(geometry.radius),
            float(view_center[0]),
            float(view_center[1]),
            float(content_fov_deg),
            None if observer_lat_deg is None else float(observer_lat_deg),
            None if observer_lon_deg is None else float(observer_lon_deg),
            float(observer_height_m),
            bool(draw_sky_disc),
            bool(show_guidelines),
            float(terrain_horizon_opacity),
            float(earth_guide_opacity),
            float(night_light_opacity),
            float(ridge_glow_opacity),
            None if night_light_sun_alt_deg is None else round(float(night_light_sun_alt_deg), 3),
            None
            if sun_altaz is None
            else tuple(round(float(value), 3) for value in sun_altaz),
            float(never_rises_opacity),
            bool(fast_mode),
            night_light_key,
            glow_key,
            None
            if ground_reset_rgba is None
            else tuple(int(np.clip(c, 0, 255)) for c in ground_reset_rgba),
            float(GLOW_MASK_SCALE),
            tuple(int(c) for c in GLOW_MASK_TINT_RGB),
            str(sky_disc_altaz_rings),
            edge_glow_key,
            None
            if theme is None
            else (
                tuple(int(c) for c in theme.window_background.base_rgb),
                tuple(int(c) for c in theme.window_background.delta_rgb),
                int(theme.window_background.outer_alpha),
                int(theme.window_background.edge_alpha),
                bool(theme.window_background.flat_background),
            ),
        )

        if self._composite_key != comp_key or self._composited_img is None:
            def _scaled(qimg: QImage | None) -> QImage | None:
                if qimg is None:
                    return None
                if qimg.width() == w and qimg.height() == h:
                    return qimg
                return _scale_qimage_preserving_aspect(qimg, w, h)

            def _black_disc_image() -> QImage:
                img = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
                img.fill(Qt.GlobalColor.transparent)
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

            sky_s = _scaled(sky_img) if draw_sky_disc else None
            if draw_sky_disc and sky_s is None:
                sky_s = _black_disc_image()
            if draw_sky_disc:
                assert sky_s is not None
            if draw_sky_disc and sky_disc_altaz_rings == "dimalt" and sky_s is not None:
                sky_s = apply_altitude_ring_highlights(
                    sky_s,
                    geometry,
                    sky_projection,
                    theme=theme,
                    altaz_rings_mode="dimalt",
                )
            if not draw_sky_disc:
                composited = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
                composited.fill(Qt.GlobalColor.transparent)
            else:
                assert sky_s is not None
                composited = sky_s
            if draw_sky_disc:
                # Clip the combined sky layers only after their low-resolution
                # surfaces have been enlarged, so the terrain edge is smooth.
                if molecular_cloud_overlay is not None:
                    composited = _additive_rgb_overlay(composited, molecular_cloud_overlay)
                composited = _clip_below_terrain_horizon(
                    composited,
                    geometry=geometry,
                    projection=replace(
                        sky_projection,
                        content_fov_deg=content_fov_deg + SKY_DISC_OVERSCAN_DEG,
                    ),
                    terrain_profile_altaz=terrain_profile_altaz,
                )
                composited = _clip_sky_image_to_disc(
                    composited,
                    geometry=geometry,
                    viewport_offset=(x, y),
                    edge_fov_deg=edge_fov_deg,
                    content_fov_deg=content_fov_deg + SKY_DISC_OVERSCAN_DEG,
                )
            earth_viewer_data = viewer_data
            composited = _overlay_earth_guide(
                composited,
                geometry=geometry,
                viewer_data=earth_viewer_data,
                terrain_profile_altaz=terrain_profile_altaz,
                earth_guide_opacity=earth_guide_opacity,
                visibility_boost=earth_guide_visibility_boost,
                fast_mode=fast_mode,
                theme=theme,
                sun_altaz=sun_altaz,
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
                edge_glow_image = _glow_mask_to_qimage(edge_glow_mask, EDGE_GLOW_MASK_TINT_RGB)
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
                    projection=sky_projection,
                    observer_lat_deg=observer_lat_deg,
                    never_rises_opacity=never_rises_opacity,
                )

            self._composited_img = composited
            self._composite_key = comp_key

        painter.save()
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
        painter.drawImage(QRect(x, y, w, h), self._composited_img)
        painter.restore()
