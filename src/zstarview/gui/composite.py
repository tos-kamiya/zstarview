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
from dataclasses import dataclass, replace
from typing import cast

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
    NIGHT_LIGHT_DEFAULT_OPACITY,
    PALETTE_NEVER_RISES_GUIDE_RGB,
    RIDGE_GLOW_DEFAULT_OPACITY,
    CloudLayerStyle,
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
from ..render.sky_disc import SKY_DISC_OVERSCAN_DEG
from ..types import ScreenGeometry, ViewerData, ViewProjection
from .cloud_render import (
    CLOUD_DAY_RGB,
    CLOUD_NIGHT_BOOST,
    _cloud_tint_rgb_for_theme,
    _mask_cloud_alpha_by_missing_rgba,
    _render_alpha_scaled_cloud_stripes_rgba_from_altaz_grid,
    _render_halftone_cloud_rgba_from_altaz_grid,
    _render_variable_width_cloud_stripes_rgba_from_altaz_grid,
    _scale_qimage_preserving_aspect,
    _smooth_cloud_amount_grid,
    _sunset_cloud_tint_rgb,
    compose_cloud_over_sky,
    mask_cloud_alpha_by_missing,  # noqa: F401 - public compatibility export
)


@dataclass(frozen=True)
class GlowMask:
    """Low-resolution alpha-only glow buffer in normalized screen space."""

    alpha: np.ndarray
    scale: float

NEVER_RISES_GUIDE_WIDTH_SCALE = 4.5
NEVER_RISES_GUIDE_ALPHA_SCALE = 0.5
ALT_RING_DIMALT_SAMPLE_AZ_STEP_DEG = 30.0
HALFTONE_MIN_GRID_DELTA_PX = 22.0
HALFTONE_LEVEL_DIAMETER_BASE_SCALE = 0.9


def _cloud_shell_grids(grid: CloudAltAzGrid) -> tuple[CloudAltAzGrid, ...]:
    """Return shell-specific grids, or the legacy single grid when absent."""
    if not grid.shell_amounts:
        return (grid,)
    return tuple(
        replace(grid, amount=amount, shell_amounts=None)
        for amount in grid.shell_amounts
    )


def _render_cloud_grid_rgba(
    cache: SkyCompositorCache,
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    *,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry,
    projection: ViewProjection,
    target_stripes: int,
    width_factor: float,
    density_reference_size: tuple[int, int] | None,
) -> np.ndarray:
    """Render one cloud grid using the selected cloud stripe mode."""
    if cache._cloud_stripe_mode == "alpha":
        return _render_alpha_scaled_cloud_stripes_rgba_from_altaz_grid(
            grid,
            width,
            height,
            hatch_cfg,
            geometry=geometry,
            projection=projection,
            target_stripes=target_stripes,
            width_factor=width_factor,
            density_reference_size=density_reference_size,
        )
    if cache._cloud_stripe_mode == "halftone":
        return _render_halftone_cloud_rgba_from_altaz_grid(
            grid,
            width,
            height,
            hatch_cfg,
            geometry=geometry,
            projection=projection,
            target_stripes=target_stripes,
            width_factor=width_factor,
            density_reference_size=density_reference_size,
        )
    return _render_variable_width_cloud_stripes_rgba_from_altaz_grid(
        grid,
        width,
        height,
        hatch_cfg,
        geometry=geometry,
        projection=projection,
        target_stripes=target_stripes,
        width_factor=width_factor,
        density_reference_size=density_reference_size,
    )


def _combine_cloud_shell_rgba(
    layers: list[tuple[np.ndarray, tuple[int, int, int]]],
) -> np.ndarray | None:
    """Composite shell images from far/high to near/low."""
    if not layers:
        return None
    shape = layers[0][0].shape
    out_premultiplied = np.zeros(shape[:2] + (3,), dtype=np.float32)
    out_alpha = np.zeros(shape[:2], dtype=np.float32)
    for image, tint_rgb in reversed(layers):
        src_alpha = image[..., 3].astype(np.float32) / 255.0
        if not np.any(src_alpha > 0.0):
            continue
        src_rgb = np.asarray(tint_rgb, dtype=np.float32)[None, None, :]
        out_premultiplied = (
            src_rgb * src_alpha[..., None]
            + out_premultiplied * (1.0 - src_alpha[..., None])
        )
        out_alpha = src_alpha + out_alpha * (1.0 - src_alpha)
    out = np.zeros(shape, dtype=np.uint8)
    out_rgb = np.zeros_like(out_premultiplied)
    np.divide(
        out_premultiplied,
        np.maximum(out_alpha[..., None], 1.0e-6),
        out=out_rgb,
        where=out_alpha[..., None] > 0.0,
    )
    out[..., :3] = np.clip(np.round(out_rgb), 0, 255).astype(np.uint8)
    out[..., 3] = np.clip(np.round(out_alpha * 255.0), 0, 255).astype(np.uint8)
    return out


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
    clipped.fill(Qt.transparent)
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
GLOW_MASK_NOISE_VARIATION = 0.16
GLOW_MASK_NIGHT_LIGHT_HEIGHT_DEG = 30.0
GLOW_MASK_NIGHT_LIGHT_DECAY_RATE = 2.4
GLOW_MASK_NIGHT_LIGHT_ALTITUDE_CROP_ALPHA_THRESHOLD = 1.0e-4
GLOW_MASK_NIGHT_LIGHT_ALTITUDE_CROP_PAD_ROWS = 1





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
    view_center: tuple[float, float],
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
    tint_rgba: tuple[int, int, int, int] = CLOUD_MISSING_TINT_RGBA,
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




def _inverse_project_disc(
    width: int,
    height: int,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
    content_fov_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project square composited pixels up to the requested content FOV."""
    return _shared_inverse_project_disc(
        width,
        height,
        geometry,
        view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
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
    view_center: tuple[float, float],
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
    view_center: tuple[float, float],
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
    theme: ThemeStyle | None = None,
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
            layer_style=None if theme is None else theme.overlays.earth_guide,
        )
    finally:
        painter.end()
    return out


def _additive_rgb_overlay(base_img: QImage, overlay: np.ndarray) -> QImage:
    """Add an already-opacity-scaled RGB overlay to an image."""
    base = qimage_to_np_rgba(
        base_img
        if base_img.format() == QImage.Format_RGBA8888
        else base_img.convertToFormat(QImage.Format_RGBA8888)
    )
    if overlay.shape != base.shape:
        return base_img
    rgb = base[..., :3].astype(np.uint16) + overlay[..., :3].astype(np.uint16)
    base[..., :3] = np.minimum(rgb, 255).astype(np.uint8)
    return np_rgba_to_qimage(base)


class SkyCompositorCache:
    """Manage compositing and reuse the last composited image via a cache key."""

    def __init__(
        self,
        *,
        hatch_cfg: HatchConfig = CLOUD_HATCH_DEFAULT,
        gray_mix: float = 1.0,
        cloud_target_stripes: int = 30,
        cloud_stripe_width_factor: float = 1.7,
        cloud_stripe_mode: str = "halftone",
        missing_tint_rgba: tuple[int, int, int, int] = CLOUD_MISSING_TINT_RGBA,
        ground_tint_opacity: float = 0.025,
    ) -> None:
        self._hatch_cfg = hatch_cfg
        self._gray_mix = gray_mix
        self._cloud_target_stripes = max(1, int(cloud_target_stripes))
        self._cloud_stripe_width_factor = max(0.01, float(cloud_stripe_width_factor))
        mode = str(cloud_stripe_mode).strip().lower()
        self._cloud_stripe_mode = mode if mode in ("alpha", "halftone") else "width"
        self._missing_tint_rgba: tuple[int, int, int, int] = (
            int(np.clip(missing_tint_rgba[0], 0, 255)),
            int(np.clip(missing_tint_rgba[1], 0, 255)),
            int(np.clip(missing_tint_rgba[2], 0, 255)),
            int(np.clip(missing_tint_rgba[3], 0, 255)),
        )
        self._ground_tint_opacity = float(np.clip(ground_tint_opacity, 0.0, 1.0))
        self._composited_img: QImage | None = None
        self._composite_key: tuple | None = None
        self._glow_mask_cache_stamp: tuple | None = None
        self._glow_mask_cache: GlowMask | None = None
        self._edge_glow_mask_cache_stamp: tuple | None = None
        self._edge_glow_mask_cache: GlowMask | None = None
        self._atlas_cloud_cache_key: tuple[object, ...] | None = None
        self._atlas_cloud_cache_images: tuple[QImage, QImage | None] | None = None

    def invalidate(self) -> None:
        self._composite_key = None
        self._composited_img = None
        self._glow_mask_cache_stamp = None
        self._glow_mask_cache = None
        self._edge_glow_mask_cache_stamp = None
        self._edge_glow_mask_cache = None
        self._atlas_cloud_cache_key = None
        self._atlas_cloud_cache_images = None
    @property
    def cloud_target_stripes(self) -> int:
        return self._cloud_target_stripes

    @property
    def cloud_stripe_width_factor(self) -> float:
        return self._cloud_stripe_width_factor

    def render_atlas_cloud_layer(
        self,
        *,
        width: int,
        height: int,
        geometry: ScreenGeometry,
        projection: ViewProjection,
        grid: CloudAltAzGrid,
        missing_mask: np.ndarray | None,
        target_stripes: int,
        width_factor: float,
        opacity: float,
        style: CloudLayerStyle,
    ) -> tuple[QImage, QImage | None]:
        """Return cached Atlas cloud and missing-data images."""
        w = max(1, int(width))
        h = max(1, int(height))
        cloud_opacity = float(np.clip(opacity * style.alpha_scale, 0.0, 1.0))
        tint = tuple(int(np.clip(value, 0, 255)) for value in style.missing_rgba)
        key = (
            id(grid),
            id(missing_mask) if missing_mask is not None else 0,
            w,
            h,
            tuple(round(float(value), 3) for value in geometry.center),
            round(float(geometry.radius), 3),
            tuple(round(float(value), 3) for value in projection.view_center),
            round(float(projection.edge_fov_deg), 3),
            round(float(projection.content_fov_deg), 3),
            max(1, int(target_stripes)),
            round(float(width_factor * style.width_scale), 5),
            round(cloud_opacity, 5),
            tuple(int(value) for value in style.rgb),
            tint,
        )
        if self._atlas_cloud_cache_key == key and self._atlas_cloud_cache_images is not None:
            return self._atlas_cloud_cache_images

        hatch_cfg = HatchConfig(
            self._hatch_cfg.tile_w_px,
            self._hatch_cfg.tile_h_px,
            self._hatch_cfg.line_px,
            int(round(float(self._hatch_cfg.strength) * cloud_opacity)),
        )
        cloud_rgba = _render_variable_width_cloud_stripes_rgba_from_altaz_grid(
            grid,
            w,
            h,
            hatch_cfg,
            geometry=geometry,
            projection=projection,
            target_stripes=target_stripes,
            width_factor=width_factor * style.width_scale,
            stripe_rgb=style.rgb,
        )
        cloud_rgba = _mask_cloud_alpha_by_missing_rgba(cloud_rgba, missing_mask) if missing_mask is not None else cloud_rgba
        cloud_image = np_rgba_to_qimage(cloud_rgba)

        missing_image: QImage | None = None
        if missing_mask is not None and tint[3] > 0:
            missing = np.asarray(missing_mask)
            if missing.shape != (h, w):
                y_idx = np.rint(np.linspace(0, missing.shape[0] - 1, h)).astype(np.int32)
                x_idx = np.rint(np.linspace(0, missing.shape[1] - 1, w)).astype(np.int32)
                missing = missing[y_idx][:, x_idx]
            missing_rgba = np.zeros((h, w, 4), dtype=np.uint8)
            missing_rgba[..., :3] = np.asarray(tint[:3], dtype=np.uint8)
            missing_rgba[..., 3] = np.clip(
                np.round(missing.astype(np.float32) * (float(tint[3]) / 255.0)),
                0,
                255,
            ).astype(np.uint8)
            missing_image = np_rgba_to_qimage(missing_rgba)

        result = (cloud_image, missing_image)
        self._atlas_cloud_cache_key = key
        self._atlas_cloud_cache_images = result
        return result

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
        cloud_alpha: float,
        density_reference_size: tuple[int, int] | None = None,
        view_center: tuple[float, float] = (0.0, 0.0),
        observer_lat_deg: float | None = None,
        observer_lon_deg: float | None = None,
        observer_height_m: float = 0.0,
        cloud_altaz_grid: CloudAltAzGrid | None = None,
        missing_mask: np.ndarray | None = None,
        show_guidelines: bool = True,
        terrain_profile_altaz: list[tuple[float, float]] | None = None,
        terrain_profile_distances_m: list[float] | None = None,
        terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
        terrain_secondary_ridges_distances_m_layers: list[list[float]] | None = None,
        night_light_glow_profile: NightLightGlowProfile | None = None,
        terrain_horizon_opacity: float = 0.003,
        earth_guide_opacity: float = 0.028,
        earth_guide_visibility_boost: float = 1.0,
        night_light_opacity: float = NIGHT_LIGHT_DEFAULT_OPACITY,
        ridge_glow_opacity: float = RIDGE_GLOW_DEFAULT_OPACITY,
        night_light_sun_alt_deg: float | None = None,
        sun_altaz: tuple[float, float] | None = None,
        aerosol_optical_depth: float | None = None,
        molecular_cloud_overlay: np.ndarray | None = None,
        never_rises_opacity: float = 0.2,
        ground_reset_rgba: tuple[int, int, int, int] | None = None,
        theme: ThemeStyle | None = None,
        edge_fov_deg: float = 90.0,
        content_fov_deg: float,
        fast_mode: bool = False,
        draw_sky_disc: bool = True,
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
        molecular_cloud_ck = id(molecular_cloud_overlay) if molecular_cloud_overlay is not None else 0
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
                id(cloud_altaz_grid.shell_amounts),
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
        effective_cloud_alpha = float(np.clip(cloud_alpha, 0.0, 1.0))
        if night_light_sun_alt_deg is not None:
            effective_cloud_alpha *= 1.0 + (
                CLOUD_NIGHT_BOOST
                * float(night_light_strength_factor(night_light_sun_alt_deg))
            )
            effective_cloud_alpha = float(np.clip(effective_cloud_alpha, 0.0, 1.0))
        hatch_key = (
            self._hatch_cfg.tile_w_px,
            self._hatch_cfg.tile_h_px,
            self._hatch_cfg.line_px,
            self._hatch_cfg.strength,
        )
        comp_key = (
            "comp",
            sky_ck,
            molecular_cloud_ck,
            missing_ck,
            terrain_key,
            terrain_distance_key,
            terrain_secondary_key,
            terrain_secondary_distance_key,
            w,
            h,
            tuple(geometry.center),
            int(geometry.radius),
            effective_cloud_alpha,
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
            None
            if aerosol_optical_depth is None
            else round(float(aerosol_optical_depth), 5),
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

            sky_s = _scaled(sky_img) if draw_sky_disc else None
            if draw_sky_disc and sky_s is None:
                sky_s = _black_disc_image()
            missing_s = missing_mask
            cloud_s: np.ndarray | None = None
            cloud_tint_rgb = _cloud_tint_rgb_for_theme(
                theme, night_light_sun_alt_deg
            )

            if draw_sky_disc and effective_cloud_alpha > 0.0:
                if cloud_altaz_grid is not None:
                    if cloud_altaz_grid.shell_amounts:
                        layers: list[tuple[np.ndarray, tuple[int, int, int]]] = []
                        for shell_index, shell_grid in enumerate(
                            _cloud_shell_grids(cloud_altaz_grid)
                        ):
                            shell_image = _render_cloud_grid_rgba(
                                self,
                                shell_grid,
                                w,
                                h,
                                hatch_cfg=self._hatch_cfg,
                                geometry=geometry,
                                projection=cloud_projection,
                                target_stripes=self._cloud_target_stripes,
                                width_factor=self._cloud_stripe_width_factor,
                                density_reference_size=density_reference_size,
                            )
                            if missing_s is not None:
                                shell_image = _mask_cloud_alpha_by_missing_rgba(
                                    shell_image, missing_s
                                )
                            layers.append(
                                (
                                    shell_image,
                                    _sunset_cloud_tint_rgb(
                                        sun_altaz,
                                        base_rgb=cloud_tint_rgb,
                                        shell_index=shell_index,
                                        aerosol_optical_depth=aerosol_optical_depth,
                                    ),
                                )
                            )
                        cloud_s = _combine_cloud_shell_rgba(layers)
                        cloud_tint_rgb = CLOUD_DAY_RGB
                    else:
                        cloud_s = _render_cloud_grid_rgba(
                            self,
                            cloud_altaz_grid,
                            w,
                            h,
                            geometry=geometry,
                            projection=cloud_projection,
                            hatch_cfg=self._hatch_cfg,
                            target_stripes=self._cloud_target_stripes,
                            width_factor=self._cloud_stripe_width_factor,
                            density_reference_size=density_reference_size,
                        )
                if (
                    missing_s is not None
                    and cloud_s is not None
                    and not cloud_altaz_grid.shell_amounts
                ):
                    cloud_s = _mask_cloud_alpha_by_missing_rgba(cloud_s, missing_s)
            if draw_sky_disc and sky_disc_altaz_rings == "dimalt" and sky_s is not None:
                sky_s = apply_altitude_ring_highlights(
                    sky_s,
                    geometry,
                    view_center,
                    theme=theme,
                    edge_fov_deg=edge_fov_deg,
                    altaz_rings_mode="dimalt",
                )
                if cloud_s is None or effective_cloud_alpha <= 0.0:
                    composited = sky_s
                else:
                    composited = compose_cloud_over_sky(
                        sky_img=sky_s,
                        cloud_img_rgba=cloud_s,
                        dest_rect=QRect(0, 0, w, h),
                        geometry=geometry,
                        cloud_opacity=effective_cloud_alpha,
                        gray_mix=self._gray_mix,
                        edge_fov_deg=edge_fov_deg,
                        content_fov_deg=content_fov_deg,
                        sun_alt_deg=night_light_sun_alt_deg,
                        cloud_tint_rgb=cloud_tint_rgb,
                        transparent_sky_rgb=None
                        if theme is None
                        else tuple(int(c) for c in theme.window_background.inner_rgba[:3]),
                    )
            if not draw_sky_disc:
                composited = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
                composited.fill(Qt.transparent)
            elif cloud_s is None or effective_cloud_alpha <= 0.0:
                composited = sky_s
            else:
                composited = compose_cloud_over_sky(
                    sky_img=sky_s,
                    cloud_img_rgba=cloud_s,
                    dest_rect=QRect(0, 0, w, h),
                    geometry=geometry,
                    cloud_opacity=effective_cloud_alpha,
                    gray_mix=self._gray_mix,
                    edge_fov_deg=edge_fov_deg,
                    content_fov_deg=content_fov_deg,
                    sun_alt_deg=night_light_sun_alt_deg,
                    cloud_tint_rgb=cloud_tint_rgb,
                    transparent_sky_rgb=None
                    if theme is None
                    else tuple(int(c) for c in theme.window_background.inner_rgba[:3]),
                )
            if draw_sky_disc:
                composited = _apply_ground_reset(
                    composited,
                    geometry=geometry,
                    view_center=view_center,
                    terrain_profile_altaz=terrain_profile_altaz,
                    ground_reset_rgba=ground_reset_rgba,
                    edge_fov_deg=edge_fov_deg,
                    content_fov_deg=content_fov_deg + SKY_DISC_OVERSCAN_DEG,
                )
                if molecular_cloud_overlay is not None:
                    composited = _additive_rgb_overlay(composited, molecular_cloud_overlay)
                composited = _clip_sky_image_to_disc(
                    composited,
                    geometry=geometry,
                    viewport_offset=(x, y),
                    edge_fov_deg=edge_fov_deg,
                    content_fov_deg=content_fov_deg + SKY_DISC_OVERSCAN_DEG,
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
                theme=theme,
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
