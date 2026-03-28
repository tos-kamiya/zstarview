# -*- coding: utf-8 -*-
"""
Sky and cloud compositing utilities and cache.

This module provides:
- Pure helpers to generate hatch tiles and apply them to cloud alpha.
- A function to composite sky and cloud layers without relying on any global state.
- A small cache class that handles scaling, hatching, compositing, and reuse.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional, Tuple, cast

import numpy as np
from PySide6.QtCore import QRect, Qt
from PySide6.QtGui import QImage, QPainter

from ..paths import CLOUD_HATCH_DEFAULT, CLOUD_MISSING_TINT_RGBA, HatchConfig
from ..render.draw_sky_disc import GROUND_TINT_RGB, NEVER_RISES_TINT_RGB, NEVER_RISES_TINT_STRENGTH
from ..types import ScreenGeometry
from ..render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


@dataclass(frozen=True)
class StripeDensityField:
    """Compact cloud-density field in normalized 45-degree (u, v) space."""

    density: np.ndarray  # float32 in [0,1], shape=(bins_u, bins_v)
    u_min: float
    u_max: float
    v_min: float
    v_max: float
    source_cache_key: int


def cloud_with_hatched_alpha(cloud_img: QImage, hatch_cfg: HatchConfig) -> QImage:
    """Apply continuous 45-degree hatch directly to alpha (no tile seams)."""
    out = cloud_img if cloud_img.format() == QImage.Format_RGBA8888 else cloud_img.convertToFormat(QImage.Format_RGBA8888)
    cloud = qimage_to_np_rgba(out)

    h, w = cloud.shape[:2]
    xs = np.arange(w, dtype=np.int32)[None, :]
    ys = np.arange(h, dtype=np.int32)[:, None]

    period = max(2, int(round(math.hypot(hatch_cfg.tile_w_px, hatch_cfg.tile_h_px) * 0.5)))
    band = max(1, int(round(hatch_cfg.line_px)))

    # Exact 45-degree lines: x - y = const, repeated by period.
    u = xs - ys
    u_mod = np.mod(u, period)
    dist = np.minimum(u_mod, period - u_mod)
    line_mask = dist <= (band / 2.0)
    keep_mask = ~line_mask

    # Convert the "cloud band" into a cloud-amount indicator:
    # for each diagonal band (same stripe index), flatten alpha to the
    # band-average value so cloud contours are not visible inside the band.
    # Only include on-disc pixels (RGB>0 for this cloud layer) in the average.
    stripe_id = np.floor_divide(u, period)
    inside_disc = cloud[..., 0] > 0
    avg_mask = keep_mask & inside_disc
    if np.any(avg_mask):
        # Work in normalized alpha [0.0, 1.0].
        src_alpha01 = cloud[..., 3].astype(np.float32) / 255.0
        stripe_idx = stripe_id[avg_mask].astype(np.int64, copy=False)
        alpha_vals = src_alpha01[avg_mask].astype(np.float64, copy=False)
        sid_min = int(stripe_idx.min())
        sid_norm = stripe_idx - sid_min
        sums = np.bincount(sid_norm, weights=alpha_vals)
        counts = np.bincount(sid_norm)
        means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)
        src_alpha01[avg_mask] = means[sid_norm].astype(np.float32, copy=False)
        cloud[..., 3] = np.clip(np.round(src_alpha01 * 255.0), 0, 255).astype(np.uint8)

    erase = float(np.clip(hatch_cfg.strength, 0, 255)) / 255.0
    keep = 1.0 - erase

    # Match DestinationOut semantics: attenuate both alpha and premultiplied color.
    rgb = cloud[..., :3].astype(np.float32)
    alpha = cloud[..., 3].astype(np.float32)
    rgb[line_mask] = rgb[line_mask] * keep
    alpha[line_mask] = alpha[line_mask] * keep

    cloud[..., :3] = np.clip(np.round(rgb), 0, 255).astype(np.uint8)
    cloud[..., 3] = np.clip(np.round(alpha), 0, 255).astype(np.uint8)
    return np_rgba_to_qimage(cloud)


def build_stripe_density_field(cloud_img: QImage, *, bins: int = 192) -> StripeDensityField:
    """Build a compact density field from a cloud image in normalized (u, v) space."""
    cloud = qimage_to_np_rgba(cloud_img if cloud_img.format() == QImage.Format_RGBA8888 else cloud_img.convertToFormat(QImage.Format_RGBA8888))
    h, w = cloud.shape[:2]

    alpha01 = cloud[..., 3].astype(np.float32) / 255.0
    inside = alpha01 > 0.0
    if not np.any(inside):
        density = np.zeros((bins, bins), dtype=np.float32)
        return StripeDensityField(density=density, u_min=-2.0, u_max=2.0, v_min=-2.0, v_max=2.0, source_cache_key=int(cloud_img.cacheKey()))

    cy, cx = (h - 1) * 0.5, (w - 1) * 0.5
    r = max(1.0, min(cx, cy))
    y, x = np.ogrid[:h, :w]
    xn = (x - cx) / r
    yn = (y - cy) / r
    u = xn - yn
    v = xn + yn

    u_min, u_max = -2.0, 2.0
    v_min, v_max = -2.0, 2.0
    bins_u = max(32, int(bins))
    bins_v = bins_u

    u_idx = np.clip(((u - u_min) * (bins_u / (u_max - u_min))).astype(np.int32), 0, bins_u - 1)
    v_idx = np.clip(((v - v_min) * (bins_v / (v_max - v_min))).astype(np.int32), 0, bins_v - 1)

    ids = (u_idx[inside] * bins_v + v_idx[inside]).astype(np.int64, copy=False)
    vals = alpha01[inside].astype(np.float64, copy=False)
    size = bins_u * bins_v
    sums = np.bincount(ids, weights=vals, minlength=size)
    counts = np.bincount(ids, minlength=size)
    means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0).astype(np.float32, copy=False)
    density = cast(np.ndarray, np.asarray(means, dtype=np.float32).reshape((bins_u, bins_v)))

    return StripeDensityField(
        density=density,
        u_min=u_min,
        u_max=u_max,
        v_min=v_min,
        v_max=v_max,
        source_cache_key=int(cloud_img.cacheKey()),
    )


def render_hatched_cloud_from_density(
    density: StripeDensityField,
    width: int,
    height: int,
    hatch_cfg: HatchConfig,
    geometry: ScreenGeometry | None = None,
    *,
    target_stripes: int = 50,
    width_factor: float = 0.2,
    content_fov_deg: float = 90.0,
) -> QImage:
    """Render a sharp hatch cloud image from density field at the target window size."""
    w = max(1, int(width))
    h = max(1, int(height))

    diameter_px = float(min(w, h))
    stripes = max(1, int(target_stripes))
    wf = max(0.01, float(width_factor))
    period = int(np.clip(round(diameter_px / stripes), 14, 64))
    band = max(1, int(round(period * wf)))

    xs = np.arange(w, dtype=np.int32)[None, :]
    ys = np.arange(h, dtype=np.int32)[:, None]
    u_pix = xs - ys
    u_mod = np.mod(u_pix, period)
    dist = np.minimum(u_mod, period - u_mod)
    line_mask = dist <= (band / 2.0)

    if geometry is None:
        cx = (w - 1) * 0.5
        cy = (h - 1) * 0.5
        rr = max(1.0, min(cx, cy))
    else:
        cx = float(geometry.center[0])
        cy = float(geometry.center[1])
        rr = max(1.0, float(geometry.radius))
    max_r = max(0.0, float(content_fov_deg) / 90.0)
    y, x = np.ogrid[:h, :w]
    inside_disc = ((x - cx) ** 2 + (y - cy) ** 2) <= ((rr * max_r) + 0.25) ** 2
    draw_mask = line_mask & inside_disc

    out = np.zeros((h, w, 4), dtype=np.uint8)
    if not np.any(draw_mask):
        return np_rgba_to_qimage(out)

    xn = (x - cx) / rr
    yn = (y - cy) / rr
    u = xn - yn
    v = xn + yn

    bins_u, bins_v = density.density.shape
    u_idx = np.clip(((u - density.u_min) * (bins_u / (density.u_max - density.u_min))).astype(np.int32), 0, bins_u - 1)
    v_idx = np.clip(((v - density.v_min) * (bins_v / (density.v_max - density.v_min))).astype(np.int32), 0, bins_v - 1)

    sampled = density.density[u_idx, v_idx]
    strength = float(np.clip(hatch_cfg.strength, 0, 255)) / 255.0
    alpha = np.zeros((h, w), dtype=np.float32)
    alpha[draw_mask] = sampled[draw_mask] * 255.0 * strength

    out[..., :3][draw_mask] = 255
    out[..., 3] = np.clip(np.round(alpha), 0, 255).astype(np.uint8)
    return np_rgba_to_qimage(out)


def compose_cloud_over_sky(
    sky_img: QImage,
    cloud_img_rgba: QImage,
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
    if cloud_img_rgba.width() != w or cloud_img_rgba.height() != h:
        cloud_img_rgba = cloud_img_rgba.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

    sky_np = qimage_to_np_rgba(sky_img)
    cloud_np = qimage_to_np_rgba(cloud_img_rgba)

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
    missing_mask_img: QImage,
    *,
    tint_rgba: Tuple[int, int, int, int] = CLOUD_MISSING_TINT_RGBA,
) -> QImage:
    """Overlay missing-data regions with a faint yellow solid tint."""
    w, h = base_img.width(), base_img.height()
    if missing_mask_img.width() != w or missing_mask_img.height() != h:
        missing_mask_img = missing_mask_img.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

    out = qimage_to_np_rgba(base_img if base_img.format() == QImage.Format_RGBA8888 else base_img.convertToFormat(QImage.Format_RGBA8888))
    mask_rgba = qimage_to_np_rgba(
        missing_mask_img
        if missing_mask_img.format() == QImage.Format_RGBA8888
        else missing_mask_img.convertToFormat(QImage.Format_RGBA8888)
    )

    # Treat missing-mask alpha as coverage indicator in [0,1].
    m = mask_rgba[..., 3].astype(np.float32) / 255.0
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


def mask_cloud_alpha_by_missing(
    cloud_img: QImage,
    missing_mask_img: QImage,
) -> QImage:
    """Set cloud alpha to zero where missing-mask alpha is present."""
    w, h = cloud_img.width(), cloud_img.height()
    if missing_mask_img.width() != w or missing_mask_img.height() != h:
        missing_mask_img = missing_mask_img.scaled(w, h, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)
    cloud = qimage_to_np_rgba(cloud_img if cloud_img.format() == QImage.Format_RGBA8888 else cloud_img.convertToFormat(QImage.Format_RGBA8888))
    missing = qimage_to_np_rgba(
        missing_mask_img
        if missing_mask_img.format() == QImage.Format_RGBA8888
        else missing_mask_img.convertToFormat(QImage.Format_RGBA8888)
    )
    cut = missing[..., 3] > 0
    if np.any(cut):
        cloud[..., 3][cut] = 0
        cloud[..., :3][cut] = 0
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


def apply_ground_tint(
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
        cloud_opacity_scale: float = 0.7,
        gray_mix: float = 1.0,
        cloud_target_stripes: int = 50,
        cloud_stripe_width_factor: float = 0.2,
        missing_tint_rgba: Tuple[int, int, int, int] = CLOUD_MISSING_TINT_RGBA,
        ground_tint_opacity: float = 1.0,
    ) -> None:
        self._hatch_cfg = hatch_cfg
        self._cloud_opacity_scale = cloud_opacity_scale
        self._gray_mix = gray_mix
        self._cloud_target_stripes = max(1, int(cloud_target_stripes))
        self._cloud_stripe_width_factor = max(0.01, float(cloud_stripe_width_factor))
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
        cloud_img: Optional[QImage],
        *,
        cloud_alpha: float,
        view_center: Tuple[float, float] = (0.0, 0.0),
        observer_lat_deg: float | None = None,
        stripe_density: Optional[StripeDensityField] = None,
        missing_mask: Optional[QImage] = None,
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
        cloud_ck = int(cloud_img.cacheKey()) if cloud_img else 0
        density_ck = int(stripe_density.source_cache_key) if stripe_density is not None else 0
        missing_ck = int(missing_mask.cacheKey()) if missing_mask is not None else 0
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
            density_ck,
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
            self._cloud_opacity_scale,
            self._gray_mix,
            self._cloud_target_stripes,
            self._cloud_stripe_width_factor,
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
            cloud_s = _scaled(cloud_img)
            missing_s = _scaled(missing_mask)

            if cloud_s is not None and cloud_alpha > 0.0:
                if stripe_density is not None:
                    cloud_s = render_hatched_cloud_from_density(
                        stripe_density,
                        w,
                        h,
                        self._hatch_cfg,
                        geometry=geometry,
                        target_stripes=self._cloud_target_stripes,
                        width_factor=self._cloud_stripe_width_factor,
                        content_fov_deg=content_fov_deg,
                    )
                else:
                    cloud_s = cloud_with_hatched_alpha(cloud_s, self._hatch_cfg)
                if missing_s is not None:
                    cloud_s = mask_cloud_alpha_by_missing(cloud_s, missing_s)

            if cloud_s is None or cloud_alpha <= 0.0:
                composited = sky_s
            else:
                composited = compose_cloud_over_sky(
                    sky_img=sky_s,
                    cloud_img_rgba=cloud_s,
                    dest_rect=QRect(0, 0, w, h),
                    geometry=geometry,
                    cloud_opacity=cloud_alpha * self._cloud_opacity_scale,
                    gray_mix=self._gray_mix,
                    content_fov_deg=content_fov_deg,
                )
            composited = apply_ground_tint(
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
