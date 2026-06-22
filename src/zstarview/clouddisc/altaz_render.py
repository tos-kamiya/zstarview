# -*- coding: utf-8 -*-
"""Render a `CloudAltAzGrid` into a screen-space RGBA image.

The MVP renderer draws soft white circles whose radius and opacity scale
with the cloud amount stored in each alt/az cell.
"""

from __future__ import annotations

import logging

import numpy as np

from .altaz_constants import (
    ALT_AZ_CIRCLE_AMOUNT_THRESHOLD,
    ALT_AZ_CIRCLE_BASE_RADIUS_PX,
    ALT_AZ_CIRCLE_MAX_RADIUS_PX,
    ALT_AZ_CIRCLE_OPACITY_SCALE,
)
from .altaz_grid import CloudAltAzGrid
from .altaz_projection import altaz_to_screen_coords

logger = logging.getLogger(__name__)


def _gaussian_stamp(radius_px: float) -> np.ndarray:
    """Return a 2D gaussian stamp of radius ``radius_px`` (sigma = radius/2)."""
    r = max(1, int(np.ceil(radius_px)))
    y, x = np.mgrid[-r : r + 1, -r : r + 1]
    dist_sq = x * x + y * y
    sigma = max(0.5, radius_px / 2.0)
    return np.exp(-dist_sq / (2.0 * sigma * sigma))


def _bin_centers(bins: int, min_val: float, max_val: float) -> np.ndarray:
    """Return center coordinates for each bin in a uniform grid."""
    edges = np.linspace(min_val, max_val, bins + 1)
    return (edges[:-1] + edges[1:]) * 0.5


def render_altaz_grid_circles(
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    *,
    center_alt_deg: float,
    center_az_deg: float,
    edge_fov_deg: float,
    mask_fov_deg: float = 90.0,
    base_radius_px: float = ALT_AZ_CIRCLE_BASE_RADIUS_PX,
    max_radius_px: float = ALT_AZ_CIRCLE_MAX_RADIUS_PX,
    opacity_scale: float = ALT_AZ_CIRCLE_OPACITY_SCALE,
) -> np.ndarray:
    """Render cloud cells as white circles into a (H, W, 4) uint8 RGBA image.

    Args:
        grid: the camera-independent alt/az cloud grid.
        width / height: output image size.
        center_alt_deg / center_az_deg: view center direction.
        edge_fov_deg: angular radius of the output disc in degrees.
        mask_fov_deg: field-of-view used to clip invisible directions.
        base_radius_px: minimum circle radius in pixels.
        max_radius_px: maximum circle radius at cloud amount = 1.
        opacity_scale: global alpha multiplier.

    Returns:
        uint8 RGBA array of shape ``(height, width, 4)``.
    """
    w = max(1, int(width))
    h = max(1, int(height))

    alt_centers = _bin_centers(
        grid.amount.shape[0], grid.alt_min_deg, grid.alt_max_deg
    )
    az_centers = _bin_centers(
        grid.amount.shape[1], grid.az_min_deg, grid.az_max_deg
    )
    alt_grid, az_grid = np.meshgrid(alt_centers, az_centers, indexing="ij")

    amount = grid.amount.astype(np.float32, copy=False)
    active = amount > ALT_AZ_CIRCLE_AMOUNT_THRESHOLD
    if not np.any(active):
        return np.zeros((h, w, 4), dtype=np.uint8)

    active_alt = alt_grid[active]
    active_az = az_grid[active]
    active_amount = amount[active]

    x_px, y_px = altaz_to_screen_coords(
        active_alt,
        active_az,
        width=w,
        height=h,
        center_alt_deg=center_alt_deg,
        center_az_deg=center_az_deg,
        edge_fov_deg=edge_fov_deg,
        mask_fov_deg=mask_fov_deg,
        observer_lat_deg=grid.observer_lat,
        observer_lon_deg=grid.observer_lon,
    )

    valid = np.isfinite(x_px) & np.isfinite(y_px)
    if not np.any(valid):
        return np.zeros((h, w, 4), dtype=np.uint8)

    x_px = x_px[valid]
    y_px = y_px[valid]
    active_amount = active_amount[valid]

    alpha_buffer = np.zeros((h, w), dtype=np.float32)
    base_r = max(0.5, float(base_radius_px))
    max_r = max(base_r + 0.5, float(max_radius_px))

    for ix in range(x_px.size):
        cx = int(round(float(x_px[ix])))
        cy = int(round(float(y_px[ix])))
        amount_v = float(active_amount[ix])
        if amount_v <= 0.0 or not (0 <= cx < w and 0 <= cy < h):
            continue

        radius = base_r + amount_v * (max_r - base_r)
        stamp = _gaussian_stamp(radius)
        alpha_max = amount_v * opacity_scale

        stamp_h, stamp_w = stamp.shape
        y0 = max(0, cy - stamp_h // 2)
        y1 = min(h, cy - stamp_h // 2 + stamp_h)
        x0 = max(0, cx - stamp_w // 2)
        x1 = min(w, cx - stamp_w // 2 + stamp_w)

        stamp_y0 = y0 - (cy - stamp_h // 2)
        stamp_y1 = stamp_y0 + (y1 - y0)
        stamp_x0 = x0 - (cx - stamp_w // 2)
        stamp_x1 = stamp_x0 + (x1 - x0)

        alpha_buffer[y0:y1, x0:x1] += (
            stamp[stamp_y0:stamp_y1, stamp_x0:stamp_x1] * alpha_max
        )

    np.clip(alpha_buffer, 0.0, 1.0, out=alpha_buffer)
    alpha_u8 = (alpha_buffer * 255.0).astype(np.uint8)

    out = np.zeros((h, w, 4), dtype=np.uint8)
    positive = alpha_u8 > 0
    out[..., :3][positive] = 255
    out[..., 3] = alpha_u8
    return out


def render_altaz_missing_mask(
    grid: CloudAltAzGrid,
    width: int,
    height: int,
    *,
    center_alt_deg: float,
    center_az_deg: float,
    edge_fov_deg: float,
    mask_fov_deg: float = 90.0,
    stamp_radius_px: float = 2.0,
) -> np.ndarray:
    """Project the alt/az missing-data mask to a screen-space uint8 alpha image.

    Each missing alt/az cell is forward-projected to screen space and drawn as a
    small solid disc.  This is much faster than an inverse nearest-neighbour
    search over the full pixel grid.

    Returns:
        uint8 array of shape ``(height, width)`` with 0 / 255 values.
    """
    w = max(1, int(width))
    h = max(1, int(height))
    out = np.zeros((h, w), dtype=np.uint8)

    missing_cells = grid.missing_mask > 0
    if not np.any(missing_cells):
        return out

    alt_centers = _bin_centers(
        grid.amount.shape[0], grid.alt_min_deg, grid.alt_max_deg
    )
    az_centers = _bin_centers(
        grid.amount.shape[1], grid.az_min_deg, grid.az_max_deg
    )
    alt_grid, az_grid = np.meshgrid(alt_centers, az_centers, indexing="ij")

    x_px, y_px = altaz_to_screen_coords(
        alt_grid[missing_cells],
        az_grid[missing_cells],
        width=w,
        height=h,
        center_alt_deg=center_alt_deg,
        center_az_deg=center_az_deg,
        edge_fov_deg=edge_fov_deg,
        mask_fov_deg=mask_fov_deg,
        observer_lat_deg=grid.observer_lat,
        observer_lon_deg=grid.observer_lon,
    )

    valid = np.isfinite(x_px) & np.isfinite(y_px)
    x_px = x_px[valid]
    y_px = y_px[valid]

    radius = max(1, int(round(stamp_radius_px)))
    y_idx, x_idx = np.mgrid[-radius : radius + 1, -radius : radius + 1]
    disc = x_idx * x_idx + y_idx * y_idx <= radius * radius

    for ix in range(x_px.size):
        cx = int(round(float(x_px[ix])))
        cy = int(round(float(y_px[ix])))
        y0 = max(0, cy - radius)
        y1 = min(h, cy + radius + 1)
        x0 = max(0, cx - radius)
        x1 = min(w, cx + radius + 1)
        dy0 = y0 - (cy - radius)
        dy1 = dy0 + (y1 - y0)
        dx0 = x0 - (cx - radius)
        dx1 = dx0 + (x1 - x0)
        out[y0:y1, x0:x1] = np.maximum(
            out[y0:y1, x0:x1], disc[dy0:dy1, dx0:dx1] * np.uint8(255)
        )

    return out
