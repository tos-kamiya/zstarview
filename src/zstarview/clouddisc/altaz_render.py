"""Render a `CloudAltAzGrid` into a screen-space RGBA image.

The MVP renderer draws white discs whose radius and opacity scale with the
cloud amount stored in each alt/az cell.
"""

from __future__ import annotations

import logging

import numpy as np

from .._numba_kernels import stamp_cloud_cells_numba
from .altaz_constants import (
    ALT_AZ_CIRCLE_AMOUNT_THRESHOLD,
    ALT_AZ_CIRCLE_BASE_RADIUS_PX,
    ALT_AZ_CIRCLE_MAX_RADIUS_PX,
    ALT_AZ_CIRCLE_OPACITY_SCALE,
)
from .altaz_grid import CloudAltAzGrid
from .altaz_projection import altaz_to_screen_coords

logger = logging.getLogger(__name__)


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
    """Render cloud cells as white discs with mild 8-neighbor softening.

    Projection remains vectorized in NumPy. The compiled kernel performs the
    per-cell disc rasterization, alpha accumulation, and a small local alpha
    filter without allocating a separate stamp array for every cell.
    """
    w = max(1, int(width))
    h = max(1, int(height))

    amount = grid.amount.astype(np.float32, copy=False)
    active_alt_idx, active_az_idx = np.nonzero(amount > ALT_AZ_CIRCLE_AMOUNT_THRESHOLD)
    if active_alt_idx.size == 0:
        return np.zeros((h, w, 4), dtype=np.uint8)

    alt_centers = _bin_centers(
        grid.amount.shape[0], grid.alt_min_deg, grid.alt_max_deg
    )
    az_centers = _bin_centers(
        grid.amount.shape[1], grid.az_min_deg, grid.az_max_deg
    )
    active_alt = alt_centers[active_alt_idx]
    active_az = az_centers[active_az_idx]
    active_amount = amount[active_alt_idx, active_az_idx]

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

    base_r = max(0.5, float(base_radius_px))
    max_r = max(base_r + 0.5, float(max_radius_px))
    alpha_buffer = stamp_cloud_cells_numba(
        np.ascontiguousarray(x_px[valid], dtype=np.float64),
        np.ascontiguousarray(y_px[valid], dtype=np.float64),
        np.ascontiguousarray(active_amount[valid], dtype=np.float32),
        h,
        w,
        base_r,
        max_r,
        float(opacity_scale),
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

    disc_size = disc.shape[0]
    for ix in range(x_px.size):
        cx = int(round(float(x_px[ix])))
        cy = int(round(float(y_px[ix])))
        y0 = max(0, cy - radius)
        y1 = min(h, cy + radius + 1)
        x0 = max(0, cx - radius)
        x1 = min(w, cx + radius + 1)
        if y0 >= y1 or x0 >= x1:
            continue
        sy0 = y0 - (cy - radius)
        sy1 = sy0 + (y1 - y0)
        sx0 = x0 - (cx - radius)
        sx1 = sx0 + (x1 - x0)
        # Clamp source slice to the actual disc footprint.
        if sy0 < 0:
            y0 -= sy0
            sy0 = 0
        if sx0 < 0:
            x0 -= sx0
            sx0 = 0
        if sy1 > disc_size:
            y1 -= sy1 - disc_size
            sy1 = disc_size
        if sx1 > disc_size:
            x1 -= sx1 - disc_size
            sx1 = disc_size
        if y0 >= y1 or x0 >= x1 or sy0 >= sy1 or sx0 >= sx1:
            continue
        out[y0:y1, x0:x1] = np.maximum(
            out[y0:y1, x0:x1], disc[sy0:sy1, sx0:sx1] * np.uint8(255)
        )

    return out
