# -*- coding: utf-8 -*-
"""Coordinate conversions for the alt/az cloud grid.

This module contains the geometric helpers needed both during ingestion
(geographic + height → altitude/azimuth) and during rendering
(altitude/azimuth → screen pixels).
"""

from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

from .projectors.az import altaz_to_dir_ecef, enu_basis, geodetic_to_ecef

logger = logging.getLogger(__name__)


def geodetic_to_altaz(
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> Tuple[float, float]:
    """Convert a point on a cloud shell to (altitude, azimuth) for an observer.

    Args:
        lat_deg: point latitude in degrees.
        lon_deg: point longitude in degrees.
        radius_km: distance from Earth's center in km (e.g. 6374 km for 3 km
            above a spherical Earth).
        observer_lat_deg: observer latitude in degrees.
        observer_lon_deg: observer longitude in degrees.

    Returns:
        (altitude_deg, azimuth_deg) where azimuth is 0° at north and increases
        eastward to 360°.
    """
    observer_pos = geodetic_to_ecef(observer_lat_deg, observer_lon_deg)
    point_pos = geodetic_to_ecef(lat_deg, lon_deg, r_km=radius_km)
    vec = point_pos - observer_pos
    dist = float(np.linalg.norm(vec))
    if dist < 1e-9:
        return 90.0, 0.0

    east, north, up = enu_basis(observer_lat_deg, observer_lon_deg)
    enu = np.array([np.dot(vec, east), np.dot(vec, north), np.dot(vec, up)])
    alt_rad = np.arcsin(np.clip(enu[2] / dist, -1.0, 1.0))
    az_rad = np.arctan2(enu[0], enu[1])

    alt_deg = float(np.degrees(alt_rad))
    az_deg = float(np.degrees(az_rad) % 360.0)
    return alt_deg, az_deg


def altaz_to_screen_coords(
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    *,
    width: int,
    height: int,
    center_alt_deg: float,
    center_az_deg: float,
    edge_fov_deg: float,
    mask_fov_deg: float = 90.0,
    observer_lat_deg: float = 0.0,
    observer_lon_deg: float = 0.0,
    radius_px: float | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Project arrays of alt/az directions to screen pixel coordinates.

    Uses an azimuthal-equidistant projection centred on
    ``(center_alt_deg, center_az_deg)``.  Pixel coordinates follow image
    conventions: x increases right, y increases down.

    Args:
        alt_deg / az_deg: arrays of directions in degrees.
        width / height: output image dimensions.
        center_alt_deg / center_az_deg: centre of the view in degrees.
        edge_fov_deg: angular radius represented by the disc edge.
        mask_fov_deg: directions farther than this from the centre are clipped.
        observer_lat_deg / observer_lon_deg: observer location used to build
            the local ENU frame.
        radius_px: optional pixel radius; defaults to ``min(width, height)/2``.

    Returns:
        (x_px, y_px) arrays with the same shape as the inputs; points outside
        the visible cone are set to NaN.
    """
    if radius_px is None:
        radius_px = min(width, height) / 2.0

    cx = (width - 1) * 0.5
    cy = (height - 1) * 0.5

    # Direction vectors for the view centre and each input direction.
    center_dir = altaz_to_dir_ecef(center_alt_deg, center_az_deg, observer_lat_deg, observer_lon_deg)
    dir_vecs = altaz_to_dir_ecef_array(alt_deg, az_deg, observer_lat_deg, observer_lon_deg)

    # Build an orthonormal tangent basis at the view centre.
    _, _, up = enu_basis(observer_lat_deg, observer_lon_deg)
    ty_unnormalized = up - np.dot(up, center_dir) * center_dir
    ty_norm = float(np.linalg.norm(ty_unnormalized)) or 1.0
    ty = ty_unnormalized / ty_norm
    tx = np.cross(center_dir, ty)

    # Angular distance from the centre for each direction.
    cos_rho = np.clip(np.sum(dir_vecs * center_dir, axis=-1), -1.0, 1.0)
    rho_rad = np.arccos(cos_rho)
    rho_deg = np.degrees(rho_rad)

    # Project each direction onto the tangent plane to recover the bearing.
    sin_rho = np.maximum(np.sin(rho_rad), 1e-12)
    tangent = dir_vecs - cos_rho[..., None] * center_dir
    u = np.sum(tangent * tx, axis=-1) / sin_rho
    v = np.sum(tangent * ty, axis=-1) / sin_rho
    psi = np.arctan2(v, u)

    # Map angular distance to pixel radius.
    pixel_radius = rho_deg * (radius_px / max(1e-6, float(edge_fov_deg)))
    x_px = cx + pixel_radius * np.cos(psi)
    y_px = cy - pixel_radius * np.sin(psi)

    # Mask out directions behind the view plane or outside the mask FOV.
    valid = (cos_rho > 0.0) & (rho_deg <= mask_fov_deg + 1e-6)
    x_px = np.where(valid, x_px, np.nan)
    y_px = np.where(valid, y_px, np.nan)

    return x_px, y_px


def altaz_to_dir_ecef_array(
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    lat0_deg: float,
    lon0_deg: float,
) -> np.ndarray:
    """Vectorised `altaz_to_dir_ecef`.

    Returns an array of shape ``alt_deg.shape + (3,)`` containing normalised
    ECEF direction vectors.
    """
    east, north, up = enu_basis(lat0_deg, lon0_deg)
    alt = np.radians(alt_deg)
    az = np.radians(az_deg)
    cos_alt, sin_alt = np.cos(alt), np.sin(alt)
    sin_az, cos_az = np.sin(az), np.cos(az)

    direction = (
        sin_alt[..., None] * up
        + cos_alt[..., None] * (sin_az[..., None] * east + cos_az[..., None] * north)
    )
    norm = np.linalg.norm(direction, axis=-1, keepdims=True)
    norm = np.maximum(norm, 1e-12)
    return direction / norm


def altaz_to_bin_indices(
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    *,
    alt_bins: int,
    az_bins: int,
    alt_min_deg: float = 0.0,
    alt_max_deg: float = 90.0,
    az_min_deg: float = 0.0,
    az_max_deg: float = 360.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Map altitude/azimuth arrays to integer grid indices.

    Out-of-range altitudes are clipped to the nearest edge bin. Azimuths are
    wrapped to [az_min, az_max).

    Returns:
        (alt_idx, az_idx) integer arrays of the same shape as the inputs.
    """
    alt_span = max(1e-9, alt_max_deg - alt_min_deg)
    az_span = max(1e-9, az_max_deg - az_min_deg)

    az_normalized = (az_deg - az_min_deg) % az_span
    az_idx = (az_normalized / az_span * az_bins).astype(np.int64)
    az_idx = np.clip(az_idx, 0, az_bins - 1)

    alt_idx = ((alt_deg - alt_min_deg) / alt_span * alt_bins).astype(np.int64)
    alt_idx = np.clip(alt_idx, 0, alt_bins - 1)

    return alt_idx, az_idx


def _geodetic_to_ecef_array(
    lat_deg: np.ndarray,
    lon_deg: np.ndarray,
    radius_km: float,
) -> np.ndarray:
    """Vectorized ECEF conversion with shape (..., 3)."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    cos_lat, sin_lat = np.cos(lat), np.sin(lat)
    cos_lon, sin_lon = np.cos(lon), np.sin(lon)
    x = radius_km * cos_lat * cos_lon
    y = radius_km * cos_lat * sin_lon
    z = radius_km * sin_lat
    return np.stack([x, y, z], axis=-1)


def geodetic_to_altaz_array(
    lat_deg: np.ndarray,
    lon_deg: np.ndarray,
    radius_km: float,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorized version of `geodetic_to_altaz`.

    Returns altitude and azimuth arrays of the same shape as the inputs.
    """
    observer_pos = geodetic_to_ecef(observer_lat_deg, observer_lon_deg)
    point_pos = _geodetic_to_ecef_array(lat_deg, lon_deg, radius_km)
    vec = point_pos - observer_pos
    dist = np.linalg.norm(vec, axis=-1)

    east, north, up = enu_basis(observer_lat_deg, observer_lon_deg)
    enu_east = vec[..., 0] * east[0] + vec[..., 1] * east[1] + vec[..., 2] * east[2]
    enu_north = vec[..., 0] * north[0] + vec[..., 1] * north[1] + vec[..., 2] * north[2]
    enu_up = vec[..., 0] * up[0] + vec[..., 1] * up[1] + vec[..., 2] * up[2]

    alt_rad = np.arcsin(np.clip(enu_up / np.maximum(dist, 1e-9), -1.0, 1.0))
    az_rad = np.arctan2(enu_east, enu_north)

    alt_deg = np.degrees(alt_rad)
    az_deg = np.degrees(az_rad) % 360.0
    return alt_deg, az_deg
