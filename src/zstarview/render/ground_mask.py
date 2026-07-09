"""Pixel masks for the region below the terrain horizon."""

from __future__ import annotations

import math

import numpy as np

from ..types import ScreenGeometry


def inverse_project_disc(
    width: int,
    height: int,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    *,
    edge_fov_deg: float,
    content_fov_deg: float,
    origin: tuple[float, float] = (0.0, 0.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse-project image pixels into altitude/azimuth coordinates."""
    cx = float(geometry.center[0])
    cy = float(geometry.center[1])
    origin_x, origin_y = origin
    radius = max(1.0, float(geometry.radius))
    ys = (np.arange(height, dtype=np.float32) + origin_y - cy) / radius
    xs = (np.arange(width, dtype=np.float32) + origin_x - cx) / radius
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
    y = np.sin(psi) * sin_theta * cos_phi1
    x = cos_theta - sin_phi1 * sin_phi2
    lam2 = lam1 + np.arctan2(y, x)
    alt = np.degrees(phi2).astype(np.float32)
    az = (np.degrees(lam2) + 360.0) % 360.0
    return alt, az.astype(np.float32), inside


def interpolate_terrain_horizon_altitude(
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


def build_ground_mask(
    width: int,
    height: int,
    geometry: ScreenGeometry,
    view_center: tuple[float, float],
    terrain_profile_altaz: list[tuple[float, float]] | None,
    *,
    edge_fov_deg: float,
    content_fov_deg: float,
    origin: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Return a pixel mask for points below the terrain horizon."""
    alt, az, inside = inverse_project_disc(
        width,
        height,
        geometry,
        view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
        origin=origin,
    )
    mask = np.zeros((height, width), dtype=bool)
    if alt.size == 0:
        return mask
    horizon_alt = interpolate_terrain_horizon_altitude(az, terrain_profile_altaz)
    mask[inside] = alt < horizon_alt
    return mask
