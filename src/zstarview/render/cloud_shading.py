"""Cloud cell illumination and compositing for the rendering prototype."""

from __future__ import annotations

import numpy as np

from ..clouddisc.altaz_grid import (
    CloudAltAzGrid,
    _altaz_grid_centers,
    _intersect_altaz_rays_with_shell,
)
from ..clouddisc.altaz_projection import altaz_to_bin_indices, altaz_to_dir_ecef_array
from ..clouddisc.projectors.az import geodetic_to_ecef

EARTH_RADIUS_KM = 6371.0
SUN_EXTINCTION = 1.8
VIEW_EXTINCTION = 2.4
AMBIENT_FLOOR = 0.18


def _ecef_to_altaz(
    points: np.ndarray, observer_lat_deg: float, observer_lon_deg: float
) -> tuple[np.ndarray, np.ndarray]:
    """Convert ECEF points to observer-relative altitude and azimuth."""
    lat = np.radians(observer_lat_deg)
    lon = np.radians(observer_lon_deg)
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north = np.array(
        [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)]
    )
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    observer = geodetic_to_ecef(observer_lat_deg, observer_lon_deg)
    relative = points - observer
    e = relative @ east
    n = relative @ north
    u = relative @ up
    alt = np.degrees(np.arctan2(u, np.hypot(e, n)))
    az = np.mod(np.degrees(np.arctan2(e, n)), 360.0)
    return alt, az


def _sample_grid(
    grid: CloudAltAzGrid,
    values: np.ndarray,
    alt: np.ndarray,
    az: np.ndarray,
) -> np.ndarray:
    ai, zi = altaz_to_bin_indices(
        alt,
        az,
        alt_bins=grid.amount.shape[0],
        az_bins=grid.amount.shape[1],
        alt_min_deg=grid.alt_min_deg,
        alt_max_deg=grid.alt_max_deg,
        az_min_deg=grid.az_min_deg,
        az_max_deg=grid.az_max_deg,
    )
    in_range = (
        (alt >= grid.alt_min_deg)
        & (alt <= grid.alt_max_deg)
        & (az >= grid.az_min_deg)
        & (az <= grid.az_max_deg)
    )
    result = np.zeros(alt.shape, dtype=np.float32)
    result[in_range] = values[ai[in_range], zi[in_range]]
    missing = np.zeros(alt.shape, dtype=bool)
    missing[in_range] = grid.missing_mask[ai[in_range], zi[in_range]] > 0
    result[missing] = 0.0
    return result


def _sun_transmission(
    grid: CloudAltAzGrid, sun_alt_deg: float, sun_az_deg: float
) -> tuple[np.ndarray, ...]:
    """Estimate direct sunlight through upper cloud groups at each cell."""
    groups = grid.physical_shell_amounts or grid.shell_amounts
    if not groups:
        groups = (grid.amount,)
    alt_centers, az_centers = _altaz_grid_centers(
        *grid.amount.shape,
        alt_min_deg=grid.alt_min_deg,
        alt_max_deg=grid.alt_max_deg,
        az_min_deg=grid.az_min_deg,
        az_max_deg=grid.az_max_deg,
    )
    alt_centers = np.broadcast_to(alt_centers, grid.amount.shape)
    az_centers = np.broadcast_to(az_centers, grid.amount.shape)
    shells = grid.shells_km
    if len(groups) == 3 and len(shells) >= 9:
        radii = tuple(float(np.mean(shells[i : i + 3])) for i in (0, 3, 6))
    elif len(shells) == len(groups):
        radii = tuple(float(value) for value in shells)
    else:
        radii = tuple(float(v) for v in shells[: len(groups)])
    if len(radii) != len(groups):
        radii = tuple(EARTH_RADIUS_KM + 5.0 + 2.0 * i for i in range(len(groups)))

    sun_dirs = altaz_to_dir_ecef_array(
        np.full(grid.amount.shape, sun_alt_deg, dtype=np.float32),
        np.full(grid.amount.shape, sun_az_deg, dtype=np.float32),
        grid.observer_lat,
        grid.observer_lon,
    )
    transmissions: list[np.ndarray] = []
    for group_index, radius in enumerate(radii):
        lon, lat, valid = _intersect_altaz_rays_with_shell(
            alt_centers,
            az_centers,
            observer_lat=grid.observer_lat,
            observer_lon=grid.observer_lon,
            shell_km=radius,
        )
        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)
        points = np.stack(
            (
                radius * np.cos(lat_rad) * np.cos(lon_rad),
                radius * np.cos(lat_rad) * np.sin(lon_rad),
                radius * np.sin(lat_rad),
            ),
            axis=-1,
        )
        # Sunlight reaching a cell is reduced by cloud intersections in higher groups.
        tau = np.zeros(grid.amount.shape, dtype=np.float32)
        for upper_index in range(group_index + 1, len(groups)):
            upper_radius = radii[upper_index]
            direction = sun_dirs
            b = np.sum(points * direction, axis=-1)
            c = np.sum(points * points, axis=-1) - upper_radius * upper_radius
            discriminant = b * b - c
            t = -b + np.sqrt(np.maximum(discriminant, 0.0))
            intersects = valid & (discriminant >= 0.0) & (t > 0.0)
            upper_points = points + t[..., None] * direction
            upper_alt, upper_az = _ecef_to_altaz(
                upper_points, grid.observer_lat, grid.observer_lon
            )
            amount = _sample_grid(grid, groups[upper_index], upper_alt, upper_az)
            radial = points / np.maximum(
                np.linalg.norm(points, axis=-1, keepdims=True), 1.0e-6
            )
            incidence = np.maximum(0.18, np.abs(np.sum(radial * direction, axis=-1)))
            tau += np.where(intersects, amount * SUN_EXTINCTION / incidence, 0.0)
        # The Earth blocks direct sunlight when the sun is below the local horizon.
        earth_b = np.sum(points * sun_dirs, axis=-1)
        earth_c = np.sum(points * points, axis=-1) - EARTH_RADIUS_KM**2
        earth_disc = earth_b * earth_b - earth_c
        earth_near = -earth_b - np.sqrt(np.maximum(earth_disc, 0.0))
        earth_blocked = (earth_disc >= 0.0) & (earth_near > 1.0e-3)
        transmissions.append(
            np.where(
                earth_blocked,
                AMBIENT_FLOOR,
                AMBIENT_FLOOR + (1.0 - AMBIENT_FLOOR) * np.exp(-tau),
            )
        )
    return tuple(transmissions)


def shade_cloud_cells(
    grid: CloudAltAzGrid,
    base_rgb: np.ndarray,
    altitudes_deg: np.ndarray,
    azimuths_deg: np.ndarray,
    sun_alt_deg: float,
    sun_az_deg: float,
    *,
    opacity: float = 0.85,
) -> np.ndarray:
    """Return RGB cloud compositing over the base sky, one row per input ray."""
    groups = grid.physical_shell_amounts or grid.shell_amounts or (grid.amount,)
    transmissions = _sun_transmission(grid, sun_alt_deg, sun_az_deg)
    ai, zi = altaz_to_bin_indices(
        altitudes_deg,
        azimuths_deg,
        alt_bins=grid.amount.shape[0],
        az_bins=grid.amount.shape[1],
        alt_min_deg=grid.alt_min_deg,
        alt_max_deg=grid.alt_max_deg,
        az_min_deg=grid.az_min_deg,
        az_max_deg=grid.az_max_deg,
    )
    in_range = (
        (altitudes_deg >= grid.alt_min_deg)
        & (altitudes_deg <= grid.alt_max_deg)
        & (azimuths_deg >= grid.az_min_deg)
        & (azimuths_deg <= grid.az_max_deg)
    )
    rgb = np.asarray(base_rgb, dtype=np.float32).reshape((-1, 3)).copy()
    for index in reversed(range(len(groups))):
        amount = np.zeros(altitudes_deg.shape, dtype=np.float32)
        amount[in_range] = groups[index][ai[in_range], zi[in_range]]
        missing = np.zeros(altitudes_deg.shape, dtype=bool)
        missing[in_range] = grid.missing_mask[ai[in_range], zi[in_range]] > 0
        amount[missing] = 0.0
        sun = np.ones(altitudes_deg.shape, dtype=np.float32)
        sun[in_range] = transmissions[min(index, len(transmissions) - 1)][
            ai[in_range], zi[in_range]
        ]
        # Saturating white cloud radiance with a slight warm daylight tint.
        density = 1.0 - np.exp(-3.0 * amount)
        cloud = (
            density[:, None]
            * sun[:, None]
            * np.array([0.96, 0.975, 1.0], dtype=np.float32)
        )
        transmittance = np.exp(
            -VIEW_EXTINCTION * amount * float(np.clip(opacity, 0.0, 1.0))
        )
        rgb = rgb * transmittance[:, None] + cloud * (1.0 - transmittance[:, None])
    return np.clip(rgb, 0.0, 1.0)
