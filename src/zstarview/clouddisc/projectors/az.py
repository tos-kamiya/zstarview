# -*- coding: utf-8 -*-
"""
Azimuthal projection and coordinate transformation utilities.

This module provides functions for projecting a celestial sphere view, as seen
by an observer, onto a 2D grid. It includes helpers for converting between
Geodetic, ECEF (Earth-Centered, Earth-Fixed), and local ENU (East, North, Up)
coordinate systems, which are essential for accurately mapping sky coordinates
to screen pixels.
"""

import math
from typing import Tuple

import numpy as np

# Standard Earth radius in kilometers, assuming a spherical Earth.
EARTH_R_KM = 6371.0


def deg2rad(degrees: float) -> float:
    """Converts degrees to radians."""
    return degrees * math.pi / 180.0


def enu_basis(lat_deg: float, lon_deg: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculates the East, North, Up (ENU) basis vectors for a given geodetic location.

    These vectors form a local Cartesian coordinate system at the observer's position.

    Args:
        lat_deg: Latitude in degrees.
        lon_deg: Longitude in degrees.

    Returns:
        A tuple containing the East, North, and Up basis vectors as 3D numpy arrays.
    """
    lat, lon = deg2rad(lat_deg), deg2rad(lon_deg)
    cos_lat, sin_lat = math.cos(lat), math.sin(lat)
    cos_lon, sin_lon = math.cos(lon), math.sin(lon)

    up = np.array([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])
    east = np.array([-sin_lon, cos_lon, 0.0])
    north = np.array([-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat])

    return east, north, up


def geodetic_to_ecef(lat_deg: float, lon_deg: float, r_km: float = EARTH_R_KM) -> np.ndarray:
    """
    Converts geodetic coordinates to Earth-Centered, Earth-Fixed (ECEF) coordinates.

    Args:
        lat_deg: Latitude in degrees.
        lon_deg: Longitude in degrees.
        r_km: The radius from the Earth's center in kilometers.

    Returns:
        The ECEF coordinates [x, y, z] as a numpy array.
    """
    lat, lon = deg2rad(lat_deg), deg2rad(lon_deg)
    cos_lat, sin_lat = np.cos(lat), np.sin(lat)
    cos_lon, sin_lon = np.cos(lon), np.sin(lon)

    x = r_km * cos_lat * cos_lon
    y = r_km * cos_lat * sin_lon
    z = r_km * sin_lat

    return np.array([x, y, z], dtype=float).T


def azalt_to_dir_ecef(az_deg: float, alt_deg: float, lat0_deg: float, lon0_deg: float) -> np.ndarray:
    """
    Converts a local Azimuth/Altitude direction to a normalized ECEF direction vector.

    Args:
        az_deg: Azimuth in degrees (0=N, 90=E).
        alt_deg: Altitude (elevation) in degrees.
        lat0_deg: Observer's latitude in degrees.
        lon0_deg: Observer's longitude in degrees.

    Returns:
        A normalized ECEF direction vector [x, y, z].
    """
    east, north, up = enu_basis(lat0_deg, lon0_deg)
    az, alt = deg2rad(az_deg), deg2rad(alt_deg)
    cos_alt, sin_alt = math.cos(alt), math.sin(alt)
    sin_az, cos_az = math.sin(az), math.cos(az)

    # Convert Az/Alt to a vector in the local ENU frame, then transform to ECEF.
    direction_vector = sin_alt * up + cos_alt * (sin_az * east + cos_az * north)

    return direction_vector / (np.linalg.norm(direction_vector) or 1.0)


def az_project_lonlat_grid(
    lat0_deg: float,
    lon0_deg: float,
    alt0_deg: float,
    az0_deg: float,
    radius_px: int,
    cloud_shell_km: float,
    alt_min_deg: float = -2.0,
    mask_fov_deg: float = 90.0,
    edge_fov_deg: float = 90.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Projects a 2D image grid onto a spherical shell to find the corresponding longitudes and latitudes.

    This is the core inverse projection function. For each pixel in a square output
    image, it calculates a ray from the observer and finds where that ray intersects
    a spherical "cloud shell". It then returns the geo-coordinates of that intersection.

    Args:
        lat0_deg: Observer's latitude in degrees.
        lon0_deg: Observer's longitude in degrees.
        alt0_deg: Observer's viewing altitude in degrees.
        az0_deg: Observer's viewing azimuth in degrees.
        radius_px: Radius of the output disc in pixels (image is 2R x 2R).
        cloud_shell_km: Radius of the cloud sphere in kilometers.
        alt_min_deg: Minimum altitude to be considered visible (horizon cutoff).
        mask_fov_deg: Field-of-view angle for the visibility mask.
        edge_fov_deg: Field-of-view angle that the disc radius represents (projection scale).

    Returns:
        A tuple (lon_grid, lat_grid, mask_inside) where grids are 2D numpy arrays.
    """
    image_size = 2 * radius_px

    # --- Parameter Validation ---
    fov_limit = edge_fov_deg * math.sqrt(2.0)
    if mask_fov_deg > fov_limit + 1e-6:
        raise ValueError(f"mask_fov_deg ({mask_fov_deg}°) exceeds geometric limit ({fov_limit:.2f}°)")
    if not (0 < edge_fov_deg <= 180):
        raise ValueError("edge_fov_deg must be in (0, 180]")

    # --- Step 1: Create a grid of pixel coordinates and project them to angular distances ---
    pixel_coords = (np.arange(image_size) + 0.5) - radius_px
    x, y = np.meshgrid(pixel_coords, -pixel_coords)  # Y is inverted for image coordinates
    pixel_radius = np.hypot(x, y)
    # rho is the angular distance from the view center for each pixel.
    rho_deg = (edge_fov_deg * pixel_radius / radius_px).astype(np.float32)

    # --- Step 2: Define the observer's reference frame in ECEF coordinates ---
    observer_pos_ecef = geodetic_to_ecef(lat0_deg, lon0_deg)
    _, _, up_vec = enu_basis(lat0_deg, lon0_deg)
    center_view_dir = azalt_to_dir_ecef(az0_deg, alt0_deg, lat0_deg, lon0_deg)

    # Create an orthonormal basis for the view plane (tangent plane).
    ty_unnormalized = up_vec - np.dot(up_vec, center_view_dir) * center_view_dir
    ty = ty_unnormalized / (np.linalg.norm(ty_unnormalized) or 1.0)
    tx = np.cross(center_view_dir, ty)

    # --- Step 3: Calculate the ECEF direction vector for each pixel's ray ---
    psi = np.arctan2(y, x)  # Angle of each pixel on the view plane
    cos_rho, sin_rho = np.cos(np.deg2rad(rho_deg)), np.sin(np.deg2rad(rho_deg))
    b = np.cos(psi)[..., None] * tx + np.sin(psi)[..., None] * ty  # Direction on tangent plane
    d = cos_rho[..., None] * center_view_dir + sin_rho[..., None] * b  # Final direction vector
    d /= np.linalg.norm(d, axis=2, keepdims=True) + 1e-12  # Normalize

    # --- Step 4: Create visibility masks ---
    alt_rad = np.arcsin(np.dot(d, up_vec))
    visible_mask = np.degrees(alt_rad) >= alt_min_deg  # Above the horizon
    fov_mask = rho_deg <= mask_fov_deg + 1e-6  # Within the specified field of view
    mask_inside = fov_mask & visible_mask

    # --- Step 5: Solve for the intersection of each ray with the cloud sphere ---
    # This is a standard line-sphere intersection problem, solving a quadratic equation for t.
    b_quad = 2.0 * np.sum(observer_pos_ecef * d, axis=2)
    c_quad = np.dot(observer_pos_ecef, observer_pos_ecef) - cloud_shell_km * cloud_shell_km
    discriminant = b_quad * b_quad - 4.0 * c_quad
    valid_intersection = discriminant >= 0

    sqrt_disc = np.sqrt(np.maximum(discriminant, 0.0))
    t1 = (-b_quad - sqrt_disc) / 2.0
    t2 = (-b_quad + sqrt_disc) / 2.0
    # Find the smallest positive intersection distance `t`.
    t = np.where(t1 > 1e-6, t1, np.where(t2 > 1e-6, t2, np.nan))

    # --- Step 6: Convert intersection points back to longitude and latitude ---
    lon_grid = np.full((image_size, image_size), np.nan, dtype=np.float32)
    lat_grid = np.full((image_size, image_size), np.nan, dtype=np.float32)

    # Calculate ECEF coordinates of the intersection points.
    P = observer_pos_ecef + d * t[..., None]
    x_int, y_int, z_int = P[..., 0], P[..., 1], P[..., 2]

    # Convert ECEF back to geodetic coordinates.
    lon_grid[valid_intersection] = np.degrees(np.arctan2(y_int, x_int))[valid_intersection]
    hyp = np.hypot(x_int, y_int)
    lat_grid[valid_intersection] = np.degrees(np.arctan2(z_int, hyp))[valid_intersection]

    # Apply the final visibility mask.
    lon_grid[~mask_inside] = np.nan
    lat_grid[~mask_inside] = np.nan

    return lon_grid, lat_grid, mask_inside
