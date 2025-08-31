"""
Azimuthal projection and coordinate transformation utilities.

This module provides functions for converting between different coordinate systems
(Geodetic, ECEF) and for projecting a celestial sphere view onto a 2D grid.
"""

import math
from typing import Optional, Tuple

import numpy as np

# Standard Earth radius in kilometers
EARTH_R_KM = 6371.0


def deg2rad(degrees: float) -> float:
    """Converts degrees to radians."""
    return degrees * math.pi / 180.0


# def rad2deg(radians: float) -> float:
#     """Converts radians to degrees."""
#     return radians * 180.0 / math.pi


def enu_basis(lat_deg: float, lon_deg: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculates the East, North, Up (ENU) basis vectors for a given geodetic location.

    Args:
        lat_deg: Latitude in degrees.
        lon_deg: Longitude in degrees.

    Returns:
        A tuple containing the East, North, and Up basis vectors as numpy arrays.
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
    Converts geodetic coordinates (latitude, longitude) to Earth-Centered,
    Earth-Fixed (ECEF) Cartesian coordinates.

    Can handle scalar or numpy array inputs.

    Args:
        lat_deg: Latitude in degrees.
        lon_deg: Longitude in degrees.
        r_km: The radius of the Earth in kilometers.

    Returns:
        The ECEF coordinates [x, y, z] as a numpy array.
        Shape is (3,) for scalar inputs, (N, 3) for array inputs of length N.
    """
    lat, lon = deg2rad(lat_deg), deg2rad(lon_deg)
    cos_lat, sin_lat = np.cos(lat), np.sin(lat)
    cos_lon, sin_lon = np.cos(lon), np.sin(lon)

    x = r_km * cos_lat * cos_lon
    y = r_km * cos_lat * sin_lon
    z = r_km * sin_lat

    # For scalar inputs, this produces a (3,) array.
    # For array inputs (N,), this produces a (N, 3) array.
    return np.array([x, y, z], dtype=float).T


def azalt_to_dir_ecef(az_deg: float, alt_deg: float, lat0_deg: float, lon0_deg: float) -> np.ndarray:
    """
    Converts Azimuth/Altitude at a given location to a normalized ECEF direction vector.

    Args:
        az_deg: Azimuth in degrees.
        alt_deg: Altitude in degrees.
        lat0_deg: Observer's latitude in degrees.
        lon0_deg: Observer's longitude in degrees.

    Returns:
        A normalized ECEF direction vector [x, y, z] as a numpy array.
    """
    east, north, up = enu_basis(lat0_deg, lon0_deg)
    az, alt = deg2rad(az_deg), deg2rad(alt_deg)
    cos_alt, sin_alt = math.cos(alt), math.sin(alt)
    sin_az, cos_az = math.sin(az), math.cos(az)

    # Direction vector in local ENU coordinates, then converted to ECEF
    direction_vector = sin_alt * up + cos_alt * (sin_az * east + cos_az * north)

    return direction_vector / (np.linalg.norm(direction_vector) or 1.0)


# def line_sphere_intersection(origin: np.ndarray, direction: np.ndarray, radius_km: float) -> Optional[float]:
#     """
#     Finds the intersection of a line (ray) and a sphere.

#     Solves the equation |origin + t * direction| = radius_km for t.

#     Args:
#         origin: The origin of the line (e.g., observer's ECEF position).
#         direction: The normalized direction vector of the line.
#         radius_km: The radius of the sphere.

#     Returns:
#         The smallest positive distance (t) to an intersection, or None if no intersection.
#     """
#     # Quadratic equation coefficients for t: a*t^2 + b*t + c = 0
#     # Since direction is normalized, a = 1.
#     b = 2.0 * np.dot(origin, direction)
#     c = np.dot(origin, origin) - radius_km * radius_km

#     discriminant = b * b - 4.0 * c
#     if discriminant < 0:
#         return None

#     sqrt_discriminant = math.sqrt(discriminant)
#     t1 = (-b - sqrt_discriminant) / 2.0
#     t2 = (-b + sqrt_discriminant) / 2.0

#     # Find the smallest positive intersection distance
#     positive_intersections = [t for t in (t1, t2) if t > 1e-6]
#     return min(positive_intersections) if positive_intersections else None


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
    Projects a grid of image pixels onto a spherical shell (the "cloud shell") and
    returns the corresponding longitude, latitude, and a visibility mask.

    Scale vs. Mask:
        - Projection scale is determined by `edge_fov_deg`:
          the image disc of radius `radius_px` corresponds to a field-of-view
          of `edge_fov_deg` (center → rim = edge_fov_deg).
        - Visibility masking is controlled separately by `mask_fov_deg`:
          pixels are considered inside if their angular distance from the center
          is <= mask_fov_deg and altitude >= alt_min_deg.
        - Setting mask_fov_deg > edge_fov_deg allows drawing beyond the disc
          into the square corners of the image.

    Geometric limit:
        - The farthest angular distance covered by the image *corners* is
          (edge_fov_deg) * √2.
        - Thus the maximum meaningful mask_fov_deg is approximately:
              edge_fov_deg * √2
        - A ValueError is raised if mask_fov_deg exceeds this limit.

    Args:
        lat0_deg: Observer's latitude in degrees.
        lon0_deg: Observer's longitude in degrees.
        alt0_deg: Observer's viewing altitude (elevation) in degrees.
        az0_deg: Observer's viewing azimuth (bearing) in degrees.
        radius_px: Radius of the output disc in pixels (image is 2R × 2R).
        cloud_shell_km: Radius of the cloud sphere in kilometers.
        alt_min_deg: Minimum altitude to be considered visible (horizon cutoff).
        mask_fov_deg: Field-of-view angle for the visibility mask.
        edge_fov_deg: Field-of-view angle that the *disc* represents
                      (projection scaling).

    Returns:
        A tuple containing:
        - lon_grid: 2D numpy array of longitudes for each pixel (NaN if masked).
        - lat_grid: 2D numpy array of latitudes for each pixel (NaN if masked).
        - mask_inside: 2D boolean array, True if within the mask FOV and above horizon.

    Raises:
        ValueError: If mask_fov_deg > edge_fov_deg * sqrt(2).
    """
    image_size = 2 * radius_px

    # --- Validate parameters ---
    fov_limit = edge_fov_deg * math.sqrt(2.0)
    if mask_fov_deg > fov_limit + 1e-6:
        raise ValueError(f"`mask_fov_deg` ({mask_fov_deg}°) exceeds geometric limit for edge_fov_deg={edge_fov_deg}°. " f"Maximum allowed ≈ {fov_limit:.2f}°.")

    if edge_fov_deg <= 0 or edge_fov_deg > 180:
        raise ValueError("`edge_fov_deg` must be in (0, 180].")
    if mask_fov_deg <= 0:
        raise ValueError("`mask_fov_deg` must be > 0.")

    # --- Pixel coordinates centered at (0,0); y inverted for image coords ---
    pixel_coords = (np.arange(image_size) + 0.5) - radius_px
    x, y = np.meshgrid(pixel_coords, -pixel_coords)
    pixel_radius = np.hypot(x, y)

    # --- Projection scale: disc rim (r = R) → rho = edge_fov_deg/2 ---
    rho_deg = (edge_fov_deg * pixel_radius / radius_px).astype(np.float32)

    # --- Observer position and view basis (ECEF) ---
    observer_pos_ecef = geodetic_to_ecef(lat0_deg, lon0_deg)
    _, _, up_vec = enu_basis(lat0_deg, lon0_deg)
    center_view_dir = azalt_to_dir_ecef(az0_deg, alt0_deg, lat0_deg, lon0_deg)

    ty_unnormalized = up_vec - np.dot(up_vec, center_view_dir) * center_view_dir
    ty = ty_unnormalized / (np.linalg.norm(ty_unnormalized) or 1.0)
    tx = np.cross(center_view_dir, ty)

    psi = np.arctan2(y, x)  # angle on the view plane
    cos_rho, sin_rho = np.cos(np.deg2rad(rho_deg)), np.sin(np.deg2rad(rho_deg))

    # Unit vector on the view plane
    b = np.cos(psi)[..., None] * tx + np.sin(psi)[..., None] * ty
    # Direction vector for each pixel
    d = cos_rho[..., None] * center_view_dir + sin_rho[..., None] * b
    d /= np.linalg.norm(d, axis=2, keepdims=True) + 1e-12

    # --- Visibility masks ---
    alt_rad = np.arcsin(np.dot(d, up_vec))
    visible_mask = np.degrees(alt_rad) >= alt_min_deg
    fov_mask = rho_deg <= mask_fov_deg + 1e-6
    mask_inside = fov_mask & visible_mask

    # --- Intersect with the cloud shell ---
    lon_grid = np.full((image_size, image_size), np.nan, dtype=np.float32)
    lat_grid = np.full((image_size, image_size), np.nan, dtype=np.float32)

    b_quad = 2.0 * np.sum(observer_pos_ecef * d, axis=2)
    c_quad = np.dot(observer_pos_ecef, observer_pos_ecef) - cloud_shell_km * cloud_shell_km
    discriminant = b_quad * b_quad - 4.0 * c_quad
    valid = discriminant >= 0

    sqrt_disc = np.sqrt(np.maximum(discriminant, 0.0))
    t1 = (-b_quad - sqrt_disc) / 2.0
    t2 = (-b_quad + sqrt_disc) / 2.0
    t = np.where(t1 > 1e-6, t1, np.where(t2 > 1e-6, t2, np.nan))

    P = observer_pos_ecef + d * t[..., None]
    x_int, y_int, z_int = P[..., 0], P[..., 1], P[..., 2]

    lon_grid[valid] = np.degrees(np.arctan2(y_int, x_int))[valid]
    hyp = np.hypot(x_int, y_int)
    lat_grid[valid] = np.degrees(np.arctan2(z_int, hyp))[valid]

    # Apply mask
    lon_grid[~mask_inside] = np.nan
    lat_grid[~mask_inside] = np.nan

    return lon_grid, lat_grid, mask_inside
