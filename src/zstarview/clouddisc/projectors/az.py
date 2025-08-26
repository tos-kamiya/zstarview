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


def rad2deg(radians: float) -> float:
    """Converts radians to degrees."""
    return radians * 180.0 / math.pi


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

    Args:
        lat_deg: Latitude in degrees.
        lon_deg: Longitude in degrees.
        r_km: The radius of the Earth in kilometers.

    Returns:
        The ECEF coordinates [x, y, z] as a numpy array.
    """
    lat, lon = deg2rad(lat_deg), deg2rad(lon_deg)
    cos_lat, sin_lat = math.cos(lat), math.sin(lat)
    cos_lon, sin_lon = math.cos(lon), math.sin(lon)

    x = r_km * cos_lat * cos_lon
    y = r_km * cos_lat * sin_lon
    z = r_km * sin_lat

    return np.array([x, y, z], dtype=float)


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


def line_sphere_intersection(origin: np.ndarray, direction: np.ndarray, radius_km: float) -> Optional[float]:
    """
    Finds the intersection of a line (ray) and a sphere.

    Solves the equation |origin + t * direction| = radius_km for t.

    Args:
        origin: The origin of the line (e.g., observer's ECEF position).
        direction: The normalized direction vector of the line.
        radius_km: The radius of the sphere.

    Returns:
        The smallest positive distance (t) to an intersection, or None if no intersection.
    """
    # Quadratic equation coefficients for t: a*t^2 + b*t + c = 0
    # Since direction is normalized, a = 1.
    b = 2.0 * np.dot(origin, direction)
    c = np.dot(origin, origin) - radius_km * radius_km

    discriminant = b * b - 4.0 * c
    if discriminant < 0:
        return None

    sqrt_discriminant = math.sqrt(discriminant)
    t1 = (-b - sqrt_discriminant) / 2.0
    t2 = (-b + sqrt_discriminant) / 2.0

    # Find the smallest positive intersection distance
    positive_intersections = [t for t in (t1, t2) if t > 1e-6]
    return min(positive_intersections) if positive_intersections else None


def az_project_lonlat_grid(
    lat0_deg: float, lon0_deg: float, alt0_deg: float, az0_deg: float, radius_px: int, cloud_shell_km: float, alt_min_deg: float = -2.0, fov_deg: float = 90.0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Projects a grid of pixels onto a spherical shell (the "cloud shell") and
    returns the corresponding longitude and latitude for each pixel.

    This function simulates a camera view from a specific location on Earth,
    looking in a specific direction (azimuth, altitude).

    Args:
        lat0_deg: Observer's latitude in degrees.
        lon0_deg: Observer's longitude in degrees.
        alt0_deg: Observer's viewing altitude (vertical angle) in degrees.
        az0_deg: Observer's viewing azimuth (horizontal angle) in degrees.
        radius_px: The radius of the output image disc in pixels.
        cloud_shell_km: The radius of the cloud sphere in kilometers.
        alt_min_deg: The minimum altitude to be considered visible (horizon).
        fov_deg: The field of view, determining the mapping from pixels to angular distance.

    Returns:
        A tuple containing:
        - lon_grid: A 2D numpy array of longitudes for each pixel.
        - lat_grid: A 2D numpy array of latitudes for each pixel.
        - mask_inside: A 2D boolean array, True for pixels within the FOV and above the horizon.
    """
    image_size = 2 * radius_px

    # Create a grid of pixel coordinates (X, Y) centered at (0,0)
    pixel_coords = (np.arange(image_size) + 0.5) - radius_px
    x, y = np.meshgrid(pixel_coords, -pixel_coords)  # Y is inverted for image coordinates
    pixel_radius = np.hypot(x, y)

    # Convert pixel radius to angular distance (rho) from the center of view
    # This maps the image radius to half the field of view.
    rho_deg = (90.0 / 2.0 * pixel_radius / radius_px).astype(np.float32)

    # --- Calculate direction vectors for each pixel ---

    # Observer's position and orientation in ECEF
    observer_pos_ecef = geodetic_to_ecef(lat0_deg, lon0_deg)
    _, _, up_vec = enu_basis(lat0_deg, lon0_deg)
    center_view_dir = azalt_to_dir_ecef(az0_deg, alt0_deg, lat0_deg, lon0_deg)

    # Create an orthonormal basis for the view plane
    # tY is the "up" direction in the view plane, tX is the "right"
    ty_unnormalized = up_vec - np.dot(up_vec, center_view_dir) * center_view_dir
    ty = ty_unnormalized / (np.linalg.norm(ty_unnormalized) or 1.0)
    tx = np.cross(center_view_dir, ty)

    # Calculate the direction vector `d` for each pixel using spherical coordinates
    # on the view plane.
    psi = np.arctan2(y, x)  # Angle on the view plane
    cos_rho, sin_rho = np.cos(np.deg2rad(rho_deg)), np.sin(np.deg2rad(rho_deg))

    # `b` is the direction on the unit circle of the view plane
    b = np.cos(psi)[..., None] * tx + np.sin(psi)[..., None] * ty
    # `d` is the final direction vector for each pixel in ECEF
    d = cos_rho[..., None] * center_view_dir + sin_rho[..., None] * b
    d /= np.linalg.norm(d, axis=2, keepdims=True) + 1e-12

    # --- Masking and Intersection ---

    # Create a mask for pixels below the horizon
    alt_rad = np.arcsin(np.dot(d, up_vec))
    visible_mask = np.degrees(alt_rad) >= alt_min_deg

    # Create a mask for pixels within the specified field of view
    fov_mask = rho_deg <= fov_deg / 2.0 + 1e-6
    mask_inside = fov_mask & visible_mask

    # --- Find intersection with the cloud shell ---

    lon_grid = np.full((image_size, image_size), np.nan, dtype=np.float32)
    lat_grid = np.full((image_size, image_size), np.nan, dtype=np.float32)

    # Vectorized calculation of line-sphere intersection for all pixels
    b_quad = 2.0 * np.sum(observer_pos_ecef * d, axis=2)
    c_quad = np.dot(observer_pos_ecef, observer_pos_ecef) - cloud_shell_km * cloud_shell_km
    discriminant = b_quad * b_quad - 4.0 * c_quad

    valid_intersection = discriminant >= 0

    sqrt_discriminant = np.sqrt(np.maximum(discriminant, 0.0))
    t1 = (-b_quad - sqrt_discriminant) / 2.0
    t2 = (-b_quad + sqrt_discriminant) / 2.0

    # Choose the smallest positive intersection distance
    t = np.where(t1 > 1e-6, t1, np.where(t2 > 1e-6, t2, np.nan))

    # Calculate intersection points P = O + t*d
    intersection_points = observer_pos_ecef + d * t[..., None]

    # Convert intersection points back to geodetic coordinates
    x_intersect, y_intersect, z_intersect = intersection_points[..., 0], intersection_points[..., 1], intersection_points[..., 2]
    lon_grid[valid_intersection] = np.degrees(np.arctan2(y_intersect, x_intersect))[valid_intersection]
    hyp = np.hypot(x_intersect, y_intersect)
    lat_grid[valid_intersection] = np.degrees(np.arctan2(z_intersect, hyp))[valid_intersection]

    # Apply visibility mask
    lon_grid[~mask_inside] = np.nan
    lat_grid[~mask_inside] = np.nan

    return lon_grid, lat_grid, mask_inside
