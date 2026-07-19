"""Small, display-only approximations for short-term star motion."""

from __future__ import annotations

import math

import numpy as np


_SIDEREAL_ROTATION_DEG_PER_SECOND = 360.0 / 86164.0905


def _screen_to_direction(
    x: np.ndarray,
    y: np.ndarray,
    *,
    width_px: float,
    height_px: float,
    geometry_center: tuple[float, float],
    geometry_radius: float,
    view_center_altaz_deg: tuple[float, float],
    edge_fov_deg: float,
) -> np.ndarray:
    center_alt_deg, center_az_deg = view_center_altaz_deg
    center_alt = math.radians(float(center_alt_deg))
    center_az = math.radians(float(center_az_deg))
    center_x, center_y = geometry_center
    radius = max(1.0, float(geometry_radius))
    nx = (x - float(center_x)) / radius
    ny = (y - float(center_y)) / radius
    rho = np.hypot(nx, ny)
    theta = np.radians(float(edge_fov_deg)) * rho

    center = np.array(
        [
            math.cos(center_alt) * math.cos(center_az),
            math.cos(center_alt) * math.sin(center_az),
            math.sin(center_alt),
        ],
        dtype=float,
    )
    east = np.array([-math.sin(center_az), math.cos(center_az), 0.0], dtype=float)
    up = np.array(
        [
            -math.sin(center_alt) * math.cos(center_az),
            -math.sin(center_alt) * math.sin(center_az),
            math.cos(center_alt),
        ],
        dtype=float,
    )

    safe_rho = np.where(rho > 1.0e-12, rho, 1.0)
    tangent = (
        nx[:, None] / safe_rho[:, None] * east[None, :]
        - ny[:, None] / safe_rho[:, None] * up[None, :]
    )
    tangent[rho <= 1.0e-12] = 0.0
    return np.cos(theta)[:, None] * center[None, :] + np.sin(theta)[:, None] * tangent


def _rotate_about_axis(vectors: np.ndarray, axis: np.ndarray, angle_rad: float) -> np.ndarray:
    axis = axis / max(1.0e-12, float(np.linalg.norm(axis)))
    cross = np.cross(axis[None, :], vectors)
    dot = vectors @ axis
    cosine = math.cos(angle_rad)
    sine = math.sin(angle_rad)
    return (
        vectors * cosine
        + cross * sine
        + dot[:, None] * axis[None, :] * (1.0 - cosine)
    )


def _direction_to_screen(
    vectors: np.ndarray,
    *,
    width_px: float,
    height_px: float,
    geometry_center: tuple[float, float],
    geometry_radius: float,
    view_center_altaz_deg: tuple[float, float],
    edge_fov_deg: float,
) -> np.ndarray:
    center_alt_deg, center_az_deg = view_center_altaz_deg
    center_alt = math.radians(float(center_alt_deg))
    center_az = math.radians(float(center_az_deg))
    center_x, center_y = geometry_center
    radius = max(1.0, float(geometry_radius))
    center = np.array(
        [
            math.cos(center_alt) * math.cos(center_az),
            math.cos(center_alt) * math.sin(center_az),
            math.sin(center_alt),
        ],
        dtype=float,
    )
    east = np.array([-math.sin(center_az), math.cos(center_az), 0.0], dtype=float)
    up = np.array(
        [
            -math.sin(center_alt) * math.cos(center_az),
            -math.sin(center_alt) * math.sin(center_az),
            math.cos(center_alt),
        ],
        dtype=float,
    )
    cosine = np.clip(vectors @ center, -1.0, 1.0)
    theta = np.arccos(cosine)
    east_component = vectors @ east
    up_component = vectors @ up
    tangent_length = np.hypot(east_component, up_component)
    safe_length = np.where(tangent_length > 1.0e-12, tangent_length, 1.0)
    radial_scale = theta / math.radians(max(1.0e-6, float(edge_fov_deg)))
    nx = radial_scale * east_component / safe_length
    ny = -radial_scale * up_component / safe_length
    nx[tangent_length <= 1.0e-12] = 0.0
    ny[tangent_length <= 1.0e-12] = 0.0
    return np.column_stack(
        [
            float(center_x) + nx * radius,
            float(center_y) + ny * radius,
        ]
    )


def _fit_homography(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    rows: list[list[float]] = []
    for (x, y), (u, v) in zip(source, target):
        rows.append([-x, -y, -1.0, 0.0, 0.0, 0.0, x * u, y * u, u])
        rows.append([0.0, 0.0, 0.0, -x, -y, -1.0, x * v, y * v, v])
    _, _, vh = np.linalg.svd(np.asarray(rows, dtype=float))
    matrix = vh[-1].reshape(3, 3)
    scale = matrix[2, 2]
    if abs(float(scale)) <= 1.0e-12:
        return np.eye(3, dtype=float)
    return matrix / scale


def apply_homography(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    homogeneous = np.column_stack([points, np.ones(len(points), dtype=float)])
    mapped = homogeneous @ matrix.T
    denominator = mapped[:, 2]
    safe = np.where(np.abs(denominator) > 1.0e-12, denominator, 1.0)
    return mapped[:, :2] / safe[:, None]


def build_star_interpolation_homography(
    *,
    width_px: int,
    height_px: int,
    geometry_center: tuple[float, float],
    geometry_radius: float,
    view_center_altaz_deg: tuple[float, float],
    observer_lat_deg: float,
    edge_fov_deg: float,
    elapsed_seconds: float,
    sample_margin: float = 0.04,
) -> np.ndarray:
    """Approximate short-term sidereal motion in screen coordinates.

    The returned matrix maps coordinates from the snapshot-time projected
    surface to the current projected surface. It is intentionally a display
    approximation; the expensive celestial snapshot remains unchanged.
    """
    width = max(1.0, float(width_px))
    height = max(1.0, float(height_px))
    elapsed = float(elapsed_seconds)
    if abs(elapsed) <= 1.0e-9:
        return np.eye(3, dtype=float)

    margin = min(0.25, max(0.0, float(sample_margin)))
    xs = np.array([margin, 0.5, 1.0 - margin], dtype=float) * width
    ys = np.array([margin, 0.5, 1.0 - margin], dtype=float) * height
    grid_x, grid_y = np.meshgrid(xs, ys)
    source = np.column_stack([grid_x.ravel(), grid_y.ravel()])

    directions = _screen_to_direction(
        source[:, 0],
        source[:, 1],
        width_px=width,
        height_px=height,
        geometry_center=geometry_center,
        geometry_radius=geometry_radius,
        view_center_altaz_deg=view_center_altaz_deg,
        edge_fov_deg=edge_fov_deg,
    )
    latitude = math.radians(float(observer_lat_deg))
    pole_axis = np.array([math.cos(latitude), 0.0, math.sin(latitude)], dtype=float)
    angle = math.radians(_SIDEREAL_ROTATION_DEG_PER_SECOND * elapsed)
    rotated = _rotate_about_axis(directions, pole_axis, angle)
    target = _direction_to_screen(
        rotated,
        width_px=width,
        height_px=height,
        geometry_center=geometry_center,
        geometry_radius=geometry_radius,
        view_center_altaz_deg=view_center_altaz_deg,
        edge_fov_deg=edge_fov_deg,
    )
    return _fit_homography(source, target)
