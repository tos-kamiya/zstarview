"""Small, display-only scintillation helpers for the scenic viewer."""

from __future__ import annotations

import math
import random

import numpy as np

from ..astro import is_in_fov_vectorized
from ..types import StarsTable

SCINTILLATION_MAX_DISTANCE_DEG = 3.0
SCINTILLATION_TARGET_COUNT = 10
SCINTILLATION_BRIGHT_STAR_EXCLUSIVE_VMAG = 4.0
SCINTILLATION_BASE_ALPHA = 0.6
SCINTILLATION_DOT_SIZE_MULTIPLIER = 1.0


def scintillation_strength(alt_deg: float) -> float:
    """Return the normalized altitude strength, or zero below the horizon."""
    altitude = float(alt_deg)
    if not math.isfinite(altitude) or altitude <= 0.0:
        return 0.0
    effective_altitude = max(altitude, 10.0)
    return math.sin(math.radians(10.0)) / math.sin(math.radians(effective_altitude))


def scintillation_alpha(alt_deg: float) -> float:
    """Return the black-mask alpha for a star at the given altitude."""
    return float(np.clip(SCINTILLATION_BASE_ALPHA * scintillation_strength(alt_deg), 0.0, 1.0))


def spherical_distance_deg(
    alt_deg: float,
    az_deg: float,
    target_alt_deg: float,
    target_az_deg: float,
) -> float:
    """Return angular separation between two alt/az directions in degrees."""
    alt1 = math.radians(float(alt_deg))
    alt2 = math.radians(float(target_alt_deg))
    delta_az = math.radians(float(az_deg) - float(target_az_deg))
    cosine = (
        math.sin(alt1) * math.sin(alt2)
        + math.cos(alt1) * math.cos(alt2) * math.cos(delta_az)
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def nearest_scintillation_star_index(
    stars: StarsTable,
    *,
    target_alt_deg: float,
    target_az_deg: float,
    view_center: tuple[float, float],
    content_fov_deg: float,
    vmag_limit: float,
    max_distance_deg: float = SCINTILLATION_MAX_DISTANCE_DEG,
) -> int | None:
    """Find the nearest eligible faint star within the angular search radius."""
    if stars["alt"].size == 0 or float(target_alt_deg) <= 0.0:
        return None
    candidate_mask = (
        (stars["alt"] > 0.0)
        & (stars["vmag"] > SCINTILLATION_BRIGHT_STAR_EXCLUSIVE_VMAG)
        & (stars["vmag"] <= float(vmag_limit))
        & is_in_fov_vectorized(
            stars["alt"],
            stars["az"],
            view_center,
            fov_deg=float(content_fov_deg),
        )
    )
    candidate_indices = np.flatnonzero(candidate_mask)
    if candidate_indices.size == 0:
        return None
    alt = stars["alt"][candidate_indices]
    az = stars["az"][candidate_indices]
    alt1 = np.radians(alt)
    alt2 = math.radians(float(target_alt_deg))
    delta_az = np.radians(az - float(target_az_deg))
    cosine = (
        np.sin(alt1) * math.sin(alt2)
        + np.cos(alt1) * math.cos(alt2) * np.cos(delta_az)
    )
    distances = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    nearest = int(np.argmin(distances))
    if float(distances[nearest]) > float(max_distance_deg):
        return None
    return int(stars["star_index"][candidate_indices[nearest]])


def sample_scintillation_direction(
    view_center: tuple[float, float],
    content_fov_deg: float,
    *,
    rng: random.Random,
) -> tuple[float, float]:
    """Sample an approximate uniformly distributed direction in the view."""
    center_alt, center_az = map(float, view_center)
    radius = math.radians(max(0.0, float(content_fov_deg)))
    # Uniform sampling over a spherical cap: cos(rho) is uniform in the cap's
    # solid angle, unlike a planar sqrt-area approximation at wide FOVs.
    rho = math.acos(1.0 - rng.random() * (1.0 - math.cos(radius)))
    bearing = rng.random() * math.tau
    center = np.array(
        [
            math.cos(math.radians(center_alt)) * math.cos(math.radians(center_az)),
            math.cos(math.radians(center_alt)) * math.sin(math.radians(center_az)),
            math.sin(math.radians(center_alt)),
        ]
    )
    east = np.array([-math.sin(math.radians(center_az)), math.cos(math.radians(center_az)), 0.0])
    up = np.array(
        [
            -math.sin(math.radians(center_alt)) * math.cos(math.radians(center_az)),
            -math.sin(math.radians(center_alt)) * math.sin(math.radians(center_az)),
            math.cos(math.radians(center_alt)),
        ]
    )
    direction = (
        math.cos(rho) * center
        + math.sin(rho) * (math.cos(bearing) * east + math.sin(bearing) * up)
    )
    altitude = math.degrees(math.asin(float(np.clip(direction[2], -1.0, 1.0))))
    azimuth = math.degrees(math.atan2(float(direction[1]), float(direction[0]))) % 360.0
    return altitude, azimuth
