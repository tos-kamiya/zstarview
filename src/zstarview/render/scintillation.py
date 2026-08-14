"""Small, display-only scintillation helpers for the scenic viewer."""

from __future__ import annotations

import math
import random

import numpy as np

from ..astro import is_in_fov_vectorized
from ..types import StarsTable

SCINTILLATION_MAX_DISTANCE_DEG = 3.0
SCINTILLATION_TARGET_COUNT = 15
SCINTILLATION_BRIGHT_STAR_EXCLUSIVE_VMAG = 4.0
SCINTILLATION_BASE_ALPHA = 0.5
SCINTILLATION_MIN_ALPHA = 0.0
SCINTILLATION_DOT_SIZE_MULTIPLIER = 1.0


def scintillation_strength(alt_deg: float) -> float:
    """Return the linear altitude strength above the horizon."""
    altitude = float(alt_deg)
    if not math.isfinite(altitude) or altitude < 0.0:
        return 0.0
    effective_altitude = min(max(altitude, 10.0), 90.0)
    alpha_range = SCINTILLATION_BASE_ALPHA - SCINTILLATION_MIN_ALPHA
    return 1.0 - alpha_range * (effective_altitude - 10.0) / (
        SCINTILLATION_BASE_ALPHA * 80.0
    )


def scintillation_alpha(alt_deg: float) -> float:
    """Return the black-mask alpha for a star at the given altitude."""
    if float(alt_deg) < 0.0:
        return 0.0
    return float(
        np.clip(SCINTILLATION_BASE_ALPHA * scintillation_strength(alt_deg), 0.0, 1.0)
    )


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
    *,
    rng: random.Random,
) -> tuple[float, float]:
    """Sample an alt/az direction using the requested altitude distribution."""
    altitude = 10.0 + 80.0 * rng.random() ** 2.5
    azimuth = rng.random() * 360.0
    return altitude, azimuth
