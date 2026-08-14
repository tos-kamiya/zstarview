"""Small, display-only twinkle helpers for the scenic viewer."""

from __future__ import annotations

import math
import random

import numpy as np

from ..astro import is_in_fov_vectorized
from ..types import StarsTable

TWINKLE_MAX_DISTANCE_DEG = 3.0
TWINKLE_TARGET_COUNT = 30
TWINKLE_BRIGHT_STAR_EXCLUSIVE_VMAG = 4.0
TWINKLE_BASE_ALPHA = 0.5
TWINKLE_MIN_ALPHA = 0.0
TWINKLE_DOT_SIZE_MULTIPLIER = 1.0


def twinkle_strength(alt_deg: float) -> float:
    """Return the linear altitude strength above the horizon."""
    altitude = float(alt_deg)
    if not math.isfinite(altitude) or altitude < 0.0:
        return 0.0
    effective_altitude = min(max(altitude, 10.0), 90.0)
    alpha_range = TWINKLE_BASE_ALPHA - TWINKLE_MIN_ALPHA
    return 1.0 - alpha_range * (effective_altitude - 10.0) / (
        TWINKLE_BASE_ALPHA * 80.0
    )


def twinkle_alpha(alt_deg: float) -> float:
    """Return the black-mask alpha for a star at the given altitude."""
    if float(alt_deg) < 0.0:
        return 0.0
    return float(
        np.clip(TWINKLE_BASE_ALPHA * twinkle_strength(alt_deg), 0.0, 1.0)
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


def nearest_twinkle_star_index(
    stars: StarsTable,
    *,
    target_alt_deg: float,
    target_az_deg: float,
    view_center: tuple[float, float],
    content_fov_deg: float,
    vmag_limit: float,
    max_distance_deg: float = TWINKLE_MAX_DISTANCE_DEG,
) -> int | None:
    """Find the nearest eligible faint star within the angular search radius."""
    rows = nearest_twinkle_star_rows(
        stars,
        target_alt_deg=np.asarray([target_alt_deg], dtype=float),
        target_az_deg=np.asarray([target_az_deg], dtype=float),
        view_center=view_center,
        content_fov_deg=content_fov_deg,
        vmag_limit=vmag_limit,
        max_distance_deg=max_distance_deg,
    )
    row = int(rows[0])
    if row < 0:
        return None
    return int(stars["star_index"][row])


def nearest_twinkle_star_rows(
    stars: StarsTable,
    *,
    target_alt_deg: np.ndarray,
    target_az_deg: np.ndarray,
    view_center: tuple[float, float],
    content_fov_deg: float,
    vmag_limit: float,
    max_distance_deg: float = TWINKLE_MAX_DISTANCE_DEG,
) -> np.ndarray:
    """Return source rows nearest to all target directions, or -1 when absent."""
    target_alt = np.asarray(target_alt_deg, dtype=float).reshape(-1)
    target_az = np.asarray(target_az_deg, dtype=float).reshape(-1)
    if target_alt.shape != target_az.shape:
        raise ValueError("target altitude and azimuth arrays must have equal shape")
    result = np.full(target_alt.shape, -1, dtype=np.int32)
    if stars["alt"].size == 0 or target_alt.size == 0:
        return result

    candidate_mask = (
        (stars["alt"] > 0.0)
        & (stars["vmag"] > TWINKLE_BRIGHT_STAR_EXCLUSIVE_VMAG)
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
        return result

    candidate_alt_rad = np.radians(stars["alt"][candidate_indices])
    candidate_az_rad = np.radians(stars["az"][candidate_indices])
    candidate_cos_alt = np.cos(candidate_alt_rad)
    candidate_vectors = np.column_stack(
        (
            candidate_cos_alt * np.cos(candidate_az_rad),
            candidate_cos_alt * np.sin(candidate_az_rad),
            np.sin(candidate_alt_rad),
        )
    )

    valid_targets = np.isfinite(target_alt) & np.isfinite(target_az) & (target_alt > 0.0)
    if not np.any(valid_targets):
        return result
    valid_alt_rad = np.radians(target_alt[valid_targets])
    valid_az_rad = np.radians(target_az[valid_targets])
    valid_cos_alt = np.cos(valid_alt_rad)
    target_vectors = np.column_stack(
        (
            valid_cos_alt * np.cos(valid_az_rad),
            valid_cos_alt * np.sin(valid_az_rad),
            np.sin(valid_alt_rad),
        )
    )
    cosine = target_vectors @ candidate_vectors.T
    nearest_candidate = np.argmax(cosine, axis=1)
    nearest_cosine = cosine[np.arange(nearest_candidate.size), nearest_candidate]
    within_radius = nearest_cosine >= math.cos(math.radians(float(max_distance_deg)))
    valid_result = np.full(nearest_candidate.shape, -1, dtype=np.int32)
    valid_result[within_radius] = candidate_indices[nearest_candidate[within_radius]]
    result[valid_targets] = valid_result
    return result


def sample_twinkle_direction(
    *,
    rng: random.Random,
) -> tuple[float, float]:
    """Sample an alt/az direction using the requested altitude distribution."""
    altitude = 10.0 + 80.0 * rng.random() ** 2.5
    azimuth = rng.random() * 360.0
    return altitude, azimuth
