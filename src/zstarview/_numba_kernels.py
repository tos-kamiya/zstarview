"""Small numerical kernels used by the rendering/data-preparation paths."""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def accumulate_glow_strengths_numba(
    source_azimuths: np.ndarray,
    source_altitudes: np.ndarray,
    source_strengths: np.ndarray,
    target_azimuths: np.ndarray,
    target_altitudes: np.ndarray,
    az_lut: np.ndarray,
    alt_lut: np.ndarray,
    lut_step_deg: float,
    max_az_delta_deg: float,
) -> np.ndarray:
    result = np.zeros(target_azimuths.shape[0], dtype=np.float64)
    for source_index in range(source_strengths.shape[0]):
        source_az = source_azimuths[source_index]
        source_alt = source_altitudes[source_index]
        strength = source_strengths[source_index]
        for target_index in range(target_azimuths.shape[0]):
            delta_az = abs((source_az - target_azimuths[target_index] + 180.0) % 360.0 - 180.0)
            if delta_az > max_az_delta_deg:
                continue
            az_lut_index = int(np.rint(delta_az / lut_step_deg))
            if az_lut_index >= az_lut.shape[0]:
                az_lut_index = az_lut.shape[0] - 1
            delta_alt = abs(source_alt - target_altitudes[target_index])
            alt_lut_index = int(np.rint(delta_alt / lut_step_deg))
            if alt_lut_index >= alt_lut.shape[0]:
                alt_lut_index = alt_lut.shape[0] - 1
            result[target_index] += strength * az_lut[az_lut_index] * alt_lut[alt_lut_index]
    return result


@njit(cache=True)
def accumulate_glow_field_numba(
    source_azimuths: np.ndarray,
    source_altitudes: np.ndarray,
    source_strengths: np.ndarray,
    target_azimuths: np.ndarray,
    target_altitudes: np.ndarray,
    az_lut: np.ndarray,
    alt_lut: np.ndarray,
    lut_step_deg: float,
    max_az_delta_deg: float,
) -> np.ndarray:
    result = np.zeros((target_altitudes.shape[0], target_azimuths.shape[0]), dtype=np.float64)
    for source_index in range(source_strengths.shape[0]):
        source_az = source_azimuths[source_index]
        source_alt = source_altitudes[source_index]
        strength = source_strengths[source_index]
        if strength < 0.0:
            strength = 0.0
        for az_index in range(target_azimuths.shape[0]):
            delta_az = abs((source_az - target_azimuths[az_index] + 180.0) % 360.0 - 180.0)
            if delta_az > max_az_delta_deg:
                continue
            az_lut_index = int(np.rint(delta_az / lut_step_deg))
            if az_lut_index >= az_lut.shape[0]:
                az_lut_index = az_lut.shape[0] - 1
            az_weight = az_lut[az_lut_index]
            for alt_index in range(target_altitudes.shape[0]):
                alt_lut_index = int(np.rint(abs(source_alt - target_altitudes[alt_index]) / lut_step_deg))
                if alt_lut_index >= alt_lut.shape[0]:
                    alt_lut_index = alt_lut.shape[0] - 1
                result[alt_index, az_index] += strength * az_weight * alt_lut[alt_lut_index]
    return result


@njit(cache=True)
def accumulate_glow_strengths_weighted_numba(
    source_altitudes: np.ndarray,
    source_strengths: np.ndarray,
    target_altitudes: np.ndarray,
    azimuth_weights: np.ndarray,
    alt_lut: np.ndarray,
    lut_step_deg: float,
) -> np.ndarray:
    result = np.zeros(target_altitudes.shape[0], dtype=np.float64)
    for source_index in range(source_strengths.shape[0]):
        for target_index in range(target_altitudes.shape[0]):
            alt_lut_index = int(np.rint(abs(source_altitudes[source_index] - target_altitudes[target_index]) / lut_step_deg))
            if alt_lut_index >= alt_lut.shape[0]:
                alt_lut_index = alt_lut.shape[0] - 1
            result[target_index] += source_strengths[source_index] * azimuth_weights[source_index, target_index] * alt_lut[alt_lut_index]
    return result


@njit(cache=True)
def accumulate_glow_field_weighted_numba(
    source_altitudes: np.ndarray,
    source_strengths: np.ndarray,
    target_azimuths: np.ndarray,
    target_altitudes: np.ndarray,
    azimuth_weights: np.ndarray,
    alt_lut: np.ndarray,
    lut_step_deg: float,
) -> np.ndarray:
    result = np.zeros((target_altitudes.shape[0], target_azimuths.shape[0]), dtype=np.float64)
    for source_index in range(source_strengths.shape[0]):
        strength = max(0.0, source_strengths[source_index])
        for az_index in range(target_azimuths.shape[0]):
            az_weight = azimuth_weights[source_index, az_index]
            for alt_index in range(target_altitudes.shape[0]):
                alt_lut_index = int(np.rint(abs(source_altitudes[source_index] - target_altitudes[alt_index]) / lut_step_deg))
                if alt_lut_index >= alt_lut.shape[0]:
                    alt_lut_index = alt_lut.shape[0] - 1
                result[alt_index, az_index] += strength * az_weight * alt_lut[alt_lut_index]
    return result


@njit(cache=True)
def scatter_cloud_samples_numba(
    flat_indices: np.ndarray,
    finite_mask: np.ndarray,
    valid_amount_mask: np.ndarray,
    shell_amounts: np.ndarray,
    amount_values: np.ndarray,
    delta_values: np.ndarray,
    sample_count: np.ndarray,
    delta_sum: np.ndarray,
    delta_count: np.ndarray,
) -> None:
    for index in range(flat_indices.shape[0]):
        flat_index = flat_indices[index]
        if finite_mask[index]:
            sample_count[flat_index] += 1
        if valid_amount_mask[index] and amount_values[index] > shell_amounts[flat_index]:
            shell_amounts[flat_index] = amount_values[index]
        if finite_mask[index] and np.isfinite(delta_values[index]):
            delta_sum[flat_index] += delta_values[index]
            delta_count[flat_index] += 1
