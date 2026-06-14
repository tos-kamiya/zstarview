# -*- coding: utf-8 -*-
"""
Heuristics for estimating warm and cold brightness temperature thresholds.

These functions analyze the brightness temperature (BT) data to find suitable
"white point" (coldest temperature) and "black point" (warmest temperature)
values. This dynamic range adjustment is crucial for creating a well-contrasted
grayscale image of the clouds.
"""

from typing import Sequence, Tuple

import numpy as np
import xarray as xr

# --- Helper Functions ---


def _safe_int(x: float) -> int:
    """Safely converts a float to an integer, clipping to the 32-bit signed int range."""
    return int(np.clip(x, -2_147_000_000, 2_147_000_000))


def _percentile_ignore_nan(vals: Sequence[float], p: float) -> float:
    """Calculates the percentile of a sequence, ignoring any NaN values."""
    finite_vals = np.asarray(vals, dtype=np.float64)
    finite_vals = finite_vals[np.isfinite(finite_vals)]
    return float(np.percentile(finite_vals, p)) if finite_vals.size > 0 else np.nan


def _window_mean(da: xr.DataArray, cx: float, cy: float, half: int = 5) -> float:
    """Calculates the mean value within a small square window of a DataArray."""
    x0 = max(0, _safe_int(cx - half))
    x1 = min(da.sizes["x"], _safe_int(cx + half + 1))
    y0 = max(0, _safe_int(cy - half))
    y1 = min(da.sizes["y"], _safe_int(cy + half + 1))
    if x1 <= x0 or y1 <= y0:
        return np.nan

    block = da.isel(x=slice(x0, x1), y=slice(y0, y1)).values
    if not np.isfinite(block).any():
        return np.nan
    return float(np.nanmean(block))


# --- Estimation Functions ---


def estimate_bt_warm_from_equator_band(
    da_b13: xr.DataArray,
    lon_center_deg: float,
    delta_lon: float = 60.0,
    equator_lat: float = 0.0,
    step_deg: float = 1.0,
    half: int = 3,
    warm_p: float = 97.0,
    equator_lat_half_band_deg: float = 5.0,
) -> Tuple[float, np.ndarray]:
    """
    Estimates the 'warm' brightness temperature by sampling a band along the equator.

    This strategy provides a stable reference for the warmest temperature (clear sky
    or ground) by sampling a narrow tropical belt around the equator, which is less
    likely to have large, cold cloud systems than temperate or polar regions.

    Args:
        da_b13: The input DataArray for band 13.
        lon_center_deg: The central longitude for the sampling band.
        delta_lon: The half-width of the longitude band to sample.
        equator_lat: The center latitude for the reference belt (defaults to 0.0).
        step_deg: The longitude step for sampling.
        half: The half-width of the averaging window at each sample point.
        warm_p: The percentile to use for the final estimate.
        equator_lat_half_band_deg: Half-width of the latitude band to sample around
            `equator_lat`. The default `5.0` samples from `-5°` to `+5°`.

    Returns:
        A tuple containing the estimated warm BT and the array of samples taken.
    """
    area = da_b13.attrs.get("area")
    if area is None:
        return 310.0, np.array([], dtype=np.float32)

    lons = np.arange(lon_center_deg - delta_lon, lon_center_deg + delta_lon + 1, step_deg)
    lat_half_band = max(0.0, float(equator_lat_half_band_deg))
    lats = np.arange(
        float(equator_lat) - lat_half_band,
        float(equator_lat) + lat_half_band + float(step_deg),
        float(step_deg),
    )
    sample_arr = []
    for lat in lats:
        for lon in lons:
            try:
                x, y = area.get_xy_from_lonlat(lon, lat)
                v = _window_mean(da_b13, x, y, half=half)
                if np.isfinite(v):
                    sample_arr.append(v)
            except Exception:
                continue

    if not sample_arr:
        return 310.0, np.array(sample_arr, dtype=np.float32)

    bt_warm = _percentile_ignore_nan(sample_arr, warm_p)
    return float(np.clip(bt_warm, 180.0, 315.0)), np.array(sample_arr, dtype=np.float32)


def estimate_bt_cold_hybrid(
    bt_view: np.ndarray,
    mask_inside: np.ndarray,
    eq_samples: np.ndarray,
    bt_warm: float,
    cold_local_p: float = 5.0,
    cold_eq_p: float = 3.0,
    beta_max: float = 0.7,
    beta_min: float = 0.15,
    clear_std_thresh: float = 2.5,
    guard_range_min: float = 20.0,
    guard: Tuple[float, float] = (180.0, 315.0),
) -> float:
    """
    Estimates the 'cold' brightness temperature using a hybrid approach.

    This function blends a "local" estimate (from the current view) and a "global"
    estimate (from the equatorial band) to get a robust cold temperature value.
    The blending factor `beta` is dynamic: it gives more weight to the local view
    if the view contains significant cloud cover (high variance), and more weight
    to the global equatorial data if the view is mostly clear (low variance).

    Args:
        bt_view: The projected brightness temperature data for the current view.
        mask_inside: A boolean mask for the visible area.
        eq_samples: An array of BT samples taken from the equatorial band.
        bt_warm: The pre-estimated warm temperature threshold.
        cold_local_p: The percentile to use for the local cold estimate.
        cold_eq_p: The percentile to use for the equatorial cold estimate.
        beta_max: The maximum weight to give to the local estimate (when cloudy).
        beta_min: The minimum weight to give to the local estimate (when clear).
        clear_std_thresh: The standard deviation threshold to determine if the sky is clear.
        guard_range_min: The minimum enforced dynamic range (bt_warm - bt_cold).
        guard: A tuple (min, max) for clipping the final result.

    Returns:
        The estimated cold brightness temperature in Kelvin.
    """
    # 1. Calculate local statistics from the visible pixels in the current view.
    inside_mask = (mask_inside.astype(bool)) & np.isfinite(bt_view)
    if inside_mask.sum() < 50:
        bt_cold_local, loc_std = np.nan, np.nan  # Not enough data for local stats
    else:
        vals = bt_view[inside_mask].astype(np.float64)
        bt_cold_local = _percentile_ignore_nan(vals, cold_local_p)
        loc_std = float(np.nanstd(vals))

    # 2. Calculate the cold temperature from the global equatorial samples.
    bt_cold_eq = _percentile_ignore_nan(eq_samples.tolist(), cold_eq_p)

    # 3. Determine the blending factor `beta` based on local variance.
    #    If the local standard deviation is low, the scene is likely clear, so we
    #    rely more on the global estimate (low beta). If std is high, it's cloudy,
    #    so we rely more on the local estimate (high beta).
    if np.isfinite(loc_std):
        t = np.clip((loc_std - 1.0) / (2 * clear_std_thresh), 0.0, 1.0)
        beta = beta_min * (1.0 - t) + beta_max * t
    else:
        beta = 0.4  # Default weight if local stats are unavailable

    # 4. Perform the hybrid blend using available data.
    parts, weights = [], []
    if np.isfinite(bt_cold_eq):
        parts.append(bt_cold_eq)
        weights.append(1.0 - beta)
    if np.isfinite(bt_cold_local):
        parts.append(bt_cold_local)
        weights.append(beta)

    if not parts:
        bt_cold = 190.0  # Fallback if no data is available
    else:
        w_norm = np.array(weights, dtype=np.float64) / np.sum(weights)
        bt_cold = float(np.dot(w_norm, np.array(parts, dtype=np.float64)))

    # 5. Apply guardrails to ensure the result is physically plausible and has a minimum range.
    lo_g, hi_g = guard
    bt_cold = float(np.clip(bt_cold, lo_g, hi_g))
    if bt_warm - bt_cold < guard_range_min:
        bt_cold = bt_warm - guard_range_min

    return bt_cold
