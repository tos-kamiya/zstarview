from typing import Sequence, Tuple

import numpy as np
import xarray as xr


def _safe_int(x: float) -> int:
    """
    Safely converts a float to an integer by clipping it to the standard 32-bit signed integer range.

    Args:
        x: The float value to convert.

    Returns:
        The clipped integer value.
    """
    return int(np.clip(x, -2_147_000_000, 2_147_000_000))


def _lon_norm(lon: float) -> float:
    """
    Normalizes a longitude to the range [-180, 180).

    Args:
        lon: The longitude in degrees.

    Returns:
        The normalized longitude.
    """
    return ((lon + 180.0) % 360.0) - 180.0


def _window_mean(da: xr.DataArray, cx: float, cy: float, half: int = 5) -> float:
    """
    Calculates the mean value within a small square window of a DataArray.
    Ignores pixels outside the data array's bounds. The window size is (2*half+1) x (2*half+1).

    Args:
        da: The input DataArray.
        cx: The center x-coordinate (in pixels) of the window.
        cy: The center y-coordinate (in pixels) of the window.
        half: The half-width of the window in pixels.

    Returns:
        The mean of the values in the window, or np.nan if the window is empty.
    """
    x0 = max(0, _safe_int(cx - half))
    x1 = min(da.sizes["x"], _safe_int(cx + half + 1))
    y0 = max(0, _safe_int(cy - half))
    y1 = min(da.sizes["y"], _safe_int(cy + half + 1))
    if x1 <= x0 or y1 <= y0:
        return np.nan
    block = da.isel(x=slice(x0, x1), y=slice(y0, y1)).values
    return float(np.nanmean(block))


def _percentile_ignore_nan(vals: Sequence[float], p: float) -> float:
    """
    Calculates the percentile of a sequence of values, ignoring NaNs.

    Args:
        vals: The sequence of values.
        p: The percentile to compute (0-100).

    Returns:
        The computed percentile, or np.nan if the sequence is empty after removing NaNs.
    """
    vals = np.asarray(vals, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    return float(np.percentile(vals, p))


def estimate_bt_warm_from_equator_band(
    da_b13: xr.DataArray, lon_center_deg: float, delta_lon: float = 60.0, equator_lat: float = 0.0, step_deg: float = 1.0, half: int = 3, warm_p: float = 97.0
) -> Tuple[float, np.ndarray]:
    """
    Estimates the 'warm' brightness temperature by sampling a band along the equator
    centered on the observer's longitude [lon-Δ, lon+Δ].

    Args:
        da_b13: The input DataArray for band 13.
        lon_center_deg: The central longitude for sampling.
        delta_lon: The half-width of the longitude band to sample.
        equator_lat: The latitude of the equator.
        step_deg: The longitude step for sampling.
        half: The half-width of the averaging window at each sample point in pixels.
        warm_p: The percentile to use for the final estimate.

    Returns:
        The estimated warm brightness temperature in Kelvin.
    """
    area = da_b13.attrs.get("area", None)
    if area is None:
        return 310.0

    lons = np.arange(lon_center_deg - delta_lon, lon_center_deg + delta_lon + 1, step_deg)
    sample_arr = []
    for lon in lons:
        try:
            x, y = area.get_xy_from_lonlat(lon, equator_lat)
            # Windowed average
            v = _window_mean(da_b13, x, y, half=half)
            if np.isfinite(v):
                sample_arr.append(v)
        except Exception:
            continue

    if not sample_arr:
        return 310.0

    bt_warm = np.percentile(sample_arr, warm_p)
    return float(np.clip(bt_warm, 180.0, 315.0)), sample_arr


def _percentile_nan(a, p):
    a = np.asarray(a, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0: return np.nan
    return float(np.percentile(a, p))

def estimate_bt_cold_hybrid(
    bt_view: np.ndarray,          # 視界に投影済みBT (mask_insideと同shape)
    mask_inside: np.ndarray,      # 可視円盤マスク (bool/0-1)
    eq_samples: np.ndarray,       # 赤道帯サンプル（-60°〜+60°の小窓平均を並べた1Dなど）
    bt_warm: float,               # 先に推定した bt_warm（赤道帯の上位分位）
    cold_local_p: float = 5.0,    # 視界内の下位パーセンタイル
    cold_eq_p: float    = 3.0,    # 赤道帯側の下位パーセンタイル
    beta_max: float     = 0.7,    # 雲が多いときの視界内重みの上限
    beta_min: float     = 0.15,   # 晴天時の視界内重みの下限（赤道帯重視）
    clear_std_thresh: float = 2.5,# 晴天とみなす標準偏差しきい値[K]の中点目安
    guard_range_min: float = 20.0,# 最小ダイナミックレンジ[K]
    guard=(180.0, 315.0)
):
    # 1) 局所統計
    inside = (mask_inside.astype(bool)) & np.isfinite(bt_view)
    if inside.sum() < 50:
        # データ不足 → 赤道帯のみ
        bt_cold_local = np.nan
        loc_std = np.nan
    else:
        vals = bt_view[inside].astype(np.float64)
        bt_cold_local = _percentile_nan(vals, cold_local_p)
        loc_std = float(np.nanstd(vals))

    # 2) 赤道帯の冷側
    bt_cold_eq = _percentile_nan(eq_samples, cold_eq_p)

    # 3) βを局所分散で可変（晴れほどβ↓）
    #    loc_std ≈ 1〜2K: 晴れ, 5K以上: 雲のばらつき大 という経験則を想定
    if np.isfinite(loc_std):
        t = np.clip((loc_std - 1.0) / (2*clear_std_thresh), 0.0, 1.0)  # ざっくり [1K, ~6K]を0..1に
        beta = beta_min*(1.0 - t) + beta_max*t
    else:
        beta = 0.4  # 適当な中庸

    # 4) ハイブリッド合成（利用可能なものだけで）
    parts = []
    weights = []
    if np.isfinite(bt_cold_eq):
        parts.append(bt_cold_eq); weights.append(1.0 - beta)
    if np.isfinite(bt_cold_local):
        parts.append(bt_cold_local); weights.append(beta)

    if not parts:
        bt_cold = 190.0  # フォールバック
    else:
        w = np.array(weights, dtype=np.float64)
        w = w / w.sum()
        bt_cold = float(np.dot(w, np.array(parts, dtype=np.float64)))

    # 5) 物理ガード & 最小幅
    lo_g, hi_g = guard
    bt_cold = float(np.clip(bt_cold, lo_g, hi_g))
    if bt_warm - bt_cold < guard_range_min:
        bt_cold = bt_warm - guard_range_min

    return bt_cold
