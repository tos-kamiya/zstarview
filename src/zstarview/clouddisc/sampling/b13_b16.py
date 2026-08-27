"""Diagnostic calculations for collocated B13/B16 brightness temperatures."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BT_VALID_MIN_K = 150.0
BT_VALID_MAX_K = 350.0
DELTA_BT_MIN_K = 0.0
DELTA_BT_MAX_K = 30.0
HIGH_SCORE_DELTA_LOW_K = 4.0
HIGH_SCORE_DELTA_HIGH_K = 18.0
MAX_REDISTRIBUTION_STRENGTH = 0.30
DIAGNOSTIC_PERCENTILES = (0.0, 1.0, 5.0, 25.0, 50.0, 75.0, 95.0, 99.0, 100.0)


@dataclass(frozen=True, slots=True)
class B13B16Diagnostic:
    """Collocated diagnostic fields; this does not drive cloud rendering."""

    bt13_k: np.ndarray
    bt16_k: np.ndarray
    delta_bt_k: np.ndarray
    valid_mask: np.ndarray


def collocate_b13_b16(
    bt13_k: np.ndarray,
    bt16_k: np.ndarray,
    *,
    valid_min_k: float = BT_VALID_MIN_K,
    valid_max_k: float = BT_VALID_MAX_K,
) -> B13B16Diagnostic:
    """Validate equally shaped BT samples and calculate ``BT13 - BT16``."""
    bt13 = np.asarray(bt13_k, dtype=np.float32)
    bt16 = np.asarray(bt16_k, dtype=np.float32)
    if bt13.shape != bt16.shape:
        raise ValueError(f"B13 shape {bt13.shape} != B16 shape {bt16.shape}")
    valid = (
        np.isfinite(bt13)
        & np.isfinite(bt16)
        & (bt13 >= float(valid_min_k))
        & (bt13 <= float(valid_max_k))
        & (bt16 >= float(valid_min_k))
        & (bt16 <= float(valid_max_k))
    )
    delta = np.full(bt13.shape, np.nan, dtype=np.float32)
    np.subtract(bt13, bt16, out=delta, where=valid)
    return B13B16Diagnostic(bt13, bt16, delta, valid)


def diagnostic_summary(diagnostic: B13B16Diagnostic) -> dict[str, object]:
    """Return JSON-serializable coverage and percentile statistics."""
    total = int(diagnostic.valid_mask.size)
    valid = int(np.count_nonzero(diagnostic.valid_mask))
    summary: dict[str, object] = {
        "sample_count": total,
        "valid_count": valid,
        "valid_ratio": float(valid / total) if total else 0.0,
        "delta_definition": "BT13 - BT16",
    }
    for label, values in (
        ("bt13_k", diagnostic.bt13_k[diagnostic.valid_mask]),
        ("bt16_k", diagnostic.bt16_k[diagnostic.valid_mask]),
        ("delta_bt_k", diagnostic.delta_bt_k[diagnostic.valid_mask]),
    ):
        if values.size == 0:
            summary[label] = None
            continue
        percentiles = np.percentile(values.astype(np.float64), DIAGNOSTIC_PERCENTILES)
        summary[label] = {
            f"p{percentile:g}": float(value)
            for percentile, value in zip(
                DIAGNOSTIC_PERCENTILES, percentiles, strict=True
            )
        }
    return summary


def high_cloud_score(delta_bt_k: np.ndarray) -> np.ndarray:
    """Return a provisional upper-cloud hint; smaller delta means higher score."""
    clipped = np.clip(
        np.asarray(delta_bt_k, dtype=np.float32), DELTA_BT_MIN_K, DELTA_BT_MAX_K
    )
    span = HIGH_SCORE_DELTA_HIGH_K - HIGH_SCORE_DELTA_LOW_K
    t = np.clip((clipped - HIGH_SCORE_DELTA_LOW_K) / span, 0.0, 1.0)
    smooth = t * t * (3.0 - 2.0 * t)
    return (1.0 - smooth).astype(np.float32)
