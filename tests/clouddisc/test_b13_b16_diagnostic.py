from __future__ import annotations

import numpy as np
import pytest

from zstarview.clouddisc.sampling.b13_b16 import (
    collocate_b13_b16,
    diagnostic_summary,
)


def test_collocate_b13_b16_uses_bt13_minus_bt16() -> None:
    result = collocate_b13_b16(
        np.array([[280.0, 220.0]], dtype=np.float32),
        np.array([[260.0, 219.0]], dtype=np.float32),
    )
    np.testing.assert_allclose(result.delta_bt_k, [[20.0, 1.0]])
    assert result.valid_mask.all()


def test_collocate_b13_b16_masks_missing_and_extreme_values() -> None:
    result = collocate_b13_b16(
        np.array([[280.0, np.nan, 500.0]], dtype=np.float32),
        np.array([[260.0, 250.0, 250.0]], dtype=np.float32),
    )
    np.testing.assert_array_equal(result.valid_mask, [[True, False, False]])
    assert np.isnan(result.delta_bt_k[0, 1:]).all()
    summary = diagnostic_summary(result)
    assert summary["valid_count"] == 1
    assert summary["valid_ratio"] == pytest.approx(1.0 / 3.0)


def test_collocate_b13_b16_rejects_shape_mismatch() -> None:
    with pytest.raises(ValueError, match="B13 shape"):
        collocate_b13_b16(np.zeros((2, 2)), np.zeros((2, 3)))
