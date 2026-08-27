from __future__ import annotations

import numpy as np
import pytest

from zstarview.clouddisc.sampling.b13_b16 import (
    MAX_REDISTRIBUTION_STRENGTH,
    collocate_b13_b16,
    diagnostic_summary,
    high_cloud_score,
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


def test_high_cloud_score_clips_delta_and_uses_provisional_strength() -> None:
    score = high_cloud_score(np.array([-5.0, 4.0, 18.0, 40.0]))
    assert score[0] == pytest.approx(1.0)
    assert score[1] == pytest.approx(1.0)
    assert score[2] == pytest.approx(0.0)
    assert score[3] == pytest.approx(0.0)
    assert MAX_REDISTRIBUTION_STRENGTH == pytest.approx(0.30)
