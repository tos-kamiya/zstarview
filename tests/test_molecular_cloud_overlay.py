from __future__ import annotations

import numpy as np

from zstarview.render.molecular_cloud_overlay import _apply_molecular_cloud_value_knee


def test_apply_molecular_cloud_value_knee_preserves_hue_and_saturation() -> None:
    rgb = np.array(
        [[0.2, 0.8, 0.4], [1.0, 0.0, 0.7]],
        dtype=np.float32,
    )

    result = _apply_molecular_cloud_value_knee(rgb)

    np.testing.assert_allclose(result, [[1.0 / 9.0, 4.0 / 9.0, 2.0 / 9.0], [0.5, 0.0, 0.35]])
    np.testing.assert_allclose(rgb, [[0.2, 0.8, 0.4], [1.0, 0.0, 0.7]])
