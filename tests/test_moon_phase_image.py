from __future__ import annotations

import numpy as np

from zstarview.utils.image import generate_moon_phase_rgba


def test_generate_moon_phase_rgba_returns_rgba_buffer() -> None:
    rgba = generate_moon_phase_rgba(
        16,
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
    )

    assert rgba.shape == (16, 16, 4)
    assert rgba.dtype == np.uint8
    assert int(rgba[0, 0, 3]) == 0
    assert int(rgba[8, 8, 3]) > 0


def test_generate_moon_phase_rgba_applies_tint() -> None:
    untinted = generate_moon_phase_rgba(
        24,
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
    )
    tinted = generate_moon_phase_rgba(
        24,
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
        tint_color=(255, 64, 64, 128),
    )

    center = (12, 12)
    assert int(tinted[center][0]) >= int(untinted[center][0])
    assert int(tinted[center][1]) < int(untinted[center][1])
    assert int(tinted[center][2]) < int(untinted[center][2])
