from __future__ import annotations

import numpy as np

from zstarview.utils.image import (
    generate_flat_moon_phase_rgba,
    generate_moon_phase_rgba,
)


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


def test_generate_flat_moon_phase_uses_uniform_full_moon_color() -> None:
    rgba = generate_flat_moon_phase_rgba(
        64,
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
        edge_soft_px=0.0,
    )

    interior = rgba[16:48, 16:48]
    opaque_rgb = interior[interior[..., 3] > 0, :3]
    assert np.unique(opaque_rgb, axis=0).shape[0] == 1


def test_generate_flat_moon_phase_shows_bright_and_dark_halves() -> None:
    rgba = generate_flat_moon_phase_rgba(
        64,
        np.array([1.0, 0.0, 0.0], dtype=np.float32),
        edge_soft_px=0.0,
    )

    assert int(rgba[32, 48, 0]) > int(rgba[32, 16, 0])
    assert int(rgba[32, 48, 3]) == 255
    assert int(rgba[32, 16, 3]) == 255


def test_generate_flat_moon_phase_applies_eclipse_tint() -> None:
    untinted = generate_flat_moon_phase_rgba(
        64,
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
    )
    tinted = generate_flat_moon_phase_rgba(
        64,
        np.array([0.0, 0.0, 1.0], dtype=np.float32),
        tint_color=(80, 10, 10, 180),
    )

    assert int(tinted[32, 32, 0]) < int(untinted[32, 32, 0])
    assert int(tinted[32, 32, 0]) > int(tinted[32, 32, 1])
