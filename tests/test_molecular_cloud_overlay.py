from __future__ import annotations

import numpy as np
import pytest

from zstarview.render import molecular_cloud_overlay
from zstarview.render.molecular_cloud_overlay import (
    _apply_akari_two_band_mapping,
    _apply_creative_hubble_mapping,
    _apply_molecular_cloud_value_knee,
    _upscale_molecular_cloud_overlay,
    molecular_cloud_cache_path,
)


def test_apply_molecular_cloud_value_knee_preserves_hue_and_saturation() -> None:
    rgb = np.array(
        [[0.2, 0.8, 0.4], [1.0, 0.0, 0.7]],
        dtype=np.float32,
    )

    result = _apply_molecular_cloud_value_knee(rgb)

    np.testing.assert_allclose(result, [[1.0 / 9.0, 4.0 / 9.0, 2.0 / 9.0], [0.5, 0.0, 0.35]])
    np.testing.assert_allclose(rgb, [[0.2, 0.8, 0.4], [1.0, 0.0, 0.7]])


def test_creative_hubble_mapping_matches_article_formula() -> None:
    rgb = np.array([[0.8, 0.4, 0.2]], dtype=np.float32)

    result = _apply_creative_hubble_mapping(rgb)

    np.testing.assert_allclose(result, [[0.8, 0.5, 0.2]])


def test_akari_two_band_mapping_uses_blended_green_channel() -> None:
    bands = {
        90: np.array([0.2, 0.8], dtype=np.float32),
        140: np.array([0.7, 0.3], dtype=np.float32),
    }

    result = _apply_akari_two_band_mapping(bands, np.zeros(2, dtype=np.float32))

    np.testing.assert_allclose(result, [[0.7, 0.45, 0.2], [0.3, 0.55, 0.8]])


def test_upscale_molecular_cloud_overlay_returns_target_size() -> None:
    overlay = np.zeros((2, 2, 4), dtype=np.uint8)
    overlay[..., 3] = 255
    overlay[0, 0, 0] = 255

    result = _upscale_molecular_cloud_overlay(overlay, width=8, height=6)

    assert result.shape == (6, 8, 4)
    assert 0 < int(result[2, 2, 0]) < 255
    assert int(result[3, 7, 0]) == 0


def test_molecular_cloud_cache_path_uses_explicit_source() -> None:
    assert molecular_cloud_cache_path("gaia") == (
        molecular_cloud_overlay.GAIA_MOLECULAR_CLOUD_CACHE
    )
    assert molecular_cloud_cache_path("akari") == (
        molecular_cloud_overlay.AKARI_MOLECULAR_CLOUD_CACHE
    )


def test_molecular_cloud_cache_path_rejects_unknown_source() -> None:
    with pytest.raises(ValueError, match="unsupported diffuse sky source"):
        molecular_cloud_cache_path("unknown")
