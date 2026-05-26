from __future__ import annotations

import io
from pathlib import Path

import pytest
import numpy as np
from PIL import Image

from zstarview.geosatellite.cache import compute_digest, raw_cache_stem
from zstarview.geosatellite.mask import build_common_mask, fill_masked_regions
from zstarview.geosatellite.pipeline import run_geo_satellite_pipeline
from zstarview.geosatellite.projection import render_gray_image_to_cloud_rgba
from zstarview.geosatellite.pipeline import is_within_europe_band
from zstarview.geosatellite.proxy import build_cloud_proxy


def test_is_within_europe_band() -> None:
    assert is_within_europe_band(51.5, -0.1)
    assert not is_within_europe_band(10.0, 0.0)
    assert not is_within_europe_band(51.5, 120.0)


def test_build_cloud_proxy_masks_logo_area() -> None:
    rgb = np.zeros((20, 20, 3), dtype=np.uint8)
    rgb[..., 0] = 180
    rgb[..., 1] = 175
    rgb[..., 2] = 170
    rgb[:4, :4] = (255, 255, 255)
    proxy = build_cloud_proxy(
        Image.fromarray(rgb, mode="RGB"),
        contrast_low=0.0,
        contrast_high=100.0,
    )
    proxy_arr = np.asarray(proxy, dtype=np.uint8)
    assert proxy.mode == "L"
    assert proxy_arr.shape == (20, 20)
    assert int(proxy_arr[0, 0]) == 0
    assert int(proxy_arr[10, 10]) > 0


def test_build_common_mask_keeps_gray_pixels() -> None:
    base = np.full((10, 10, 3), 120, dtype=np.uint8)
    base[0, 0] = (0, 0, 255)
    other = np.full((10, 10, 3), 125, dtype=np.uint8)
    other[0, 0] = (0, 255, 0)
    mask = build_common_mask([Image.fromarray(base, mode="RGB"), Image.fromarray(other, mode="RGB")])
    mask_arr = np.asarray(mask, dtype=np.uint8)
    assert mask.mode == "L"
    assert int(mask_arr[0, 0]) == 0
    assert int(mask_arr[5, 5]) == 255


def test_fill_masked_regions_propagates_neighbors() -> None:
    image = np.array(
        [
            [10, 10, 10],
            [10, 0, 30],
            [40, 40, 40],
        ],
        dtype=np.float32,
    )
    mask = np.array(
        [
            [False, False, False],
            [False, True, False],
            [False, False, False],
        ],
        dtype=bool,
    )
    filled = fill_masked_regions(image, mask, radius=1, max_iterations=8)
    assert filled.shape == image.shape
    assert float(filled[1, 1]) > 0.0


def test_render_gray_image_to_cloud_rgba_uses_alpha_for_cloud_amount() -> None:
    gray = np.array([[0, 64], [128, 255]], dtype=np.uint8)
    rgba = render_gray_image_to_cloud_rgba(gray)

    assert rgba.shape == (2, 2, 4)
    assert np.array_equal(rgba[..., 0], np.full((2, 2), 255, dtype=np.uint8))
    assert np.array_equal(rgba[..., 1], np.full((2, 2), 255, dtype=np.uint8))
    assert np.array_equal(rgba[..., 2], np.full((2, 2), 255, dtype=np.uint8))
    assert np.array_equal(rgba[..., 3], gray)


def test_run_geo_satellite_pipeline_separates_mask_digest(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from zstarview.geosatellite import cache as geo_cache

    cache_root = tmp_path / "geo-cache"
    monkeypatch.setattr(geo_cache, "GEOSATELLITE_CACHE_ROOT_DIR", str(cache_root))

    raw_image = Image.new("RGB", (64, 64), (128, 128, 128))
    raw_buffer = io.BytesIO()
    raw_image.save(raw_buffer, format="PNG")
    raw_png = raw_buffer.getvalue()

    mask_a = tmp_path / "mask-a.png"
    Image.new("L", (64, 64), 255).save(mask_a)
    mask_b = tmp_path / "mask-b.png"
    mask_img_b = np.zeros((64, 64), dtype=np.uint8)
    mask_img_b[16:48, 16:48] = 255
    Image.fromarray(mask_img_b, mode="L").save(mask_b)

    result_a = run_geo_satellite_pipeline(
        observer_lat=51.5,
        observer_lon=-0.1,
        alt=10.0,
        az=180.0,
        fov_deg=60.0,
        raw_png=raw_png,
        mask_path=mask_a,
        grid_npz=Path("src/zstarview/data/geosatellite/eqdc_lonlat.npz"),
    )
    result_b = run_geo_satellite_pipeline(
        observer_lat=51.5,
        observer_lon=-0.1,
        alt=10.0,
        az=180.0,
        fov_deg=60.0,
        raw_png=raw_png,
        mask_path=mask_b,
        grid_npz=Path("src/zstarview/data/geosatellite/eqdc_lonlat.npz"),
    )

    assert result_a.disc_gray.shape == result_b.disc_gray.shape
    assert result_a.intermediate.manifest is not None
    assert result_b.intermediate.manifest is not None
    assert result_a.intermediate.manifest["mask_digest"] != result_b.intermediate.manifest["mask_digest"]
    assert result_a.intermediate.proxy_path != result_b.intermediate.proxy_path
    assert result_a.intermediate.inpainted_path != result_b.intermediate.inpainted_path


def test_cache_key_uses_minute_precision() -> None:
    import datetime as dt

    stamp = dt.datetime(2026, 5, 26, 9, 30, 45, tzinfo=dt.timezone.utc)
    assert raw_cache_stem(fetched_at_utc=stamp, kind="infrared") == "20260526T0930Z_infrared"
    assert compute_digest(b"abc") == compute_digest(b"abc")
