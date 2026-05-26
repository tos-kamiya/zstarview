from __future__ import annotations

import io
import datetime as dt
from pathlib import Path

import pytest
import numpy as np
from PIL import Image

from zstarview.geosatellite.cache import compute_digest, raw_cache_stem
from zstarview.geosatellite.types import GeoSatelliteDownloadResult
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


def test_run_geo_satellite_pipeline_uses_raw_cache_hit(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    from zstarview.geosatellite import pipeline as geo_pipeline

    cache_root = tmp_path / "geo-cache"
    monkeypatch.setattr(
        geo_pipeline,
        "download_latest_image",
        lambda **_kwargs: (_ for _ in ()).throw(AssertionError("download should not run on cache hit")),
    )
    monkeypatch.setattr(
        geo_pipeline,
        "read_raw_cache",
        lambda **_kwargs: GeoSatelliteDownloadResult(
            fetched_at_utc=dt.datetime(2026, 5, 26, 9, 31, tzinfo=dt.timezone.utc),
            captured_at_utc=dt.datetime(2026, 5, 26, 9, 30, tzinfo=dt.timezone.utc),
            kind="infrared",
            source_url="cache",
            png_bytes=raw_png,
            content_type="image/png",
            cache_path=cache_root / "raw.png",
            metadata_path=cache_root / "raw.json",
        ),
    )
    monkeypatch.setattr(geo_pipeline, "read_intermediate_cache", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(geo_pipeline, "write_raw_cache", lambda result: (cache_root / "raw.png", cache_root / "raw.json"))
    monkeypatch.setattr(geo_pipeline, "write_intermediate_cache", lambda **_kwargs: (cache_root / "proxy.png", cache_root / "inpainted.png", cache_root / "manifest.json"))
    monkeypatch.setattr(geo_pipeline, "build_proxy_image", lambda image, **_kwargs: Image.fromarray(np.full((64, 64), 128, dtype=np.uint8), mode="L"))
    monkeypatch.setattr(geo_pipeline, "build_inpainted_image", lambda proxy_image, mask_image, **_kwargs: proxy_image)
    monkeypatch.setattr(geo_pipeline, "project_gray_image_to_disc", lambda *_args, **_kwargs: np.zeros((8, 8), dtype=np.uint8))

    raw_image = Image.new("RGB", (64, 64), (128, 128, 128))
    raw_buffer = io.BytesIO()
    raw_image.save(raw_buffer, format="PNG")
    raw_png = raw_buffer.getvalue()

    with caplog.at_level("INFO", logger="zstarview.geosatellite.pipeline"):
        result = run_geo_satellite_pipeline(
            observer_lat=51.5,
            observer_lon=-0.1,
            alt=10.0,
            az=180.0,
            fov_deg=60.0,
            download_time_utc=dt.datetime(2026, 5, 26, 9, 30, tzinfo=dt.timezone.utc),
            grid_npz=Path("src/zstarview/data/geosatellite/eqdc_lonlat.npz"),
        )

    assert result.download.source_url == "cache"
    assert "Geo-sat raw cache hit: infrared 20260526T0930Z" in caplog.text
    assert result.download.captured_at_utc == dt.datetime(2026, 5, 26, 9, 30, tzinfo=dt.timezone.utc)


def test_run_geo_satellite_pipeline_uses_available_latest_time(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    from zstarview.geosatellite import cache as geo_cache
    from zstarview.geosatellite import pipeline as geo_pipeline

    cache_root = tmp_path / "geo-cache"
    monkeypatch.setattr(geo_cache, "GEOSATELLITE_CACHE_ROOT_DIR", str(cache_root))
    monkeypatch.setattr(geo_pipeline, "read_latest_available_cache", lambda **_kwargs: None)
    monkeypatch.setattr(
        geo_pipeline,
        "fetch_latest_available_image_time",
        lambda **_kwargs: (
            dt.datetime(2026, 5, 26, 13, 15, tzinfo=dt.timezone.utc),
            42,
            "https://api.met.no/weatherapi/geosatellite/1.4/available.json",
        ),
    )
    written_slots: list[str] = []
    monkeypatch.setattr(
        geo_pipeline,
        "write_latest_available_cache",
        lambda **kwargs: written_slots.append(kwargs["available_time_utc"].isoformat()) or (cache_root / "latest.json"),
    )
    monkeypatch.setattr(
        geo_pipeline,
        "download_latest_image",
        lambda **_kwargs: (_ for _ in ()).throw(AssertionError("download should not run when available time is cached")),
    )
    monkeypatch.setattr(
        geo_pipeline,
        "read_raw_cache",
        lambda **kwargs: GeoSatelliteDownloadResult(
            fetched_at_utc=dt.datetime(2026, 5, 26, 13, 16, tzinfo=dt.timezone.utc),
            captured_at_utc=kwargs["image_time_utc"],
            kind="infrared",
            source_url="cache",
            png_bytes=raw_png,
            content_type="image/png",
            cache_path=cache_root / "raw.png",
            metadata_path=cache_root / "raw.json",
        ),
    )
    monkeypatch.setattr(geo_pipeline, "read_intermediate_cache", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(geo_pipeline, "write_raw_cache", lambda result: (cache_root / "raw.png", cache_root / "raw.json"))
    monkeypatch.setattr(geo_pipeline, "write_intermediate_cache", lambda **_kwargs: (cache_root / "proxy.png", cache_root / "inpainted.png", cache_root / "manifest.json"))
    monkeypatch.setattr(geo_pipeline, "build_proxy_image", lambda image, **_kwargs: Image.fromarray(np.full((64, 64), 128, dtype=np.uint8), mode="L"))
    monkeypatch.setattr(geo_pipeline, "build_inpainted_image", lambda proxy_image, mask_image, **_kwargs: proxy_image)
    monkeypatch.setattr(geo_pipeline, "project_gray_image_to_disc", lambda *_args, **_kwargs: np.zeros((8, 8), dtype=np.uint8))

    raw_image = Image.new("RGB", (64, 64), (128, 128, 128))
    raw_buffer = io.BytesIO()
    raw_image.save(raw_buffer, format="PNG")
    raw_png = raw_buffer.getvalue()

    with caplog.at_level("INFO", logger="zstarview.geosatellite.pipeline"):
        result = run_geo_satellite_pipeline(
            observer_lat=51.5,
            observer_lon=-0.1,
            alt=10.0,
            az=180.0,
            fov_deg=60.0,
            grid_npz=Path("src/zstarview/data/geosatellite/eqdc_lonlat.npz"),
        )

    assert result.download.captured_at_utc == dt.datetime(2026, 5, 26, 13, 15, tzinfo=dt.timezone.utc)
    assert written_slots == [dt.datetime(2026, 5, 26, 13, 15, tzinfo=dt.timezone.utc).isoformat()]
    assert "Geo-sat raw cache hit: infrared 20260526T1315Z" in caplog.text


def test_cache_key_uses_minute_precision() -> None:
    import datetime as dt

    stamp = dt.datetime(2026, 5, 26, 9, 30, 45, tzinfo=dt.timezone.utc)
    assert raw_cache_stem(image_time_utc=stamp, kind="infrared") == "20260526T0930Z_infrared"
    assert compute_digest(b"abc") == compute_digest(b"abc")
