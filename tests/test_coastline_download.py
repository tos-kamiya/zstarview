from __future__ import annotations

import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path

import pytest

from zstarview.cli import download_coastline as download_coastline_cli
from zstarview.coastline_download import (
    COASTLINE_DATASET_VERSION,
    WATER_MASK_ASSET_NAME,
    CoastlineDownloadError,
    _extract_checked,
    download_coastline_data,
    download_water_mask_25m,
    select_columns,
)


def _archive(*, column: int = 0, unsafe_name: str | None = None) -> bytes:
    payload = json.dumps(
        {
            "type": "FeatureCollection",
            "features": [],
        }
    ).encode("utf-8")
    compressed = io.BytesIO()
    with gzip.GzipFile(fileobj=compressed, mode="wb", mtime=0) as handle:
        handle.write(payload)
    name = unsafe_name or f"grid-32x16/y07/x{column:02d}/tile.geojson.gz"
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as archive:
        archive.writestr(name, compressed.getvalue())
    return output.getvalue()


def _manifest(archive_bytes: bytes) -> bytes:
    assets = []
    for column in range(32):
        name = f"coastline-grid-x{column:02d}-{COASTLINE_DATASET_VERSION}.zip"
        assets.append(
            {
                "x": column,
                "name": name,
                "bytes": len(archive_bytes),
                "sha256": hashlib.sha256(archive_bytes).hexdigest(),
            }
        )
    return json.dumps(
        {
            "schema": 1,
            "dataset": "coastline-vector-tiles",
            "version": COASTLINE_DATASET_VERSION,
            "grid": {"columns": 32, "rows": 16},
            "assets": assets,
        }
    ).encode("utf-8")


def test_select_columns_expands_to_grid_boundaries() -> None:
    assert select_columns(lon_min=121.5, lon_max=123.75) == (26, 27)
    assert select_columns(all_columns=True) == tuple(range(32))


def test_select_columns_rejects_invalid_selection() -> None:
    with pytest.raises(ValueError):
        select_columns(lon_min=1.0)
    with pytest.raises(ValueError):
        select_columns(lon_min=2.0, lon_max=1.0)
    with pytest.raises(ValueError):
        select_columns(all_columns=True, lon_min=0.0, lon_max=1.0)


def test_download_installs_selected_column(tmp_path: Path) -> None:
    archive_bytes = _archive()
    manifest_bytes = _manifest(archive_bytes)
    base_url = "https://example.test/release"
    responses = {
        f"{base_url}/coastline-manifest-20260725.json": manifest_bytes,
        f"{base_url}/coastline-grid-x00-20260725.zip": archive_bytes,
    }

    def read_url(url: str, *, timeout_s: float) -> bytes:
        assert timeout_s > 0
        return responses[url]

    assert download_coastline_data(
        lon_min=-180.0,
        lon_max=-170.0,
        cache_dir=tmp_path,
        base_url=base_url,
        read_url=read_url,
    ) == (0,)
    root = tmp_path / "osm-water-polygons" / "20260725" / "schema-1"
    assert (root / "READY").exists()
    assert (root / "manifest.json").exists()
    assert (root / "grid-32x16/y07/x00/tile.geojson.gz").exists()
    assert download_coastline_data(
        lon_min=-180.0,
        lon_max=-170.0,
        cache_dir=tmp_path,
        base_url=base_url,
        read_url=read_url,
    ) == (0,)


def test_extract_rejects_unsafe_zip_entry(tmp_path: Path) -> None:
    archive = tmp_path / "unsafe.zip"
    archive.write_bytes(_archive(unsafe_name="grid-32x16/y07/x00/../../escape"))
    with pytest.raises(CoastlineDownloadError):
        _extract_checked(archive, tmp_path / "out", 0)


def test_download_installs_25m_water_mask_without_extracting(tmp_path: Path) -> None:
    archive_bytes = b"test-water-mask-zip"
    base_url = "https://example.test/release"
    manifest = {
        "schema": 1,
        "data_source_date": "2026-07-25",
        "coverage": {"columns": 32, "latitude_rows": 16},
        "raster": {"resolution_m": 25},
        "assets": [
            {
                "name": WATER_MASK_ASSET_NAME,
                "bytes": len(archive_bytes),
                "sha256": hashlib.sha256(archive_bytes).hexdigest(),
            }
        ],
    }
    responses = {
        f"{base_url}/water-manifest-20260725.json": json.dumps(manifest).encode(),
        f"{base_url}/{WATER_MASK_ASSET_NAME}": archive_bytes,
    }

    def read_url(url: str, *, timeout_s: float) -> bytes:
        assert timeout_s > 0
        return responses[url]

    root = download_water_mask_25m(
        cache_dir=tmp_path,
        base_url=base_url,
        read_url=read_url,
    )
    assert root == tmp_path / "water" / "osm-water-polygons" / "20260725" / "schema-1" / "resolution-25m"
    assert (root / "READY").exists()
    assert (root / "manifest.json").exists()
    assert (root / WATER_MASK_ASSET_NAME).read_bytes() == archive_bytes


def test_cli_all_downloads_the_25m_water_mask(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []

    def fake_coastline(**kwargs: object) -> tuple[int, ...]:
        calls.append("coastline")
        assert kwargs["all_columns"] is True
        return tuple(range(32))

    def fake_water_mask(**kwargs: object) -> Path:
        calls.append("water")
        return Path("water-cache")

    monkeypatch.setattr(download_coastline_cli, "download_coastline_data", fake_coastline)
    monkeypatch.setattr(download_coastline_cli, "download_water_mask_25m", fake_water_mask)

    assert download_coastline_cli.main(["--all"]) == 0
    assert calls == ["coastline", "water"]


def test_cli_water_mask_can_be_downloaded_alone(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []

    def fail_coastline(**_kwargs: object) -> tuple[int, ...]:
        raise AssertionError("coastline download should be skipped")

    def fake_water_mask(**_kwargs: object) -> Path:
        calls.append("water")
        return Path("water-cache")

    monkeypatch.setattr(download_coastline_cli, "download_coastline_data", fail_coastline)
    monkeypatch.setattr(download_coastline_cli, "download_water_mask_25m", fake_water_mask)

    assert download_coastline_cli.main(["--water-25m"]) == 0
    assert calls == ["water"]
