from __future__ import annotations

import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path

import pytest

from zstarview.coastline_download import (
    COASTLINE_DATASET_VERSION,
    CoastlineDownloadError,
    _extract_checked,
    download_coastline_data,
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
