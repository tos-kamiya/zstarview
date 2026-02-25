from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import io
import sys

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

botocore_exceptions = pytest.importorskip("botocore.exceptions")
ConnectTimeoutError = botocore_exceptions.ConnectTimeoutError
ReadTimeoutError = botocore_exceptions.ReadTimeoutError

from zstarview.clouddisc.providers._s3_io import download_s3_object, list_s3_keys
from zstarview.clouddisc.types import DownloadError, TimeoutError


class _Paginator:
    def __init__(self, pages=None, exc: Exception | None = None) -> None:
        self._pages = pages or []
        self._exc = exc

    def paginate(self, **kwargs):
        if self._exc is not None:
            raise self._exc
        return self._pages


class _S3Client:
    def __init__(self, *, paginator: _Paginator | None = None, download_exc: Exception | None = None) -> None:
        self._paginator = paginator or _Paginator()
        self._download_exc = download_exc
        self.download_calls = 0

    def get_paginator(self, name: str):
        assert name == "list_objects_v2"
        return self._paginator

    def download_fileobj(self, bucket: str, key: str, fileobj: io.BufferedWriter) -> None:
        self.download_calls += 1
        if self._download_exc is not None:
            raise self._download_exc
        fileobj.write(b"ok")


def test_list_s3_keys_success() -> None:
    s3 = _S3Client(paginator=_Paginator(pages=[{"Contents": [{"Key": "a"}, {"Key": "b"}]}, {"Contents": []}]))
    keys = list_s3_keys(
        s3_client=s3,
        bucket="bucket",
        prefix="prefix/",
        satellite="HIMAWARI",
        product="HSD/ISatSS-B13",
        time_utc=datetime.now(timezone.utc),
    )
    assert keys == ["a", "b"]


def test_list_s3_keys_timeout_maps_to_custom_error() -> None:
    s3 = _S3Client(paginator=_Paginator(exc=ConnectTimeoutError(endpoint_url="https://example.com")))
    t0 = datetime(2026, 1, 1, tzinfo=timezone.utc)
    with pytest.raises(TimeoutError) as err:
        list_s3_keys(
            s3_client=s3,
            bucket="bucket",
            prefix="prefix/",
            satellite="G16",
            product="CMIPF-C13",
            time_utc=t0,
            uri_label="S3 bucket s3://bucket/prefix/",
        )
    assert "Timeout while listing" in str(err.value)
    assert err.value.meta is not None
    assert err.value.meta.satellite == "G16"
    assert err.value.meta.product == "CMIPF-C13"
    assert err.value.meta.time_utc == t0


def test_list_s3_keys_generic_error_maps_to_download_error() -> None:
    s3 = _S3Client(paginator=_Paginator(exc=RuntimeError("boom")))
    with pytest.raises(DownloadError):
        list_s3_keys(
            s3_client=s3,
            bucket="bucket",
            prefix="prefix/",
            satellite="HIMAWARI",
            product="HSD/ISatSS-B13",
            time_utc=datetime.now(timezone.utc),
        )


def test_download_s3_object_skips_when_destination_exists(tmp_path: Path) -> None:
    dst = tmp_path / "cached.bin"
    dst.write_bytes(b"cached")
    s3 = _S3Client(download_exc=RuntimeError("must not be called"))
    out = download_s3_object(
        s3_client=s3,
        bucket="bucket",
        key="key",
        dst=dst,
        satellite="G18",
        product="CMIPF-C13",
        time_utc=datetime.now(timezone.utc),
    )
    assert out == dst
    assert dst.read_bytes() == b"cached"
    assert s3.download_calls == 0


def test_download_s3_object_timeout_maps_to_custom_error(tmp_path: Path) -> None:
    dst = tmp_path / "new.bin"
    s3 = _S3Client(download_exc=ReadTimeoutError(endpoint_url="https://example.com", error="timeout"))
    with pytest.raises(TimeoutError) as err:
        download_s3_object(
            s3_client=s3,
            bucket="bucket",
            key="k",
            dst=dst,
            satellite="HIMAWARI",
            product="HSD/ISatSS-B13",
            time_utc=datetime.now(timezone.utc),
        )
    assert "Timeout while downloading" in str(err.value)
    assert err.value.meta is not None
    assert err.value.meta.satellite == "HIMAWARI"
    assert not dst.exists()

