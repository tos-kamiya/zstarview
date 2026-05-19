from __future__ import annotations

import io
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import URLError

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.clouddisc.providers._s3_io import (  # noqa: E402
    download_s3_object,
    list_s3_keys,
)
from zstarview.clouddisc.types import (  # noqa: E402
    DownloadCancelledError,
    DownloadError,
)
from zstarview.clouddisc.types import TimeoutError as CloudTimeoutError  # noqa: E402


def _list_xml(*keys: str, is_truncated: bool = False, next_token: str | None = None) -> bytes:
    contents = "".join(
        f"<Contents><Key>{key}</Key></Contents>"
        for key in keys
    )
    token_xml = f"<NextContinuationToken>{next_token}</NextContinuationToken>" if next_token is not None else ""
    return (
        "<ListBucketResult>"
        f"{contents}"
        f"<IsTruncated>{'true' if is_truncated else 'false'}</IsTruncated>"
        f"{token_xml}"
        "</ListBucketResult>"
    ).encode("utf-8")


class _FakeResponse:
    def __init__(self, payload: bytes) -> None:
        self._buffer = io.BytesIO(payload)

    def __enter__(self) -> "_FakeResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        return None

    def read(self, size: int = -1) -> bytes:
        return self._buffer.read(size)


class _FakeUrlopen:
    def __init__(self, responses: list[bytes | Exception]) -> None:
        self._responses = list(responses)
        self.requests: list[str] = []

    def __call__(self, request, timeout=None):  # noqa: ANN001
        self.requests.append(request.full_url)
        if not self._responses:
            raise AssertionError("Unexpected urlopen call")
        response = self._responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return _FakeResponse(response)


def test_list_s3_keys_success(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_urlopen = _FakeUrlopen([_list_xml("a", "b")])
    monkeypatch.setattr("zstarview.clouddisc.providers._s3_io.urlopen", fake_urlopen)

    keys = list_s3_keys(
        bucket="bucket",
        prefix="prefix/",
        satellite="HIMAWARI",
        product="HSD/ISatSS-B13",
        time_utc=datetime.now(timezone.utc),
        timeout_s=1.0,
    )

    assert keys == ["a", "b"]
    assert "list-type=2" in fake_urlopen.requests[0]
    assert "prefix=prefix%2F" in fake_urlopen.requests[0]


def test_list_s3_keys_paginates(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_urlopen = _FakeUrlopen([
        _list_xml("a", is_truncated=True, next_token="next-token"),
        _list_xml("b"),
    ])
    monkeypatch.setattr("zstarview.clouddisc.providers._s3_io.urlopen", fake_urlopen)

    keys = list_s3_keys(
        bucket="bucket",
        prefix="prefix/",
        satellite="HIMAWARI",
        product="HSD/ISatSS-B13",
        time_utc=datetime.now(timezone.utc),
    )

    assert keys == ["a", "b"]
    assert len(fake_urlopen.requests) == 2
    assert "continuation-token=next-token" in fake_urlopen.requests[1]


def test_list_s3_keys_timeout_maps_to_custom_error(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_urlopen = _FakeUrlopen([TimeoutError("timed out")])
    monkeypatch.setattr("zstarview.clouddisc.providers._s3_io.urlopen", fake_urlopen)
    t0 = datetime(2026, 1, 1, tzinfo=timezone.utc)

    with pytest.raises(CloudTimeoutError) as err:
        list_s3_keys(
            bucket="bucket",
            prefix="prefix/",
            satellite="G19",
            product="CMIPF-C13",
            time_utc=t0,
            uri_label="S3 bucket s3://bucket/prefix/",
        )

    assert "Timeout while listing" in str(err.value)
    assert err.value.meta is not None
    assert err.value.meta.satellite == "G19"
    assert err.value.meta.product == "CMIPF-C13"
    assert err.value.meta.time_utc == t0


def test_list_s3_keys_generic_error_maps_to_download_error(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_urlopen = _FakeUrlopen([URLError("boom")])
    monkeypatch.setattr("zstarview.clouddisc.providers._s3_io.urlopen", fake_urlopen)

    with pytest.raises(DownloadError):
        list_s3_keys(
            bucket="bucket",
            prefix="prefix/",
            satellite="HIMAWARI",
            product="HSD/ISatSS-B13",
            time_utc=datetime.now(timezone.utc),
        )


def test_download_s3_object_skips_when_destination_exists(tmp_path: Path) -> None:
    dst = tmp_path / "cached.bin"
    dst.write_bytes(b"cached")

    out = download_s3_object(
        bucket="bucket",
        key="key",
        dst=dst,
        satellite="G18",
        product="CMIPF-C13",
        time_utc=datetime.now(timezone.utc),
    )
    assert out == dst
    assert dst.read_bytes() == b"cached"


def test_download_s3_object_discards_invalid_cached_file(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    dst = tmp_path / "cached.bin"
    dst.write_bytes(b"bad")
    fake_urlopen = _FakeUrlopen([b"ok"])
    seen: list[tuple[str, bytes]] = []

    def validate(path: Path) -> None:
        data = path.read_bytes()
        seen.append((path.name, data))
        if data != b"ok":
            raise ValueError("invalid payload")

    monkeypatch.setattr("zstarview.clouddisc.providers._s3_io.urlopen", fake_urlopen)

    out = download_s3_object(
        bucket="bucket",
        key="key",
        dst=dst,
        satellite="G18",
        product="CMIPF-C13",
        time_utc=datetime.now(timezone.utc),
        validate_func=validate,
    )
    assert out == dst
    assert dst.read_bytes() == b"ok"
    assert seen == [("cached.bin", b"bad"), ("cached.bin.tmp", b"ok")]


def test_download_s3_object_rejects_invalid_download(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    dst = tmp_path / "new.bin"
    fake_urlopen = _FakeUrlopen([b"bad"])

    def validate(path: Path) -> None:
        if path.read_bytes() != b"ok":
            raise ValueError("invalid payload")

    monkeypatch.setattr("zstarview.clouddisc.providers._s3_io.urlopen", fake_urlopen)

    with pytest.raises(DownloadError):
        download_s3_object(
            bucket="bucket",
            key="key",
            dst=dst,
            satellite="HIMAWARI",
            product="HSD/ISatSS-B13",
            time_utc=datetime.now(timezone.utc),
            validate_func=validate,
        )
    assert not dst.exists()


def test_download_s3_object_timeout_maps_to_custom_error(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    dst = tmp_path / "new.bin"
    fake_urlopen = _FakeUrlopen([TimeoutError("timed out")])
    monkeypatch.setattr("zstarview.clouddisc.providers._s3_io.urlopen", fake_urlopen)

    with pytest.raises(CloudTimeoutError) as err:
        download_s3_object(
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


def test_download_s3_object_cancelled_download_aborts_and_cleans_up(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    dst = tmp_path / "new.bin"
    abort_event = threading.Event()

    class _StreamingResponse:
        def __init__(self) -> None:
            self._chunks = [b"chunk1", b"chunk2"]

        def __enter__(self) -> "_StreamingResponse":
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def read(self, size: int = -1) -> bytes:
            if self._chunks:
                chunk = self._chunks.pop(0)
                if chunk == b"chunk1":
                    abort_event.set()
                return chunk
            return b""

    class _AbortAfterFirstCall:
        def __init__(self) -> None:
            self.calls = 0

        def __call__(self, request, timeout=None):  # noqa: ANN001
            self.calls += 1
            return _StreamingResponse()

    fake_urlopen = _AbortAfterFirstCall()
    monkeypatch.setattr("zstarview.clouddisc.providers._s3_io.urlopen", fake_urlopen)

    with pytest.raises(DownloadCancelledError) as err:
        download_s3_object(
            bucket="bucket",
            key="k",
            dst=dst,
            satellite="HIMAWARI",
            product="HSD/ISatSS-B13",
            time_utc=datetime.now(timezone.utc),
            abort_event=abort_event,
        )

    assert "Cancelled while downloading" in str(err.value)
    assert err.value.meta is not None
    assert err.value.meta.satellite == "HIMAWARI"
    assert not dst.exists()
    assert fake_urlopen.calls == 1
