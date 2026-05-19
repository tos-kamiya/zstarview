# -*- coding: utf-8 -*-
"""Shared S3 I/O helpers for satellite providers."""

from __future__ import annotations

import builtins
import logging
import socket
import threading
from pathlib import Path
from typing import Callable, List
from urllib.error import HTTPError, URLError
from urllib.parse import quote, unquote, urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree

from ..types import CloudMeta, DownloadCancelledError, DownloadError, TimeoutError

logger = logging.getLogger(__name__)

_S3_HOST = "s3.amazonaws.com"
_LIST_CHUNK_SIZE = 1024 * 1024


def _s3_url(bucket: str, key: str | None = None) -> str:
    if key is None:
        return f"https://{bucket}.{_S3_HOST}/"
    quoted_key = quote(key, safe="/-_.~")
    return f"https://{bucket}.{_S3_HOST}/{quoted_key}"


def _request_timeout(timeout_s: float | None) -> float | None:
    if timeout_s is None:
        return None
    timeout = float(timeout_s)
    return timeout if timeout > 0.0 else 0.0


def _list_bucket_page(
    *,
    bucket: str,
    prefix: str,
    continuation_token: str | None,
    timeout_s: float | None,
) -> tuple[list[str], bool, str | None]:
    params: dict[str, str] = {"list-type": "2", "prefix": prefix, "encoding-type": "url"}
    if continuation_token is not None:
        params["continuation-token"] = continuation_token
    url = f"{_s3_url(bucket)}?{urlencode(params)}"
    req = Request(url, method="GET")
    try:
        with urlopen(req, timeout=_request_timeout(timeout_s)) as resp:
            payload = resp.read()
    except (builtins.TimeoutError, socket.timeout) as e:
        raise TimeoutError(f"Timeout while listing s3://{bucket}/{prefix}") from e
    except HTTPError as e:
        if e.code in {301, 302, 307, 308}:
            raise DownloadError(f"Failed to list s3://{bucket}/{prefix}") from e
        if e.code == 404:
            return [], False, None
        raise DownloadError(f"Failed to list s3://{bucket}/{prefix}") from e
    except URLError as e:
        if isinstance(getattr(e, "reason", None), (builtins.TimeoutError, socket.timeout)):
            raise TimeoutError(f"Timeout while listing s3://{bucket}/{prefix}") from e
        raise DownloadError(f"Failed to list s3://{bucket}/{prefix}") from e

    try:
        root = ElementTree.fromstring(payload)
    except ElementTree.ParseError as e:
        raise DownloadError(f"Failed to parse listing for s3://{bucket}/{prefix}") from e

    def _iter_elements(node: ElementTree.Element, local_name: str):
        for child in node:
            if child.tag.rsplit("}", 1)[-1] == local_name:
                yield child
            yield from _iter_elements(child, local_name)

    keys: list[str] = []
    for key_elem in _iter_elements(root, "Key"):
        key_text = key_elem.text or ""
        keys.append(unquote(key_text))

    is_truncated_text = next((elem.text for elem in _iter_elements(root, "IsTruncated")), "false")
    next_token = next((elem.text for elem in _iter_elements(root, "NextContinuationToken")), None)
    is_truncated = str(is_truncated_text).strip().lower() == "true"
    return keys, is_truncated, next_token


def list_s3_keys(
    *,
    bucket: str,
    prefix: str,
    satellite: str,
    product: str,
    time_utc,
    uri_label: str | None = None,
    abort_event: threading.Event | None = None,
    timeout_s: float | None = None,
) -> List[str]:
    """List S3 keys under a prefix and normalize provider exceptions."""
    if uri_label is None:
        uri_label = f"s3://{bucket}/{prefix}"
    if abort_event is not None and abort_event.is_set():
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        raise DownloadCancelledError(
            f"Cancelled while listing {uri_label}",
            meta=meta,
        )

    keys: list[str] = []
    continuation_token: str | None = None
    try:
        while True:
            page_keys, is_truncated, continuation_token = _list_bucket_page(
                bucket=bucket,
                prefix=prefix,
                continuation_token=continuation_token,
                timeout_s=timeout_s,
            )
            keys.extend(page_keys)
            if abort_event is not None and abort_event.is_set():
                meta = CloudMeta(
                    satellite=satellite,
                    product=product,
                    time_utc=time_utc,
                    src_paths=[],
                )
                raise DownloadCancelledError(
                    f"Cancelled while listing {uri_label}",
                    meta=meta,
                )
            if not is_truncated:
                return keys
            if not continuation_token:
                break
    except DownloadCancelledError:
        raise
    except TimeoutError as e:
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        raise TimeoutError(f"Timeout while listing {uri_label}", meta=meta) from e
    except Exception as e:
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        raise DownloadError(f"Failed to list {uri_label}", meta=meta) from e

    meta = CloudMeta(
        satellite=satellite,
        product=product,
        time_utc=time_utc,
        src_paths=[],
    )
    raise DownloadError(f"Failed to list {uri_label}", meta=meta)


def download_s3_object(
    *,
    bucket: str,
    key: str,
    dst: Path,
    satellite: str,
    product: str,
    time_utc,
    validate_func: Callable[[Path], None] | None = None,
    abort_event: threading.Event | None = None,
    timeout_s: float | None = None,
) -> Path:
    """Download an S3 object with atomic file replacement and unified errors."""
    if dst.exists():
        if validate_func is None:
            return dst
        try:
            validate_func(dst)
            return dst
        except Exception:
            logger.warning("Discarding invalid cached file: %s", dst)
            dst.unlink(missing_ok=True)

    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = dst.with_suffix(dst.suffix + ".tmp")
    try:
        if abort_event is not None and abort_event.is_set():
            meta = CloudMeta(
                satellite=satellite,
                product=product,
                time_utc=time_utc,
                src_paths=[],
            )
            raise DownloadCancelledError(
                f"Cancelled while downloading s3://{bucket}/{key}",
                meta=meta,
            )
        req = Request(_s3_url(bucket, key), method="GET")
        with urlopen(req, timeout=_request_timeout(timeout_s)) as resp, tmp_path.open("wb") as f:
            while True:
                if abort_event is not None and abort_event.is_set():
                    raise KeyboardInterrupt()
                chunk = resp.read(_LIST_CHUNK_SIZE)
                if not chunk:
                    break
                f.write(chunk)
        if validate_func is not None:
            try:
                validate_func(tmp_path)
            except Exception as e:
                logger.warning("Discarding invalid downloaded file: %s", tmp_path)
                meta = CloudMeta(
                    satellite=satellite,
                    product=product,
                    time_utc=time_utc,
                    src_paths=[],
                )
                raise DownloadError(
                    f"Downloaded s3://{bucket}/{key} failed validation",
                    meta=meta,
                ) from e
        tmp_path.replace(dst)
    except (builtins.TimeoutError, socket.timeout) as e:
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        raise TimeoutError(
            f"Timeout while downloading s3://{bucket}/{key}",
            meta=meta,
        ) from e
    except TimeoutError as e:
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        raise TimeoutError(
            f"Timeout while downloading s3://{bucket}/{key}",
            meta=meta,
        ) from e
    except HTTPError as e:
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        if e.code in {301, 302, 307, 308}:
            raise DownloadError(
                f"Failed to download s3://{bucket}/{key}",
                meta=meta,
            ) from e
        raise DownloadError(
            f"Failed to download s3://{bucket}/{key}",
            meta=meta,
        ) from e
    except URLError as e:
        if isinstance(getattr(e, "reason", None), (builtins.TimeoutError, socket.timeout)):
            meta = CloudMeta(
                satellite=satellite,
                product=product,
                time_utc=time_utc,
                src_paths=[],
            )
            raise TimeoutError(
                f"Timeout while downloading s3://{bucket}/{key}",
                meta=meta,
            ) from e
        raise DownloadError(
            f"Failed to download s3://{bucket}/{key}",
            meta=CloudMeta(
                satellite=satellite,
                product=product,
                time_utc=time_utc,
                src_paths=[],
            ),
        ) from e
    except KeyboardInterrupt as e:
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        raise DownloadCancelledError(
            f"Cancelled while downloading s3://{bucket}/{key}",
            meta=meta,
        ) from e
    except Exception as e:
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        raise DownloadError(
            f"Failed to download s3://{bucket}/{key}",
            meta=meta,
        ) from e
    finally:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    return dst
