# -*- coding: utf-8 -*-
"""
Shared S3 I/O helpers for satellite providers.
"""

from __future__ import annotations

import logging
import threading
from pathlib import Path
from typing import Callable, List

from botocore.exceptions import ConnectTimeoutError, ReadTimeoutError

from ..types import CloudMeta, DownloadCancelledError, DownloadError, TimeoutError

logger = logging.getLogger(__name__)


def list_s3_keys(
    *,
    s3_client,
    bucket: str,
    prefix: str,
    satellite: str,
    product: str,
    time_utc,
    uri_label: str | None = None,
    abort_event: threading.Event | None = None,
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

    try:
        paginator = s3_client.get_paginator("list_objects_v2")
        keys = [
            obj["Key"]
            for page in paginator.paginate(Bucket=bucket, Prefix=prefix)
            for obj in page.get("Contents", []) or []
        ]
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
        return keys
    except (ConnectTimeoutError, ReadTimeoutError) as e:
        meta = CloudMeta(
            satellite=satellite,
            product=product,
            time_utc=time_utc,
            src_paths=[],
        )
        raise TimeoutError(
            f"Timeout while listing {uri_label}",
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
            f"Failed to list {uri_label}",
            meta=meta,
        ) from e


def download_s3_object(
    *,
    s3_client,
    bucket: str,
    key: str,
    dst: Path,
    satellite: str,
    product: str,
    time_utc,
    validate_func: Callable[[Path], None] | None = None,
    abort_event: threading.Event | None = None,
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
        with tmp_path.open("wb") as f:
            callback = None
            if abort_event is not None:
                def callback(_bytes_transferred: int) -> None:
                    if abort_event.is_set():
                        raise KeyboardInterrupt()

            s3_client.download_fileobj(bucket, key, f, Callback=callback)
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
    except (ConnectTimeoutError, ReadTimeoutError) as e:
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
