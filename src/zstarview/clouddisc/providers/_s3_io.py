# -*- coding: utf-8 -*-
"""
Shared S3 I/O helpers for satellite providers.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from botocore.exceptions import ConnectTimeoutError, ReadTimeoutError

from ..types import CloudMeta, DownloadError, TimeoutError


def list_s3_keys(
    *,
    s3_client,
    bucket: str,
    prefix: str,
    satellite: str,
    product: str,
    time_utc,
    uri_label: str | None = None,
) -> List[str]:
    """List S3 keys under a prefix and normalize provider exceptions."""
    if uri_label is None:
        uri_label = f"s3://{bucket}/{prefix}"

    try:
        paginator = s3_client.get_paginator("list_objects_v2")
        return [
            obj["Key"]
            for page in paginator.paginate(Bucket=bucket, Prefix=prefix)
            for obj in page.get("Contents", []) or []
        ]
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
) -> Path:
    """Download an S3 object with atomic file replacement and unified errors."""
    if dst.exists():
        return dst

    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = dst.with_suffix(dst.suffix + ".tmp")
    try:
        with tmp_path.open("wb") as f:
            s3_client.download_fileobj(bucket, key, f)
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
