"""
Data provider for Himawari series satellites (Himawari-8 and Himawari-9).

This module fetches data from the NOAA Open Data Dissemination (NODD) program
on AWS S3. It prioritizes Himawari Standard Data (HSD) and falls back to
ISatSS L2 data if HSD is not available. It specifically loads band 13.
"""

from __future__ import annotations

import datetime as dt
import logging
from pathlib import Path
import re
from typing import List, Optional, Tuple

import boto3
import numpy as np
import xarray as xr
from botocore import UNSIGNED
from botocore.config import Config
from satpy import Scene

from ..config import CloudDiscConfig
from ..types import DataNotFoundError, DownloadError, RenderError, CloudMeta

# Default S3 client configuration for anonymous access with retries.
# 短いタイムアウト・最小リトライで「早く失敗」
_S3CFG = Config(
    signature_version=UNSIGNED,
    retries={"max_attempts": 1, "mode": "standard"},
    connect_timeout=5,
    read_timeout=30,
)

# S3 buckets for Himawari data, ordered by priority (9 then 8).
_HIMA_BUCKETS = ["noaa-himawari9", "noaa-himawari8"]

# Regex for HSD Band 13 segments: AHI-L1b-FLDK/YYYY/MM/DD/HHMM/HS_H0[89]_YYYYMMDD_HHMM_B13_FLDK_R.._S####.DAT(.bz2)
_RE_HSD = re.compile(r"HS_H0[89]_\d{8}_\d{4}_B13_FLDK_R\d{2}_S\d{4}\.DAT(\.bz2)?$", re.IGNORECASE)

# Regex for ISatSS L2 Band 13 NetCDF files: AHI-L2-FLDK-ISatSS/YYYY/MM/DD/HHMM/....nc
_RE_IS = re.compile(r".*B13.*\.nc$", re.IGNORECASE)

DEBUG_HIMA = False
logger = logging.getLogger(__name__)


def _format_dt_for_s3_prefix(t: dt.datetime) -> Tuple[int, int, int, str]:
    """
    Formats a datetime object into year, month, day, and a 'HHM0' string.
    The time is rounded down to the nearest 10 minutes.
    """
    t_utc = t.astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
    hm_str = f"{t_utc.hour:02d}{(t_utc.minute // 10) * 10:02d}"
    return t_utc.year, t_utc.month, t_utc.day, hm_str


class HimaProvider:
    """
    A provider for fetching Himawari data from AWS S3.
    """

    def __init__(self, cfg: CloudDiscConfig):
        """
        Initializes the HimaProvider.

        Args:
            cfg: The application configuration.
        """
        self.cfg = cfg
        self.root_hsd = cfg.cache_root() / "hima_hsd"
        self.root_is = cfg.cache_root() / "hima_isatss"
        self.root_hsd.mkdir(parents=True, exist_ok=True)
        self.root_is.mkdir(parents=True, exist_ok=True)
        self.s3 = boto3.client("s3", config=_S3CFG)

    def _download(self, bucket: str, key: str, root: Path) -> Path:
        """
        Downloads a file from S3, caching it locally.

        Args:
            bucket: The S3 bucket name.
            key: The S3 object key.
            root: The root directory for the specific data type's cache.

        Returns:
            The local path to the downloaded file.
        """
        dst = root / bucket / key
        if dst.exists():
            return dst

        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = dst.with_suffix(dst.suffix + ".tmp")

        logger.info("Downloading s3://%s/%s", bucket, key, extra={"sat": "HIMAWARI", "bucket": bucket})
        try:
            with tmp_path.open("wb") as f:
                self.s3.download_fileobj(bucket, key, f)
        except Exception as e:
            meta = CloudMeta(satellite="HIMAWARI", product="HSD-B13" if "L1b" in key or "FLDK" in key else "ISatSS-B13",
                             time_utc=dt.datetime.utcnow().replace(tzinfo=dt.timezone.utc), src_paths=[])
            raise DownloadError(f"Failed to download s3://{bucket}/{key}", transient=True, meta=meta) from e
        tmp_path.replace(dst)
        return dst

    def _list(self, bucket: str, prefix: str) -> List[str]:
        """
        Lists all object keys under a given S3 prefix.

        Args:
            bucket: The S3 bucket name.
            prefix: The S3 prefix to list.

        Returns:
            A list of S3 object keys.
        """
        try:
            paginator = self.s3.get_paginator("list_objects_v2")
            pages = paginator.paginate(Bucket=bucket, Prefix=prefix)
            return [obj["Key"] for page in pages for obj in page.get("Contents", [])]
        except Exception as e:
            meta = CloudMeta(satellite="HIMAWARI", product="HSD-B13", 
                             time_utc=dt.datetime.utcnow().replace(tzinfo=dt.timezone.utc), src_paths=[])
            raise DownloadError(f"Failed to list s3://{bucket}/{prefix}", transient=True, meta=meta) from e

    def _find_hsd(self, when_utc: dt.datetime) -> Tuple[Optional[str], Optional[List[str]], Optional[dt.datetime]]:
        """
        Finds Himawari Standard Data (HSD) files for a given time.

        It searches backwards from `when_utc` within the configured time window.

        Returns:
            A tuple (bucket, [keys], timestamp) if found, otherwise (None, None, None).
        """
        for slot in range(0, self.cfg.search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=slot)
            y, m, d, hm = _format_dt_for_s3_prefix(search_time)
            for bucket in _HIMA_BUCKETS:
                prefix = f"AHI-L1b-FLDK/{y:04d}/{m:02d}/{d:02d}/{hm}/"
                keys = [k for k in self._list(bucket, prefix) if _RE_HSD.search(Path(k).name)]
                if keys:
                    keys.sort()
                    if DEBUG_HIMA:
                        logger.debug("Found %d HSD segments in s3://%s/%s", len(keys), bucket, prefix,
                                     extra={"sat": "HIMAWARI", "bucket": bucket, "prefix": prefix})
                    return bucket, keys, search_time
        return None, None, None

    def _find_isatss(self, when_utc: dt.datetime) -> Tuple[Optional[str], Optional[str], Optional[dt.datetime]]:
        """
        Finds ISatSS L2 NetCDF files for a given time.

        It searches backwards from `when_utc` within the configured time window.

        Returns:
            A tuple (bucket, key, timestamp) if found, otherwise (None, None, None).
        """
        for slot in range(0, self.cfg.search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=slot)
            y, m, d, hm = _format_dt_for_s3_prefix(search_time)
            for bucket in _HIMA_BUCKETS:
                prefix = f"AHI-L2-FLDK-ISatSS/{y:04d}/{m:02d}/{d:02d}/{hm}/"
                keys = [k for k in self._list(bucket, prefix) if _RE_IS.search(Path(k).name)]
                if keys:
                    keys.sort()
                    if DEBUG_HIMA:
                        logger.debug("Found ISatSS file in s3://%s/%s", bucket, prefix,
                                     extra={"sat": "HIMAWARI", "bucket": bucket, "prefix": prefix})
                    return bucket, keys[-1], search_time  # Return the latest one
        return None, None, None

    def fetch_bt_c13(self, when_utc: dt.datetime) -> Tuple[xr.DataArray, dt.datetime, List[Path], float]:
        """
        Fetches Himawari band 13 brightness temperature data.

        It first attempts to find and load the higher-resolution HSD data. If not found,
        it falls back to the lower-resolution, pre-processed ISatSS data.

        Args:
            when_utc: The target UTC time.

        Returns:
            A tuple containing (DataArray[K], used_time, [paths], sub-satellite_longitude).

        Raises:
            DataNotFoundError: If no data is found within the search window.
        """
        # --- 1. Try to find HSD data (higher quality) ---
        logger.info("Searching Himawari HSD B13 ...")
        hsd_bucket, hsd_keys, hsd_time = self._find_hsd(when_utc)
        if hsd_bucket and hsd_keys and hsd_time:
            try:
                paths = [self._download(hsd_bucket, k, self.root_hsd) for k in hsd_keys]
                scn = Scene(reader="ahi_hsd", filenames=[str(p) for p in paths])
                scn.load(["B13"], calibration="brightness_temperature")
                da = scn["B13"].astype(np.float32).compute()
            except DownloadError:
                raise
            except Exception as e:
                meta = CloudMeta(satellite="HIMAWARI", product="HSD-B13",
                                 time_utc=hsd_time.replace(tzinfo=dt.timezone.utc), src_paths=[])
                raise RenderError("Failed to decode Himawari HSD B13", meta=meta) from e
            # Himawari sub-satellite longitude is ~140.7 E
            used_time = hsd_time.replace(minute=(hsd_time.minute // 10) * 10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
            return da, used_time, paths, 140.7

        # --- 2. Fallback to ISatSS data (lower quality) ---
        logger.info("HSD not found, fallback to ISatSS ...")
        is_bucket, is_key, is_time = self._find_isatss(when_utc)
        if is_bucket and is_key and is_time:
            try:
                path = self._download(is_bucket, is_key, self.root_is)
                scn = Scene(reader="ahi_l2_nc", filenames=[str(path)])
                scn.load(["B13"])
                da = scn["B13"].astype(np.float32).compute()
            except DownloadError:
                raise
            except Exception as e:
                meta = CloudMeta(satellite="HIMAWARI", product="ISatSS-B13",
                                 time_utc=is_time.replace(tzinfo=dt.timezone.utc), src_paths=[])
                raise RenderError("Failed to decode Himawari ISatSS B13", meta=meta) from e
            used_time = is_time.replace(minute=(is_time.minute // 10) * 10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
            return da, used_time, [path], 140.7

        meta = CloudMeta(satellite="HIMAWARI", product="HSD/ISatSS-B13",
                         time_utc=when_utc.replace(tzinfo=dt.timezone.utc), src_paths=[])
        raise DataNotFoundError("Himawari B13 not found in search window", meta=meta)
