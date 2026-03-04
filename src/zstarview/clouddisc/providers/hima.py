# -*- coding: utf-8 -*-
"""
Data provider for Himawari series satellites (Himawari-8 and Himawari-9).

This module fetches data from the NOAA Open Data Dissemination (NODD) program
on AWS S3. It prioritizes the high-resolution Himawari Standard Data (HSD) and
falls back to the pre-processed ISatSS L2 data if HSD is not available. It
specifically loads band 13 (longwave infrared) for brightness temperatures.
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
from ._s3_io import download_s3_object, list_s3_keys

# --- Constants ---
_HIMA_BUCKETS = ["noaa-himawari9", "noaa-himawari8"]  # Priority: H9 over H8
_RE_HSD = re.compile(r"HS_H0[89]_\d{8}_\d{4}_B13_FLDK_R\d{2}_S\d{4}\.DAT(\.bz2)?$", re.IGNORECASE)
_RE_IS = re.compile(r".*B13.*\.nc$", re.IGNORECASE)

logger = logging.getLogger(__name__)


def _format_dt_for_s3_prefix(t: dt.datetime) -> Tuple[int, int, int, str]:
    """Formats a datetime for S3 prefixes, rounding down to the nearest 10 minutes."""
    t_utc = t.astimezone(dt.timezone.utc)
    hm_str = f"{t_utc.hour:02d}{(t_utc.minute // 10) * 10:02d}"
    return t_utc.year, t_utc.month, t_utc.day, hm_str


class HimaProvider:
    """
    A provider for fetching Himawari data from AWS S3.

    Handles the logic for finding and loading both HSD and ISatSS data formats,
    with a preference for the higher-quality HSD data.
    """

    def __init__(self, cfg: CloudDiscConfig):
        self.cfg = cfg
        self.root_hsd = cfg.cache_root() / "hima_hsd"
        self.root_is = cfg.cache_root() / "hima_isatss"
        self.root_hsd.mkdir(parents=True, exist_ok=True)
        self.root_is.mkdir(parents=True, exist_ok=True)

    def _s3(self) -> boto3.client:
        """Creates an anonymous boto3 S3 client."""
        s3_cfg = Config(
            signature_version=UNSIGNED,
            retries={"max_attempts": 1},
            connect_timeout=self.cfg.connect_timeout,
            read_timeout=self.cfg.read_timeout,
        )
        return boto3.client("s3", config=s3_cfg)

    def _download(self, bucket: str, key: str, root: Path) -> Path:
        """Downloads a file from S3, caching it locally using an atomic write."""
        dst = root / bucket / key
        if dst.exists():
            logger.debug("Using cached file: %s", dst)
            return dst

        logger.info("Downloading s3://%s/%s", bucket, key)
        return download_s3_object(
            s3_client=self._s3(),
            bucket=bucket,
            key=key,
            dst=dst,
            satellite="HIMAWARI",
            product="HSD/ISatSS-B13",
            time_utc=dt.datetime.now(dt.timezone.utc),
        )

    def _list(self, bucket: str, prefix: str) -> List[str]:
        """Lists all object keys under a given S3 prefix."""
        return list_s3_keys(
            s3_client=self._s3(),
            bucket=bucket,
            prefix=prefix,
            satellite="HIMAWARI",
            product="HSD/ISatSS-B13",
            time_utc=dt.datetime.now(dt.timezone.utc),
        )

    def _find_hsd(self, when_utc: dt.datetime) -> Tuple[Optional[str], Optional[List[str]], Optional[dt.datetime]]:
        """Finds Himawari Standard Data (HSD) files for a given time."""
        for slot in range(0, self.cfg.search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=slot)
            y, m, d, hm = _format_dt_for_s3_prefix(search_time)
            for bucket in _HIMA_BUCKETS:
                prefix = f"AHI-L1b-FLDK/{y:04d}/{m:02d}/{d:02d}/{hm}/"
                keys = sorted([k for k in self._list(bucket, prefix) if _RE_HSD.search(Path(k).name)])
                # HSD data is split into 10 segments, so we need all of them.
                if len(keys) == 10:
                    logger.debug("Found 10 HSD segments in s3://%s/%s", bucket, prefix)
                    return bucket, keys, search_time
        return None, None, None

    def _find_isatss(self, when_utc: dt.datetime) -> Tuple[Optional[str], Optional[str], Optional[dt.datetime]]:
        """Finds ISatSS L2 NetCDF files for a given time."""
        for slot in range(0, self.cfg.search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=slot)
            y, m, d, hm = _format_dt_for_s3_prefix(search_time)
            for bucket in _HIMA_BUCKETS:
                prefix = f"AHI-L2-FLDK-ISatSS/{y:04d}/{m:02d}/{d:02d}/{hm}/"
                keys = sorted([k for k in self._list(bucket, prefix) if _RE_IS.search(Path(k).name)])
                if keys:
                    logger.debug("Found ISatSS file in s3://%s/%s", bucket, prefix)
                    return bucket, keys[-1], search_time  # Return the latest one
        return None, None, None

    def fetch_bt_c13(self, when_utc: dt.datetime) -> Tuple[xr.DataArray, dt.datetime, List[Path]]:
        """
        Fetches Himawari band 13 brightness temperature data with a fallback strategy.

        It first attempts to find and load the higher-resolution HSD data. If not found,
        it falls back to the lower-resolution, pre-processed ISatSS data.

        Args:
            when_utc: The target UTC time.

        Returns:
            A tuple containing (DataArray[K], used_time, [paths]).

        Raises:
            DataNotFoundError: If no data is found within the search window.
            RenderError: If the data cannot be decoded by Satpy.
        """
        # --- 1. Try to find HSD data (higher quality) ---
        logger.info("Searching for Himawari HSD B13 data...")
        hsd_bucket, hsd_keys, hsd_time = self._find_hsd(when_utc)
        if hsd_bucket and hsd_keys and hsd_time:
            try:
                paths = [self._download(hsd_bucket, k, self.root_hsd) for k in hsd_keys]
                scn = Scene(reader="ahi_hsd", filenames=[str(p) for p in paths])
                scn.load(["B13"], calibration="brightness_temperature")
                da = scn["B13"].astype(np.float32).compute()
                used_time = hsd_time.replace(minute=(hsd_time.minute // 10) * 10, tzinfo=dt.timezone.utc)
                return da, used_time, paths
            except DownloadError:
                raise  # Propagate download errors immediately
            except Exception as e:
                meta = CloudMeta(satellite="HIMAWARI", product="HSD-B13", time_utc=hsd_time, src_paths=[])
                raise RenderError("Failed to decode Himawari HSD B13 data", meta=meta) from e

        # --- 2. Fallback to ISatSS data (lower quality) ---
        logger.info("HSD data not found, falling back to ISatSS...")
        is_bucket, is_key, is_time = self._find_isatss(when_utc)
        if is_bucket and is_key and is_time:
            try:
                path = self._download(is_bucket, is_key, self.root_is)
                scn = Scene(reader="ahi_l2_nc", filenames=[str(path)])
                scn.load(["B13"])
                da = scn["B13"].astype(np.float32).compute()
                used_time = is_time.replace(minute=(is_time.minute // 10) * 10, tzinfo=dt.timezone.utc)
                return da, used_time, [path]
            except DownloadError:
                raise
            except Exception as e:
                meta = CloudMeta(satellite="HIMAWARI", product="ISatSS-B13", time_utc=is_time, src_paths=[])
                raise RenderError("Failed to decode Himawari ISatSS B13 data", meta=meta) from e

        meta = CloudMeta(satellite="HIMAWARI", product="HSD/ISatSS-B13", time_utc=when_utc, src_paths=[])
        raise DataNotFoundError("Himawari B13 data not found in search window", meta=meta)
