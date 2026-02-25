# -*- coding: utf-8 -*-
"""
Data provider for GOES-R series satellites (GOES-16 and GOES-18).

This module fetches Cloud and Moisture Imagery Product (CMIPF) data from the
NOAA Open Data Dissemination (NODD) program on AWS S3. It specifically targets
channel 13 (longwave infrared) to get brightness temperatures for cloud rendering.
"""

import datetime as dt
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import boto3
import xarray as xr
from botocore import UNSIGNED
from botocore.config import Config
from satpy import Scene

from ..config import CloudDiscConfig
from ..types import DataNotFoundError, CloudMeta
from ._s3_io import download_s3_object, list_s3_keys

# --- Constants ---
_GOES_BUCKET = {"G16": "noaa-goes16", "G18": "noaa-goes18"}
_GOES_REGION = {"noaa-goes16": "us-east-1", "noaa-goes18": "us-west-2"}

logger = logging.getLogger(__name__)


def _doy(dt_utc: dt.datetime) -> int:
    """Calculates the day of the year (1-366)."""
    return dt_utc.timetuple().tm_yday


class GoesProvider:
    """
    A provider for fetching GOES-R ABI L2 CMIPF data from AWS S3.

    This class handles the logic for finding the correct data files on S3,
    downloading them to a local cache, and loading them into an xarray DataArray.
    """

    def __init__(self, cfg: CloudDiscConfig):
        self.cfg = cfg
        self.root = cfg.cache_root() / "goes_cmipf"
        self.root.mkdir(parents=True, exist_ok=True)
        self._list_cache: Dict[Tuple, Any] = {}

    def _s3(self, bucket: str) -> boto3.client:
        """Creates an anonymous boto3 S3 client for the specified bucket."""
        cfg = Config(
            signature_version=UNSIGNED,  # No credentials needed for public bucket
            retries={"max_attempts": 1, "mode": "standard"},
            connect_timeout=self.cfg.connect_timeout,
            read_timeout=self.cfg.read_timeout,
        )
        return boto3.client("s3", region_name=_GOES_REGION[bucket], config=cfg)

    def _list_hour(self, bucket: str, t: dt.datetime) -> List[str]:
        """
        Lists all object keys for a given hour in the S3 bucket, with in-memory caching.
        The S3 path is structured as `ABI-L2-CMIPF/YYYY/DOY/HH/`.
        """
        cache_key = (bucket, t.year, _doy(t), t.hour)
        if cache_key in self._list_cache:
            return self._list_cache[cache_key]

        prefix = f"ABI-L2-CMIPF/{t.year:04d}/{_doy(t):03d}/{t.hour:02d}/"
        s3 = self._s3(bucket)
        logger.debug("Listing s3://%s/%s", bucket, prefix)

        keys = list_s3_keys(
            s3_client=s3,
            bucket=bucket,
            prefix=prefix,
            satellite="G16" if "16" in bucket else "G18",
            product="CMIPF-C13",
            time_utc=t,
            uri_label=f"S3 bucket s3://{bucket}/{prefix}",
        )

        logger.debug("Found %d objects under %s", len(keys), prefix)
        self._list_cache[cache_key] = keys
        return keys

    def _download(self, bucket: str, key: str) -> Path:
        """Downloads a file from S3, caching it locally using an atomic write."""
        dst = self.root / bucket / key
        if dst.exists():
            logger.info("Using cached file: %s", dst)
            return dst

        s3 = self._s3(bucket)

        logger.info("Downloading s3://%s/%s", bucket, key)
        return download_s3_object(
            s3_client=s3,
            bucket=bucket,
            key=key,
            dst=dst,
            satellite="G16" if "16" in bucket else "G18",
            product="CMIPF-C13",
            time_utc=dt.datetime.now(dt.timezone.utc),
        )

    def _fetch_bt_c13_once(self, sat: str, when_utc: dt.datetime, search_back_minutes: int) -> Optional[Tuple[xr.DataArray, dt.datetime, List[Path]]]:
        """
        Searches for and loads a single C13 brightness temp file for a given satellite and time.
        """
        bucket = _GOES_BUCKET[sat]

        # Iterate backwards from the target time to find the most recent available file.
        for mback in range(0, search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=mback)
            keys = self._list_hour(bucket, search_time)
            if not keys:
                continue

            # Filter for Channel 13 (C13) files and get the latest one.
            keys_c13 = sorted([k for k in keys if "-M6C13_" in k or "-C13_" in k])
            if not keys_c13:
                continue
            key = keys_c13[-1]

            path = self._download(bucket, key)

            try:
                # Use Satpy to load the NetCDF file, as it correctly handles projection info.
                scn = Scene(reader="abi_l2_nc", filenames=[str(path)])
                scn.load(["C13"])
                da = scn["C13"].astype("float32").compute()
                used_time = search_time.replace(minute=(search_time.minute // 10) * 10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
                return da, used_time, [path]
            except Exception as e:
                logger.warning("Satpy load failed for %s (%s), trying xarray fallback", Path(key).name, e)
                # If Satpy fails, try a direct xarray read as a fallback.
                try:
                    with xr.open_dataset(path, engine="netcdf4") as ds:
                        if "CMI" in ds.variables:
                            da = ds["CMI"].astype("float32").compute()
                            return da, used_time, [path]
                except Exception as ex:
                    logger.error("xarray fallback failed for %s: %s", Path(key).name, ex)
        return None

    def fetch_bt_c13_with_failover(
        self, sat: str, when_utc: dt.datetime, extra_back_minutes: int = 30
    ) -> Tuple[Tuple[xr.DataArray, dt.datetime, List[Path]], str]:
        """
        Fetches C13 data with a two-pass failover strategy.

        Pass 1: Tries the primary satellite, then the secondary, with the standard search window.
        Pass 2: If Pass 1 fails, it widens the search window and retries both satellites.

        Args:
            sat: The primary satellite to try first ("G16" or "G18").
            when_utc: The target UTC time.
            extra_back_minutes: Extra minutes to add to the search window for the second pass.

        Returns:
            A tuple: ((data_array, used_time, [path]), used_satellite_name).

        Raises:
            DataNotFoundError: If no data is found after all attempts.
        """
        primary, secondary = sat, "G18" if sat == "G16" else "G16"

        # --- Pass 1: Standard search window ---
        logger.info("Searching GOES (primary=%s, window=%dmin)", primary, self.cfg.search_back_minutes)
        if res := self._fetch_bt_c13_once(primary, when_utc, self.cfg.search_back_minutes):
            return res, primary

        logger.info("No data from %s, trying failover satellite %s", primary, secondary)
        if res := self._fetch_bt_c13_once(secondary, when_utc, self.cfg.search_back_minutes):
            return res, secondary

        # --- Pass 2: Widened search window ---
        widen_minutes = self.cfg.search_back_minutes + extra_back_minutes
        logger.info("Widening search window to %d minutes and retrying both satellites", widen_minutes)
        if res := self._fetch_bt_c13_once(primary, when_utc, widen_minutes):
            return res, primary
        if res := self._fetch_bt_c13_once(secondary, when_utc, widen_minutes):
            return res, secondary

        meta = CloudMeta(satellite=primary, product="CMIPF-C13", time_utc=when_utc, src_paths=[])
        raise DataNotFoundError("GOES CMIPF C13 data not found after all attempts", meta=meta)
