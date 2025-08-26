"""
Data provider for GOES-R series satellites (GOES-16 and GOES-18).

This module fetches Cloud and Moisture Imagery Product (CMIPF) data from the
NOAA Open Data Dissemination (NODD) program on AWS S3.
It specifically targets channel 13 (longwave infrared) to get brightness temperatures.

Note: This module uses print() for logging. In a production environment,
it would be better to use Python's `logging` module.
"""

import datetime as dt
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import boto3
import xarray as xr
from botocore import UNSIGNED
from botocore.config import Config
from satpy import Scene

from ..config import CloudDiscConfig
from ..types import DataNotFoundError

# S3 bucket names for GOES satellites
_GOES_BUCKET = {"G16": "noaa-goes16", "G18": "noaa-goes18"}
# Corresponding AWS regions for the buckets
_GOES_REGION = {"noaa-goes16": "us-east-1", "noaa-goes18": "us-west-2"}


def _doy(dt_utc: dt.datetime) -> int:
    """Calculates the day of the year (1-366)."""
    return dt_utc.timetuple().tm_yday


class GoesProvider:
    """
    A provider for fetching GOES-R ABI L2 CMIPF data from AWS S3.
    """

    def __init__(self, cfg: CloudDiscConfig):
        """
        Initializes the GoesProvider.

        Args:
            cfg: The application configuration.
        """
        self.cfg = cfg
        self.root = cfg.cache_root() / "goes_cmipf"
        self.root.mkdir(parents=True, exist_ok=True)
        self._list_cache: Dict[Tuple, Any] = {}
        # TODO: Implement a cleanup mechanism for _list_cache (e.g., every 24 hours).

    def _s3(self, bucket: str) -> boto3.client:
        """
        Creates an anonymous S3 client for the specified bucket.

        Args:
            bucket: The name of the S3 bucket.

        Returns:
            A boto3 S3 client instance.
        """
        return boto3.client("s3", region_name=_GOES_REGION[bucket], config=Config(signature_version=UNSIGNED, retries={"max_attempts": 3}))

    def _list_hour(self, bucket: str, t: dt.datetime) -> List[str]:
        """
        Lists all object keys for a given hour in the S3 bucket.
        Results are cached in memory.

        Args:
            bucket: The S3 bucket name.
            t: The UTC datetime.

        Returns:
            A list of S3 object keys.
        """
        cache_key = (bucket, t.year, _doy(t), t.hour)
        if cache_key in self._list_cache:
            return self._list_cache[cache_key]

        prefix = f"ABI-L2-CMIPF/{t.year:04d}/{_doy(t):03d}/{t.hour:02d}/"
        s3 = self._s3(bucket)
        print(f"[DEBUG] Listing s3://{bucket}/{prefix} (region: {s3.meta.region_name})")

        paginator = s3.get_paginator("list_objects_v2")
        pages = paginator.paginate(Bucket=bucket, Prefix=prefix)

        keys = [obj["Key"] for page in pages for obj in page.get("Contents", []) or []]

        print(f"[DEBUG] Found {len(keys)} objects.")
        self._list_cache[cache_key] = keys
        return keys

    def _download(self, bucket: str, key: str) -> Path:
        """
        Downloads a file from S3, caching it locally.

        Args:
            bucket: The S3 bucket name.
            key: The S3 object key.

        Returns:
            The local path to the downloaded file.
        """
        dst = self.root / bucket / key
        if dst.exists():
            return dst

        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = dst.with_suffix(dst.suffix + ".tmp")
        s3 = self._s3(bucket)

        print(f"[INFO] Downloading s3://{bucket}/{key}")
        with tmp_path.open("wb") as f:
            s3.download_fileobj(bucket, key, f)
        tmp_path.replace(dst)
        return dst

    def _fetch_bt_c13_once(self, sat: str, when_utc: dt.datetime, search_back_minutes: int) -> Optional[Tuple[xr.DataArray, dt.datetime, List[Path]]]:
        """
        Searches for and loads a single C13 brightness temperature file for a given satellite and time window.

        Args:
            sat: The satellite to use ("G16" or "G18").
            when_utc: The target UTC time.
            search_back_minutes: The number of minutes to search backward from the target time.

        Returns:
            A tuple (data_array, used_time, [path]) if found, otherwise None.
        """
        bucket = _GOES_BUCKET[sat]

        # Iterate backwards from the target time in 10-minute intervals
        for mback in range(0, search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=mback)
            keys = self._list_hour(bucket, search_time)
            if not keys:
                continue

            # Filter for Channel 13 (C13) files
            keys_c13 = [k for k in keys if ("-M6C13_" in k or "-C13_" in k)]
            if not keys_c13:
                continue

            # The latest file in the hour is the most recent one
            keys_c13.sort()
            key = keys_c13[-1]
            print(f"[DEBUG] Found candidate C13: {Path(key).name}")

            path = self._download(bucket, key)

            try:
                # Try loading with Satpy, which handles projection info correctly
                scn = Scene(reader="abi_l2_nc", filenames=[str(path)])
                scn.load(["C13"])
                da = scn["C13"].astype("float32").compute()
                used_time = search_time.replace(minute=(search_time.minute // 10) * 10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
                return da, used_time, [path]
            except Exception as e:
                print(f"[WARN] Satpy load failed for {Path(key).name}: {e}. Falling back to xarray.")
                # Fallback: Try opening with xarray directly. The variable name is often 'CMI'.
                try:
                    with xr.open_dataset(path, engine="netcdf4") as ds:
                        if "CMI" in ds.variables:
                            da = ds["CMI"].astype("float32").compute()
                            used_time = search_time.replace(minute=(search_time.minute // 10) * 10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
                            return da, used_time, [path]
                except Exception as ex:
                    print(f"[ERROR] xarray fallback also failed for {Path(key).name}: {ex}")
                    # If both fail, continue to the next candidate time

        return None  # Not found in the given time window

    def fetch_bt_c13_with_failover(
        self, sat: str, when_utc: dt.datetime, extra_back_minutes: int = 30
    ) -> Tuple[Tuple[xr.DataArray, dt.datetime, List[Path]], str]:
        """
        Fetches C13 data with failover and retry logic.

        It first tries the primary satellite. If that fails, it tries the secondary one.
        If both fail, it widens the search window by `extra_back_minutes` and retries both.

        Args:
            sat: The primary satellite to try first ("G16" or "G18").
            when_utc: The target UTC time.
            extra_back_minutes: Extra minutes to add to the search window on the second pass.

        Returns:
            A tuple containing ((data_array, used_time, [path]), used_satellite_name).

        Raises:
            DataNotFoundError: If no data can be found after all attempts.
        """
        primary = sat
        secondary = "G18" if sat == "G16" else "G16"

        # 1st pass: Use the standard search window
        print(f"[INFO] Searching for GOES data (primary: {primary}, window: {self.cfg.search_back_minutes} min)")
        res = self._fetch_bt_c13_once(primary, when_utc, self.cfg.search_back_minutes)
        if res:
            return res, primary

        print(f"[INFO] No data from {primary}, trying failover satellite {secondary}.")
        res = self._fetch_bt_c13_once(secondary, when_utc, self.cfg.search_back_minutes)
        if res:
            return res, secondary

        # 2nd pass: Widen the search window and retry
        widen_minutes = self.cfg.search_back_minutes + extra_back_minutes
        print(f"[INFO] Widening search window to {widen_minutes} minutes and retrying.")

        res = self._fetch_bt_c13_once(primary, when_utc, widen_minutes)
        if res:
            return res, primary

        res = self._fetch_bt_c13_once(secondary, when_utc, widen_minutes)
        if res:
            return res, secondary

        raise DataNotFoundError("GOES CMIPF C13 not found (after failover and widened window)")
