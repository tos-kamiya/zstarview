# -*- coding: utf-8 -*-
"""
Data provider for GOES-R series satellites (GOES-19 and GOES-18).

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

from ..config import CloudDiscConfig
from ..types import DataNotFoundError, CloudMeta
from ._goes_abi import load_cmi_with_area
from ._s3_io import download_s3_object, list_s3_keys
from .select import GOES_SATELLITES, normalize_satellite_name

# --- Constants ---
_GOES_BUCKET = {"G19": "noaa-goes19", "G18": "noaa-goes18"}
_GOES_REGION = {"noaa-goes19": "us-east-1", "noaa-goes18": "us-west-2"}
_GOES_BUCKET_TO_SATELLITE = {bucket: sat for sat, bucket in _GOES_BUCKET.items()}

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
            satellite=_GOES_BUCKET_TO_SATELLITE[bucket],
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
            logger.debug("Using cached file: %s", dst)
            return dst

        s3 = self._s3(bucket)

        logger.info("Downloading s3://%s/%s", bucket, key)
        return download_s3_object(
            s3_client=s3,
            bucket=bucket,
            key=key,
            dst=dst,
            satellite=_GOES_BUCKET_TO_SATELLITE[bucket],
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
                da = load_cmi_with_area(path)
                used_time = search_time.replace(minute=(search_time.minute // 10) * 10, second=0, microsecond=0, tzinfo=dt.timezone.utc)
                return da, used_time, [path]
            except Exception as e:
                logger.error("Direct GOES CMI load failed for %s: %s", Path(key).name, e)
        return None

    def fetch_bt_c13_with_failover(
        self,
        sat: str,
        when_utc: dt.datetime,
        extra_back_minutes: int = 30,
        allowed_sats: Optional[Tuple[str, ...]] = None,
    ) -> Tuple[Tuple[xr.DataArray, dt.datetime, List[Path]], str]:
        """
        Fetches C13 data with a two-pass failover strategy.

        Pass 1: Tries the primary satellite, then the secondary, with the standard search window.
        Pass 2: If Pass 1 fails, it widens the search window and retries both satellites.

        Args:
            sat: The primary satellite to try first ("G19" or "G18").
            when_utc: The target UTC time.
            extra_back_minutes: Extra minutes to add to the search window for the second pass.

        Returns:
            A tuple: ((data_array, used_time, [path]), used_satellite_name).

        Raises:
            DataNotFoundError: If no data is found after all attempts.
        """
        sat = normalize_satellite_name(sat)
        if allowed_sats is None:
            allowed: Tuple[str, ...] = GOES_SATELLITES
        else:
            allowed = tuple(
                normalized
                for normalized in (normalize_satellite_name(s) for s in allowed_sats)
                if normalized in GOES_SATELLITES
            )

        if sat in allowed:
            order = [sat] + [s for s in allowed if s != sat]
        else:
            order = list(allowed)
        if not order:
            meta = CloudMeta(satellite=sat, product="CMIPF-C13", time_utc=when_utc, src_paths=[])
            raise DataNotFoundError("No allowed GOES satellites for this location", meta=meta)

        # --- Pass 1: Standard search window ---
        logger.info("Searching GOES (order=%s, window=%dmin)", ",".join(order), self.cfg.search_back_minutes)
        for sat_name in order:
            if res := self._fetch_bt_c13_once(sat_name, when_utc, self.cfg.search_back_minutes):
                return res, sat_name

        # --- Pass 2: Widened search window ---
        widen_minutes = self.cfg.search_back_minutes + extra_back_minutes
        logger.info("Widening search window to %d minutes and retrying order=%s", widen_minutes, ",".join(order))
        for sat_name in order:
            if res := self._fetch_bt_c13_once(sat_name, when_utc, widen_minutes):
                return res, sat_name

        meta = CloudMeta(satellite=order[0], product="CMIPF-C13", time_utc=when_utc, src_paths=[])
        raise DataNotFoundError("GOES CMIPF C13 data not found after all attempts", meta=meta)
