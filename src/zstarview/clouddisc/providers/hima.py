# -*- coding: utf-8 -*-
"""Himawari provider using direct ISatSS tile stitching instead of Satpy."""

from __future__ import annotations

import datetime as dt
import logging
from pathlib import Path
from typing import List, Optional, Tuple

import boto3
import numpy as np
import xarray as xr
from botocore import UNSIGNED
from botocore.config import Config

from ..config import CloudDiscConfig
from ..types import CloudMeta, DataNotFoundError, DownloadError
from ._hima_isatss import (
    DATA_VAR,
    attach_area_to_dataset,
    extract_tile_token,
    find_matching_keys,
    format_prefix,
    load_template_from_tile,
    select_equator_tiles,
    select_needed_tiles,
    stitch_tiles_from_paths,
)
from ._s3_io import download_s3_object


logger = logging.getLogger(__name__)


class HimaProvider:
    """Fetch Himawari ISatSS `M1C13` tiles and stitch them into one BT field."""

    def __init__(self, cfg: CloudDiscConfig):
        self.cfg = cfg
        self.root_is = cfg.cache_root() / "hima_isatss"
        self.root_is.mkdir(parents=True, exist_ok=True)

    def _s3(self):
        s3_cfg = Config(
            signature_version=UNSIGNED,
            retries={"max_attempts": 1},
            connect_timeout=self.cfg.connect_timeout,
            read_timeout=self.cfg.read_timeout,
        )
        return boto3.client("s3", config=s3_cfg)

    def _download(self, bucket: str, key: str) -> Path:
        dst = self.root_is / bucket / key
        if dst.exists():
            logger.debug("Using cached file: %s", dst)
            return dst
        logger.debug("Downloading s3://%s/%s", bucket, key)
        return download_s3_object(
            s3_client=self._s3(),
            bucket=bucket,
            key=key,
            dst=dst,
            satellite="HIMAWARI",
            product="ISatSS-B13",
            time_utc=dt.datetime.now(dt.timezone.utc),
        )

    def _find_isatss(self, when_utc: dt.datetime) -> Tuple[Optional[str], Optional[List[str]], Optional[dt.datetime]]:
        """Find available ISatSS C13 tile keys, searching backwards by slot."""
        s3_client = self._s3()
        for slot in range(0, self.cfg.search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=slot)
            try:
                bucket, keys = find_matching_keys(
                    s3_client,
                    search_time,
                    satellite="HIMAWARI",
                    product="ISatSS-B13",
                )
            except FileNotFoundError:
                continue
            template_path = self._download(bucket, keys[0])
            expected_tile_count = load_template_from_tile(template_path, bucket=bucket).tile_count
            if len(keys) < expected_tile_count:
                logger.info(
                    "Skipping incomplete %s slot under %s: found %d/%d ISatSS M1C13 tiles",
                    bucket,
                    format_prefix(search_time),
                    len(keys),
                    expected_tile_count,
                )
                continue
            logger.info(
                "Checked %s and found %d ISatSS M1C13 tiles under %s",
                bucket,
                len(keys),
                format_prefix(search_time),
            )
            return bucket, keys, search_time
        return None, None, None

    def _select_keys_for_observer(
        self,
        *,
        bucket: str,
        keys: List[str],
        when_utc: dt.datetime,
        observer_lat: float,
        observer_lon: float,
        cloud_shell_km: float,
        azimuth_samples: int,
        margin_tiles: int,
        equator_margin_tiles: int = 0,
    ) -> List[str]:
        template_path = self._download(bucket, keys[0])
        meta = load_template_from_tile(template_path, bucket=bucket)
        render_tiles, _poly_x, _poly_y = select_needed_tiles(
            lat_deg=observer_lat,
            lon_deg=observer_lon,
            meta=meta,
            cloud_shell_km=cloud_shell_km,
            azimuth_samples=azimuth_samples,
            margin_tiles=margin_tiles,
        )
        equator_tiles, _eq_x, _eq_y = select_equator_tiles(
            lon_center_deg=observer_lon,
            meta=meta,
            delta_lon=60.0,
            equator_lat=0.0,
            step_deg=1.0,
            margin_tiles=equator_margin_tiles,
        )
        selected_tokens = {record.token for record in render_tiles} | {record.token for record in equator_tiles}
        key_map = {extract_tile_token(Path(key).name): key for key in keys}
        selected_keys = [key_map[token] for token in sorted(selected_tokens) if token in key_map]
        if not selected_keys:
            raise DataNotFoundError(
                "No Himawari ISatSS tiles selected for observer footprint",
                meta=CloudMeta(satellite="HIMAWARI", product="ISatSS-B13", time_utc=when_utc, src_paths=[]),
            )
        logger.info(
            "Selected %d/%d Himawari tiles for observer lat=%.3f lon=%.3f (render=%d equator=%d)",
            len(selected_keys),
            len(keys),
            observer_lat,
            observer_lon,
            len(render_tiles),
            len(equator_tiles),
        )
        return selected_keys

    def _stitch_local_paths(
        self,
        paths: List[Path],
        *,
        source_label: str,
        observer_lat: float | None,
        observer_lon: float | None,
    ) -> xr.DataArray:
        stitched = stitch_tiles_from_paths(paths, source_label=source_label)
        prepared = attach_area_to_dataset(stitched)
        da = prepared[DATA_VAR].astype(np.float32)
        da.attrs = dict(da.attrs)
        da.attrs["source_key_count"] = len(paths)
        if observer_lat is not None and observer_lon is not None:
            da.attrs["observer_lat"] = float(observer_lat)
            da.attrs["observer_lon"] = float(observer_lon)
        return da

    def fetch_bt_c13(
        self,
        when_utc: dt.datetime,
        *,
        observer_lat: float | None = None,
        observer_lon: float | None = None,
        cloud_shell_km: float = 6376.0,
        azimuth_samples: int = 1440,
        margin_tiles: int = 1,
        equator_margin_tiles: int = 0,
    ) -> Tuple[xr.DataArray, dt.datetime, List[Path]]:
        """
        Fetch Himawari ISatSS C13 brightness temperature data.

        If `observer_lat/lon` are provided, only tiles overlapping the observer-visible
        horizon footprint are downloaded and stitched. Otherwise all tiles for the
        chosen timeslot are used.
        """
        logger.info("Searching for Himawari ISatSS M1C13 data...")
        bucket, keys, used_time = self._find_isatss(when_utc)
        if not bucket or not keys or used_time is None:
            meta = CloudMeta(satellite="HIMAWARI", product="ISatSS-B13", time_utc=when_utc, src_paths=[])
            raise DataNotFoundError("Himawari ISatSS B13 data not found in search window", meta=meta)

        if (observer_lat is None) ^ (observer_lon is None):
            raise ValueError("observer_lat and observer_lon must be provided together")

        selected_keys = keys
        if observer_lat is not None and observer_lon is not None:
            selected_keys = self._select_keys_for_observer(
                bucket=bucket,
                keys=keys,
                when_utc=used_time,
                observer_lat=float(observer_lat),
                observer_lon=float(observer_lon),
                cloud_shell_km=float(cloud_shell_km),
                azimuth_samples=int(azimuth_samples),
                margin_tiles=int(margin_tiles),
                equator_margin_tiles=int(equator_margin_tiles),
            )

        try:
            paths = [self._download(bucket, key) for key in selected_keys]
            da = self._stitch_local_paths(
                paths,
                source_label=f"s3://{bucket}/{format_prefix(used_time)}",
                observer_lat=observer_lat,
                observer_lon=observer_lon,
            )
            da.attrs["source_bucket"] = bucket
            rounded = used_time.replace(minute=(used_time.minute // 10) * 10, tzinfo=dt.timezone.utc)
            return da, rounded, paths
        except DownloadError:
            raise
        except Exception as e:
            meta = CloudMeta(satellite="HIMAWARI", product="ISatSS-B13", time_utc=used_time, src_paths=[])
            raise DownloadError("Failed to fetch or stitch Himawari ISatSS B13 data", meta=meta) from e
