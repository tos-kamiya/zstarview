"""Himawari provider using direct ISatSS tile stitching instead of Satpy."""

from __future__ import annotations

import datetime as dt
import logging
import threading
from pathlib import Path

import numpy as np
import xarray as xr

from ..config import CloudDiscConfig
from ..types import CloudMeta, DataNotFoundError, DownloadError
from ._hima_isatss import (
    DATA_VAR,
    ObserverTileSelection,
    TileRecord,
    attach_area_to_dataset,
    extract_tile_token,
    find_matching_keys,
    format_prefix,
    load_template_from_tile,
    select_equator_band_tiles,
    select_needed_tiles,
    stitch_tiles_from_paths,
    tile_distance_km,
)
from ._s3_io import download_s3_object

logger = logging.getLogger(__name__)
FAR_CLEAR_SKY_DISTANCE_KM = 50.0
FAR_CLEAR_SKY_BT_K = 315.0


class HimaProvider:
    """Fetch Himawari ISatSS `M1C13` tiles and stitch them into one BT field."""

    def __init__(self, cfg: CloudDiscConfig):
        self.cfg = cfg
        self.root_is = cfg.cache_root() / "hima_isatss"
        self.root_is.mkdir(parents=True, exist_ok=True)

    def _download(self, bucket: str, key: str, *, abort_event: threading.Event | None = None) -> Path:
        dst = self.root_is / bucket / key
        logger.debug("Downloading s3://%s/%s", bucket, key)
        return download_s3_object(
            bucket=bucket,
            key=key,
            dst=dst,
            satellite="HIMAWARI",
            product="ISatSS-B13",
            time_utc=dt.datetime.now(dt.timezone.utc),
            validate_func=lambda path: load_template_from_tile(path, bucket=bucket),
            abort_event=abort_event,
            timeout_s=max(self.cfg.connect_timeout, self.cfg.read_timeout),
        )

    def _find_isatss(
        self,
        when_utc: dt.datetime,
        *,
        abort_event: threading.Event | None = None,
    ) -> tuple[str | None, list[str] | None, dt.datetime | None, int | None]:
        """Find available ISatSS C13 tile keys, searching backwards by slot."""
        for slot in range(0, self.cfg.search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=slot)
            try:
                bucket, keys = find_matching_keys(
                    search_time,
                    satellite="HIMAWARI",
                    product="ISatSS-B13",
                    timeout_s=max(self.cfg.connect_timeout, self.cfg.read_timeout),
                )
            except FileNotFoundError:
                continue
            template_path = self._download(bucket, keys[0], abort_event=abort_event)
            expected_tile_count = load_template_from_tile(template_path, bucket=bucket).tile_count
            logger.info(
                "Checked %s and found %d/%d ISatSS M1C13 tiles under %s",
                bucket,
                len(keys),
                expected_tile_count,
                format_prefix(search_time),
            )
            return bucket, keys, search_time, expected_tile_count
        return None, None, None, None

    def _build_observer_selection(
        self,
        *,
        bucket: str,
        keys: list[str],
        when_utc: dt.datetime,
        observer_lat: float,
        observer_lon: float,
        cloud_shell_km: float,
        azimuth_samples: int,
        margin_tiles: int,
        equator_margin_tiles: int = 0,
        abort_event: threading.Event | None = None,
    ) -> ObserverTileSelection:
        template_path = self._download(bucket, keys[0], abort_event=abort_event)
        meta = load_template_from_tile(template_path, bucket=bucket)
        render_tiles, _poly_x, _poly_y = select_needed_tiles(
            lat_deg=observer_lat,
            lon_deg=observer_lon,
            meta=meta,
            cloud_shell_km=cloud_shell_km,
            azimuth_samples=azimuth_samples,
            margin_tiles=margin_tiles,
        )
        equator_tiles, _eq_x, _eq_y = select_equator_band_tiles(
            lon_center_deg=observer_lon,
            meta=meta,
            delta_lon=60.0,
            equator_lat=0.0,
            equator_lat_half_band_deg=5.0,
            step_deg=1.0,
            margin_tiles=equator_margin_tiles,
        )
        key_map = {extract_tile_token(Path(key).name): key for key in keys}
        missing_render_tiles = [record for record in render_tiles if record.token not in key_map]
        far_missing_render_tiles = [
            record
            for record in missing_render_tiles
            if tile_distance_km(record, meta, observer_lat=observer_lat, observer_lon=observer_lon) >= FAR_CLEAR_SKY_DISTANCE_KM
        ]
        near_missing_render_tiles = [
            record
            for record in missing_render_tiles
            if tile_distance_km(record, meta, observer_lat=observer_lat, observer_lon=observer_lon) < FAR_CLEAR_SKY_DISTANCE_KM
        ]
        if near_missing_render_tiles:
            logger.info(
                "Skipping %s slot under %s: missing %d near-field ISatSS tiles for observer selection",
                bucket,
                format_prefix(when_utc),
                len(near_missing_render_tiles),
            )
            raise DataNotFoundError(
                "Required Himawari ISatSS tiles are missing for observer footprint",
                meta=CloudMeta(satellite="HIMAWARI", product="ISatSS-B13", time_utc=when_utc, src_paths=[]),
            )
        missing_equator_tiles = [record for record in equator_tiles if record.token not in key_map]
        if missing_equator_tiles:
            logger.info(
                "Equator-band warm-threshold tiles missing for %s under %s: %d tiles will reuse prior bt_warm",
                bucket,
                format_prefix(when_utc),
                len(missing_equator_tiles),
            )
        selected_tokens = {record.token for record in render_tiles if record.token in key_map} | {
            record.token for record in equator_tiles if record.token in key_map
        }
        selected_keys = [key_map[token] for token in sorted(selected_tokens)]
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
        return ObserverTileSelection(
            selected_keys=selected_keys,
            render_tiles=render_tiles,
            equator_tiles=equator_tiles,
            far_missing_render_tiles=far_missing_render_tiles,
            near_missing_render_tiles=near_missing_render_tiles,
            missing_equator_tiles=missing_equator_tiles,
        )

    def _select_keys_for_observer(
        self,
        *,
        bucket: str,
        keys: list[str],
        when_utc: dt.datetime,
        observer_lat: float,
        observer_lon: float,
        cloud_shell_km: float,
        azimuth_samples: int,
        margin_tiles: int,
        equator_margin_tiles: int = 0,
    ) -> list[str]:
        selection = self._build_observer_selection(
            bucket=bucket,
            keys=keys,
            when_utc=when_utc,
            observer_lat=observer_lat,
            observer_lon=observer_lon,
            cloud_shell_km=cloud_shell_km,
            azimuth_samples=azimuth_samples,
            margin_tiles=margin_tiles,
            equator_margin_tiles=equator_margin_tiles,
            abort_event=None,
        )
        return selection.selected_keys

    def _stitch_local_paths(
        self,
        paths: list[Path],
        *,
        source_label: str,
        observer_lat: float | None,
        observer_lon: float | None,
        far_missing_render_tiles: list[TileRecord] | None = None,
    ) -> xr.DataArray:
        stitched = stitch_tiles_from_paths(paths, source_label=source_label)
        if far_missing_render_tiles:
            tile_height = int(stitched.attrs["stitched_tile_height"])
            tile_width = int(stitched.attrs["stitched_tile_width"])
            data = stitched[DATA_VAR].values
            coverage = stitched["tile_coverage_mask"].values
            filled_tiles = 0
            for tile in far_missing_render_tiles:
                row0 = int(tile.row_offset)
                col0 = int(tile.col_offset)
                row1 = row0 + tile_height
                col1 = col0 + tile_width
                data[row0:row1, col0:col1] = FAR_CLEAR_SKY_BT_K
                coverage[row0:row1, col0:col1] = 1
                filled_tiles += 1
            stitched["tile_coverage_mask"] = (("y", "x"), coverage)
            stitched[DATA_VAR].attrs["far_clear_fill_tile_count"] = filled_tiles
            stitched[DATA_VAR].attrs["far_clear_fill_value_k"] = FAR_CLEAR_SKY_BT_K
            stitched.attrs["far_clear_fill_tile_count"] = filled_tiles
            stitched.attrs["far_clear_fill_value_k"] = FAR_CLEAR_SKY_BT_K
            stitched[DATA_VAR].attrs["coverage_fraction"] = float(coverage.mean())
        prepared = attach_area_to_dataset(stitched)
        da = prepared[DATA_VAR].astype(np.float32)
        da.attrs = dict(da.attrs)
        da.attrs["source_key_count"] = len(paths)
        da.attrs["coverage_fraction"] = 1.0 if far_missing_render_tiles else float(prepared["tile_coverage_mask"].values.mean())
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
        abort_event: threading.Event | None = None,
    ) -> tuple[xr.DataArray, dt.datetime, list[Path]]:
        """
        Fetch Himawari ISatSS C13 brightness temperature data.

        If `observer_lat/lon` are provided, only tiles overlapping the observer-visible
        horizon footprint are downloaded and stitched. Otherwise all tiles for the
        chosen timeslot are used.
        """
        logger.info("Searching for Himawari ISatSS M1C13 data...")
        found = self._find_isatss(when_utc, abort_event=abort_event)
        if len(found) == 4:
            bucket, keys, used_time, expected_tile_count = found
        else:
            bucket, keys, used_time = found
            expected_tile_count = len(keys) if keys is not None else None
        if not bucket or not keys or used_time is None or expected_tile_count is None:
            meta = CloudMeta(satellite="HIMAWARI", product="ISatSS-B13", time_utc=when_utc, src_paths=[])
            raise DataNotFoundError("Himawari ISatSS B13 data not found in search window", meta=meta)

        if (observer_lat is None) ^ (observer_lon is None):
            raise ValueError("observer_lat and observer_lon must be provided together")

        selection = ObserverTileSelection(selected_keys=keys, render_tiles=[], equator_tiles=[], far_missing_render_tiles=[], near_missing_render_tiles=[], missing_equator_tiles=[])
        if observer_lat is not None and observer_lon is not None:
            logger.info("Himawari tile download started...")
            selection = self._build_observer_selection(
                bucket=bucket,
                keys=keys,
                when_utc=used_time,
                observer_lat=float(observer_lat),
                observer_lon=float(observer_lon),
                cloud_shell_km=float(cloud_shell_km),
                azimuth_samples=int(azimuth_samples),
                margin_tiles=int(margin_tiles),
                equator_margin_tiles=int(equator_margin_tiles),
                abort_event=abort_event,
            )

        try:
            paths = [self._download(bucket, key, abort_event=abort_event) for key in selection.selected_keys]
            da = self._stitch_local_paths(
                paths,
                source_label=f"s3://{bucket}/{format_prefix(used_time)}",
                observer_lat=observer_lat,
                observer_lon=observer_lon,
                far_missing_render_tiles=selection.far_missing_render_tiles,
            )
            if selection.missing_equator_tiles:
                da.attrs["equator_band_missing"] = True
            da.attrs["source_bucket"] = bucket
            da.attrs["source_expected_count"] = int(expected_tile_count)
            da.attrs["source_available_count"] = len(keys)
            da.attrs["source_completeness_ratio"] = float(len(keys)) / float(expected_tile_count)
            rounded = used_time.replace(minute=(used_time.minute // 10) * 10, tzinfo=dt.timezone.utc)
            logger.info("Himawari tile download ready.")
            return da, rounded, paths
        except DownloadError:
            raise
        except Exception as e:
            meta = CloudMeta(satellite="HIMAWARI", product="ISatSS-B13", time_utc=used_time, src_paths=[])
            raise DownloadError("Failed to fetch or stitch Himawari ISatSS B13 data", meta=meta) from e
