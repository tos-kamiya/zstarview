"""
Data provider for GOES-R series satellites (GOES-19 and GOES-18).

This module fetches Cloud and Moisture Imagery Product (CMIPF) data from the
NOAA Open Data Dissemination (NODD) program on AWS S3. It specifically targets
channel 13 (longwave infrared) to get brightness temperatures for cloud rendering.
"""

import datetime as dt
import logging
import re
import threading
from pathlib import Path

import xarray as xr

from ..config import CloudDiscConfig
from ..diagnostics import DiagnosticSink, emit_diagnostic
from ..types import CloudMeta, DataNotFoundError
from ._goes_abi import load_cmi_with_area
from ._s3_io import download_s3_object, list_s3_keys
from .select import GOES_SATELLITES, normalize_satellite_name

# --- Constants ---
_GOES_BUCKET = {"G19": "noaa-goes19", "G18": "noaa-goes18"}
_GOES_REGION = {"noaa-goes19": "us-east-1", "noaa-goes18": "us-west-2"}
_GOES_BUCKET_TO_SATELLITE = {bucket: sat for sat, bucket in _GOES_BUCKET.items()}

logger = logging.getLogger(__name__)
_SCAN_START_RE = re.compile(r"_s(\d{13,14})_")


def _scan_start_token(path_or_key: str | Path) -> str:
    """Return the ABI scan-start token embedded in an object name."""
    match = _SCAN_START_RE.search(Path(path_or_key).name)
    if match is None:
        raise ValueError(f"GOES filename has no scan-start token: {path_or_key}")
    return match.group(1)


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
        self._list_cache: dict[tuple[int, int, int, str], list[str]] = {}

    def _list_hour(
        self,
        bucket: str,
        t: dt.datetime,
        *,
        abort_event: threading.Event | None = None,
        diagnostic_sink: DiagnosticSink | None = None,
    ) -> list[str]:
        """
        Lists all object keys for a given hour in the S3 bucket, with in-memory caching.
        The S3 path is structured as `ABI-L2-CMIPF/YYYY/DOY/HH/`.
        """
        t_utc = t if t.tzinfo is not None else t.replace(tzinfo=dt.timezone.utc)
        t_utc = t_utc.astimezone(dt.timezone.utc)
        cache_key = (t_utc.year, _doy(t_utc), t_utc.hour, bucket)
        if self._should_cache_hour_listing(t_utc) and cache_key in self._list_cache:
            logger.debug(
                "GOES list cache hit: bucket=%s year=%04d doy=%03d hour=%02d",
                bucket,
                t_utc.year,
                _doy(t_utc),
                t_utc.hour,
            )
            return self._list_cache[cache_key]

        prefix = f"ABI-L2-CMIPF/{t_utc.year:04d}/{_doy(t_utc):03d}/{t_utc.hour:02d}/"
        logger.debug("GOES list start: s3://%s/%s", bucket, prefix)
        emit_diagnostic(
            diagnostic_sink,
            "list_s3_prefix",
            "start",
            "Listing GOES S3 prefix",
            bucket=bucket,
            prefix=prefix,
            time_utc=t_utc,
        )

        keys = list_s3_keys(
            bucket=bucket,
            prefix=prefix,
            satellite=_GOES_BUCKET_TO_SATELLITE[bucket],
            product="CMIPF-C13",
            time_utc=t_utc,
            uri_label=f"S3 bucket s3://{bucket}/{prefix}",
            abort_event=abort_event,
            timeout_s=max(self.cfg.connect_timeout, self.cfg.read_timeout),
        )

        logger.debug("GOES list done: s3://%s/%s keys=%d", bucket, prefix, len(keys))
        emit_diagnostic(
            diagnostic_sink,
            "list_s3_prefix",
            "ok",
            "GOES S3 prefix listed",
            bucket=bucket,
            prefix=prefix,
            key_count=len(keys),
        )
        if self._should_cache_hour_listing(t_utc):
            self._list_cache[cache_key] = keys
        return keys

    @staticmethod
    def _should_cache_hour_listing(t_utc: dt.datetime) -> bool:
        """
        Cache only closed UTC hours.

        The current hour can receive newly published GOES files while the app is
        running, so reusing an in-memory listing for that hour can pin the
        provider to an older scene across multiple timer refreshes.
        """
        now_utc = dt.datetime.now(dt.timezone.utc)
        return (
            t_utc.year,
            _doy(t_utc),
            t_utc.hour,
        ) != (
            now_utc.year,
            _doy(now_utc),
            now_utc.hour,
        )

    def _download(
        self,
        bucket: str,
        key: str,
        *,
        channel: int = 13,
        abort_event: threading.Event | None = None,
        diagnostic_sink: DiagnosticSink | None = None,
    ) -> Path:
        """Downloads a file from S3, caching it locally using an atomic write."""
        dst = self.root / bucket / key

        def validate(path: Path) -> None:
            load_cmi_with_area(path, diagnostic_sink=diagnostic_sink)

        logger.debug("GOES download start: s3://%s/%s", bucket, key)
        emit_diagnostic(
            diagnostic_sink,
            "download_source",
            "start",
            "Downloading GOES source file",
            bucket=bucket,
            key=key,
            destination=dst,
        )
        path = download_s3_object(
            bucket=bucket,
            key=key,
            dst=dst,
            satellite=_GOES_BUCKET_TO_SATELLITE[bucket],
            product=f"CMIPF-C{channel:02d}",
            time_utc=dt.datetime.now(dt.timezone.utc),
            validate_func=validate,
            abort_event=abort_event,
            timeout_s=max(self.cfg.connect_timeout, self.cfg.read_timeout),
        )
        emit_diagnostic(
            diagnostic_sink,
            "download_source",
            "ok",
            "GOES source file ready",
            bucket=bucket,
            key=key,
            path=path,
        )
        return path

    def fetch_bt_c16_for_c13(
        self,
        sat: str,
        c13_time_utc: dt.datetime,
        c13_path: Path,
        *,
        abort_event: threading.Event | None = None,
        diagnostic_sink: DiagnosticSink | None = None,
    ) -> tuple[xr.DataArray, dt.datetime, list[Path]]:
        """Fetch C16 from exactly the same ABI scan as an existing C13 file.

        This diagnostic-stage API deliberately performs no backward search.  A
        missing companion must not be silently paired with another scan.
        """
        sat = normalize_satellite_name(sat)
        bucket = _GOES_BUCKET[sat]
        scan_start = _scan_start_token(c13_path)
        keys = self._list_hour(
            bucket,
            c13_time_utc,
            abort_event=abort_event,
            diagnostic_sink=diagnostic_sink,
        )
        candidates = [
            key
            for key in keys
            if ("-M6C16_" in key or "-C16_" in key)
            and _SCAN_START_RE.search(Path(key).name)
            and _scan_start_token(key) == scan_start
        ]
        if not candidates:
            meta = CloudMeta(
                satellite=sat,
                product="CMIPF-C16",
                time_utc=c13_time_utc,
                src_paths=[],
            )
            raise DataNotFoundError(
                "GOES CMIPF C16 companion not found for the selected C13 scan",
                meta=meta,
            )
        key = sorted(candidates)[-1]
        path = self._download(
            bucket,
            key,
            channel=16,
            abort_event=abort_event,
            diagnostic_sink=diagnostic_sink,
        )
        return load_cmi_with_area(path, diagnostic_sink=diagnostic_sink), c13_time_utc, [path]

    def _fetch_bt_c13_once(
        self,
        sat: str,
        when_utc: dt.datetime,
        search_back_minutes: int,
        *,
        abort_event: threading.Event | None = None,
        diagnostic_sink: DiagnosticSink | None = None,
    ) -> tuple[xr.DataArray, dt.datetime, list[Path]] | None:
        """
        Searches for and loads a single C13 brightness temp file for a given satellite and time.
        """
        bucket = _GOES_BUCKET[sat]
        logger.debug(
            "GOES fetch attempt start: sat=%s bucket=%s when=%s search_back=%dmin",
            sat,
            bucket,
            when_utc.isoformat(),
            search_back_minutes,
        )
        emit_diagnostic(
            diagnostic_sink,
            "select_product",
            "start",
            "Searching GOES C13 product",
            satellite=sat,
            bucket=bucket,
            when_utc=when_utc,
            search_back_minutes=search_back_minutes,
        )

        # Iterate backwards from the target time to find the most recent available file.
        for mback in range(0, search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=mback)
            keys = self._list_hour(
                bucket,
                search_time,
                abort_event=abort_event,
                diagnostic_sink=diagnostic_sink,
            )
            if not keys:
                continue

            # Filter for Channel 13 (C13) files and get the latest one.
            keys_c13 = sorted([k for k in keys if "-M6C13_" in k or "-C13_" in k])
            if not keys_c13:
                continue
            key = keys_c13[-1]
            logger.debug(
                "GOES candidate selected: sat=%s bucket=%s key=%s search_time=%s",
                sat,
                bucket,
                Path(key).name,
                search_time.isoformat(),
            )
            emit_diagnostic(
                diagnostic_sink,
                "select_product",
                "ok",
                "Selected GOES C13 product",
                satellite=sat,
                bucket=bucket,
                key=key,
                search_time_utc=search_time,
            )

            path = self._download(
                bucket,
                key,
                abort_event=abort_event,
                diagnostic_sink=diagnostic_sink,
            )

            try:
                logger.debug("GOES load start: %s", path)
                da = load_cmi_with_area(path, diagnostic_sink=diagnostic_sink)
                logger.debug("GOES load done: %s", path)
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
        allowed_sats: tuple[str, ...] | None = None,
        abort_event: threading.Event | None = None,
        diagnostic_sink: DiagnosticSink | None = None,
    ) -> tuple[tuple[xr.DataArray, dt.datetime, list[Path]], str]:
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
            allowed: tuple[str, ...] = GOES_SATELLITES
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
        emit_diagnostic(
            diagnostic_sink,
            "resolve_source",
            "info",
            "Searching GOES source",
            order=order,
            window_minutes=self.cfg.search_back_minutes,
        )
        for sat_name in order:
            logger.debug("GOES pass 1 satellite start: %s", sat_name)
            if res := self._fetch_bt_c13_once(
                sat_name,
                when_utc,
                self.cfg.search_back_minutes,
                abort_event=abort_event,
                diagnostic_sink=diagnostic_sink,
            ):
                return res, sat_name

        # --- Pass 2: Widened search window ---
        widen_minutes = self.cfg.search_back_minutes + extra_back_minutes
        logger.info("Widening search window to %d minutes and retrying order=%s", widen_minutes, ",".join(order))
        for sat_name in order:
            logger.debug("GOES pass 2 satellite start: %s", sat_name)
            if res := self._fetch_bt_c13_once(
                sat_name,
                when_utc,
                widen_minutes,
                abort_event=abort_event,
                diagnostic_sink=diagnostic_sink,
            ):
                return res, sat_name

        meta = CloudMeta(satellite=order[0], product="CMIPF-C13", time_utc=when_utc, src_paths=[])
        emit_diagnostic(
            diagnostic_sink,
            "select_product",
            "failed",
            "GOES C13 data not found after all attempts",
            order=order,
            when_utc=when_utc,
        )
        raise DataNotFoundError("GOES CMIPF C13 data not found after all attempts", meta=meta)
