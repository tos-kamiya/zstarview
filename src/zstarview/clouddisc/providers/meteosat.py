# -*- coding: utf-8 -*-
"""
Data provider scaffold for Meteosat/SEVIRI brightness temperature.

This provider is intentionally isolated from runtime routing until integration
is completed in a follow-up change. The class exposes the same contract used by
other providers:

    fetch_bt_c13(when_utc) -> (DataArray[K], used_time_utc, [src_paths])

The actual source lookup/downloading backend is left behind a private hook so
it can be implemented without changing call sites.
"""

from __future__ import annotations

import datetime as dt
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import xarray as xr

from ..config import CloudDiscConfig
from ..types import CloudMeta, DataNotFoundError, DownloadError, RenderError

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class MeteosatFetchResult:
    """A normalized result container for a single fetch attempt."""

    da: xr.DataArray
    used_time: dt.datetime
    src_paths: List[Path]
    product: str


class MeteosatProvider:
    """
    Provider entry point for Meteosat IR BT fetching.

    Notes:
    - This PR2 implementation establishes interface/flow and exception behavior.
    - Backend discovery/download is intentionally left as a hook (`_fetch_slot`)
      and returns no data by default.
    """

    def __init__(self, cfg: CloudDiscConfig):
        self.cfg = cfg
        self.root = cfg.cache_root() / "meteosat"
        self.root.mkdir(parents=True, exist_ok=True)

    def _fetch_slot(self, search_time_utc: dt.datetime) -> Optional[MeteosatFetchResult]:
        """
        Backend hook for fetching one slot of Meteosat BT data.

        Return:
            MeteosatFetchResult if data was found and decoded.
            None if data is not available for this slot.

        Raises:
            DownloadError / RenderError for hard failures.
        """
        del search_time_utc
        return None

    def fetch_bt_c13(self, when_utc: dt.datetime) -> Tuple[xr.DataArray, dt.datetime, List[Path]]:
        """
        Fetch Meteosat IR BT (B13-equivalent) with search-back behavior.

        This mirrors the behavior style used by existing providers.
        """
        for mback in range(0, self.cfg.search_back_minutes + 1, 10):
            search_time = when_utc - dt.timedelta(minutes=mback)
            try:
                result = self._fetch_slot(search_time)
            except (DownloadError, RenderError):
                raise
            except Exception as e:
                meta = CloudMeta(
                    satellite="METEOSAT",
                    product="SEVIRI-IR",
                    time_utc=search_time,
                    src_paths=[],
                )
                raise RenderError("Unexpected Meteosat decode/fetch failure", meta=meta) from e

            if result is None:
                continue

            if result.used_time.tzinfo is None:
                used_time = result.used_time.replace(tzinfo=dt.timezone.utc)
            else:
                used_time = result.used_time.astimezone(dt.timezone.utc)
            return result.da, used_time, result.src_paths

        meta = CloudMeta(
            satellite="METEOSAT",
            product="SEVIRI-IR",
            time_utc=when_utc,
            src_paths=[],
        )
        raise DataNotFoundError("Meteosat IR data not found in search window", meta=meta)
