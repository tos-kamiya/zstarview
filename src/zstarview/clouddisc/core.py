# -*- coding: utf-8 -*-
"""
Core module for CloudDisc rendering.

This module contains the main `CloudDisc` class, which orchestrates the entire
process of fetching satellite data, processing it, and rendering a cloud image
from a specific observer's perspective.
"""

import datetime as dt
import logging
import threading
from dataclasses import replace
from typing import Any, Optional, Sequence


from .altaz_grid import build_altaz_grid
from .config import CloudDiscConfig
from .diagnostics import DiagnosticSink
from .providers.goes import GoesProvider
from .providers.hima import HimaProvider
from .providers.select import (
    GOES_SATELLITES,
    SUPPORTED_SATELLITES,
    pick_satellite,
    visible_satellites,
)
from .workers.cloud_source import CloudSourceFetchRequest, fetch_cloud_source
from .workers.constants import DEFAULT_CLOUD_SHELLS_KM
from .types import (
    CloudSourceData,
    SourceKey,
    VisibilityError,
    round_down_utc_to_slot,
)

logger = logging.getLogger(__name__)


class CloudDisc:
    """
    A class to render cloud images from satellite data.

    This class brings together all the components of the clouddisc library:
    - Data providers (GOES, Himawari)
    - Projection logic
    - Sampling and interpolation
    - Image rendering
    """

    def __init__(self, cfg: CloudDiscConfig, **kwargs):
        """
        Initializes the CloudDisc renderer.

        Args:
            cfg: An instance of CloudDiscConfig containing the configuration.
            **kwargs: Additional configuration options to override `cfg`.
        """
        if kwargs:
            cfg = replace(cfg, **kwargs)

        self.cfg: CloudDiscConfig = cfg
        self.goes: GoesProvider = GoesProvider(cfg)
        self.hima: HimaProvider = HimaProvider(cfg)

    def _now_rounded(self) -> dt.datetime:
        """
        Gets the current UTC time rounded down to the nearest 10-minute interval.
        Satellite data is typically published in 10-minute cycles.

        Returns:
            The rounded datetime object.
        """
        t = dt.datetime.now(dt.timezone.utc).replace(second=0, microsecond=0)
        return t.replace(minute=(t.minute // 10) * 10)

    def _select_satellite(self, lat: float, lon: float) -> str:
        supported_visible = tuple(visible_satellites(lat, lon, SUPPORTED_SATELLITES))
        if not supported_visible:
            raise VisibilityError("No supported satellite for this region")
        sat = pick_satellite(
            lat,
            lon,
            priority=self.cfg.sat_priority,
        )
        logger.debug(
            "Selected satellite=%s for observer at (lat=%.2f, lon=%.2f)", sat, lat, lon
        )
        return sat

    def make_source_key(
        self,
        *,
        lat: float,
        lon: float,
        when_utc: Optional[dt.datetime] = None,
    ) -> SourceKey:
        """Build a source key for cloud data fetch/cache lookup."""
        when = (
            self._now_rounded()
            if when_utc is None
            else round_down_utc_to_slot(when_utc)
        )
        sat = self._select_satellite(lat, lon)
        provider = "GOES" if sat in GOES_SATELLITES else "HIMAWARI"
        return SourceKey(
            satellite=sat,
            provider=provider,
            timeslot_utc=when,
            sat_priority=self.cfg.sat_priority,
        )

    def fetch_source(
        self,
        *,
        lat: float,
        lon: float,
        when_utc: Optional[dt.datetime] = None,
        cloud_shells_km: Sequence[float] = DEFAULT_CLOUD_SHELLS_KM,
        abort_event: threading.Event | None = None,
        diagnostic_sink: DiagnosticSink | None = None,
    ) -> CloudSourceData:
        """Fetch cloud source data independently from camera-dependent rendering."""
        request = CloudSourceFetchRequest(
            lat=lat,
            lon=lon,
            when_utc=when_utc,
            cloud_shells_km=tuple(float(v) for v in cloud_shells_km),
        )
        return fetch_cloud_source(
            self,
            request,
            abort_event=abort_event,
            diagnostic_sink=diagnostic_sink,
        )

    def build_altaz_grid_from_source(
        self,
        *,
        source: CloudSourceData,
        lat: float,
        lon: float,
        cloud_shells_km: Sequence[float] = DEFAULT_CLOUD_SHELLS_KM,
    ) -> Any:
        """Build a camera-independent (altitude, azimuth) grid from source data."""
        return build_altaz_grid(
            source,
            lat,
            lon,
            shells_km=cloud_shells_km,
        )
