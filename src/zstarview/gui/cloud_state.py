"""
Cloud image state holder.

Keeps the latest rendered cloud disc image, metadata, and status banner text
in one place. This allows the window to delegate state bookkeeping and keep UI
code focused on orchestration.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Any

import numpy as np

from ..clouddisc.types import RenderKey, SourceKey


@dataclass
class CloudImageState:
    image: np.ndarray | None = None
    missing_mask: np.ndarray | None = None
    cloud_amount_field: Any | None = None
    altaz_grid: Any | None = None
    meta: Any | None = None
    banner_text: str | None = None
    current_satellite: str | None = None
    source_refreshed_at_utc: datetime | None = None
    last_az: float | None = None
    last_time_utc: datetime | None = None
    source_key: SourceKey | None = None
    render_key: RenderKey | None = None
    request_id: int | None = None
    coverage_ratio: float | None = None
    missing_mask_key: int | None = None
    source_expected_count: int | None = None
    source_available_count: int | None = None
    source_completeness_ratio: float | None = None

    def set_result(
        self,
        image: np.ndarray,
        meta: Any | None,
        *,
        az: float,
        time_utc: datetime,
        cloud_amount_field: Any | None = None,
        altaz_grid: Any | None = None,
        missing_mask: np.ndarray | None = None,
        source_key: SourceKey | None = None,
        render_key: RenderKey | None = None,
        request_id: int | None = None,
        coverage_ratio: float | None = None,
        missing_mask_key: int | None = None,
        source_expected_count: int | None = None,
        source_available_count: int | None = None,
        source_completeness_ratio: float | None = None,
    ) -> None:
        self.image = image
        self.missing_mask = missing_mask
        self.cloud_amount_field = cloud_amount_field
        self.altaz_grid = altaz_grid
        self.meta = meta
        sat = getattr(meta, "satellite", None) if meta is not None else None
        if sat:
            self.current_satellite = str(sat)
        self.last_az = az
        self.last_time_utc = time_utc
        self.source_key = source_key
        self.render_key = render_key
        self.request_id = request_id
        self.coverage_ratio = coverage_ratio
        self.missing_mask_key = missing_mask_key
        self.source_expected_count = source_expected_count
        self.source_available_count = source_available_count
        self.source_completeness_ratio = source_completeness_ratio
        self.banner_text = None

    def set_source_ready(
        self,
        *,
        refreshed_at_utc: datetime,
        satellite: str | None = None,
        source_key: SourceKey | None = None,
        banner_text: str | None = None,
        altaz_grid: Any | None = None,
    ) -> None:
        self.source_refreshed_at_utc = refreshed_at_utc
        if satellite:
            self.current_satellite = str(satellite)
        if source_key is not None:
            self.source_key = source_key
        if altaz_grid is not None:
            self.altaz_grid = altaz_grid
        if banner_text is not None:
            self.banner_text = banner_text
        else:
            self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
