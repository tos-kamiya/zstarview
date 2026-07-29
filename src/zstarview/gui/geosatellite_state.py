from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Any

import numpy as np

from ..clouddisc.types import RenderKey, SourceKey


@dataclass
class GeoSatelliteState:
    image: np.ndarray | None = None
    missing_mask: np.ndarray | None = None
    cloud_amount_field: Any | None = None
    altaz_grid: Any | None = None
    meta: Any | None = None
    banner_text: str | None = None
    current_source: str | None = None
    source_refreshed_at_utc: datetime | None = None
    captured_at_utc: datetime | None = None
    fetched_at_utc: datetime | None = None
    last_az: float | None = None
    last_time_utc: datetime | None = None
    source_key: SourceKey | None = None
    render_key: RenderKey | None = None
    request_id: int | None = None
    coverage_ratio: float | None = None
    missing_mask_key: int | None = None

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
        current_source: str | None = None,
        captured_at_utc: datetime | None = None,
        fetched_at_utc: datetime | None = None,
    ) -> None:
        self.image = image
        self.missing_mask = missing_mask
        self.cloud_amount_field = cloud_amount_field
        self.altaz_grid = altaz_grid
        self.meta = meta
        self.current_source = current_source or "Geo-sat"
        self.captured_at_utc = captured_at_utc
        self.fetched_at_utc = fetched_at_utc
        self.last_az = az
        self.last_time_utc = time_utc
        self.source_key = source_key
        self.render_key = render_key
        self.request_id = request_id
        self.coverage_ratio = coverage_ratio
        self.missing_mask_key = missing_mask_key
        self.banner_text = None

    def set_source_ready(
        self,
        *,
        refreshed_at_utc: datetime,
        banner_text: str | None = None,
        current_source: str | None = None,
    ) -> None:
        self.source_refreshed_at_utc = refreshed_at_utc
        if current_source:
            self.current_source = current_source
        if banner_text is not None:
            self.banner_text = banner_text
        else:
            self.banner_text = None

    def set_banner(self, text: str) -> None:
        self.banner_text = text

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
