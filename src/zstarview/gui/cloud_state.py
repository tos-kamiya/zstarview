# -*- coding: utf-8 -*-
"""
Cloud image state holder.

Keeps the latest rendered cloud disc image, metadata, and status banner text
in one place. This allows the window to delegate state bookkeeping and keep UI
code focused on orchestration.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Any, Optional

import numpy as np

from ..clouddisc.types import RenderKey, SourceKey


@dataclass
class CloudImageState:
    image: Optional[np.ndarray] = None
    missing_mask: Optional[np.ndarray] = None
    cloud_amount_field: Optional[Any] = None
    meta: Optional[Any] = None
    banner_text: Optional[str] = None
    current_satellite: Optional[str] = None
    source_refreshed_at_utc: Optional[datetime] = None
    last_az: Optional[float] = None
    last_time_utc: Optional[datetime] = None
    source_key: Optional[SourceKey] = None
    render_key: Optional[RenderKey] = None
    request_id: Optional[int] = None
    coverage_ratio: Optional[float] = None
    missing_mask_key: Optional[int] = None

    def set_result(
        self,
        image: np.ndarray,
        meta: Optional[Any],
        *,
        az: float,
        time_utc: datetime,
        cloud_amount_field: Optional[Any] = None,
        missing_mask: Optional[np.ndarray] = None,
        source_key: Optional[SourceKey] = None,
        render_key: Optional[RenderKey] = None,
        request_id: Optional[int] = None,
        coverage_ratio: Optional[float] = None,
        missing_mask_key: Optional[int] = None,
    ) -> None:
        self.image = image
        self.missing_mask = missing_mask
        self.cloud_amount_field = cloud_amount_field
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
        self.banner_text = None

    def set_source_ready(
        self,
        *,
        refreshed_at_utc: datetime,
        satellite: str | None = None,
        source_key: Optional[SourceKey] = None,
        banner_text: str | None = None,
    ) -> None:
        self.source_refreshed_at_utc = refreshed_at_utc
        if satellite:
            self.current_satellite = str(satellite)
        if source_key is not None:
            self.source_key = source_key
        if banner_text is not None:
            self.banner_text = banner_text
        else:
            self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
