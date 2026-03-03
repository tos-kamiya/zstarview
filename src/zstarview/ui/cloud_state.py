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

from PySide6.QtGui import QImage


@dataclass
class CloudImageState:
    image: Optional[QImage] = None
    stripe_density: Optional[Any] = None
    meta: Optional[Any] = None
    banner_text: Optional[str] = None
    current_satellite: Optional[str] = None
    last_az: Optional[float] = None
    last_time_utc: Optional[datetime] = None

    def set_result(
        self,
        image: QImage,
        meta: Optional[Any],
        *,
        az: float,
        time_utc: datetime,
        stripe_density: Optional[Any] = None,
    ) -> None:
        self.image = image
        self.stripe_density = stripe_density
        self.meta = meta
        sat = getattr(meta, "satellite", None) if meta is not None else None
        if sat:
            self.current_satellite = str(sat)
        self.last_az = az
        self.last_time_utc = time_utc
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
