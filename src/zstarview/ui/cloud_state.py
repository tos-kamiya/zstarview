# -*- coding: utf-8 -*-
"""
Cloud image state holder.

Keeps the latest rendered cloud disc image, metadata, status banner text,
and simple cleanup ticking logic in one place. This allows the window to
delegate state bookkeeping and keep UI code focused on orchestration.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Optional

from PySide6.QtGui import QImage


@dataclass
class CloudImageState:
    image: Optional[QImage] = None
    stripe_density: Optional[Any] = None
    meta: Optional[dict] = None
    banner_text: Optional[str] = None
    last_az: Optional[float] = None
    last_time_utc: Optional[datetime] = None

    cleanup_interval: int = 10
    cleanup_counter: int = field(default=0, init=False)

    def set_result(
        self,
        image: QImage,
        meta: Optional[dict],
        *,
        az: float,
        time_utc: datetime,
        stripe_density: Optional[Any] = None,
    ) -> None:
        self.image = image
        self.stripe_density = stripe_density
        self.meta = meta
        self.last_az = az
        self.last_time_utc = time_utc
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text

    def tick_cleanup(self) -> bool:
        """Increment the cleanup counter; return True when it's time to clean."""
        run = (self.cleanup_counter % self.cleanup_interval) == 0
        self.cleanup_counter += 1
        return run
