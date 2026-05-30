from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Optional

from ..tropical_cyclones.models import TropicalCycloneSnapshot


@dataclass
class TropicalCycloneState:
    snapshot: Optional[TropicalCycloneSnapshot] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    cached_at_utc: Optional[datetime] = None
    last_checked_utc: Optional[datetime] = None
    next_check_utc: Optional[datetime] = None
    next_refresh_utc: Optional[datetime] = None
    projection_next_refresh_utc: Optional[datetime] = None
    source_url: Optional[str] = None
    current_source: Optional[str] = None

    def set_result(
        self,
        snapshot: TropicalCycloneSnapshot,
        *,
        cached_at_utc: datetime | None = None,
        last_checked_utc: datetime | None = None,
        next_check_utc: datetime | None = None,
        next_refresh_utc: datetime | None = None,
        banner_text: str | None = None,
    ) -> None:
        self.snapshot = snapshot
        self.cached_at_utc = cached_at_utc
        self.last_checked_utc = last_checked_utc
        self.next_check_utc = next_check_utc
        self.next_refresh_utc = next_refresh_utc
        self.projection_next_refresh_utc = None
        self.source_url = snapshot.source_url or None
        self.current_source = snapshot.service_name or None
        self.failed_this_session = False
        self.banner_text = banner_text

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True
