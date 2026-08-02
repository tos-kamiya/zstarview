from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

from ..tropical_cyclones.models import (
    TropicalCycloneSnapshot,
    TropicalCycloneSnapshotCollection,
)


@dataclass
class TropicalCycloneState:
    snapshots: tuple[TropicalCycloneSnapshot, ...] = ()
    snapshot_collection: TropicalCycloneSnapshotCollection | None = None
    banner_text: str | None = None
    failed_this_session: bool = False
    cached_at_utc: datetime | None = None
    last_checked_utc: datetime | None = None
    next_check_utc: datetime | None = None
    next_refresh_utc: datetime | None = None
    projection_next_refresh_utc: datetime | None = None
    source_url: str | None = None
    current_source: str | None = None

    def set_result(
        self,
        snapshots: tuple[TropicalCycloneSnapshot, ...],
        *,
        cached_at_utc: datetime | None = None,
        last_checked_utc: datetime | None = None,
        next_check_utc: datetime | None = None,
        next_refresh_utc: datetime | None = None,
        banner_text: str | None = None,
        source_url: str | None = None,
        current_source: str | None = None,
    ) -> None:
        self.snapshots = tuple(snapshots)
        if source_url is None and self.snapshots:
            source_url = self.snapshots[0].source_url or None
        if current_source is None and self.snapshots:
            current_source = self.snapshots[0].service_name or None
        self.snapshot_collection = TropicalCycloneSnapshotCollection(
            snapshots=self.snapshots,
            source_url=source_url or "",
            service_name=current_source or "",
            refreshed_at_utc=cached_at_utc,
        )
        self.cached_at_utc = cached_at_utc
        self.last_checked_utc = last_checked_utc
        self.next_check_utc = next_check_utc
        self.next_refresh_utc = next_refresh_utc
        self.projection_next_refresh_utc = None
        self.source_url = source_url
        self.current_source = current_source
        self.failed_this_session = False
        self.banner_text = banner_text

    def set_error_banner(self, text: str) -> None:
        self.snapshots = ()
        self.snapshot_collection = None
        self.cached_at_utc = None
        self.last_checked_utc = None
        self.next_check_utc = None
        self.next_refresh_utc = None
        self.projection_next_refresh_utc = None
        self.source_url = None
        self.current_source = None
        self.banner_text = text
        self.failed_this_session = True
