from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Optional

from ..aircraft.opensky import AircraftBoundingBox
from ..aircraft.types import AircraftSnapshot


@dataclass
class AircraftState:
    snapshots: Optional[list[AircraftSnapshot]] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    last_success_utc: Optional[datetime] = None
    last_bbox: Optional[AircraftBoundingBox] = None

    def set_result(
        self,
        snapshots: list[AircraftSnapshot],
        *,
        bbox: AircraftBoundingBox | None = None,
        refreshed_at_utc: datetime | None = None,
    ) -> None:
        self.snapshots = snapshots
        self.last_bbox = bbox
        self.last_success_utc = refreshed_at_utc
        self.failed_this_session = False
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def set_banner(self, text: str) -> None:
        self.banner_text = text
