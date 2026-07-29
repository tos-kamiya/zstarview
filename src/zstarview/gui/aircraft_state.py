from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

from ..aircraft.opensky import AircraftBoundingBox
from ..aircraft.types import AircraftSnapshot


@dataclass
class AircraftState:
    snapshots: list[AircraftSnapshot] | None = None
    banner_text: str | None = None
    failed_this_session: bool = False
    last_success_utc: datetime | None = None
    last_bbox: AircraftBoundingBox | None = None

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
