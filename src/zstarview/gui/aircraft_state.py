from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Optional

from ..aircraft.types import AircraftOverlayPoint, AircraftSnapshot
from ..aircraft.opensky import AircraftBoundingBox


@dataclass
class AircraftState:
    snapshots: Optional[list[AircraftSnapshot]] = None
    overlay_points: Optional[list[AircraftOverlayPoint]] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    last_success_utc: Optional[datetime] = None
    last_bbox: Optional[AircraftBoundingBox] = None

    def set_result(
        self,
        snapshots: list[AircraftSnapshot],
        *,
        overlay_points: list[AircraftOverlayPoint] | None = None,
        bbox: AircraftBoundingBox | None = None,
        refreshed_at_utc: datetime | None = None,
    ) -> None:
        self.snapshots = snapshots
        self.overlay_points = overlay_points
        self.last_bbox = bbox
        self.last_success_utc = refreshed_at_utc
        self.failed_this_session = False
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def set_banner(self, text: str) -> None:
        self.banner_text = text
