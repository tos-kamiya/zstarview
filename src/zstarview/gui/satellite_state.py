from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional

from ..satellites.types import SatelliteOmmRecord, SatelliteOverlayPoint


@dataclass
class SatelliteState:
    records_by_group: dict[str, list[SatelliteOmmRecord]] = field(default_factory=dict)
    overlay_points: Optional[list[SatelliteOverlayPoint]] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    element_epoch_utc: Optional[datetime] = None

    def set_result(
        self,
        records_by_group: dict[str, list[SatelliteOmmRecord]],
        *,
        overlay_points: list[SatelliteOverlayPoint] | None = None,
        element_epoch_utc: datetime | None = None,
    ) -> None:
        self.records_by_group = dict(records_by_group)
        self.overlay_points = overlay_points
        self.element_epoch_utc = element_epoch_utc
        self.failed_this_session = False
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def set_banner(self, text: str) -> None:
        self.banner_text = text
