from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..water_overlay import WaterOverlayPoint


@dataclass
class WaterOverlayState:
    points: Optional[list[WaterOverlayPoint]] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    current_source: Optional[str] = None

    def set_result(
        self,
        points: list[WaterOverlayPoint] | None,
        *,
        source: str,
    ) -> None:
        self.points = points
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def clear_points(self) -> None:
        self.points = None
