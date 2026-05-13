from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..water_overlay import WaterOverlayPoint


@dataclass
class WaterOverlayState:
    sea_level_points: Optional[list[WaterOverlayPoint]] = None
    dem_points: Optional[list[WaterOverlayPoint]] = None
    points: Optional[list[WaterOverlayPoint]] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    current_source: Optional[str] = None
    current_mode: Optional[str] = None

    def set_sea_level_result(
        self,
        points: list[WaterOverlayPoint] | None,
        *,
        source: str,
    ) -> None:
        self.sea_level_points = points
        self.current_mode = "sea"
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def set_dem_result(
        self,
        points: list[WaterOverlayPoint] | None,
        *,
        source: str,
    ) -> None:
        self.dem_points = points
        self.current_mode = "dem"
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def select_active_points(self, *, use_dem: bool) -> None:
        if use_dem and self.dem_points is not None:
            self.points = self.dem_points
            self.current_mode = "dem"
        elif use_dem and self.sea_level_points is not None:
            self.points = self.sea_level_points
            self.current_mode = "sea"
        elif not use_dem and self.sea_level_points is not None:
            self.points = self.sea_level_points
            self.current_mode = "sea"
        else:
            self.points = None
            self.current_mode = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def clear_points(self) -> None:
        self.sea_level_points = None
        self.dem_points = None
        self.points = None
