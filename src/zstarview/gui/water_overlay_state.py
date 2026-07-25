from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..water_overlay import WaterOverlayPoint, WaterOverlayPolyline


@dataclass
class WaterOverlayState:
    sea_level_dots: Optional[list[WaterOverlayPoint]] = None
    inland_dots: Optional[list[WaterOverlayPoint]] = None
    dem_dots: Optional[list[WaterOverlayPoint]] = None
    dots: Optional[list[WaterOverlayPoint]] = None
    polylines: Optional[list[WaterOverlayPolyline]] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    current_source: Optional[str] = None
    current_mode: Optional[str] = None

    def set_sea_level_dots_result(
        self,
        dots: list[WaterOverlayPoint] | None,
        *,
        source: str,
    ) -> None:
        self.sea_level_dots = dots
        self.inland_dots = None
        self.current_mode = "sea"
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def set_dem_dots_result(
        self,
        dots: list[WaterOverlayPoint] | None,
        *,
        source: str,
    ) -> None:
        self.dem_dots = dots
        self.inland_dots = None
        self.current_mode = "dem"
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def select_active_dots(self, *, use_dem: bool) -> None:
        if use_dem and self.dem_dots is not None:
            self.dots = self.dem_dots
            self.current_mode = "dem"
        elif use_dem and self.sea_level_dots is not None:
            self.dots = self._combined_dots()
            self.current_mode = "sea"
        elif not use_dem:
            self.dots = self._combined_dots()
            self.current_mode = "sea"
        else:
            self.dots = None
            self.current_mode = None

    def set_polylines(self, polylines: list[WaterOverlayPolyline] | None) -> None:
        self.polylines = polylines

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def _combined_dots(self) -> Optional[list[WaterOverlayPoint]]:
        if self.sea_level_dots is None and self.inland_dots is None:
            return None
        combined: list[WaterOverlayPoint] = []
        if self.sea_level_dots is not None:
            combined.extend(self.sea_level_dots)
        if self.inland_dots is not None:
            combined.extend(self.inland_dots)
        return combined
