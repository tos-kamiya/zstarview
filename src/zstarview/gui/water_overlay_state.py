from __future__ import annotations

from dataclasses import dataclass

from ..water_overlay import WaterOverlayPoint, WaterOverlayPolyline


@dataclass
class WaterOverlayState:
    sea_level_dots: list[WaterOverlayPoint] | None = None
    inland_dots: list[WaterOverlayPoint] | None = None
    dem_dots: list[WaterOverlayPoint] | None = None
    dots: list[WaterOverlayPoint] | None = None
    polylines: list[WaterOverlayPolyline] | None = None
    banner_text: str | None = None
    failed_this_session: bool = False
    current_source: str | None = None
    current_mode: str | None = None

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
        elif use_dem and self.sea_level_dots is not None or not use_dem:
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

    def _combined_dots(self) -> list[WaterOverlayPoint] | None:
        if self.sea_level_dots is None and self.inland_dots is None:
            return None
        combined: list[WaterOverlayPoint] = []
        if self.sea_level_dots is not None:
            combined.extend(self.sea_level_dots)
        if self.inland_dots is not None:
            combined.extend(self.inland_dots)
        return combined
