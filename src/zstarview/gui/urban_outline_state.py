from __future__ import annotations

from dataclasses import dataclass

from ..types import UrbanOutlinePolyline


@dataclass
class UrbanOutlineState:
    outlines: list[UrbanOutlinePolyline] | None = None
    base_outline_count: int | None = None
    skyscraper_outline_count: int | None = None
    banner_text: str | None = None
    failed_this_session: bool = False
    current_source: str | None = None

    def set_result(
        self,
        outlines: list[UrbanOutlinePolyline] | None,
        *,
        source: str,
        base_outline_count: int | None = None,
        skyscraper_outline_count: int | None = None,
    ) -> None:
        self.outlines = outlines
        self.base_outline_count = base_outline_count
        self.skyscraper_outline_count = skyscraper_outline_count
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def clear_outlines(self) -> None:
        self.outlines = None
        self.base_outline_count = None
        self.skyscraper_outline_count = None
