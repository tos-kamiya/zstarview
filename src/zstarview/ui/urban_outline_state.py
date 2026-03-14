from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class UrbanOutlineState:
    outlines: Optional[list[list[tuple[float, float]]]] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    current_source: Optional[str] = None

    def set_result(
        self,
        outlines: list[list[tuple[float, float]]] | None,
        *,
        source: str,
    ) -> None:
        self.outlines = outlines
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def clear_outlines(self) -> None:
        self.outlines = None
