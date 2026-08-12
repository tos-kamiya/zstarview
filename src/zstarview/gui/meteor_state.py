from __future__ import annotations

from dataclasses import dataclass

from ..meteors.types import MeteorWindowResult


@dataclass
class MeteorState:
    result: MeteorWindowResult | None = None
    banner_text: str | None = None

    def set_result(self, result: MeteorWindowResult) -> None:
        self.result = result
        self.banner_text = None

    def set_banner(self, text: str) -> None:
        self.banner_text = text
