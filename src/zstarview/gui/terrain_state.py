from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class TerrainHorizonState:
    profile_altaz: Optional[list[tuple[float, float]]] = None
    profile_distances_m: Optional[list[float]] = None
    secondary_profile_altaz_layers: Optional[list[list[tuple[float, float]]]] = None
    secondary_profile_distances_m_layers: Optional[list[list[float]]] = None
    banner_text: Optional[str] = None
    failed_this_session: bool = False
    current_source: Optional[str] = None

    def set_result(
        self,
        profile_altaz: list[tuple[float, float]],
        *,
        profile_distances_m: Optional[list[float]] = None,
        secondary_profile_altaz_layers: Optional[list[list[tuple[float, float]]]] = None,
        secondary_profile_distances_m_layers: Optional[list[list[float]]] = None,
        source: str,
    ) -> None:
        self.profile_altaz = profile_altaz
        self.profile_distances_m = profile_distances_m
        self.secondary_profile_altaz_layers = secondary_profile_altaz_layers
        self.secondary_profile_distances_m_layers = secondary_profile_distances_m_layers
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def clear_profile(self) -> None:
        self.profile_altaz = None
        self.profile_distances_m = None
        self.secondary_profile_altaz_layers = None
        self.secondary_profile_distances_m_layers = None
