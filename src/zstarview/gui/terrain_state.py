from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class TerrainHorizonState:
    profile_altaz: list[tuple[float, float]] | None = None
    profile_distances_m: list[float] | None = None
    secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None
    secondary_ridges_distances_m_layers: list[list[float]] | None = None
    sample_distances_m: np.ndarray | None = None
    sample_terrain_elevation_m: np.ndarray | None = None
    # Retained ground elevation from the latest successful terrain update.
    ground_elevation_m: float | None = None
    banner_text: str | None = None
    failed_this_session: bool = False
    current_source: str | None = None

    def set_result(
        self,
        profile_altaz: list[tuple[float, float]],
        *,
        profile_distances_m: list[float] | None = None,
        secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None,
        secondary_ridges_distances_m_layers: list[list[float]] | None = None,
        sample_distances_m: np.ndarray | None = None,
        sample_terrain_elevation_m: np.ndarray | None = None,
        source: str,
    ) -> None:
        self.profile_altaz = profile_altaz
        self.profile_distances_m = profile_distances_m
        self.secondary_ridges_altaz_layers = secondary_ridges_altaz_layers
        self.secondary_ridges_distances_m_layers = secondary_ridges_distances_m_layers
        self.sample_distances_m = sample_distances_m
        self.sample_terrain_elevation_m = sample_terrain_elevation_m
        self.current_source = source
        self.failed_this_session = False
        self.banner_text = None

    def set_error_banner(self, text: str) -> None:
        self.banner_text = text
        self.failed_this_session = True

    def clear_profile(self) -> None:
        self.profile_altaz = None
        self.profile_distances_m = None
        self.secondary_ridges_altaz_layers = None
        self.secondary_ridges_distances_m_layers = None
        self.sample_distances_m = None
        self.sample_terrain_elevation_m = None
