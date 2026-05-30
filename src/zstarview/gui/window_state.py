from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Optional

from PySide6.QtCore import QPoint
from PySide6.QtGui import QImage

from ..night_lights import NightLightGlowProfile
from ..search.models import SearchJumpTarget
from ..types import CelestialData, StarsTable, UrbanOutlinePolyline, ViewCenterAltAz
from ..water_overlay import WaterOverlayPoint


@dataclass
class SkyWindowState:
    """Mutable UI/runtime state kept separate from the window shell."""

    render_view_center: ViewCenterAltAz
    rotation_step: float = 5.0
    interaction_idle_ms: int = 300
    interaction_mode: bool = False
    viewport_interaction_idle_ms: int = 700
    viewport_interaction_mode: bool = False
    viewport_interaction_stars: Optional[StarsTable] = None
    mouse_pos: Optional[QPoint] = None
    overlay_info_bottom_left: bool = False
    jump_highlight_name: Optional[str] = None
    jump_highlight_altaz: Optional[ViewCenterAltAz] = None
    jump_highlight_until_ms: float = 0.0
    persistent_search_target: Optional[SearchJumpTarget] = None
    persistent_search_reference_time_utc: Optional[datetime] = None
    persistent_search_next_refresh_utc: Optional[datetime] = None
    persistent_search_last_refresh_utc: Optional[datetime] = None
    persistent_search_last_error: Optional[str] = None
    sky_next_refresh_utc: Optional[datetime] = None
    cloud_next_refresh_utc: Optional[datetime] = None
    cloud_projection_next_refresh_utc: Optional[datetime] = None
    satellite_next_refresh_utc: Optional[datetime] = None
    aircraft_next_refresh_utc: Optional[datetime] = None
    satellite_projection_next_refresh_utc: Optional[datetime] = None
    aircraft_projection_next_refresh_utc: Optional[datetime] = None
    tropical_cyclone_projection_next_refresh_utc: Optional[datetime] = None
    sky_update_pending: bool = False
    pending_star_vmag_limit: Optional[float] = None
    cloud_repaint_deferred: bool = False
    viewport_interaction_release_pending: bool = False
    last_star_render_stats: Optional[tuple[int, int, int, int]] = None
    celestial_data: Optional[CelestialData] = None
    sky_disc_base_size: int = 1024
    sky_disc_image: Optional[QImage] = None
    cloud_base_size: int = 256
    terrain_horizon_profile: Optional[list[tuple[float, float]]] = None
    terrain_horizon_profile_distances_m: Optional[list[float]] = None
    terrain_secondary_ridges_altaz_layers: Optional[list[list[tuple[float, float]]]] = None
    terrain_secondary_ridges_distances_m_layers: Optional[list[list[float]]] = None
    urban_outlines: Optional[list[UrbanOutlinePolyline]] = None
    water_overlay_dots: Optional[list[WaterOverlayPoint]] = None
    night_light_glow_profile: Optional[NightLightGlowProfile] = None
