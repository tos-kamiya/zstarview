from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

import astropy.time
from PySide6.QtCore import QPoint
from PySide6.QtGui import QImage

from ..moon_hover import MoonHoverImage
from ..night_lights import NightLightGlowProfile
from ..precipitation import PrecipitationRenderItem
from ..road_night_lights import RoadNightLightPolyline
from ..search.models import SearchJumpTarget
from ..types import (
    CelestialData,
    PlanetBody,
    StarsTable,
    UrbanOutlinePolyline,
    ViewCenterAltAz,
)
from ..water_overlay import WaterOverlayPoint


@dataclass
class SkyWindowState:
    """Mutable UI/runtime state kept separate from the window shell."""

    render_view_center: ViewCenterAltAz
    rotation_step: float = 5.0
    interaction_idle_ms: int = 300
    interaction_mode: bool = False
    viewport_interaction_idle_ms: int = 1000
    viewport_interaction_mode: bool = False
    viewport_interaction_stars: StarsTable | None = None
    twinkle_bucket: int | None = None
    twinkle_targets: tuple[tuple[int, float], ...] = ()
    simplified_view_enabled: bool = False
    simplified_view_labels_enabled: bool = True
    default_display_mode: str = "normal"
    current_display_mode: str = "normal"
    mouse_pos: QPoint | None = None
    overlay_info_bottom_left: bool = False
    jump_highlight_name: str | None = None
    jump_highlight_altaz: ViewCenterAltAz | None = None
    jump_highlight_until_ms: float = 0.0
    persistent_search_target: SearchJumpTarget | None = None
    persistent_search_reference_time_utc: datetime | None = None
    persistent_search_next_refresh_utc: datetime | None = None
    persistent_search_last_refresh_utc: datetime | None = None
    persistent_search_last_error: str | None = None
    sky_next_refresh_utc: datetime | None = None
    sky_disc_next_refresh_utc: datetime | None = None
    cloud_next_refresh_utc: datetime | None = None
    cloud_projection_next_refresh_utc: datetime | None = None
    satellite_next_refresh_utc: datetime | None = None
    aircraft_next_refresh_utc: datetime | None = None
    meteor_next_refresh_utc: datetime | None = None
    satellite_projection_next_refresh_utc: datetime | None = None
    aircraft_projection_next_refresh_utc: datetime | None = None
    tropical_cyclone_projection_next_refresh_utc: datetime | None = None
    sky_update_pending: bool = False
    pending_star_vmag_limit: float | None = None
    cloud_repaint_deferred: bool = False
    viewport_interaction_release_pending: bool = False
    viewport_interaction_completion_reason: str | None = None
    last_star_render_stats: tuple[int, int, int, int] | None = None
    celestial_data: CelestialData | None = None
    dynamic_display_time: astropy.time.Time | None = None
    dynamic_display_second: int | None = None
    dynamic_planets: list[PlanetBody] | None = None
    dynamic_planet_bucket: int | None = None
    prepared_dynamic_planets: list[PlanetBody] | None = None
    prepared_dynamic_planet_bucket: int | None = None
    dynamic_planet_requested_bucket: int | None = None
    moon_hover_image_key: datetime | None = None
    moon_hover_image: MoonHoverImage | None = None
    sky_disc_base_size: int = 1024
    sky_disc_image: QImage | None = None
    cloud_base_size: int = 256
    terrain_horizon_profile: list[tuple[float, float]] | None = None
    terrain_horizon_profile_distances_m: list[float] | None = None
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None = None
    terrain_secondary_ridges_distances_m_layers: list[list[float]] | None = None
    urban_outlines: list[UrbanOutlinePolyline] | None = None
    water_overlay_dots: list[WaterOverlayPoint] | None = None
    road_night_light_polylines: list[RoadNightLightPolyline] | None = None
    precipitation_columns: list[PrecipitationRenderItem] | None = None
    precipitation_next_refresh_utc: datetime | None = None
    night_light_glow_profile: NightLightGlowProfile | None = None
