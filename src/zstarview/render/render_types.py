"""Shared data structures used by render presentations."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

import astropy.time
import numpy as np
from PySide6.QtCore import QPoint, QRect
from PySide6.QtGui import QFont, QImage

from ..aircraft.types import AircraftSnapshot
from ..clouddisc.altaz_grid import CloudAltAzGrid
from ..meteors.types import MeteorTrail
from ..night_lights import NightLightGlowProfile
from ..paths import (
    NIGHT_LIGHT_DEFAULT_OPACITY,
    RIDGE_GLOW_DEFAULT_OPACITY,
    ROAD_LIGHT_DEFAULT_OPACITY,
    THEME_STYLES_BY_PRESET,
    ThemeStyle,
)
from ..precipitation import PrecipitationRenderItem
from ..road_night_lights import RoadNightLightPolyline
from ..satellites.types import SatelliteOmmRecord
from ..tropical_cyclones.models import TropicalCycloneSnapshot
from ..types import (
    CelestialData,
    PlanetBody,
    ScreenGeometry,
    StarsTable,
    UrbanOutlinePolyline,
    ViewerData,
)
from ..water_overlay import WaterOverlayPoint, WaterOverlayPolyline


@dataclass(frozen=True)
class FrameContext:
    viewer: ViewerData
    time_obj: astropy.time.Time | None
    geometry: ScreenGeometry
    viewport_rect: QRect
    sky_update_interval: int = 60


@dataclass(frozen=True)
class RenderSceneData:
    viewer: ViewerData
    celestial_data: CelestialData
    sky_disc_image: QImage | None
    cloud_missing_mask: np.ndarray | None
    cloud_altaz_grid: CloudAltAzGrid | None
    terrain_horizon_profile: list[tuple[float, float]] | None
    terrain_horizon_profile_distances_m: list[float] | None
    terrain_secondary_ridges_altaz_layers: list[list[tuple[float, float]]] | None
    terrain_secondary_ridges_distances_m_layers: list[list[float]] | None
    urban_outlines: list[UrbanOutlinePolyline] | None
    satellite_element_epoch_utc: datetime | None = None
    satellite_records_by_group: dict[str, list[SatelliteOmmRecord]] | None = None
    aircraft_snapshots: list[AircraftSnapshot] | None = None
    meteor_trails: tuple[MeteorTrail, ...] | None = None
    meteor_window_end_utc: datetime | None = None
    time_obj: astropy.time.Time | None = None
    night_light_glow_profile: NightLightGlowProfile | None = None
    water_overlay_dots: list[WaterOverlayPoint] | None = None
    water_overlay_polylines: list[WaterOverlayPolyline] | None = None
    road_night_light_polylines: list[RoadNightLightPolyline] | None = None
    precipitation_columns: list[PrecipitationRenderItem] | None = None
    tropical_cyclone_snapshots: tuple[TropicalCycloneSnapshot, ...] | None = None
    dynamic_planets: list[PlanetBody] | None = None


@dataclass(frozen=True)
class RenderStyle:
    visual_preset: str
    text_font: QFont
    status_line_font: QFont
    show_background_gradient: bool
    show_custom_window_frame: bool
    show_observation_info: bool
    show_dso: bool
    show_asterisms: bool
    show_guidelines: bool
    enlarge_moon: bool
    bright_bodies_mode: str
    star_base_radius: float
    star_visibility_boost: float
    asterism_visibility_boost: float
    earth_guide_visibility_boost: float
    vmag_limit: float
    sky_disc_altaz_rings: str
    sky_disc_altaz_rings_hover: str
    cloud_disc_alpha: float
    satellite_opacity: float
    terrain_horizon_opacity: float
    earth_guide_opacity: float
    night_light_opacity: float = NIGHT_LIGHT_DEFAULT_OPACITY
    akari_ir_bands_opacity: float = 0.10
    ridge_glow_opacity: float = RIDGE_GLOW_DEFAULT_OPACITY
    urban_outline_opacity: float = 0.2
    show_urban_outline_layer: bool = True
    inverted_city_enabled: bool = False
    water_overlay_opacity: float = 0.4
    road_night_lights_opacity: float = ROAD_LIGHT_DEFAULT_OPACITY
    precipitation_opacity: float = 0.0
    aircraft_opacity: float = 0.0
    meteor_opacity: float = 0.0
    tropical_cyclone_opacity: float = 0.4
    show_tropical_cyclone_overlay: bool = True
    star_render_expected_width: int = 600
    ground_tint_opacity: float = 0.025
    theme: ThemeStyle = THEME_STYLES_BY_PRESET["night"]
    light_background_star_outline: bool = False
    sky_disc_alpha: float = 0.15
    presentation_id: str = "scenic"


@dataclass(frozen=True)
class RenderHudState:
    mouse_pos: QPoint | None
    overlay_info_bottom_left: bool
    viewport_interaction_mode: bool
    viewport_interaction_stars: StarsTable | None
    status_message: str | None
    # Exported images do not draw the observation HUD, so their marker
    # position may be chosen independently of the HUD collision-avoidance
    # position used by the window presentations.
    time_of_day_marker_bottom_left: bool | None = None
    simplified_view_enabled: bool = False
    simplified_view_labels_enabled: bool = True
