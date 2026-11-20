"""Attribute contract shared by the SkyWindow mixins.

The SkyWindow behavior is split across several mixins that read and write
state owned by ``SkyWindowCoreMixin.__init__``. This class only declares the
types of that shared state so each mixin can be type-checked on its own.
It intentionally holds annotations only: no values are assigned here, and the
runtime behavior of the composed window is unchanged.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .window_shared_methods import SkyWindowSharedMethods

if TYPE_CHECKING:
    from PySide6.QtGui import QAction

    from ..types import ViewerData
    from .aircraft_state import AircraftState
    from .cloud_state import CloudImageState
    from .geosatellite_state import GeoSatelliteState
    from .satellite_state import SatelliteState
    from .terrain_state import TerrainHorizonState
    from .tropical_cyclone_state import TropicalCycloneState
    from .urban_outline_state import UrbanOutlineState
    from .water_overlay_state import WaterOverlayState
    from .window_state import SkyWindowState


class SkyWindowSharedState(SkyWindowSharedMethods):
    """Type declarations for state shared across SkyWindow mixins."""

    viewer_data: ViewerData
    state: SkyWindowState

    # Layer visibility toggles.
    show_asterisms: bool
    show_guidelines: bool
    show_tropical_cyclone_overlay: bool
    twinkle_enabled: bool

    # Layer opacities (0.0 disables the layer).
    sky_disc_alpha: float
    cloud_disc_alpha: float
    satellite_opacity: float
    aircraft_opacity: float
    meteor_opacity: float
    precipitation_opacity: float
    tropical_cyclone_opacity: float
    terrain_horizon_opacity: float
    earth_guide_opacity: float
    water_overlay_opacity: float
    night_light_opacity: float
    road_night_lights_opacity: float
    akari_ir_bands_opacity: float
    urban_outline_opacity: float

    # Feature availability and startup flags.
    _geo_satellite_enabled: bool
    _cloud_toggle_supported: bool
    _startup_initial_load_started: bool
    _startup_initial_data_loaded: bool
    _client_geometry_sync_done: bool

    # Landscape-mode annotation fade.
    _landscape_fade_was_active: bool
    _landscape_annotation_opacity: float

    # Menu actions. Each is None until the window menu builds it, and an
    # optional layer's action stays None when the layer is unavailable.
    _action_moon_option: QAction | None
    _action_raise_view: QAction | None
    _action_lower_view: QAction | None
    _action_square_window: QAction | None
    _action_toggle_akari_ir_bands: QAction | None
    _action_toggle_aircraft: QAction | None
    _action_toggle_asterisms: QAction | None
    _action_toggle_clouds: QAction | None
    _action_toggle_dso: QAction | None
    _action_toggle_earth_guide: QAction | None
    _action_toggle_geo_satellite: QAction | None
    _action_toggle_guidelines: QAction | None
    _action_toggle_meteors: QAction | None
    _action_toggle_night_lights: QAction | None
    _action_toggle_observation_info: QAction | None
    _action_toggle_precipitation: QAction | None
    _action_toggle_road_lights: QAction | None
    _action_toggle_satellites: QAction | None
    _action_toggle_sky_disc: QAction | None
    _action_toggle_terrain_horizon: QAction | None
    _action_toggle_tropical_cyclone: QAction | None
    _action_toggle_twinkle: QAction | None
    _action_toggle_urban_outline: QAction | None
    _action_toggle_water_overlay: QAction | None

    # Layer data states.
    aircraft_state: AircraftState
    cloud_state: CloudImageState
    geosatellite_state: GeoSatelliteState
    satellite_state: SatelliteState
    terrain_horizon_state: TerrainHorizonState
    tropical_cyclone_state: TropicalCycloneState
    urban_outline_state: UrbanOutlineState
    water_overlay_state: WaterOverlayState
