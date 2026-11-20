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
    from datetime import timedelta
    from typing import Any

    import numpy as np
    from PySide6.QtCore import QTimer, SignalInstance
    from PySide6.QtGui import QAction

    from ..astro import DeepSkyCatalogArrays, StarCatalogArrays
    from ..clouddisc.core import CloudDisc
    from ..paths import ThemeStyle
    from ..types import StarCatalogMeta, ViewerData
    from .aircraft_controller import AircraftController
    from .aircraft_state import AircraftState
    from .application_services import ApplicationServices
    from .cloud_controller import CloudController
    from .cloud_state import CloudImageState
    from .composite import SkyCompositorCache
    from .geosatellite_controller import GeoSatelliteController
    from .geosatellite_state import GeoSatelliteState
    from .jpl_small_body_controller import JplSmallBodyController
    from .satellite_controller import SatelliteController
    from .satellite_state import SatelliteState
    from .sky_worker import SkyDataWorker
    from .terrain_controller import TerrainHorizonController
    from .terrain_state import TerrainHorizonState
    from .tropical_cyclone_controller import TropicalCycloneController
    from .tropical_cyclone_state import TropicalCycloneState
    from .urban_outline_controller import UrbanOutlineController
    from .urban_outline_state import UrbanOutlineState
    from .water_overlay_controller import WaterOverlayController
    from .water_overlay_state import WaterOverlayState
    from .window_inputs import SkyWindowRuntimeOptions
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
    _startup_initial_terrain_loaded: bool
    _startup_initial_urban_loaded: bool
    _aircraft_gui_allowed: bool
    _aircraft_toggle_supported: bool
    _akari_ir_bands_toggle_supported: bool
    _cloud_gui_allowed: bool
    _earth_guide_gui_allowed: bool
    _night_light_toggle_supported: bool
    _precipitation_toggle_supported: bool
    _road_light_toggle_supported: bool
    _satellite_gui_allowed: bool
    _satellite_toggle_supported: bool
    _sky_disc_gui_allowed: bool
    _terrain_horizon_gui_allowed: bool
    _tropical_cyclone_toggle_supported: bool
    _urban_outline_gui_allowed: bool

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

    # Background data controllers.
    _aircraft_controller: AircraftController | None
    _cloud_controller: CloudController | None
    _geosatellite_controller: GeoSatelliteController | None
    _jpl_small_body_controller: JplSmallBodyController | None
    _satellite_controller: SatelliteController | None
    _terrain_horizon_controller: TerrainHorizonController | None
    _tropical_cyclone_controller: TropicalCycloneController | None
    _urban_outline_controller: UrbanOutlineController | None
    _water_overlay_controller: WaterOverlayController | None

    # Other window services, catalogs, and settings.
    _is_shutting_down: bool
    _sky_worker: SkyDataWorker
    sky_update_interval: int
    vmag_limit: float
    _compositor: SkyCompositorCache
    _viewport_rotation_keys_down: set[int]
    star_catalog_lod6_indices: np.ndarray
    star_catalog_np: StarCatalogArrays
    _disc_generation: int
    _ephemeris: Any
    _services: ApplicationServices
    delta_t: timedelta
    dso_catalog_np: DeepSkyCatalogArrays | None
    star_catalog_meta: StarCatalogMeta
    theme: ThemeStyle
    _clouddisc: CloudDisc | None
    content_fov_deg: float
    runtime_options: SkyWindowRuntimeOptions
    _frameless_window: bool
    observation_info_mode: str
    observation_info_pinned: bool

    # Timers.
    _hover_repaint_timer: QTimer
    _interaction_idle_timer: QTimer
    _viewport_interaction_idle_timer: QTimer

    # Qt signals declared on the window class.
    initial_data_loaded: SignalInstance
    cloud_repaint_requested: SignalInstance

    # Opacities restored when a layer is toggled back on.
    _aircraft_opacity_when_enabled: float
    _akari_ir_bands_opacity_when_enabled: float
    _cloud_alpha_when_enabled: float
    _earth_guide_opacity_when_enabled: float
    _meteor_opacity_when_enabled: float
    _night_light_opacity_when_enabled: float
    _precipitation_opacity_when_enabled: float
    _road_light_opacity_when_enabled: float
    _satellite_opacity_when_enabled: float
    _sky_disc_alpha_when_enabled: float
    _terrain_horizon_opacity_when_enabled: float
    _urban_outline_opacity_when_enabled: float
    _water_overlay_opacity_when_enabled: float
