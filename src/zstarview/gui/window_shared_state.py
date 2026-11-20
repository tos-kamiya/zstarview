"""Attribute contract shared by the SkyWindow mixins.

The SkyWindow behavior is split across several mixins that read and write
state owned by ``SkyWindowCoreMixin.__init__``. This class only declares the
types of that shared state so each mixin can be type-checked on its own.
It intentionally holds annotations only: no values are assigned here, and the
runtime behavior of the composed window is unchanged.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..types import ViewerData
    from .window_state import SkyWindowState


class SkyWindowSharedState:
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
    _startup_initial_data_loaded: bool
