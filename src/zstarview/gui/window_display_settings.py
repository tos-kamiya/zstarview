"""Mutable display settings owned by the sky window."""

from __future__ import annotations

from dataclasses import dataclass

from ..paths import THEME_STYLES_BY_PRESET, ThemeStyle


@dataclass
class SkyWindowDisplaySettings:
    """User-facing display values that can change during a GUI session.

    The window exposes compatibility properties for the old flat attribute
    names while this object remains the single storage location.
    """

    visual_preset: str = "night"
    presentation_id: str = "scenic"
    theme: ThemeStyle = THEME_STYLES_BY_PRESET["night"]
    display_tone_curve: tuple[float, ...] | None = None
    show_dso: bool = False
    show_asterisms: bool = True
    show_guidelines: bool = True
    show_observation_info: bool = True
    moon_style: str = "marker"
    moon_scale: int = 1
    bright_bodies_mode: str = "standard"
    star_base_radius: float = 1.0
    star_visibility_boost: float = 1.0
    asterism_visibility_boost: float = 1.0
    earth_guide_visibility_boost: float = 1.0
    vmag_limit: float = 6.0
    sky_disc_altaz_rings: str = "off"
    sky_disc_altaz_rings_hover: str = "off"
    cloud_disc_alpha: float = 0.0
    satellite_opacity: float = 0.0
    terrain_horizon_opacity: float = 0.0
    earth_guide_opacity: float = 0.0
    night_light_opacity: float = 0.0
    diffuse_sky_source: str = "gaia"
    akari_ir_bands_opacity: float = 0.0
    ridge_glow_opacity: float = 0.0
    urban_outline_opacity: float = 0.0
    show_urban_outline_layer: bool = False
    inverted_city_enabled: bool = False
    water_overlay_opacity: float = 0.0
    road_night_lights_opacity: float = 0.0
    precipitation_opacity: float = 0.0
    aircraft_opacity: float = 0.0
    meteor_opacity: float = 0.0
    tropical_cyclone_opacity: float = 0.0
    show_tropical_cyclone_overlay: bool = False
    star_render_expected_width: int = 600
    light_background_star_outline: bool = False
    asterism_opacity: float | None = None
    sky_disc_alpha: float = 0.15


DISPLAY_SETTING_NAMES = frozenset(SkyWindowDisplaySettings.__dataclass_fields__)
