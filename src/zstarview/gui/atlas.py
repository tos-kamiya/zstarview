"""White-background Atlas application entry point."""

from __future__ import annotations

from ..paths import APP_DISPLAY_NAME, ATLAS_THEME_PRESET
from . import viewer


ATLAS_DISPLAY_NAME = f"{APP_DISPLAY_NAME} Atlas"

_ATLAS_DEFAULTS: dict[str, object] = {
    "presentation_id": "instrument",
    "star_data_policy": "positional_static",
    "theme": ATLAS_THEME_PRESET,
    "sky_opacity": 0.0,
    "sky_disc_style": "smooth",
    "sky_disc_altaz_rings": "off",
    "sky_disc_altaz_rings_hover": "off",
    "edge_fov_deg": 90.0,
    "content_fov_deg": 110.0,
    "cloud_opacity": 0.15,
    "cloud_stripe": ("width", 50, 0.85),
    "vmag_limit": 4.0,
    "show_dso_initial": True,
    "show_asterisms_initial": True,
    "show_guidelines_initial": True,
    "ground_tint_opacity": 0.08,
    "night_light_opacity": 0.0,
    "ridge_glow_opacity": 0.0,
    "observation_info": "bottom",
    "light_background_star_outline": True,
}


def apply_atlas_profile(args: object) -> None:
    """Apply Atlas defaults without overriding explicit CLI options."""
    explicit_options = set(getattr(args, "_explicit_options", set()) or set())
    if "presentation_id" not in explicit_options:
        setattr(args, "presentation_id", "instrument")
    if "star_data_policy" not in explicit_options:
        setattr(args, "star_data_policy", "positional_static")
    for key, value in _ATLAS_DEFAULTS.items():
        if key in explicit_options:
            continue
        setattr(args, key, value)


def main() -> None:
    """Run the white-background Atlas viewer."""
    viewer.main(
        apply_app_profile=apply_atlas_profile,
        app_name=ATLAS_DISPLAY_NAME,
        parser_default_overrides=_ATLAS_DEFAULTS,
        parser_description="White-background Atlas sky map",
        include_scenic_arguments=False,
        vmag_limit_max=6.0,
    )
