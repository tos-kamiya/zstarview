"""White-background object-viewer application entry point."""

from __future__ import annotations

from ..paths import APP_DISPLAY_NAME, OBJECT_VIEWER_THEME_PRESET
from . import viewer


OBJECT_VIEWER_DISPLAY_NAME = f"{APP_DISPLAY_NAME} Object Viewer"

_OBJECT_VIEWER_DEFAULTS: dict[str, object] = {
    "presentation_id": "instrument",
    "star_data_policy": "positional_static",
    "theme": OBJECT_VIEWER_THEME_PRESET,
    "sky_opacity": 0.0,
    "sky_disc_style": "smooth",
    "sky_disc_altaz_rings": "off",
    "sky_disc_altaz_rings_hover": "off",
    "vmag_limit": 4.0,
    "show_dso_initial": False,
    "show_asterisms_initial": False,
    "ground_tint_opacity": 0.0,
    "night_light_opacity": 0.0,
    "ridge_glow_opacity": 0.0,
    "observation_info": "off",
    "light_background_star_outline": True,
}


def apply_object_viewer_profile(args: object) -> None:
    """Apply object-viewer defaults without overriding explicit CLI options."""
    explicit_options = set(getattr(args, "_explicit_options", set()) or set())
    if "presentation_id" not in explicit_options:
        setattr(args, "presentation_id", "instrument")
    if "star_data_policy" not in explicit_options:
        setattr(args, "star_data_policy", "positional_static")
    for key, value in _OBJECT_VIEWER_DEFAULTS.items():
        if key in explicit_options:
            continue
        setattr(args, key, value)
    if "light_background_star_outline" not in explicit_options:
        setattr(args, "light_background_star_outline", True)


def main() -> None:
    """Run the white-background object viewer."""
    viewer.main(
        apply_app_profile=apply_object_viewer_profile,
        app_name=OBJECT_VIEWER_DISPLAY_NAME,
        parser_default_overrides=_OBJECT_VIEWER_DEFAULTS,
        parser_description="White-background sky object viewer",
        include_scenic_arguments=False,
        vmag_limit_max=6.0,
    )
