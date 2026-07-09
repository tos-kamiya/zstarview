from zstarview.aircraft_constants import AIRCRAFT_OVERLAY_LINE_COLOR_RGB
from zstarview.paths import (
    CELESTIAL_EQUATOR_COLOR,
    HORIZON_LINE_COLOR,
    OBJECT_VIEWER_THEME_PRESET,
    PALETTE_AIRCRAFT_AND_SATELLITE_RGB,
    PALETTE_ASTERISM_RGB,
    PALETTE_ASTERISM_LABEL_RGB,
    PALETTE_EARTH_GUIDE_RGB,
    PALETTE_NEVER_RISES_GUIDE_RGB,
    PALETTE_NEVER_RISES_RGB,
    PALETTE_HORIZON_AND_LABEL_RGB,
    PALETTE_TERRAIN_HORIZON_RGB,
    TERRAIN_HORIZON_LINE_COLOR,
    THEME_STYLES_BY_PRESET,
)
from zstarview.render.deep_sky_objects import DSO_LABEL_TEXT_RGB
from zstarview.satellite_constants import SATELLITE_OVERLAY_MARKER_COLOR_RGB


def test_overlay_colors_follow_the_palette_swatches() -> None:
    assert PALETTE_EARTH_GUIDE_RGB == (112, 99, 89)
    assert PALETTE_NEVER_RISES_RGB == (192, 192, 192)
    assert PALETTE_NEVER_RISES_GUIDE_RGB == PALETTE_NEVER_RISES_RGB
    assert PALETTE_TERRAIN_HORIZON_RGB == (216, 206, 192)
    assert TERRAIN_HORIZON_LINE_COLOR == PALETTE_TERRAIN_HORIZON_RGB
    assert CELESTIAL_EQUATOR_COLOR == PALETTE_NEVER_RISES_RGB
    assert PALETTE_HORIZON_AND_LABEL_RGB == (206, 240, 122)
    assert HORIZON_LINE_COLOR == PALETTE_HORIZON_AND_LABEL_RGB
    assert PALETTE_AIRCRAFT_AND_SATELLITE_RGB == (227, 108, 240)
    assert AIRCRAFT_OVERLAY_LINE_COLOR_RGB == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert SATELLITE_OVERLAY_MARKER_COLOR_RGB == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert PALETTE_ASTERISM_RGB == (122, 226, 240)
    assert PALETTE_ASTERISM_LABEL_RGB == (111, 207, 219)
    assert DSO_LABEL_TEXT_RGB == (111, 158, 219)


def test_default_overlay_styles_preserve_existing_palette_swatches() -> None:
    overlays = THEME_STYLES_BY_PRESET["night"].overlays

    assert overlays.aircraft.rgb == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert overlays.aircraft.label_rgb == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert overlays.satellite.rgb == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert overlays.satellite.label_rgb == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert overlays.terrain_horizon.rgb == PALETTE_TERRAIN_HORIZON_RGB
    assert overlays.urban_outline.rgb == (214, 232, 255)
    assert overlays.water.rgb == (122, 218, 240)
    assert overlays.earth_guide.rgb == PALETTE_EARTH_GUIDE_RGB


def test_object_viewer_overlay_styles_use_white_background_palette() -> None:
    overlays = THEME_STYLES_BY_PRESET[OBJECT_VIEWER_THEME_PRESET].overlays

    assert overlays.aircraft.rgb == (178, 31, 107)
    assert overlays.aircraft.outline_rgba == (42, 24, 34, 80)
    assert overlays.aircraft.label_rgb == (143, 24, 88)
    assert overlays.satellite.rgb == overlays.aircraft.rgb
    assert overlays.satellite.outline_rgba == overlays.aircraft.outline_rgba
    assert overlays.satellite.label_rgb == overlays.aircraft.label_rgb
    assert overlays.terrain_horizon.rgb == (138, 119, 101)
    assert overlays.urban_outline.rgb == (110, 110, 110)
    assert overlays.water.rgb == (35, 139, 159)
    assert overlays.earth_guide.rgb == (122, 107, 92)

    guide_style = THEME_STYLES_BY_PRESET[OBJECT_VIEWER_THEME_PRESET].guide_style
    assert guide_style.simple_reference_lines is True
    assert guide_style.reference_rgb == (0, 0, 0)
    assert guide_style.horizon_rgb == (0, 0, 0)
    assert guide_style.equator_rgb == (0, 0, 0)
    assert guide_style.ecliptic_rgb == (0, 0, 0)
    assert guide_style.label_rgb == (0, 0, 0)
    assert guide_style.ecliptic_dash_pattern == (2.0, 3.0)


def test_default_guide_style_remains_non_object_viewer_style() -> None:
    guide_style = THEME_STYLES_BY_PRESET["night"].guide_style

    assert guide_style.simple_reference_lines is False
    assert guide_style.reference_rgb == (206, 240, 122)
    assert guide_style.equator_rgb == (192, 192, 192)
    assert guide_style.ecliptic_rgb == (236, 173, 2)
