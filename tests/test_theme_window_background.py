import pytest

from zstarview.paths import THEME_PRESET_NAMES, THEME_STYLES_BY_PRESET


def test_theme_presets_are_night_and_transparent_only() -> None:
    assert THEME_PRESET_NAMES[0] == "night"
    assert all(
        name == "night" or name.startswith("transparent") for name in THEME_PRESET_NAMES
    )
    for removed in ("day", "white", "black"):
        assert removed not in THEME_PRESET_NAMES
        assert removed not in THEME_STYLES_BY_PRESET


@pytest.mark.parametrize("preset", ["night", "transparent"])
def test_theme_outline_colors_are_dark(preset: str) -> None:
    theme = THEME_STYLES_BY_PRESET[preset]
    assert theme.text.outline_rgba[:3] in {(0, 0, 0), (4, 4, 6)}
    assert theme.status_text.outline_rgba[:3] in {(0, 0, 0), (4, 4, 6)}


def test_transparent_theme_uses_low_alpha_dark_background() -> None:
    transparent = THEME_STYLES_BY_PRESET["transparent"].window_background
    night = THEME_STYLES_BY_PRESET["night"].window_background

    assert transparent.base_rgb[0] <= 8
    assert transparent.base_rgb[1] <= 8
    assert transparent.base_rgb[2] <= 9
    assert transparent.flat_background is True
    assert transparent.delta_rgb == (0, 0, 0)
    assert transparent.inner_rgba[3] == transparent.outer_alpha == transparent.edge_alpha
    assert transparent.inner_rgba[3] < night.inner_rgba[3]
    assert transparent.outer_alpha < night.outer_alpha
    assert transparent.edge_alpha < night.edge_alpha


def test_theme_style_groups_window_chrome_and_sky_disc_values() -> None:
    transparent = THEME_STYLES_BY_PRESET["transparent"]
    night = THEME_STYLES_BY_PRESET["night"]

    assert transparent.window_chrome.menu_fill_rgba == (14, 14, 15, 96)
    assert transparent.star_visibility_boost == 1.0
    assert night.star_visibility_boost == 1.0
    assert transparent.sky_disc.opacity == 0.4
    assert night.window_chrome.menu_button_text_rgb == (210, 210, 210)
