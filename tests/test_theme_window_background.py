from zstarview.paths import THEME_STYLES_BY_PRESET


def test_white_and_day_window_background_use_warm_off_white_base() -> None:
    assert THEME_STYLES_BY_PRESET["white"].window_background.base_rgb == (240, 238, 237)
    assert THEME_STYLES_BY_PRESET["day"].window_background.base_rgb == (226, 223, 222)


def test_white_and_day_window_background_share_the_same_base_rgb() -> None:
    white_bg = THEME_STYLES_BY_PRESET["white"].window_background
    day_bg = THEME_STYLES_BY_PRESET["day"].window_background

    assert white_bg.base_rgb != day_bg.base_rgb
    assert white_bg.base_rgb[0] > day_bg.base_rgb[0]


def test_white_and_day_outline_colors_use_midpoint_tones() -> None:
    assert THEME_STYLES_BY_PRESET["white"].text.outline_rgba == (79, 45, 19, 166)
    assert THEME_STYLES_BY_PRESET["white"].status_text.outline_rgba == (76, 44, 19, 169)
    assert THEME_STYLES_BY_PRESET["day"].text.outline_rgba == (78, 44, 25, 163)
    assert THEME_STYLES_BY_PRESET["day"].status_text.outline_rgba == (73, 42, 24, 165)


def test_white_and_day_text_colors_are_slightly_darker() -> None:
    assert THEME_STYLES_BY_PRESET["white"].text.foreground_rgb == (209, 149, 91)
    assert THEME_STYLES_BY_PRESET["white"].status_text.foreground_rgb == (201, 142, 86)
    assert THEME_STYLES_BY_PRESET["day"].text.foreground_rgb == (214, 136, 103)
    assert THEME_STYLES_BY_PRESET["day"].status_text.foreground_rgb == (207, 131, 98)


def test_transparent_theme_uses_low_alpha_dark_background() -> None:
    transparent = THEME_STYLES_BY_PRESET["transparent"].window_background
    black = THEME_STYLES_BY_PRESET["black"].window_background

    assert transparent.base_rgb[0] <= 8
    assert transparent.base_rgb[1] <= 8
    assert transparent.base_rgb[2] <= 9
    assert transparent.flat_background is True
    assert transparent.draw_outer_border is False
    assert transparent.delta_rgb == (0, 0, 0)
    assert transparent.inner_rgba[3] == transparent.outer_alpha == transparent.edge_alpha
    assert transparent.inner_rgba[3] < black.inner_rgba[3]
    assert transparent.outer_alpha < black.outer_alpha
    assert transparent.edge_alpha < black.edge_alpha


def test_theme_style_groups_window_chrome_and_sky_disc_values() -> None:
    transparent = THEME_STYLES_BY_PRESET["transparent"]
    white = THEME_STYLES_BY_PRESET["white"]
    night = THEME_STYLES_BY_PRESET["night"]

    assert transparent.window_chrome.menu_fill_rgba == (14, 14, 15, 96)
    assert white.window_background.draw_outer_border is True
    assert THEME_STYLES_BY_PRESET["day"].window_background.draw_outer_border is True
    assert night.window_background.draw_outer_border is False
    assert transparent.star_visibility_boost == 1.0
    assert white.star_visibility_boost == 1.12
    assert THEME_STYLES_BY_PRESET["day"].star_visibility_boost == 1.05
    assert transparent.sky_disc.opacity == 0.4
    assert white.window_chrome.menu_button_text_rgb == (70, 70, 70)
    assert night.window_chrome.menu_button_text_rgb == (210, 210, 210)
