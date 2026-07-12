from zstarview.paths import ATLAS_THEME_PRESET, THEME_STYLES_BY_PRESET


def test_white_and_day_window_background_use_warm_off_white_base() -> None:
    assert THEME_STYLES_BY_PRESET["white"].window_background.base_rgb == (240, 238, 237)
    assert THEME_STYLES_BY_PRESET["day"].window_background.base_rgb == (226, 223, 222)


def test_atlas_theme_uses_white_background_and_ink_text() -> None:
    theme = THEME_STYLES_BY_PRESET[ATLAS_THEME_PRESET]

    assert theme.window_background.base_rgb == (255, 255, 255)
    assert theme.window_background.inner_rgba == (255, 255, 255, 255)
    assert theme.window_background.flat_background is True
    assert theme.sky_disc.opacity == 0.0
    assert theme.text.foreground_rgb == (24, 24, 24, 255)
    assert theme.text.outline_rgba == (255, 255, 255, 220)
    assert theme.status_text.outline_rgba == (255, 255, 255, 220)
    assert theme.label_outline_suppressed is False


def test_white_and_day_window_background_share_the_same_base_rgb() -> None:
    white_bg = THEME_STYLES_BY_PRESET["white"].window_background
    day_bg = THEME_STYLES_BY_PRESET["day"].window_background

    assert white_bg.base_rgb != day_bg.base_rgb
    assert white_bg.base_rgb[0] > day_bg.base_rgb[0]


def test_night_style_outline_colors_are_shared_by_bright_and_black_themes() -> None:
    for preset in ("night", "white", "day", "black"):
        theme = THEME_STYLES_BY_PRESET[preset]
        assert theme.text.outline_rgba == (0, 0, 0, 76)
        assert theme.status_text.outline_rgba == (0, 0, 0, 76)


def test_white_text_matches_day_with_higher_opacity() -> None:
    assert THEME_STYLES_BY_PRESET["white"].text.foreground_rgb == (214, 136, 103, 255)
    assert THEME_STYLES_BY_PRESET["day"].text.foreground_rgb == (214, 136, 103)


def test_white_status_text_matches_day_with_higher_opacity() -> None:
    assert THEME_STYLES_BY_PRESET["day"].status_text.foreground_rgb == (190, 190, 160)
    assert THEME_STYLES_BY_PRESET["white"].status_text.foreground_rgb == (190, 190, 160, 255)
    assert THEME_STYLES_BY_PRESET["white"].status_text.foreground_rgb[:3] == THEME_STYLES_BY_PRESET[
        "day"
    ].status_text.foreground_rgb


def test_transparent_theme_uses_low_alpha_dark_background() -> None:
    transparent = THEME_STYLES_BY_PRESET["transparent"].window_background
    black = THEME_STYLES_BY_PRESET["black"].window_background

    assert transparent.base_rgb[0] <= 8
    assert transparent.base_rgb[1] <= 8
    assert transparent.base_rgb[2] <= 9
    assert transparent.flat_background is True
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
    assert transparent.star_visibility_boost == 1.0
    assert white.star_visibility_boost == 1.12
    assert THEME_STYLES_BY_PRESET["day"].star_visibility_boost == 1.05
    assert transparent.sky_disc.opacity == 0.4
    assert white.window_chrome.menu_button_text_rgb == (70, 70, 70)
    assert night.window_chrome.menu_button_text_rgb == (210, 210, 210)
