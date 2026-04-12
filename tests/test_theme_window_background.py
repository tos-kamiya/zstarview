from zstarview.paths import THEME_STYLES_BY_PRESET


def test_white_and_day_window_background_use_warm_off_white_base() -> None:
    expected_base_rgb = (226, 223, 222)
    for preset in ("white", "day"):
        assert THEME_STYLES_BY_PRESET[preset].window_background.base_rgb == expected_base_rgb


def test_white_and_day_window_background_share_the_same_base_rgb() -> None:
    white_bg = THEME_STYLES_BY_PRESET["white"].window_background
    day_bg = THEME_STYLES_BY_PRESET["day"].window_background

    assert white_bg.base_rgb == day_bg.base_rgb


def test_white_and_day_outline_colors_use_midpoint_tones() -> None:
    assert THEME_STYLES_BY_PRESET["white"].text.outline == (79, 45, 19, 166)
    assert THEME_STYLES_BY_PRESET["white"].status_text.outline == (76, 44, 19, 169)
    assert THEME_STYLES_BY_PRESET["day"].text.outline == (78, 44, 25, 163)
    assert THEME_STYLES_BY_PRESET["day"].status_text.outline == (73, 42, 24, 165)
