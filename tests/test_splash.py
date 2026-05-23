from __future__ import annotations

from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.splash import _get_splash_palette


def test_splash_info_color_uses_theme_text_color() -> None:
    for preset_name in ("night", "day", "white", "black", "transparent"):
        theme = THEME_STYLES_BY_PRESET[preset_name]
        _gradient_colors, _frame_color, info_color = _get_splash_palette(theme)

        assert (
            info_color.red(),
            info_color.green(),
            info_color.blue(),
        ) == theme.splash.info_text_rgb


def test_day_and_white_splash_info_color_is_dark() -> None:
    assert THEME_STYLES_BY_PRESET["day"].splash.info_text_rgb == (32, 32, 32)
    assert THEME_STYLES_BY_PRESET["white"].splash.info_text_rgb == (32, 32, 32)
