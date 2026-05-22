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
        ) == theme.text.foreground_rgb
