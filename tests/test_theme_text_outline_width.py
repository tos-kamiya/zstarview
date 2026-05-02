from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.render.text import get_text_outline_width


def test_theme_text_outline_width_is_three_for_all_presets() -> None:
    for preset in ("night", "black", "transparent", "day", "white"):
        theme = THEME_STYLES_BY_PRESET[preset]
        assert get_text_outline_width(theme) == 3.0
        assert get_text_outline_width(theme, status_line=True) == 3.0
