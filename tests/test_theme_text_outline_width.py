from zstarview.render.text import get_text_outline_width


def test_theme_text_outline_width_is_three_for_all_presets() -> None:
    for preset in ("night", "black", "day", "white"):
        assert get_text_outline_width(preset) == 3.0
        assert get_text_outline_width(preset, status_line=True) == 3.0
