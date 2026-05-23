from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.render.pipeline import _ground_reset_rgba_for_theme


def test_transparent_theme_uses_lower_sky_disc_opacity() -> None:
    assert THEME_STYLES_BY_PRESET["transparent"].sky_disc.opacity == 0.4


def test_non_transparent_themes_keep_full_sky_disc_opacity() -> None:
    for preset in ("night", "black", "day", "white"):
        assert THEME_STYLES_BY_PRESET[preset].sky_disc.opacity == 1.0


def test_transparent_theme_uses_lower_ground_reset_opacity() -> None:
    assert _ground_reset_rgba_for_theme(THEME_STYLES_BY_PRESET["transparent"]) == (4, 4, 4, 102)
