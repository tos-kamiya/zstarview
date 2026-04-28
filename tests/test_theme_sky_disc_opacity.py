from zstarview.paths import SKY_DISC_OPACITY_BY_PRESET


def test_transparent_theme_uses_lower_sky_disc_opacity() -> None:
    assert SKY_DISC_OPACITY_BY_PRESET["transparent"] == 0.4


def test_non_transparent_themes_keep_full_sky_disc_opacity() -> None:
    for preset in ("night", "black", "day", "white"):
        assert SKY_DISC_OPACITY_BY_PRESET[preset] == 1.0
