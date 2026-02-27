from __future__ import annotations

from zstarview.render.draw import planet_disc_style_from_vmag


def test_planet_disc_style_gets_stronger_for_brighter_vmag() -> None:
    radius_faint, alpha_faint = planet_disc_style_from_vmag(5.0)
    radius_bright, alpha_bright = planet_disc_style_from_vmag(-2.0)

    assert radius_faint < radius_bright
    assert alpha_faint < alpha_bright


def test_planet_disc_style_handles_missing_vmag() -> None:
    radius, alpha = planet_disc_style_from_vmag(None)
    assert radius > 0.0
    assert alpha > 0


def test_planet_disc_style_explicitly_clips_very_bright_values() -> None:
    saturated = planet_disc_style_from_vmag(-1.5)
    brighter_input = planet_disc_style_from_vmag(-4.9)
    assert brighter_input == saturated
