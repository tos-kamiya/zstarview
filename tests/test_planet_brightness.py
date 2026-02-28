from __future__ import annotations

from zstarview.render.draw import planet_bloom_profile_from_vmag, planet_disc_style_from_vmag


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


def test_planet_bloom_profile_gets_stronger_for_very_bright_planets() -> None:
    at_clip = planet_bloom_profile_from_vmag(-1.5, 3.0)
    very_bright = planet_bloom_profile_from_vmag(-4.5, 3.0)
    assert very_bright[0] > at_clip[0]
    assert very_bright[1] >= at_clip[1]
    assert very_bright[2] >= at_clip[2]


def test_planet_bloom_profile_is_disabled_for_faint_planets() -> None:
    radius, a0, a1 = planet_bloom_profile_from_vmag(6.0, 3.0)
    assert radius == 0.0
    assert a0 == 0
    assert a1 == 0
