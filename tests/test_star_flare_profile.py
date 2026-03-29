from __future__ import annotations

from zstarview.render.photometry import compute_flare_profile, flare_strength_from_vmag


def test_flare_strength_monotonic_for_all_magnitudes() -> None:
    mags_faint_to_bright = [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0]
    strengths = [flare_strength_from_vmag(v) for v in mags_faint_to_bright]
    assert all(a <= b for a, b in zip(strengths, strengths[1:]))


def test_flare_is_disabled_when_outer_size_is_less_than_one_pixel() -> None:
    core_scale, flare_outer_px = compute_flare_profile(vmag=6.0, core_radius_px=1.5)
    assert core_scale == 1.0
    assert flare_outer_px == 0.0


def test_perceived_radius_is_non_decreasing_for_brightening_stars() -> None:
    mags_faint_to_bright = [3.0, 2.0, 1.0, 0.0, -1.0]
    base_radius_px = 4.0
    perceived = []
    for vmag in mags_faint_to_bright:
        core_scale, flare_outer_px = compute_flare_profile(vmag=vmag, core_radius_px=base_radius_px)
        core = base_radius_px * core_scale
        perceived.append(core + flare_outer_px)
    assert all(a <= b for a, b in zip(perceived, perceived[1:]))
