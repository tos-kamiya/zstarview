from zstarview.cli.args import parse_args
from zstarview.gui.object_viewer import apply_object_viewer_profile
from zstarview.paths import OBJECT_VIEWER_THEME_PRESET


def test_object_viewer_profile_applies_white_background_defaults() -> None:
    args = parse_args([])

    apply_object_viewer_profile(args)

    assert args.presentation_id == "instrument"
    assert args.star_data_policy == "positional_static"
    assert args.theme == OBJECT_VIEWER_THEME_PRESET
    assert args.sky_opacity == 0.0
    assert args.sky_disc_altaz_rings == "off"
    assert args.vmag_limit == 6.0
    assert args.ground_tint_opacity == 0.0
    assert args.night_light_opacity == 0.0
    assert args.ridge_glow_opacity == 0.0
    assert args.observation_info == "off"
    assert args.light_background_star_outline is True


def test_object_viewer_profile_preserves_explicit_cli_options() -> None:
    args = parse_args(
        [
            "--theme",
            "night",
            "--sky-opacity",
            "0.25",
            "--vmag-limit",
            "4.5",
            "--observation-info",
            "bottom",
        ]
    )

    apply_object_viewer_profile(args)

    assert args.theme == "night"
    assert args.sky_opacity == 0.25
    assert args.vmag_limit == 4.5
    assert args.observation_info == "bottom"
    assert args.light_background_star_outline is True
