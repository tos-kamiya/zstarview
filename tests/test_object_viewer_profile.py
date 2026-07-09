import pytest

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
    assert args.vmag_limit == 4.0
    assert args.show_guidelines_initial is True
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


def test_object_viewer_parser_uses_object_viewer_vmag_defaults_and_maximum() -> None:
    args = parse_args(
        [],
        default_overrides={"vmag_limit": 4.0},
        description="White-background sky object viewer",
        include_scenic_arguments=False,
        vmag_limit_max=6.0,
    )

    assert args.vmag_limit == 4.0

    args = parse_args(
        ["-V", "7"],
        default_overrides={"vmag_limit": 4.0},
        description="White-background sky object viewer",
        include_scenic_arguments=False,
        vmag_limit_max=6.0,
    )

    assert args.vmag_limit == 6.0


@pytest.mark.parametrize("option", ["--sky-opacity", "--night-light-opacity"])
def test_object_viewer_parser_rejects_scenic_only_options(option: str) -> None:
    with pytest.raises(SystemExit):
        parse_args(
            [option, "0.5"],
            default_overrides={"vmag_limit": 4.0},
            include_scenic_arguments=False,
            vmag_limit_max=6.0,
        )
