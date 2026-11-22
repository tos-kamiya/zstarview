from __future__ import annotations

import pytest

from zstarview.cli import args as cli_args

REMOVED_OPTIONS = [
    ["--urban-outline-feature-type", "building"],
    ["--urban-outline-min-building-height-m", "10"],
    ["-b", "10"],
    ["--sky-disc-style", "smooth"],
    ["--inverted-city"],
    ["--observer-height-m", "10"],
    ["--akari-ir-bands-opacity", "0.2"],
    ["--enlarge-moon"],
    ["-m"],
]


@pytest.mark.parametrize("argv", REMOVED_OPTIONS, ids=lambda argv: argv[0])
def test_viewer_rejects_removed_options(argv: list[str]) -> None:
    with pytest.raises(SystemExit):
        cli_args.parse_args(argv)


@pytest.mark.parametrize(
    "argv",
    [argv for argv in REMOVED_OPTIONS if argv[0] != "--inverted-city"],
    ids=lambda argv: argv[0],
)
def test_export_image_rejects_removed_options(argv: list[str]) -> None:
    with pytest.raises(SystemExit):
        cli_args.parse_export_image_args([*argv, "-o", "out.png"])


def test_replacement_options_still_work() -> None:
    args = cli_args.parse_args(
        [
            "--height-add-m",
            "10",
            "--diffuse-sky-source",
            "akari",
            "--diffuse-sky-opacity",
            "0.3",
            "--moon-style",
            "sphere",
            "--moon-scale",
            "5",
            "--initial-display-mode",
            "inverted",
        ]
    )

    assert args.observer_height_m == 10.0
    assert args.diffuse_sky_source == "akari"
    assert args.akari_ir_bands_opacity == 0.3
    assert (args.moon_style, args.moon_scale) == ("sphere", 5)
    assert args.initial_display_mode == "inverted-city"


def test_launch_profile_drops_removed_option_keys() -> None:
    from zstarview.gui.launch_profile import _without_ignored_keys

    profile = _without_ignored_keys(
        {
            "sky_disc_style": "smooth",
            "urban_outline_feature_type": "both",
            "urban_outline_min_height_m": 0.0,
            "inverted_city": True,
            "enlarge_moon": True,
            "theme": "night",
        }
    )

    assert profile == {"theme": "night"}
