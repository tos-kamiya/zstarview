from __future__ import annotations

import argparse
import re

import pytest

from zstarview.cli import args as cli_args
from zstarview.paths import OVERLAY_FONT_SIZE_DEFAULT


def _build_export_like_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    cli_args.add_location_arguments(parser)
    cli_args.add_dataset_query_arguments(parser)
    cli_args.add_time_arguments(parser)
    cli_args.add_render_arguments(
        parser,
        include_window_geometry=False,
        include_window_frame=False,
        include_sky_update_interval=False,
        include_startup_overlay_arguments=False,
    )
    parser.add_argument("--edge-fov-deg", type=float, default=95.0)
    return parser


def _assert_tokens_in_order(text: str, tokens: list[str]) -> None:
    positions = [text.index(token) for token in tokens]
    assert positions == sorted(positions)


def test_render_argument_helpers_can_build_parser_without_gui_only_options() -> None:
    parser = _build_export_like_parser()

    args = parser.parse_args(
        [
            "--place",
            " Matsue Station ",
            "--content-fov-deg",
            "110",
            "--earth-guide-opacity",
            "0.12",
            "--theme",
            "day",
        ]
    )
    cli_args._normalize_location_arguments(parser, args)

    assert args.place == "Matsue Station"
    assert args.content_fov_deg == 110.0
    assert args.earth_guide_opacity == 0.12
    assert args.theme == "day"
    assert not hasattr(args, "window_geometry")
    assert not hasattr(args, "window_frame")
    assert not hasattr(args, "sky_update_interval")
    assert not hasattr(args, "show_dso_initial")
    assert not hasattr(args, "show_asterisms_initial")


def test_render_argument_helpers_reject_edge_fov_larger_than_content_fov() -> None:
    parser = _build_export_like_parser()

    args = parser.parse_args(
        [
            "--place",
            "Matsue",
            "--edge-fov-deg",
            "120",
            "--content-fov-deg",
            "110",
        ]
    )

    with pytest.raises(SystemExit):
        cli_args._validate_fov_relationship(parser, args)


def test_dataset_query_validation_works_with_parser_missing_gui_only_options() -> None:
    parser = _build_export_like_parser()
    args = parser.parse_args(["--list-viewpoints", "t", "--theme", "day"])

    cli_args._normalize_dataset_query_arguments(parser, args)

    with pytest.raises(SystemExit):
        cli_args._validate_dataset_query_compatibility(parser, args)


def test_main_help_text_is_ascii_only_for_windows_consoles() -> None:
    help_text = cli_args.build_main_argument_parser().format_help()

    assert all(ord(ch) < 128 for ch in help_text)


def test_main_help_text_uses_readme_like_groups() -> None:
    help_text = cli_args.build_main_argument_parser().format_help()
    general_match = re.search(r"\nGeneral:\n(?P<section>.*)$", help_text, re.S)
    overlays_match = re.search(r"\nOverlays:\n(?P<section>.*?)(?:\n\n[A-Z][^\n]*:\n|$)", help_text, re.S)

    assert "Observing Location and Time" in help_text
    assert "Search Objects at startup" in help_text
    assert "Sky and Stars" in help_text
    assert "Overlays" in help_text
    assert "General" in help_text
    assert general_match is not None
    assert "--observation-info" in general_match.group("section")
    assert "--visibility-boost" in general_match.group("section")
    assert "--overlay-font-size" in general_match.group("section")
    assert overlays_match is not None
    assert "--visibility-boost" not in overlays_match.group("section")
    assert "--overlay-font-size" not in overlays_match.group("section")
    _assert_tokens_in_order(
        overlays_match.group("section"),
        [
            "--cloud-opacity",
            "--cloud-stripe",
            "--cloud-missing-tint-opacity",
            "--tropical-cyclone-opacity",
            "--aircraft-opacity",
            "--satellite-opacity",
            "--show-guidelines-initial",
            "--terrain-horizon-opacity",
            "--earth-guide-opacity",
            "--ground-tint-opacity",
            "--water-surface-opacity",
            "--night-light-opacity",
            "--urban-outline-opacity",
        ],
    )
    assert re.search(r"^\s+--list\s", help_text, re.M) is None


def test_main_parser_accepts_height_add_option() -> None:
    args = cli_args.parse_args(["--height-add-m", "123", "Matsue"])

    assert float(args.observer_height_m) == 123.0


def test_main_parser_sets_default_overlay_font_size() -> None:
    args = cli_args.parse_args(["Matsue"])

    assert args.overlay_font_size == OVERLAY_FONT_SIZE_DEFAULT


def test_main_parser_accepts_overlay_font_size_override() -> None:
    args = cli_args.parse_args(["--overlay-font-size", "14.5", "Matsue"])

    assert args.overlay_font_size == 14.5


def test_main_parser_accepts_geo_satellite_option() -> None:
    args = cli_args.parse_args(["--geo-satellite", "true", "Matsue"])

    assert args.geo_satellite is True


def test_main_parser_rejects_overlay_font_size_out_of_range() -> None:
    with pytest.raises(SystemExit):
        cli_args.parse_args(["--overlay-font-size", "7", "Matsue"])


def test_export_image_help_text_uses_shared_groups() -> None:
    help_text = cli_args.build_export_image_argument_parser().format_help()
    general_match = re.search(r"\nGeneral:\n(?P<section>.*)$", help_text, re.S)
    overlays_match = re.search(r"\nOverlays:\n(?P<section>.*?)(?:\n\n[A-Z][^\n]*:\n|$)", help_text, re.S)

    assert "Observing Location and Time" in help_text
    assert "Search Objects at startup" in help_text
    assert "Overlays" in help_text
    assert "Export" in help_text
    assert "General" in help_text
    assert general_match is not None
    assert "--observation-info" in general_match.group("section")
    assert "--visibility-boost" in general_match.group("section")
    assert "--overlay-font-size" in general_match.group("section")
    assert overlays_match is not None
    assert "--visibility-boost" not in overlays_match.group("section")
    assert "--overlay-font-size" not in overlays_match.group("section")
    _assert_tokens_in_order(
        overlays_match.group("section"),
        [
            "--cloud-opacity",
            "--tropical-cyclone-opacity",
            "--aircraft-opacity",
            "--satellite-opacity",
            "--terrain-horizon-opacity",
            "--earth-guide-opacity",
            "--water-surface-opacity",
            "--night-light-opacity",
            "--urban-outline-opacity",
        ],
    )
    assert "--include-direction-grid" in help_text
    assert "--window-frame" not in help_text
    assert re.search(r"^\s+--list\s", help_text, re.M) is not None


def test_main_parser_version_option_prints_package_version(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli_args.parse_args(["--version"])

    captured = capsys.readouterr()

    assert exc_info.value.code == 0
    assert captured.out.endswith(f" {cli_args.__version__}\n")
    assert captured.err == ""


def test_parse_args_clamps_vmag_limit_to_committed_catalog_max() -> None:
    args = cli_args.parse_args(["-V", "11", "Matsue"])

    assert args.vmag_limit == 10.5


def test_parse_args_defaults_vmag_limit_to_seven() -> None:
    args = cli_args.parse_args([])

    assert args.vmag_limit == 7.0


def test_parse_args_accepts_clear_long_lived_cache() -> None:
    args = cli_args.parse_args(["--clear-long-lived-cache", "Matsue"])

    assert args.clear_long_lived_cache is True


def test_parse_args_defaults_observation_info_to_bottom() -> None:
    args = cli_args.parse_args(["Matsue"])

    assert args.observation_info == "bottom"


def test_parse_args_accepts_window_frame_mode() -> None:
    args = cli_args.parse_args(["--window-frame", "window", "Matsue"])

    assert args.window_frame == "window"


def test_parse_args_marks_explicit_view_center_values() -> None:
    args = cli_args.parse_args(["-A90", "--view-center-az=180", "Matsue"])

    assert args.view_center_alt_specified is True
    assert args.view_center_az_specified is True


def test_parse_args_accepts_earth_guide_opacity() -> None:
    args = cli_args.parse_args(["--earth-guide-opacity", "0", "Matsue"])

    assert args.earth_guide_opacity == 0.0


def test_parse_args_accepts_terrain_horizon_opacity_short_option() -> None:
    args = cli_args.parse_args(["-d", "0.5", "Matsue"])

    assert args.terrain_horizon_opacity == 0.5


def test_parse_args_defaults_terrain_horizon_opacity() -> None:
    args = cli_args.parse_args(["Matsue"])

    assert args.terrain_horizon_opacity == 0.003


def test_parse_args_accepts_earth_guide_opacity_short_option() -> None:
    args = cli_args.parse_args(["-e", "0.4", "Matsue"])

    assert args.earth_guide_opacity == 0.4


def test_parse_args_defaults_night_light_opacity() -> None:
    args = cli_args.parse_args(["Matsue"])

    assert args.night_light_opacity == 0.04


def test_parse_args_defaults_ridge_glow_opacity() -> None:
    args = cli_args.parse_args(["Matsue"])

    assert args.ridge_glow_opacity == 0.3


def test_parse_args_accepts_urban_outline_opacity_short_option() -> None:
    args = cli_args.parse_args(["-u", "0.3", "Matsue"])

    assert args.urban_outline_opacity == 0.3


def test_parse_args_defaults_sky_opacity_to_0_1() -> None:
    args = cli_args.parse_args(["Matsue"])

    assert args.sky_opacity == 0.1


def test_parse_args_defaults_sky_disc_style_to_smooth() -> None:
    args = cli_args.parse_args(["Matsue"])

    assert args.sky_disc_style == "smooth"


def test_parse_args_defaults_sky_disc_altaz_rings() -> None:
    args = cli_args.parse_args(["Matsue"])

    assert args.sky_disc_altaz_rings == "dimalt"
    assert args.sky_disc_altaz_rings_hover == "altaz"


def test_parse_args_accepts_sky_disc_altaz_ring_modes() -> None:
    args = cli_args.parse_args(
        [
            "--sky-disc-altaz-rings",
            "dimalt",
            "--sky-disc-altaz-rings-hover",
            "off",
            "Matsue",
        ]
    )

    assert args.sky_disc_altaz_rings == "dimalt"
    assert args.sky_disc_altaz_rings_hover == "off"


def test_parse_args_rejects_sky_disc_grid_style() -> None:
    with pytest.raises(SystemExit):
        cli_args.parse_args(["--sky-disc-style", "grid", "Matsue"])


def test_parse_args_rejects_multiple_search_options() -> None:
    with pytest.raises(SystemExit):
        cli_args.parse_args(["--search", "Ceres", "--search", "Mars"])
