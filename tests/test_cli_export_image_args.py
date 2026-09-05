from __future__ import annotations

import pytest
from PySide6.QtGui import QFont

import zstarview.cli.export_image as mod
from zstarview.__about__ import __version__
from zstarview.cli.args import parse_export_image_args
from zstarview.cli.export_image import (
    _clamp_view_center_altitude,
    _format_search_failure_message,
    _search_view_center_for_target,
)
from zstarview.paths import OBSERVER_MAX_ALT_DEG, OBSERVER_MIN_ALT_DEG


def test_parse_export_image_args_accepts_shared_and_export_specific_options() -> None:
    args = parse_export_image_args(
        [
            "--place",
            "Matsue Station",
            "--content-fov-deg",
            "110",
            "--overlay-font-size",
            "14.5",
            "--water-surface-opacity",
            "0.12",
            "--tropical-cyclone-opacity",
            "0.25",
            "--urban-outline-max-candidates",
            "2500",
            "--road-light-max-candidates",
            "1200",
            "--urban-outline-skyscraper-radius-km",
            "48",
            "--image-size",
            "1280,720",
            "--layer-timeout-seconds",
            "12.5",
            "--allow-partial-data",
            "--geo-satellite",
            "true",
            "--include-direction-grid",
            "-o",
            "out.png",
        ]
    )

    assert args.place == "Matsue Station"
    assert args.content_fov_deg == 110.0
    assert args.overlay_font_size == 14.5
    assert args.water_surface_opacity == 0.12
    assert args.tropical_cyclone_opacity == 0.25
    assert args.urban_outline_max_candidates == 2500
    assert args.road_light_max_candidates == 1200
    assert args.urban_outline_skyscraper_radius_km == 48.0
    assert args.image_size == (1280, 720)
    assert args.layer_timeout_seconds == 12.5
    assert args.allow_partial_data is True
    assert args.geo_satellite is True
    assert args.include_direction_grid is True
    assert args.output == "out.png"
    assert args.sixel is False


def test_parse_export_image_args_defaults_layer_timeout_to_180_seconds() -> None:
    args = parse_export_image_args(["-o", "out.png"])

    assert args.layer_timeout_seconds == 180.0


def test_load_fonts_uses_fractional_overlay_font_size(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(mod.QFontDatabase, "addApplicationFont", lambda _path: 1)
    monkeypatch.setattr(
        mod.QFontDatabase,
        "applicationFontFamilies",
        lambda _font_id: ["Test Font"],
    )

    text_font, status_line_font = mod._load_fonts(14.5)

    assert isinstance(text_font, QFont)
    assert text_font.pointSizeF() == 14.5
    assert status_line_font.pointSizeF() == 8.0


def test_parse_export_image_args_supports_version_option(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--version"])

    captured = capsys.readouterr()
    assert __version__ in captured.out


def test_parse_export_image_args_accepts_sixel_without_output() -> None:
    args = parse_export_image_args(["--sixel", "Matsue"])

    assert args.city == "Matsue"
    assert args.sixel is True
    assert args.output is None


def test_parse_export_image_args_marks_explicit_view_center_values() -> None:
    args = parse_export_image_args(["-A90", "--view-center-az=180", "-o", "out.png"])

    assert args.view_center_alt_specified is True
    assert args.view_center_az_specified is True


def test_parse_export_image_args_requires_output_or_sixel() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["Matsue"])


def test_parse_export_image_args_rejects_stdout_output_with_sixel() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["-o", "-", "--sixel", "Matsue"])


def test_parse_export_image_args_rejects_window_geometry() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--window-geometry", "restore", "-o", "out.png"])


def test_parse_export_image_args_rejects_window_frame() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--window-frame", "window", "-o", "out.png"])


@pytest.mark.parametrize(
    "option",
    ["--display-tone-curve=12,247", "--calibrate-display-tone-curve"],
)
def test_parse_export_image_args_rejects_display_tone_options(option: str) -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args([option, "-o", "out.png"])


def test_parse_export_image_args_rejects_sky_update_interval() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--sky-update-interval", "30", "-o", "out.png"])


def test_parse_export_image_args_rejects_dataset_query_options() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--list-viewpoints", "t", "-o", "out.png"])


def test_parse_export_image_args_clamps_vmag_limit_to_committed_catalog_max() -> None:
    args = parse_export_image_args(["-V", "11", "-o", "out.png"])

    assert args.vmag_limit == 10.5


def test_parse_export_image_args_accepts_clear_long_lived_cache() -> None:
    args = parse_export_image_args(["--clear-long-lived-cache", "-o", "out.png"])

    assert args.clear_long_lived_cache is True


def test_parse_export_image_args_warns_for_deprecated_urban_outline_min_height(
    capsys: pytest.CaptureFixture[str],
) -> None:
    args = parse_export_image_args(
        ["--urban-outline-min-building-height-m", "10", "-o", "out.png"]
    )

    captured = capsys.readouterr()
    assert "deprecated and ignored" in captured.err
    assert args.urban_outline_min_height_m == 0.0


def test_parse_export_image_args_ignores_deprecated_urban_outline_feature_type(
    capsys: pytest.CaptureFixture[str],
) -> None:
    args = parse_export_image_args(
        ["--urban-outline-feature-type", "building", "-o", "out.png"]
    )

    captured = capsys.readouterr()
    assert "deprecated and ignored" in captured.err
    assert args.urban_outline_feature_type == "both"


def test_parse_export_image_args_accepts_print_cache_dir_without_output() -> None:
    args = parse_export_image_args(["--print-cache-dir"])

    assert args.print_cache_dir is True


def test_parse_export_image_args_rejects_list_without_search() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--list"])


def test_parse_export_image_args_rejects_multiple_search_options() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--search", "Ceres", "--search", "Mars", "-o", "out.png"])


def test_format_search_failure_message_reports_no_results() -> None:
    assert _format_search_failure_message("Jupyter", 0) == "No search candidates found for 'Jupyter'"


def test_format_search_failure_message_reports_multiple_results() -> None:
    assert _format_search_failure_message("Mars", 2) == "Search query 'Mars' matched 2 candidates"


def test_clamp_view_center_altitude_matches_gui_floor() -> None:
    assert _clamp_view_center_altitude(-50.0) == OBSERVER_MIN_ALT_DEG
    assert _clamp_view_center_altitude(-20.0) == -20.0
    assert _clamp_view_center_altitude(120.0) == OBSERVER_MAX_ALT_DEG


def test_search_view_center_for_target_honors_fixed_axes() -> None:
    assert _search_view_center_for_target(
        base_view_center=(5.0, 180.0),
        target_alt_deg=12.5,
        target_az_deg=220.0,
        fixed_alt=True,
        fixed_az=False,
    ) == (5.0, 220.0)

    assert _search_view_center_for_target(
        base_view_center=(45.0, 210.0),
        target_alt_deg=12.5,
        target_az_deg=220.0,
        fixed_alt=False,
        fixed_az=True,
    ) == (12.5, 210.0)

    assert _search_view_center_for_target(
        base_view_center=(5.0, 210.0),
        target_alt_deg=12.5,
        target_az_deg=220.0,
        fixed_alt=True,
        fixed_az=True,
    ) == (5.0, 210.0)


def test_parse_export_image_args_rejects_skyscraper_radius_smaller_than_base_radius() -> (
    None
):
    with pytest.raises(SystemExit):
        parse_export_image_args(
            [
                "--urban-outline-radius-km",
                "12",
                "--urban-outline-skyscraper-radius-km",
                "8",
                "-o",
                "out.png",
            ]
        )


def test_parse_export_image_args_accepts_zero_skyscraper_radius() -> None:
    args = parse_export_image_args(
        ["--urban-outline-skyscraper-radius-km", "0", "-o", "out.png"]
    )

    assert args.urban_outline_skyscraper_radius_km == 0.0


def test_parse_export_image_args_rejects_print_cache_dir_with_output() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--print-cache-dir", "-o", "out.png"])
