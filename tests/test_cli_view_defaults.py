from __future__ import annotations

from zstarview.cli import args as cli_args


def test_viewer_defaults_to_low_altitude_view() -> None:
    args = cli_args.parse_args([])

    assert args.view_center_alt == 5.0
    assert args.view_center_az == 180.0
    assert args.view_center_alt_specified is False


def test_export_image_defaults_to_low_altitude_view() -> None:
    args = cli_args.build_export_image_argument_parser().parse_args([])

    assert args.view_center_alt == 5.0
