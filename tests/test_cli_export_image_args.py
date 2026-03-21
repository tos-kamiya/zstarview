from __future__ import annotations

import pytest

from zstarview.cli_args import parse_export_image_args


def test_parse_export_image_args_accepts_shared_and_export_specific_options() -> None:
    args = parse_export_image_args(
        [
            "--place",
            "Matsue Station",
            "--content-fov-deg",
            "110",
            "--image-size",
            "1280,720",
            "--layer-timeout-seconds",
            "12.5",
            "--allow-partial-data",
            "-o",
            "out.png",
        ]
    )

    assert args.place == "Matsue Station"
    assert args.content_fov_deg == 110.0
    assert args.image_size == (1280, 720)
    assert args.layer_timeout_seconds == 12.5
    assert args.allow_partial_data is True
    assert args.output == "out.png"


def test_parse_export_image_args_rejects_window_geometry() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--window-geometry", "restore", "-o", "out.png"])


def test_parse_export_image_args_rejects_sky_update_interval() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--sky-update-interval", "30", "-o", "out.png"])


def test_parse_export_image_args_rejects_dataset_query_options() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--list-viewpoints", "t", "-o", "out.png"])
