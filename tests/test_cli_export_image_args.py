from __future__ import annotations

import pytest

from zstarview.cli.args import parse_export_image_args


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
    assert args.sixel is False


def test_parse_export_image_args_accepts_sixel_without_output() -> None:
    args = parse_export_image_args(["--sixel", "Matsue"])

    assert args.city == "Matsue"
    assert args.sixel is True
    assert args.output is None


def test_parse_export_image_args_requires_output_or_sixel() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["Matsue"])


def test_parse_export_image_args_rejects_stdout_output_with_sixel() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["-o", "-", "--sixel", "Matsue"])


def test_parse_export_image_args_rejects_window_geometry() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--window-geometry", "restore", "-o", "out.png"])


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


def test_parse_export_image_args_accepts_print_cache_dir_without_output() -> None:
    args = parse_export_image_args(["--print-cache-dir"])

    assert args.print_cache_dir is True


def test_parse_export_image_args_rejects_print_cache_dir_with_output() -> None:
    with pytest.raises(SystemExit):
        parse_export_image_args(["--print-cache-dir", "-o", "out.png"])
