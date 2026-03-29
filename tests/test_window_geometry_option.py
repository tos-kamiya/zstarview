import argparse

import pytest

from zstarview.cli.args import _parse_window_geometry


def test_parse_window_geometry_accepts_restore_case_insensitive() -> None:
    assert _parse_window_geometry("restore") == "restore"
    assert _parse_window_geometry(" ReStOrE ") == "restore"


def test_parse_window_geometry_accepts_csv_integers() -> None:
    assert _parse_window_geometry("100,200,800,600") == (100, 200, 800, 600)


def test_parse_window_geometry_rejects_invalid_format() -> None:
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_window_geometry("100,200,800")
