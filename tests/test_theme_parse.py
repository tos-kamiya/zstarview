from __future__ import annotations

import argparse

from zstarview.zstarview import _parse_theme


def test_parse_theme_accepts_black_case_insensitive() -> None:
    assert _parse_theme("black") == "black"
    assert _parse_theme(" BLACK ") == "black"


def test_parse_theme_rejects_invalid_theme_with_updated_message() -> None:
    try:
        _parse_theme("dark")
    except argparse.ArgumentTypeError as e:
        message = str(e)
    else:
        raise AssertionError("Expected argparse.ArgumentTypeError for invalid theme")

    assert "night, day, white, black" in message
