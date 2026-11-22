from __future__ import annotations

import argparse

import pytest

from zstarview.cli.args import _parse_theme


def test_parse_theme_accepts_night_case_insensitive() -> None:
    assert _parse_theme("night") == "night"
    assert _parse_theme(" NIGHT ") == "night"


@pytest.mark.parametrize("removed", ["day", "white", "black", "translucent"])
def test_parse_theme_rejects_removed_presets(removed: str) -> None:
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_theme(removed)


def test_parse_theme_accepts_transparent_and_alias() -> None:
    assert _parse_theme("transparent") == "transparent-40"
    assert _parse_theme("transparent-40") == "transparent-40"


def test_parse_theme_accepts_transparent_opacity_variants() -> None:
    assert _parse_theme("transparent-10") == "transparent-10"
    assert _parse_theme("transparent-90") == "transparent-90"


def test_parse_theme_rejects_invalid_theme_with_updated_message() -> None:
    try:
        _parse_theme("dark")
    except argparse.ArgumentTypeError as e:
        message = str(e)
    else:
        raise AssertionError("Expected argparse.ArgumentTypeError for invalid theme")

    assert "night, transparent" in message
    assert "white" not in message


def test_parse_theme_rejects_non_preset_transparent_opacity() -> None:
    try:
        _parse_theme("transparent-05")
    except argparse.ArgumentTypeError as e:
        message = str(e)
    else:
        raise AssertionError("Expected argparse.ArgumentTypeError for invalid theme")

    assert "transparent-10" in message
