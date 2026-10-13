from __future__ import annotations

import pytest

from zstarview.cli.args import parse_args


def test_parse_args_accepts_place_options() -> None:
    args = parse_args(["--place", "Matsue Station", "--place-countrycode", "jp"])
    assert args.place == "Matsue Station"
    assert args.place_countrycode == "jp"
    assert args.place_lang == "en"


def test_parse_args_accepts_timezone_option() -> None:
    args = parse_args(["--timezone", "Asia/Tokyo"])
    assert args.timezone == "Asia/Tokyo"


def test_parse_args_accepts_inverted_city_startup_option() -> None:
    args = parse_args(["--inverted-city"])

    assert args.inverted_city is True


@pytest.mark.parametrize("value,expected", [
    ("default", "normal"),
    ("inverted", "inverted-city"),
    ("simple", "simple-no-labels"),
    ("labels", "simple-labels"),
    ("landscape", "landscape"),
    (" SIMPLE ", "simple-no-labels"),
])
def test_initial_display_mode_aliases(value, expected) -> None:
    assert parse_args(["--initial-display-mode", value]).initial_display_mode == expected


def test_initial_display_mode_rejects_unknown_mode() -> None:
    with pytest.raises(SystemExit) as error:
        parse_args(["--initial-display-mode", "unknown"])
    assert error.value.code == 2


def test_parse_args_rejects_place_with_location_argument() -> None:
    with pytest.raises(SystemExit):
        parse_args(["Tokyo", "--place", "Matsue Station"])


def test_parse_args_rejects_place_countrycode_without_place() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--place-countrycode", "jp"])


def test_parse_args_rejects_place_lang_without_place() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--place-lang", "ja"])


def test_parse_args_rejects_empty_timezone() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--timezone", " "])
