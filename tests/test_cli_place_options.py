from __future__ import annotations

import pytest

from zstarview.zstarview import parse_args


def test_parse_args_accepts_place_options() -> None:
    args = parse_args(["--place", "Matsue Station", "--place-countrycode", "jp"])
    assert args.place == "Matsue Station"
    assert args.place_countrycode == "jp"
    assert args.place_lang == "en"


def test_parse_args_accepts_timezone_option() -> None:
    args = parse_args(["--timezone", "Asia/Tokyo"])
    assert args.timezone == "Asia/Tokyo"


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
