from __future__ import annotations

import json

import pytest

from zstarview.mountain_viewpoints import list_mountain_all_names, list_mountain_primary_names
from zstarview.zstarview import main, parse_args


def test_parse_args_accepts_list_mountains() -> None:
    args = parse_args(["--list-mountains"])
    assert args.list_mountains is True
    assert args.list_mountain_names is False
    assert args.show_mountain_json is None


def test_parse_args_rejects_location_with_mountain_listing() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-mountains", "Mount Fuji"])


def test_parse_args_rejects_render_option_with_mountain_listing() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-mountains", "--hours", "1"])


def test_list_mountain_primary_names_includes_mount_fuji() -> None:
    names = list_mountain_primary_names()
    assert "Mount Fuji" in names
    assert "Denali" in names


def test_list_mountain_primary_names_prefers_ascii_fallback() -> None:
    names = list_mountain_primary_names()
    assert "Ayrybaba" in names
    assert "Aýrybaba" not in names


def test_list_mountain_all_names_includes_localized_name() -> None:
    names = list_mountain_all_names()
    assert "Mount Fuji" in names
    assert "富士山" in names
    assert "Ayrybaba" in names


def test_main_list_mountains_prints_names_and_exits(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--list-mountains"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    lines = capsys.readouterr().out.splitlines()
    assert "Mount Fuji" in lines


def test_main_show_mountain_json_prints_json_and_exits(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-mountain-json", "富士山"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["qid"] == "Q39231"
    assert payload["name"] == "Mount Fuji"
    assert "富士山" in payload["names"]


def test_main_show_mountain_json_includes_ascii_name(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-mountain-json", "Ayrybaba"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["name"] == "Aýrybaba"
    assert payload["ascii_name"] == "Ayrybaba"


def test_main_show_mountain_json_returns_error_for_unknown_name(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-mountain-json", "No Such Mountain"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 1
    assert "No mountain found" in capsys.readouterr().err
