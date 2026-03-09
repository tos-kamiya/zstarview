from __future__ import annotations

import json

import pytest

from zstarview.mountain_viewpoints import list_mountain_all_names, list_mountain_primary_names
from zstarview.zstarview import main, parse_args


def test_parse_args_accepts_list_viewpoints_for_mountains() -> None:
    args = parse_args(["--list-viewpoints", "m"])
    assert args.list_viewpoints == "m"
    assert args.list_viewpoint_names is None
    assert args.show_viewpoint_json is None


def test_parse_args_accepts_list_viewpoint_names_for_mountains() -> None:
    args = parse_args(["--list-viewpoint-names", "m"])
    assert args.list_viewpoint_names == "m"


def test_parse_args_rejects_location_with_mountain_listing() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-viewpoints", "m", "Mount Fuji"])


def test_parse_args_rejects_render_option_with_mountain_listing() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-viewpoints", "m", "--hours", "1"])


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


def test_main_list_viewpoints_for_mountains_prints_prefixed_names_and_exits(
    capsys, monkeypatch
) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--list-viewpoints", "m"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    lines = capsys.readouterr().out.splitlines()
    assert "m/Mount Fuji" in lines


def test_main_list_viewpoint_names_for_mountains_prints_prefixed_names_and_exits(
    capsys, monkeypatch
) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--list-viewpoint-names", "m"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    lines = capsys.readouterr().out.splitlines()
    assert "m/Mount Fuji" in lines
    assert "m/富士山" in lines


def test_main_show_viewpoint_json_prints_mountain_json_and_exits(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "m/富士山"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["qid"] == "Q39231"
    assert payload["kind"] == "mountain"
    assert payload["name"] == "Mount Fuji"
    assert "富士山" in payload["names"]
    assert payload["height_m"] == 0.0
    assert payload["meta"]["elevation_m"] == 3777.24


def test_main_show_viewpoint_json_resolves_prefixless_mountain(
    capsys, monkeypatch
) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "Mount Fuji"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["kind"] == "mountain"


def test_main_show_viewpoint_json_includes_mountain_ascii_name(
    capsys, monkeypatch
) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "m/Ayrybaba"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["name"] == "Aýrybaba"
    assert payload["ascii_name"] == "Ayrybaba"


def test_main_show_viewpoint_json_returns_error_for_unknown_mountain(
    capsys, monkeypatch
) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "m/No Such Mountain"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 1
    assert "No mountain found" in capsys.readouterr().err
