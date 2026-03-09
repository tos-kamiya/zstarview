from __future__ import annotations

import json

import pytest

from zstarview.tower_viewpoints import list_tower_all_names, list_tower_primary_names
from zstarview.zstarview import main, parse_args


def test_parse_args_accepts_list_towers() -> None:
    args = parse_args(["--list-towers"])
    assert args.list_towers is True
    assert args.list_tower_names is False
    assert args.show_tower_json is None


def test_parse_args_rejects_location_with_tower_listing() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-towers", "Tokyo Skytree"])


def test_parse_args_rejects_render_option_with_tower_listing() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-towers", "--hours", "1"])


def test_list_tower_primary_names_includes_tokyo_skytree() -> None:
    names = list_tower_primary_names()
    assert "Tokyo Skytree" in names


def test_list_tower_primary_names_prefers_ascii_fallback() -> None:
    names = list_tower_primary_names()
    assert "Tsutenkaku" in names
    assert "Tsūtenkaku" not in names


def test_list_tower_primary_names_excludes_qid_placeholders() -> None:
    names = list_tower_primary_names()
    assert "Q12049950" not in names
    assert "Q137673602" not in names


def test_list_tower_all_names_includes_localized_name() -> None:
    names = list_tower_all_names()
    assert "Tokyo Skytree" in names
    assert "東京スカイツリー" in names
    assert "Tsutenkaku" in names


def test_list_tower_all_names_excludes_qid_placeholders() -> None:
    names = list_tower_all_names()
    assert "Q12049950" not in names
    assert "Q137673602" not in names


def test_main_list_towers_prints_names_and_exits(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--list-towers"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    lines = capsys.readouterr().out.splitlines()
    assert "Tokyo Skytree" in lines


def test_main_show_tower_json_prints_json_and_exits(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-tower-json", "東京スカイツリー"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["qid"] == "Q57965"
    assert payload["kind"] == "tower"
    assert payload["name"] == "Tokyo Skytree"
    assert "東京スカイツリー" in payload["names"]
    assert payload["height_m"] == 634.0
    assert payload["meta"]["wikidata_url"] == "https://www.wikidata.org/wiki/Q57965"
    assert "Q1440300" in payload["meta"]["class_qids"]


def test_main_show_tower_json_includes_ascii_name(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-tower-json", "Tsutenkaku"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["name"] == "Tsūtenkaku"
    assert payload["ascii_name"] == "Tsutenkaku"
    assert "meta" in payload


def test_main_show_tower_json_returns_error_for_unknown_name(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-tower-json", "No Such Tower"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 1
    assert "No tower found" in capsys.readouterr().err
