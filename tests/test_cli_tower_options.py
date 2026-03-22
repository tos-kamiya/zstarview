from __future__ import annotations

import json

import pytest

from zstarview.tower_viewpoints import list_tower_all_names, list_tower_primary_names
from zstarview.viewpoints import Viewpoint
from zstarview.gui.viewer import main, parse_args


def test_parse_args_accepts_list_viewpoints_for_towers() -> None:
    args = parse_args(["--list-viewpoints", "t"])
    assert args.list_viewpoints == "t"
    assert args.list_viewpoint_names is None
    assert args.show_viewpoint_json is None


def test_parse_args_rejects_invalid_list_viewpoints_kind() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-viewpoints", "tower"])


def test_parse_args_rejects_location_with_viewpoint_listing() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-viewpoints", "t", "Tokyo Skytree"])


def test_parse_args_rejects_render_option_with_viewpoint_listing() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--list-viewpoints", "t", "--hours", "1"])


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


def test_main_list_viewpoints_for_towers_prints_prefixed_names_and_exits(
    capsys, monkeypatch
) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--list-viewpoints", "t"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    lines = capsys.readouterr().out.splitlines()
    assert "t/Tokyo Skytree" in lines


def test_main_list_viewpoint_names_for_towers_prints_prefixed_names_and_exits(
    capsys, monkeypatch
) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--list-viewpoint-names", "t"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    lines = capsys.readouterr().out.splitlines()
    assert "t/Tokyo Skytree" in lines
    assert "t/東京スカイツリー" in lines


def test_main_show_viewpoint_json_prints_tower_json_and_exits(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "t/東京スカイツリー"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["qid"] == "Q57965"
    assert payload["kind"] == "tower"
    assert payload["name"] == "Tokyo Skytree"
    assert "東京スカイツリー" in payload["names"]
    assert payload["height_m"] == 634.0
    assert payload["viewpoint_height_m"] == 634.0
    assert payload["meta"]["wikidata_url"] == "https://www.wikidata.org/wiki/Q57965"
    assert "Q1440300" in payload["meta"]["class_qids"]


def test_main_show_viewpoint_json_resolves_prefixless_tower(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "Tokyo Skytree"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["kind"] == "tower"


def test_main_show_viewpoint_json_includes_tower_ascii_name(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "t/Tsutenkaku"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["name"] == "Tsūtenkaku"
    assert payload["ascii_name"] == "Tsutenkaku"


def test_main_show_viewpoint_json_returns_error_for_unknown_tower(capsys, monkeypatch) -> None:
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "t/No Such Tower"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 1
    assert "No tower found" in capsys.readouterr().err


def test_main_show_viewpoint_json_reports_ambiguous_exact_matches(
    capsys, monkeypatch
) -> None:
    ambiguous_tower = Viewpoint(
        id="wikidata:Q1",
        qid="Q1",
        kind="tower",
        name="Shared Peak",
        labels={"en": "Shared Peak"},
        names=("Shared Peak",),
        latitude_deg=0.0,
        longitude_deg=0.0,
        height_m=100.0,
        viewpoint_height_m=100.0,
        meta={},
    )
    ambiguous_mountain = Viewpoint(
        id="wikidata:Q2",
        qid="Q2",
        kind="mountain",
        name="Shared Peak",
        labels={"en": "Shared Peak"},
        names=("Shared Peak",),
        latitude_deg=1.0,
        longitude_deg=1.0,
        height_m=0.0,
        viewpoint_height_m=None,
        meta={"elevation_m": 2000.0},
    )
    monkeypatch.setattr("zstarview.gui.viewer.load_tower_viewpoints", lambda: (ambiguous_tower,))
    monkeypatch.setattr(
        "zstarview.gui.viewer.load_mountain_viewpoints",
        lambda: (ambiguous_mountain,),
    )
    monkeypatch.setattr("zstarview.gui.viewer.resolve_tower_viewpoint", lambda _name: ambiguous_tower)
    monkeypatch.setattr(
        "zstarview.gui.viewer.resolve_mountain_viewpoint",
        lambda _name: ambiguous_mountain,
    )
    monkeypatch.setattr("sys.argv", ["zstarview", "--show-viewpoint-json", "Shared Peak"])
    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 1
    err = capsys.readouterr().err
    assert "Ambiguous viewpoint name" in err
    assert "t/Shared Peak" in err
    assert "m/Shared Peak" in err
