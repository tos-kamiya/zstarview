from __future__ import annotations

import json

from zstarview.tower_viewpoints import load_tower_viewpoints, resolve_tower_viewpoint


def test_load_tower_viewpoints_includes_tokyo_skytree() -> None:
    towers = load_tower_viewpoints()
    assert any(tower.qid == "Q57965" for tower in towers)


def test_resolve_tower_viewpoint_by_english_name() -> None:
    tower = resolve_tower_viewpoint("Tokyo Skytree")
    assert tower is not None
    assert tower.qid == "Q57965"
    assert tower.height_m == 634.0


def test_resolve_tower_viewpoint_by_ascii_fallback_name() -> None:
    tower = resolve_tower_viewpoint("Tsutenkaku")
    assert tower is not None
    assert tower.name == "Tsūtenkaku"


def test_resolve_tower_viewpoint_by_short_brand_alias() -> None:
    tower = resolve_tower_viewpoint("i360")
    assert tower is not None
    assert tower.name == "Brighton i360"


def test_resolve_tower_viewpoint_by_wikidata_key() -> None:
    tower = resolve_tower_viewpoint("wikidata:Q57965")
    assert tower is not None
    assert tower.name == "Tokyo Skytree"


def test_load_tower_viewpoints_defaults_viewpoint_height_to_height() -> None:
    tower = resolve_tower_viewpoint("Tokyo Skytree")
    assert tower is not None
    assert tower.viewpoint_height_m == tower.height_m


def test_resolve_tower_viewpoint_by_hydro_quebec_name() -> None:
    tower = resolve_tower_viewpoint("Tour d’observation Hydro-Québec de la Cité de l’énergie")
    assert tower is not None
    assert tower.qid == "Q137673602"


def test_removed_qid_only_tower_is_absent() -> None:
    towers = load_tower_viewpoints()
    assert all(tower.qid != "Q12049950" for tower in towers)


def test_resolve_tower_viewpoint_does_not_match_numeric_only_partial() -> None:
    assert resolve_tower_viewpoint("138") is None


def test_load_tower_viewpoints_ignores_legacy_observer_height_key(tmp_path) -> None:
    payload = {
        "items": [
            {
                "id": "wikidata:Q1",
                "qid": "Q1",
                "name": "Legacy Tower",
                "labels": {"en": "Legacy Tower"},
                "names": ["Legacy Tower"],
                "latitude_deg": 35.0,
                "longitude_deg": 139.0,
                "height_m": 120.0,
                "observer_height_m": 45.0,
            }
        ]
    }
    path = tmp_path / "tower_viewpoints.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    towers = load_tower_viewpoints(path)

    assert len(towers) == 1
    assert towers[0].viewpoint_height_m == 120.0
