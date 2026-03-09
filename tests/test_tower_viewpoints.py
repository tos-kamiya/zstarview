from __future__ import annotations

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


def test_resolve_tower_viewpoint_by_wikidata_key() -> None:
    tower = resolve_tower_viewpoint("wikidata:Q57965")
    assert tower is not None
    assert tower.name == "Tokyo Skytree"


def test_resolve_tower_viewpoint_by_hydro_quebec_name() -> None:
    tower = resolve_tower_viewpoint("Tour d’observation Hydro-Québec de la Cité de l’énergie")
    assert tower is not None
    assert tower.qid == "Q137673602"


def test_removed_qid_only_tower_is_absent() -> None:
    towers = load_tower_viewpoints()
    assert all(tower.qid != "Q12049950" for tower in towers)
