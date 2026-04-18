from __future__ import annotations

import pytest

from zstarview.search.models import SearchJumpTarget
from zstarview.search.query import parse_search_query, search_target_matches_query
from zstarview.search.resolver import resolve_search_targets


def test_parse_search_query_supports_label_selector() -> None:
    spec = parse_search_query("label=Ceres")

    assert spec.selector == "label"
    assert spec.value == "Ceres"
    assert spec.normalized == "ceres"


def test_parse_search_query_supports_id_selector() -> None:
    spec = parse_search_query("id=2000001")

    assert spec.selector == "id"
    assert spec.value == "2000001"
    assert spec.normalized == "2000001"


def test_search_target_matches_query_uses_label_and_id() -> None:
    target = SearchJumpTarget(
        label="Ceres",
        kind="jpl_small_body",
        sort_key=(0.0, "ceres"),
        object_key="2000001",
        command="DES=2000001;",
    )

    assert search_target_matches_query(target, parse_search_query("Ceres"))
    assert search_target_matches_query(target, parse_search_query("label=Ceres"))
    assert search_target_matches_query(target, parse_search_query("id=2000001"))


def test_resolve_search_targets_prefers_local_matches() -> None:
    local_targets = [
        SearchJumpTarget(
            label="Sirius",
            kind="star",
            sort_key=(0.0, "sirius"),
            subtitle="Vmag -1.44",
        ),
        SearchJumpTarget(
            label="Ceres",
            kind="jpl_small_body",
            sort_key=(0.0, "ceres"),
            object_key="2000001",
            command="DES=2000001;",
        ),
    ]
    resolution = resolve_search_targets(
        "Ceres",
        local_targets,
        jpl_search_callback=lambda _query: [],
    )

    assert len(resolution.candidates) == 1
    assert resolution.selected_target is not None
    assert resolution.selected_target.label == "Ceres"


def test_resolve_search_targets_prefers_known_satellites() -> None:
    local_targets = []
    satellite_calls: list[str] = []
    jpl_calls: list[str] = []

    def fake_satellite_search(query: str):
        satellite_calls.append(query)
        return [
            SearchJumpTarget(
                label="JWST",
                kind="satellite",
                sort_key=(0.0, "jwst"),
                object_key="JWST",
                alt_deg=10.0,
                az_deg=20.0,
            )
        ]

    def fake_jpl_search(query: str):
        jpl_calls.append(query)
        return []

    resolution = resolve_search_targets(
        "JWST",
        local_targets,
        satellite_search_callback=fake_satellite_search,
        jpl_search_callback=fake_jpl_search,
    )

    assert satellite_calls == ["JWST"]
    assert jpl_calls == []
    assert len(resolution.candidates) == 1
    assert resolution.selected_target is not None
    assert resolution.selected_target.label == "JWST"


def test_resolve_search_targets_propagates_missing_satellite_position() -> None:
    def fake_satellite_search(_query: str):
        raise RuntimeError("Satellite position unavailable for JWST")

    with pytest.raises(RuntimeError, match="Satellite position unavailable"):
        resolve_search_targets(
            "JWST",
            [],
            satellite_search_callback=fake_satellite_search,
            jpl_search_callback=lambda _query: [],
        )
