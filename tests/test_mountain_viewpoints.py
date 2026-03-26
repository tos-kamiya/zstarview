from __future__ import annotations

import pytest

from zstarview.mountain_viewpoints import (
    list_mountain_all_names,
    load_mountain_viewpoints,
    resolve_mountain_viewpoint,
)


def test_load_mountain_viewpoints_includes_mount_fuji() -> None:
    mountains = load_mountain_viewpoints()
    assert any(mountain.qid == "Q39231" for mountain in mountains)


def test_resolve_mountain_viewpoint_by_english_name() -> None:
    mountain = resolve_mountain_viewpoint("Mount Fuji")
    assert mountain is not None
    assert mountain.qid == "Q39231"
    assert mountain.meta["elevation_m"] == 3777.24


def test_resolve_mountain_viewpoint_mount_everest() -> None:
    mountain = resolve_mountain_viewpoint("Mount Everest")
    assert mountain is not None
    assert mountain.qid == "Q513"
    assert mountain.latitude_deg == 27.988055555
    assert mountain.longitude_deg == 86.925
    assert mountain.meta["elevation_m"] == 8848.86


def test_resolve_mountain_viewpoint_by_everest_alias() -> None:
    mountain = resolve_mountain_viewpoint("Sagarmatha")
    assert mountain is not None
    assert mountain.name == "Mount Everest"


@pytest.mark.parametrize(
    ("query", "expected_name"),
    [
        ("Kilimanjaro", "Mount Kilimanjaro"),
        ("Mount Vinson", "Mount Vinson"),
        ("Kosciuszko", "Mount Kosciuszko"),
        ("Mont Blanc", "Mont Blanc"),
        ("Matterhorn", "Matterhorn"),
        ("Lhotse", "Lhotse"),
        ("Makalu", "Makalu"),
        ("Cho Oyu", "Cho Oyu"),
        ("Dhaulagiri I", "Dhaulagiri"),
        ("Manaslu", "Manaslu"),
        ("Nanga Parbat", "Nanga Parbat"),
        ("Annapurna I", "Annapurna I"),
    ],
)
def test_resolve_mountain_viewpoint_for_added_famous_peaks(
    query: str,
    expected_name: str,
) -> None:
    mountain = resolve_mountain_viewpoint(query)

    assert mountain is not None
    assert mountain.name == expected_name


def test_resolve_mountain_viewpoint_by_localized_name() -> None:
    mountain = resolve_mountain_viewpoint("富士山")
    assert mountain is not None
    assert mountain.name == "Mount Fuji"


def test_resolve_mountain_viewpoint_by_ascii_fallback_name() -> None:
    mountain = resolve_mountain_viewpoint("Ayrybaba")
    assert mountain is not None
    assert mountain.name == "Aýrybaba"


def test_resolve_mountain_viewpoint_by_lowercase_alias() -> None:
    mountain = resolve_mountain_viewpoint("pueyosa")
    assert mountain is not None
    assert mountain.name == "Mount Wilhelm"


def test_resolve_mountain_viewpoint_by_wikidata_key() -> None:
    mountain = resolve_mountain_viewpoint("wikidata:Q39231")
    assert mountain is not None
    assert mountain.name == "Mount Fuji"


def test_resolve_pico_cristobal_colon_prefers_remaining_entry() -> None:
    mountain = resolve_mountain_viewpoint("Pico Cristóbal Colón")
    assert mountain is not None
    assert mountain.qid == "Q578862"


def test_removed_duplicate_mountain_is_absent() -> None:
    mountains = load_mountain_viewpoints()
    assert all(mountain.qid != "Q338790" for mountain in mountains)


def test_list_mountain_all_names_excludes_cleaned_problem_aliases() -> None:
    names = list_mountain_all_names()
    assert "=" not in names
    assert "==" not in names
    assert "🗻" not in names
    assert "3003" not in names
    assert "Aoraki / Mount Cook" not in names
    assert "Aoraki/Mount Cook" in names
    assert "Mt. Fuji" not in names
    assert "Mt Fuji" not in names
