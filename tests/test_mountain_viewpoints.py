from __future__ import annotations

from zstarview.mountain_viewpoints import (
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


def test_resolve_mountain_viewpoint_by_localized_name() -> None:
    mountain = resolve_mountain_viewpoint("富士山")
    assert mountain is not None
    assert mountain.name == "Mount Fuji"


def test_resolve_mountain_viewpoint_by_ascii_fallback_name() -> None:
    mountain = resolve_mountain_viewpoint("Ayrybaba")
    assert mountain is not None
    assert mountain.name == "Aýrybaba"


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
