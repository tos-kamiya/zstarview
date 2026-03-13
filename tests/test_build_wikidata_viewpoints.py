from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "build_wikidata_viewpoints.py"
    spec = importlib.util.spec_from_file_location("build_wikidata_viewpoints", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_normalize_extra_rows_deduplicates_qid_and_keeps_max_height() -> None:
    mod = _load_module()
    rows = [
        {
            "item": "http://www.wikidata.org/entity/Q1",
            "itemLabel": "Example Deck",
            "coord": "Point(139.7 35.6)",
            "height": "250",
        },
        {
            "item": "http://www.wikidata.org/entity/Q1",
            "itemLabel": "Example Deck",
            "coord": "Point(139.7 35.6)",
            "height": "300",
        },
    ]

    items = mod.normalize_extra_rows(rows)

    assert len(items) == 1
    item = items[0]
    assert item["qid"] == "Q1"
    assert item["height_m"] == 300.0
    assert item["viewpoint_height_m"] == 300.0
    assert item["viewpoint_types"] == ["tower_or_high_observation"]


def test_normalize_extra_rows_drops_zero_height_candidates() -> None:
    mod = _load_module()
    rows = [
        {
            "item": "http://www.wikidata.org/entity/Q1",
            "itemLabel": "Zero Height Deck",
            "coord": "Point(139.7 35.6)",
            "height": "0",
        },
        {
            "item": "http://www.wikidata.org/entity/Q2",
            "itemLabel": "Missing Height Deck",
            "coord": "Point(139.8 35.7)",
        },
        {
            "item": "http://www.wikidata.org/entity/Q3",
            "itemLabel": "Valid Deck",
            "coord": "Point(139.9 35.8)",
            "height": "123",
        },
    ]

    items = mod.normalize_extra_rows(rows)

    assert [item["qid"] for item in items] == ["Q3"]


def test_merge_items_adds_observer_height_to_existing_and_preserves_known_item() -> None:
    mod = _load_module()
    base_items = [
        {
            "id": "wikidata:Q57965",
            "qid": "Q57965",
            "name": "Tokyo Skytree",
            "names": ["Tokyo Skytree"],
            "labels": {"en": "Tokyo Skytree"},
            "latitude_deg": 35.710055555,
            "longitude_deg": 139.810722222,
            "height_m": 634.0,
            "viewpoint_height_m": 634.0,
        }
    ]
    extra_items = [
        {
            "id": "wikidata:Q57965",
            "qid": "Q57965",
            "name": "東京スカイツリー",
            "names": ["東京スカイツリー"],
            "labels": {},
            "latitude_deg": 35.710055555,
            "longitude_deg": 139.810722222,
            "height_m": 620.0,
            "viewpoint_height_m": 620.0,
            "viewpoint_types": ["tower_or_high_observation"],
        },
        {
            "id": "wikidata:Q16318627",
            "qid": "Q16318627",
            "name": "あべのハルカス",
            "names": ["あべのハルカス"],
            "labels": {},
            "latitude_deg": 34.645947222,
            "longitude_deg": 135.514266666,
            "height_m": 300.0,
            "viewpoint_height_m": 300.0,
            "viewpoint_types": ["tower_or_high_observation"],
        },
    ]

    merged = mod.merge_items(base_items, extra_items)

    assert [item["qid"] for item in merged] == ["Q57965", "Q16318627"]
    tokyo_skytree = merged[0]
    assert tokyo_skytree["viewpoint_height_m"] == 634.0
    assert tokyo_skytree["names"] == ["Tokyo Skytree", "東京スカイツリー"]
    abeno = merged[1]
    assert abeno["viewpoint_height_m"] == 300.0


def test_merge_entity_labels_prefers_fetched_english_name() -> None:
    mod = _load_module()
    items = [
        {
            "id": "wikidata:Q16318627",
            "qid": "Q16318627",
            "name": "あべのハルカス",
            "names": ["あべのハルカス"],
            "labels": {},
            "latitude_deg": 34.645947222,
            "longitude_deg": 135.514266666,
            "height_m": 300.0,
            "viewpoint_height_m": 300.0,
            "viewpoint_types": ["tower_or_high_observation"],
            "slug": "あべのハルカス",
        }
    ]

    got = mod.merge_entity_labels(
        items,
        {"Q16318627": {"en": "Abeno Harukas", "ja": "あべのハルカス"}},
    )

    assert len(got) == 1
    item = got[0]
    assert item["name"] == "Abeno Harukas"
    assert item["labels"] == {"en": "Abeno Harukas", "ja": "あべのハルカス"}
    assert item["names"] == ["Abeno Harukas", "あべのハルカス"]
