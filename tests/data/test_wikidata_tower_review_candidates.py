from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[2]
    mod_path = root / "dev-samples" / "build_wikidata_tower_review_candidates.py"
    spec = importlib.util.spec_from_file_location(
        "wikidata_tower_review_candidates",
        mod_path,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_normalize_candidate_rows_unions_rules_and_keeps_max_height() -> None:
    mod = _load_module()
    rows_by_input = [
        (
            Path("raw-data/query-1.json"),
            [
                {
                    "item": "http://www.wikidata.org/entity/Q1",
                    "itemLabel": "Example Tower",
                    "coord": "Point(139.7 35.6)",
                    "height": "250",
                }
            ],
        ),
        (
            Path("raw-data/query-2.json"),
            [
                {
                    "item": "http://www.wikidata.org/entity/Q1",
                    "itemLabel": "Example Tower",
                    "coord": "Point(139.7 35.6)",
                    "height": "300",
                }
            ],
        ),
    ]

    items = mod.normalize_candidate_rows(rows_by_input, min_height_m=100.0)

    assert len(items) == 1
    item = items[0]
    assert item["qid"] == "Q1"
    assert item["height_m"] == 300.0
    assert item["viewpoint_height_m"] == 300.0
    assert item["matched_rules"] == [
        "has_use_observation_deck",
        "instance_of_observation_deck",
        "instance_of_observation_tower",
    ]
    assert item["source_files"] == ["raw-data/query-1.json", "raw-data/query-2.json"]


def test_normalize_candidate_rows_respects_explicit_rule_override() -> None:
    mod = _load_module()
    rows_by_input = [
        (
            Path("raw-data/query-custom.json"),
            [
                {
                    "item": "http://www.wikidata.org/entity/Q2",
                    "itemLabel": "Custom Tower",
                    "coord": "Point(139.8 35.7)",
                    "height": "180",
                }
            ],
        ),
    ]

    items = mod.normalize_candidate_rows(
        rows_by_input,
        min_height_m=100.0,
        input_rule_overrides={
            "raw-data/query-custom.json": (
                "has_part_observation_deck",
                "member_of_all_japan_tower_association",
            )
        },
    )

    assert len(items) == 1
    item = items[0]
    assert item["matched_rules"] == [
        "has_part_observation_deck",
        "member_of_all_japan_tower_association",
    ]


def test_normalize_candidate_rows_prefers_real_name_over_qid_fallback() -> None:
    mod = _load_module()
    rows_by_input = [
        (
            Path("raw-data/query-3.json"),
            [
                {
                    "item": "http://www.wikidata.org/entity/Q3",
                    "coord": "Point(139.9 35.8)",
                    "height": "200",
                }
            ],
        ),
        (
            Path("raw-data/query-1.json"),
            [
                {
                    "item": "http://www.wikidata.org/entity/Q3",
                    "itemLabel": "Named Tower",
                    "coord": "Point(139.9 35.8)",
                    "height": "180",
                }
            ],
        ),
    ]

    items = mod.normalize_candidate_rows(rows_by_input, min_height_m=100.0)

    assert len(items) == 1
    item = items[0]
    assert item["name"] == "Named Tower"
    assert item["names"] == ["Named Tower"]


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
            "matched_rules": ["has_use_observation_deck"],
            "source_files": ["raw-data/query-2.json"],
            "viewpoint_types": ["tower_or_high_observation"],
            "wikidata_url": "https://www.wikidata.org/wiki/Q16318627",
            "location_arg": "34.645947;135.514267",
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
