from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "dev-samples" / "build_wikidata_tower_viewpoints.py"
    spec = importlib.util.spec_from_file_location("wikidata_tower_viewpoints", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_normalize_rows_deduplicates_and_keeps_max_height() -> None:
    mod = _load_module()
    rows = [
        {
            "item": "http://www.wikidata.org/entity/Q57965",
            "itemLabel": "東京スカイツリー",
            "itemLabel_en": "Tokyo Skytree",
            "class": "http://www.wikidata.org/entity/Q1068623",
            "classLabel": "電波塔",
            "coord": "Point(139.810722222 35.710055555)",
            "height": "634",
        },
        {
            "item": "http://www.wikidata.org/entity/Q57965",
            "itemLabel": "東京スカイツリー",
            "itemLabel_en": "Tokyo Skytree",
            "class": "http://www.wikidata.org/entity/Q1440300",
            "classLabel": "展望塔",
            "coord": "Point(139.810722222 35.710055555)",
            "height": "620",
        },
    ]

    towers = mod.normalize_rows(rows, min_height_m=100.0, include_lookouts=False)

    assert len(towers) == 1
    tower = towers[0]
    assert tower["qid"] == "Q57965"
    assert tower["height_m"] == 634.0
    assert tower["name"] == "Tokyo Skytree"
    assert tower["names"] == ["Tokyo Skytree", "東京スカイツリー"]
    assert tower["labels"] == {"en": "Tokyo Skytree"}
    assert tower["location_arg"] == "35.710056;139.810722"
    assert tower["classes"] == ["展望塔", "電波塔"]
    assert tower["class_qids"] == ["Q1068623", "Q1440300"]


def test_normalize_rows_filters_fire_lookouts_by_default() -> None:
    mod = _load_module()
    rows = [
        {
            "item": "http://www.wikidata.org/entity/Q1",
            "itemLabel": "Example Fire Lookout Tower",
            "class": "http://www.wikidata.org/entity/Q1440300",
            "classLabel": "展望塔",
            "coord": "Point(10.0 20.0)",
            "height": "120",
        },
        {
            "item": "http://www.wikidata.org/entity/Q2",
            "itemLabel": "Example City Tower",
            "class": "http://www.wikidata.org/entity/Q11166728",
            "classLabel": "テレビ塔",
            "coord": "Point(30.0 40.0)",
            "height": "200",
        },
    ]

    towers = mod.normalize_rows(rows, min_height_m=100.0, include_lookouts=False)

    assert towers == []


def test_normalize_rows_requires_observation_tower_by_default() -> None:
    mod = _load_module()
    rows = [
        {
            "item": "http://www.wikidata.org/entity/Q1",
            "itemLabel": "Only TV Tower",
            "class": "http://www.wikidata.org/entity/Q11166728",
            "classLabel": "テレビ塔",
            "coord": "Point(10.0 20.0)",
            "height": "300",
        },
        {
            "item": "http://www.wikidata.org/entity/Q2",
            "itemLabel": "Observation Tower",
            "class": "http://www.wikidata.org/entity/Q1440300",
            "classLabel": "展望塔",
            "coord": "Point(30.0 40.0)",
            "height": "200",
        },
    ]

    towers = mod.normalize_rows(rows, min_height_m=100.0, include_lookouts=False)

    assert [tower["qid"] for tower in towers] == ["Q2"]


def test_normalize_rows_keeps_only_english_and_japanese_labels() -> None:
    mod = _load_module()
    rows = [
        {
            "item": "http://www.wikidata.org/entity/Q57965",
            "itemLabel": "東京スカイツリー",
            "itemLabel_en": "Tokyo Skytree",
            "itemLabel_ja": "東京スカイツリー",
            "itemLabel_fr": "Tokyo Skytree FR",
            "class": "http://www.wikidata.org/entity/Q1440300",
            "classLabel": "展望塔",
            "coord": "Point(139.810722222 35.710055555)",
            "height": "634",
        },
    ]

    towers = mod.normalize_rows(rows, min_height_m=100.0, include_lookouts=False)

    assert len(towers) == 1
    tower = towers[0]
    assert tower["name"] == "Tokyo Skytree"
    assert tower["labels"] == {"en": "Tokyo Skytree", "ja": "東京スカイツリー"}
    assert tower["names"] == ["Tokyo Skytree", "東京スカイツリー"]


def test_merge_entity_labels_prefers_fetched_english_label() -> None:
    mod = _load_module()
    towers = [
        {
            "id": "wikidata:Q57965",
            "name": "東京スカイツリー",
            "names": ["東京スカイツリー"],
            "labels": {},
            "qid": "Q57965",
            "latitude_deg": 35.710055555,
            "longitude_deg": 139.810722222,
            "height_m": 634.0,
            "classes": ["展望塔"],
            "class_qids": ["Q1440300"],
            "wikidata_url": "https://www.wikidata.org/wiki/Q57965",
            "location_arg": "35.710056;139.810722",
            "slug": "東京スカイツリー",
        }
    ]

    enriched = mod.merge_entity_labels(
        towers,
        labels_by_qid={"Q57965": {"en": "Tokyo Skytree", "ja": "東京スカイツリー"}},
        include_lookouts=False,
    )

    assert len(enriched) == 1
    tower = enriched[0]
    assert tower["name"] == "Tokyo Skytree"
    assert tower["labels"] == {"en": "Tokyo Skytree", "ja": "東京スカイツリー"}
    assert tower["names"] == ["Tokyo Skytree", "東京スカイツリー"]
