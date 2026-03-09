from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "dev-samples" / "build_wikidata_mountain_candidates.py"
    spec = importlib.util.spec_from_file_location("wikidata_mountain_candidates", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_raw_rows_accepts_plain_list() -> None:
    mod = _load_module()
    rows = [{"country": "http://www.wikidata.org/entity/Q17"}]
    assert mod.load_raw_rows(rows) == rows


def test_normalize_rows_collapses_duplicates_and_marks_reviews() -> None:
    mod = _load_module()
    rows = [
        {
            "country": "http://www.wikidata.org/entity/Q17",
            "countryLabel": "日本",
            "item": "http://www.wikidata.org/entity/Q39231",
            "itemLabel": "富士山",
            "coord": "Point(138.7275 35.360555555)",
            "elevation": "3777.24",
        },
        {
            "country": "http://www.wikidata.org/entity/Q17",
            "countryLabel": "日本",
            "item": "http://www.wikidata.org/entity/Q39231",
            "itemLabel": "富士山",
            "coord": "Point(138.7275 35.360555555)",
            "elevation": "3777.24",
        },
        {
            "country": "http://www.wikidata.org/entity/Q423",
            "countryLabel": "朝鮮民主主義人民共和国",
            "item": "http://www.wikidata.org/entity/Q107635",
            "itemLabel": "白頭山",
            "coord": "Point(126.83 40.94)",
            "elevation": "2744",
        },
        {
            "country": "http://www.wikidata.org/entity/Q423",
            "countryLabel": "朝鮮民主主義人民共和国",
            "item": "http://www.wikidata.org/entity/Q107635",
            "itemLabel": "白頭山",
            "coord": "Point(128.077222222 41.992777777)",
            "elevation": "2744",
        },
        {
            "country": "http://www.wikidata.org/entity/Q148",
            "countryLabel": "中華人民共和国",
            "item": "http://www.wikidata.org/entity/Q513",
            "itemLabel": "エベレスト",
            "coord": "Point(86.925 27.988055555)",
            "elevation": "8848.86",
        },
        {
            "country": "http://www.wikidata.org/entity/Q837",
            "countryLabel": "ネパール",
            "item": "http://www.wikidata.org/entity/Q513",
            "itemLabel": "エベレスト",
            "coord": "Point(86.925 27.988055555)",
            "elevation": "8848.86",
        },
    ]

    items = mod.normalize_rows(rows)

    assert len(items) == 4

    fuji = next(item for item in items if item["item_qid"] == "Q39231")
    assert fuji["row_count"] == 2
    assert fuji["coord_count"] == 1
    assert fuji["elevation_count"] == 1
    assert fuji["needs_coord_review"] is False
    assert fuji["needs_elevation_review"] is False

    paektu = next(item for item in items if item["item_qid"] == "Q107635")
    assert paektu["coord_count"] == 2
    assert paektu["needs_coord_review"] is True

    everest_cn = next(
        item for item in items if item["item_qid"] == "Q513" and item["country_qid"] == "Q148"
    )
    assert everest_cn["is_shared_peak"] is True
    assert everest_cn["shared_with_country_qids"] == ["Q837"]


def test_merge_entity_labels_adds_item_and_country_labels() -> None:
    mod = _load_module()
    items = [
        {
            "country_qid": "Q17",
            "country_label": "日本",
            "item_qid": "Q39231",
            "item_label": "富士山",
            "row_count": 1,
            "coordinates": [{"latitude_deg": 35.3605, "longitude_deg": 138.7275}],
            "elevations_m": [3777.24],
            "latitude_deg": 35.3605,
            "longitude_deg": 138.7275,
            "location_arg": "35.360500;138.727500",
            "max_elevation_m": 3777.24,
            "coord_count": 1,
            "elevation_count": 1,
            "needs_coord_review": False,
            "needs_elevation_review": False,
            "is_shared_peak": False,
            "shared_with_country_qids": [],
            "wikidata_url": "https://www.wikidata.org/wiki/Q39231",
        }
    ]

    enriched = mod.merge_entity_labels(
        items,
        item_labels_by_qid={"Q39231": {"en": "Mount Fuji", "ja": "富士山"}},
        country_labels_by_qid={"Q17": {"en": "Japan", "ja": "日本"}},
    )

    assert enriched[0]["item_labels"] == {"en": "Mount Fuji", "ja": "富士山"}
    assert enriched[0]["country_labels"] == {"en": "Japan", "ja": "日本"}
