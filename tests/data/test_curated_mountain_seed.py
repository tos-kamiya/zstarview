from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[2]
    mod_path = root / "dev-samples" / "build_curated_mountain_seed.py"
    spec = importlib.util.spec_from_file_location("curated_mountain_seed", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_should_exclude_respects_review_flags_and_name_tokens() -> None:
    mod = _load_module()
    assert mod.should_exclude({"needs_coord_review": True}) is True
    assert mod.should_exclude({"needs_elevation_review": True}) is True
    assert mod.should_exclude({"is_shared_peak": True}) is True
    assert mod.should_exclude({"item_label": "Example Hill"}) is True
    assert mod.should_exclude({"item_labels": {"en": "Example Ridge"}}) is True
    assert mod.should_exclude({"item_label": "Mount Fuji"}) is False


def test_build_seed_keeps_only_approved_candidates() -> None:
    mod = _load_module()
    payload = {
        "source_query_result": "dev-samples/raw.json",
        "items": [
            {
                "country_label": "日本",
                "country_labels": {"en": "Japan", "ja": "日本"},
                "item_qid": "Q39231",
                "item_label": "富士山",
                "item_labels": {"en": "Mount Fuji", "ja": "富士山"},
                "latitude_deg": 35.3605,
                "longitude_deg": 138.7275,
                "max_elevation_m": 3777.24,
                "needs_coord_review": False,
                "needs_elevation_review": False,
                "is_shared_peak": False,
            },
            {
                "country_label": "日本",
                "item_qid": "Q1",
                "item_label": "Example Hill",
                "latitude_deg": 1.0,
                "longitude_deg": 2.0,
                "max_elevation_m": 10.0,
                "needs_coord_review": False,
                "needs_elevation_review": False,
                "is_shared_peak": False,
            },
        ],
    }

    output = mod.build_seed(payload)

    assert output["source_review_result"] == "dev-samples/raw.json"
    assert len(output["candidates"]) == 1
    candidate = output["candidates"][0]
    assert candidate["qid"] == "Q39231"
    assert candidate["name"] == "Mount Fuji"
    assert candidate["names"] == ["Mount Fuji", "富士山"]
    assert candidate["region_hint"] == "Japan"
