from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "dev-samples" / "build_curated_mountain_viewpoints.py"
    spec = importlib.util.spec_from_file_location("curated_mountain_viewpoints", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_normalize_seed_candidates_deduplicates_qids() -> None:
    mod = _load_module()
    payload = {
        "candidates": [
            {
                "qid": "Q39261",
                "name": "Mount Fuji",
                "names": ["Fuji"],
                "labels": {"ja": "富士山"},
                "country_codes": ["jp"],
            },
            {
                "qid": "wikidata:Q39261",
                "names": ["富士山"],
                "region_hint": "Japan",
            },
        ]
    }

    candidates = mod.normalize_seed_candidates(payload)

    assert len(candidates) == 1
    candidate = candidates[0]
    assert candidate["qid"] == "Q39261"
    assert candidate["seed_name"] == "Mount Fuji"
    assert candidate["seed_names"] == ["Fuji", "Mount Fuji", "富士山"]
    assert candidate["seed_labels"] == {"ja": "富士山"}
    assert candidate["country_codes"] == ["JP"]
    assert candidate["region_hint"] == "Japan"


def test_merge_entity_data_prefers_fetched_english_label() -> None:
    mod = _load_module()
    candidates = [
        {
            "qid": "Q39261",
            "seed_name": "富士山",
            "seed_names": ["富士山"],
            "seed_labels": {"ja": "富士山"},
            "latitude_deg": 35.3606,
            "longitude_deg": 138.7274,
            "elevation_m": 3776.0,
            "wikipedia_titles": {"ja": "富士山"},
        }
    ]
    entities = {
        "Q39261": {
            "labels": {"en": "Mount Fuji", "ja": "富士山"},
            "aliases": {"Fuji-san"},
            "coordinate": (35.3606, 138.7274),
            "elevation_m": 3776.24,
            "wikipedia_urls": {"en": "https://en.wikipedia.org/wiki/Mount_Fuji"},
        }
    }

    mountains = mod.merge_entity_data(
        candidates,
        entities,
        allow_missing_elevation=False,
    )

    assert len(mountains) == 1
    mountain = mountains[0]
    assert mountain["name"] == "Mount Fuji"
    assert mountain["labels"] == {"en": "Mount Fuji", "ja": "富士山"}
    assert "Fuji-san" in mountain["names"]
    assert mountain["elevation_m"] == 3776.0
    assert mountain["wikipedia_urls"] == {
        "en": "https://en.wikipedia.org/wiki/Mount_Fuji",
        "ja": "https://ja.wikipedia.org/wiki/富士山",
    }


def test_merge_entity_data_uses_fetched_coordinate_and_elevation() -> None:
    mod = _load_module()
    candidates = [
        {
            "qid": "Q1",
            "seed_name": "Example Peak",
            "seed_names": ["Example Peak"],
            "seed_labels": {},
            "wikipedia_titles": {},
        }
    ]
    entities = {
        "Q1": {
            "labels": {"en": "Example Peak"},
            "aliases": set(),
            "coordinate": (10.5, 20.25),
            "elevation_m": 1234.5,
            "wikipedia_urls": {},
        }
    }

    mountains = mod.merge_entity_data(
        candidates,
        entities,
        allow_missing_elevation=False,
    )

    mountain = mountains[0]
    assert mountain["latitude_deg"] == 10.5
    assert mountain["longitude_deg"] == 20.25
    assert mountain["elevation_m"] == 1234.5
    assert mountain["location_arg"] == "10.500000;20.250000"


def test_merge_entity_data_cleans_problematic_aliases() -> None:
    mod = _load_module()
    candidates = [
        {
            "qid": "Q5059",
            "seed_name": "Aoraki / Mount Cook",
            "seed_names": ["Aoraki/Mount Cook", "Mount Cook", "Aoraki"],
            "seed_labels": {"en": "Aoraki / Mount Cook"},
            "latitude_deg": -43.595,
            "longitude_deg": 170.141944,
            "elevation_m": 3724.0,
            "wikipedia_titles": {},
        }
    ]
    entities = {
        "Q5059": {
            "labels": {"en": "Aoraki / Mount Cook"},
            "aliases": {
                "Jade Mtn.",
                "Mt. Jade",
                "3003",
                "🗻",
                'Desde la ruta 39, al km 47,5 pueden entrar en los senderos de "Las Cañas", y llegar hasta el Cerro Catedral. Sino entran en el km 60 y listo. Una maravilla!',
            },
            "coordinate": (-43.595, 170.141944),
            "elevation_m": 3724.0,
            "wikipedia_urls": {},
        }
    }

    mountains = mod.merge_entity_data(
        candidates,
        entities,
        allow_missing_elevation=False,
    )

    mountain = mountains[0]
    assert mountain["name"] == "Aoraki / Mount Cook"
    assert "Aoraki/Mount Cook" in mountain["names"]
    assert "Aoraki / Mount Cook" not in mountain["names"]
    assert "Jade Mountain" in mountain["names"]
    assert "Mt. Jade" not in mountain["names"]
    assert "3003" not in mountain["names"]
    assert "🗻" not in mountain["names"]
