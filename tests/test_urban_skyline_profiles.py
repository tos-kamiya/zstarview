from __future__ import annotations

import json
from pathlib import Path

from zstarview.urban_skyline_profiles import (
    load_urban_skyline_profiles,
    resolve_urban_skyline_profile_for_city_name,
)


def test_load_urban_skyline_profiles_reads_altaz_pairs(tmp_path: Path) -> None:
    path = tmp_path / "urban.json"
    path.write_text(
        json.dumps(
            {
                "wikidata:Q57965": {
                    "name": "Tokyo Skytree",
                    "profiles": [
                        {
                            "radius_km": 0.1,
                            "profile": [
                                {"az": 0.0, "alt": -1.0},
                                {"az": 0.5, "alt": -2.0},
                            ],
                        },
                        {
                            "radius_km": 1.0,
                            "profile": [
                                {"az": 0.0, "alt": -0.5},
                            ],
                        },
                    ],
                }
            }
        ),
        encoding="utf-8",
    )

    load_urban_skyline_profiles.cache_clear()
    got = load_urban_skyline_profiles(path)

    assert got == {
        "wikidata:Q57965": [
            (0.1, [(-1.0, 0.0), (-2.0, 0.5)]),
            (1.0, [(-0.5, 0.0)]),
        ]
    }


def test_resolve_urban_skyline_profile_for_tower_city_name(tmp_path: Path) -> None:
    path = tmp_path / "urban.json"
    path.write_text(
        json.dumps(
            {
                "wikidata:Q57965": {
                    "name": "Tokyo Skytree",
                    "profiles": [
                        {
                            "radius_km": 0.1,
                            "profile": [
                                {"az": 0.0, "alt": -1.0},
                            ],
                        }
                    ],
                }
            }
        ),
        encoding="utf-8",
    )

    load_urban_skyline_profiles.cache_clear()
    got = resolve_urban_skyline_profile_for_city_name("t/Tokyo Skytree", path)

    assert got == [(0.1, [(-1.0, 0.0)])]
