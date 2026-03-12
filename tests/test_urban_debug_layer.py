from __future__ import annotations

import json
from pathlib import Path

from zstarview.urban_debug_layer import (
    load_urban_debug_layers,
    resolve_urban_debug_layer_for_city_name,
)


def test_load_urban_debug_layers_reads_outline_polylines(tmp_path: Path) -> None:
    path = tmp_path / "debug.json"
    path.write_text(
        json.dumps(
            {
                "wikidata:Q57965": {
                    "name": "Tokyo Skytree",
                    "outlines": [
                        [
                            {"az": 0.0, "alt": -1.0},
                            {"az": 0.5, "alt": -2.0},
                        ],
                        [
                            {"az": 30.0, "alt": -0.5},
                            {"az": 31.0, "alt": -0.25},
                        ],
                    ],
                }
            }
        ),
        encoding="utf-8",
    )

    load_urban_debug_layers.cache_clear()
    got = load_urban_debug_layers(path)

    assert got == {
        "wikidata:Q57965": [
            [(-1.0, 0.0), (-2.0, 0.5)],
            [(-0.5, 30.0), (-0.25, 31.0)],
        ]
    }


def test_resolve_urban_debug_layer_for_tower_city_name(tmp_path: Path) -> None:
    path = tmp_path / "debug.json"
    path.write_text(
        json.dumps(
            {
                "wikidata:Q57965": {
                    "name": "Tokyo Skytree",
                    "outlines": [
                        [
                            {"az": 0.0, "alt": -1.0},
                            {"az": 0.5, "alt": -2.0},
                        ]
                    ],
                }
            }
        ),
        encoding="utf-8",
    )

    load_urban_debug_layers.cache_clear()
    got = resolve_urban_debug_layer_for_city_name("t/Tokyo Skytree", path)

    assert got == [[(-1.0, 0.0), (-2.0, 0.5)]]
