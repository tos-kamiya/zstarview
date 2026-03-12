from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from .paths import URBAN_DEBUG_LAYER_FILE
from .tower_viewpoints import resolve_tower_viewpoint
from .viewpoints import split_prefixed_viewpoint


@lru_cache(maxsize=1)
def load_urban_debug_layers(
    path: str | Path = URBAN_DEBUG_LAYER_FILE,
) -> dict[str, list[list[tuple[float, float]]]]:
    source = Path(path)
    if not source.exists():
        return {}
    payload = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return {}

    layers: dict[str, list[list[tuple[float, float]]]] = {}
    for layer_id, item in payload.items():
        if not isinstance(layer_id, str) or not isinstance(item, dict):
            continue
        raw_outlines = item.get("outlines")
        if not isinstance(raw_outlines, list):
            continue
        outlines = _parse_outlines(raw_outlines)
        if outlines:
            layers[layer_id] = outlines
    return layers


def resolve_urban_debug_layer_for_city_name(
    city_name: str,
    path: str | Path = URBAN_DEBUG_LAYER_FILE,
) -> list[list[tuple[float, float]]] | None:
    split = split_prefixed_viewpoint(city_name)
    if split is None:
        return None
    kind, name = split
    if kind != "tower":
        return None
    tower = resolve_tower_viewpoint(name)
    if tower is None:
        return None
    return load_urban_debug_layers(path).get(tower.id)


def _parse_outlines(raw_outlines: list[object]) -> list[list[tuple[float, float]]]:
    outlines: list[list[tuple[float, float]]] = []
    for raw_outline in raw_outlines:
        if not isinstance(raw_outline, list):
            continue
        points: list[tuple[float, float]] = []
        for row in raw_outline:
            if not isinstance(row, dict):
                continue
            az = row.get("az")
            alt = row.get("alt")
            if not isinstance(az, (int, float)) or not isinstance(alt, (int, float)):
                continue
            points.append((float(alt), float(az)))
        if len(points) >= 2:
            outlines.append(points)
    return outlines
