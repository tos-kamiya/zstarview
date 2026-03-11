from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from .paths import URBAN_SKYLINE_PROFILES_FILE
from .tower_viewpoints import resolve_tower_viewpoint
from .viewpoints import split_prefixed_viewpoint


@lru_cache(maxsize=1)
def load_urban_skyline_profiles(
    path: str | Path = URBAN_SKYLINE_PROFILES_FILE,
) -> dict[str, list[tuple[float, list[tuple[float, float]]]]]:
    source = Path(path)
    if not source.exists():
        return {}
    payload = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return {}

    profiles: dict[str, list[tuple[float, list[tuple[float, float]]]]] = {}
    for profile_id, item in payload.items():
        if not isinstance(profile_id, str) or not isinstance(item, dict):
            continue
        layers: list[tuple[float, list[tuple[float, float]]]] = []
        raw_profiles = item.get("profiles")
        if isinstance(raw_profiles, list):
            for raw_layer in raw_profiles:
                if not isinstance(raw_layer, dict):
                    continue
                radius_km = raw_layer.get("radius_km")
                raw_profile = raw_layer.get("profile")
                if not isinstance(radius_km, (int, float)) or not isinstance(raw_profile, list):
                    continue
                points = _parse_profile_points(raw_profile)
                if points:
                    layers.append((float(radius_km), points))
        else:
            raw_profile = item.get("profile")
            if isinstance(raw_profile, list):
                points = _parse_profile_points(raw_profile)
                if points:
                    layers.append((0.0, points))
        if layers:
            layers.sort(key=lambda pair: pair[0])
            profiles[profile_id] = layers
    return profiles


def resolve_urban_skyline_profile_for_city_name(
    city_name: str,
    path: str | Path = URBAN_SKYLINE_PROFILES_FILE,
) -> list[tuple[float, list[tuple[float, float]]]] | None:
    split = split_prefixed_viewpoint(city_name)
    if split is None:
        return None
    kind, name = split
    if kind != "tower":
        return None
    tower = resolve_tower_viewpoint(name)
    if tower is None:
        return None
    return load_urban_skyline_profiles(path).get(tower.id)


def _parse_profile_points(raw_profile: list[object]) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    for row in raw_profile:
        if not isinstance(row, dict):
            continue
        az = row.get("az")
        alt = row.get("alt")
        if not isinstance(az, (int, float)) or not isinstance(alt, (int, float)):
            continue
        points.append((float(alt), float(az)))
    return points
