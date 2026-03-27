from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from . import nominatim


@dataclass(frozen=True)
class PlaceSearchCandidate:
    name: str
    display_name: str
    latitude_deg: float
    longitude_deg: float
    category: str
    type_name: str
    importance: float
    kind: str = "place"


def _candidate_name(item: dict[str, Any]) -> str | None:
    name = item.get("namedetails")
    if isinstance(name, dict):
        for key in ("name", "int_name"):
            value = name.get(key)
            if isinstance(value, str) and value.strip():
                return value.strip()

    for key in ("name", "display_name"):
        value = item.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip().split(",")[0]
    return None


def place_search_candidate_from_nominatim(item: dict[str, Any]) -> PlaceSearchCandidate | None:
    if not isinstance(item, dict):
        return None
    try:
        latitude_deg = float(item["lat"])
        longitude_deg = float(item["lon"])
    except (KeyError, TypeError, ValueError):
        return None

    display_name = item.get("display_name")
    if not isinstance(display_name, str) or not display_name.strip():
        return None

    name = _candidate_name(item)
    if name is None:
        return None

    category = item.get("category")
    if not isinstance(category, str) or not category.strip():
        category = item.get("class")
    if not isinstance(category, str) or not category.strip():
        category = "unknown"

    type_name = item.get("type")
    if not isinstance(type_name, str) or not type_name.strip():
        type_name = "unknown"

    try:
        importance = float(item.get("importance") or 0.0)
    except (TypeError, ValueError):
        importance = 0.0

    return PlaceSearchCandidate(
        name=name,
        display_name=display_name.strip(),
        latitude_deg=latitude_deg,
        longitude_deg=longitude_deg,
        category=category.strip(),
        type_name=type_name.strip(),
        importance=importance,
    )


def normalize_place_search_candidates(items: Iterable[dict[str, Any]]) -> tuple[PlaceSearchCandidate, ...]:
    candidates = []
    for item in items:
        candidate = place_search_candidate_from_nominatim(item)
        if candidate is not None:
            candidates.append(candidate)
    candidates.sort(key=lambda candidate: candidate.importance, reverse=True)
    return tuple(candidates)


def search_place_candidates(
    query: str,
    *,
    limit: int = 5,
    countrycode: str | None = None,
    language: str = "en",
    user_agent: str = nominatim._DEFAULT_USER_AGENT,
) -> tuple[PlaceSearchCandidate, ...]:
    url = nominatim._build_url(query, limit=limit, countrycode=countrycode)
    raw_results = nominatim._fetch(url, language=language, user_agent=user_agent)
    return normalize_place_search_candidates(raw_results)


__all__ = [
    "PlaceSearchCandidate",
    "normalize_place_search_candidates",
    "place_search_candidate_from_nominatim",
    "search_place_candidates",
]
