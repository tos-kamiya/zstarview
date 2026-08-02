from __future__ import annotations

import json
import logging
import os
import tempfile
import time
import urllib.error
from collections.abc import Iterable
from contextlib import contextmanager
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from filelock import FileLock

from ..paths import NOMINATIM_CACHE_DIR
from . import nominatim
from .place_cache import (
    build_place_cache_key,
    load_place_cache,
    save_place_cache,
)

logger = logging.getLogger(__name__)
RATE_LIMIT_INTERVAL_SECONDS = 1.0
RATE_LIMIT_FILE_NAME = "last_request.json"
REQUEST_LOCK_FILE_NAME = "request.lock"


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
    cache_fetched_at_utc: datetime | None = None


class PlaceSearchNetworkError(RuntimeError):
    """Nominatim was unreachable and no matching cache entry exists."""


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


def normalize_place_search_candidates(
    items: Iterable[dict[str, Any]],
) -> tuple[PlaceSearchCandidate, ...]:
    candidates = []
    for item in items:
        candidate = place_search_candidate_from_nominatim(item)
        if candidate is not None:
            candidates.append(candidate)
    candidates.sort(key=lambda candidate: candidate.importance, reverse=True)
    return tuple(candidates)


def _candidate_to_cache_dict(candidate: PlaceSearchCandidate) -> dict[str, object]:
    return {
        "category": candidate.category,
        "display_name": candidate.display_name,
        "importance": candidate.importance,
        "lat": candidate.latitude_deg,
        "lon": candidate.longitude_deg,
        "name": candidate.name,
        "type": candidate.type_name,
    }


def _read_last_request_time(path: Path) -> float | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        value = payload.get("last_request_started_at")
        return float(value) if isinstance(payload, dict) and value is not None else None
    except (FileNotFoundError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def _write_last_request_time(path: Path, value: float) -> None:
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="ascii",
            dir=path.parent,
            prefix=path.name + ".",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_name = temporary.name
            json.dump({"last_request_started_at": value}, temporary)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_name, path)
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


@contextmanager
def _nominatim_request_slot(
    cache_root: str | Path,
    *,
    now_func=None,
    sleep_func=None,
):
    now = now_func or time.time
    sleep = sleep_func or time.sleep
    root = Path(cache_root)
    root.mkdir(parents=True, exist_ok=True)
    marker_path = root / RATE_LIMIT_FILE_NAME
    with FileLock(str(root / REQUEST_LOCK_FILE_NAME)):
        current = float(now())
        previous = _read_last_request_time(marker_path)
        if previous is not None:
            remaining = RATE_LIMIT_INTERVAL_SECONDS - (current - previous)
            if remaining > 0.0:
                sleep(min(remaining, RATE_LIMIT_INTERVAL_SECONDS))
        started_at = float(now())
        _write_last_request_time(marker_path, started_at)
        yield


def search_place_candidates(
    query: str,
    *,
    limit: int = 5,
    countrycode: str | None = None,
    language: str = "en",
    user_agent: str = nominatim._DEFAULT_USER_AGENT,
    cache_root: str | Path = NOMINATIM_CACHE_DIR,
    _now_func=None,
    _sleep_func=None,
) -> tuple[PlaceSearchCandidate, ...]:
    key = build_place_cache_key(query, countrycode, language, limit=limit)
    url = nominatim._build_url(query, limit=limit, countrycode=countrycode)
    try:
        with _nominatim_request_slot(
            cache_root, now_func=_now_func, sleep_func=_sleep_func
        ):
            raw_results = nominatim._fetch(
                url, language=language, user_agent=user_agent
            )
    except urllib.error.HTTPError:
        raise
    except OSError as exc:
        cached = load_place_cache(key, cache_root=cache_root)
        if cached is None:
            raise PlaceSearchNetworkError(
                "Place search requires a network connection, and no cache is "
                f"available for '{query}'."
            ) from exc
        candidates = normalize_place_search_candidates(cached.results)
        if not candidates:
            raise PlaceSearchNetworkError(
                f"The cached place search for '{query}' is invalid."
            ) from exc
        logger.warning(
            "Nominatim unavailable; using cached place results for '%s' fetched at %s",
            query,
            cached.fetched_at_utc.isoformat(),
        )
        return tuple(
            replace(candidate, cache_fetched_at_utc=cached.fetched_at_utc)
            for candidate in candidates
        )
    candidates = normalize_place_search_candidates(raw_results)
    if candidates:
        try:
            save_place_cache(
                key,
                tuple(_candidate_to_cache_dict(candidate) for candidate in candidates),
                original_query=query,
                fetched_at_utc=datetime.now(timezone.utc),
                cache_root=cache_root,
            )
        except OSError:
            logger.warning(
                "Failed to save Nominatim cache for '%s'", query, exc_info=True
            )
    return candidates


__all__ = [
    "PlaceSearchCandidate",
    "PlaceSearchNetworkError",
    "normalize_place_search_candidates",
    "place_search_candidate_from_nominatim",
    "search_place_candidates",
]
