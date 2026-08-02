from __future__ import annotations

import hashlib
import json
import logging
import os
import tempfile
import unicodedata
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from filelock import FileLock

from ..paths import NOMINATIM_CACHE_DIR

logger = logging.getLogger(__name__)

CACHE_SCHEMA_VERSION = 1
KEY_SCHEMA_VERSION = 1
PROVIDER = "nominatim.openstreetmap.org"


@dataclass(frozen=True)
class PlaceCacheKey:
    query: str
    countrycode: str | None
    language: str
    limit: int = 5
    provider: str = PROVIDER
    key_schema_version: int = KEY_SCHEMA_VERSION

    def to_dict(self) -> dict[str, object]:
        return {
            "countrycode": self.countrycode,
            "key_schema_version": self.key_schema_version,
            "language": self.language,
            "limit": self.limit,
            "provider": self.provider,
            "query": self.query,
        }


@dataclass(frozen=True)
class CachedPlaceResults:
    fetched_at_utc: datetime
    original_query: str
    results: tuple[dict[str, object], ...]


def normalize_place_query(value: str) -> str:
    normalized = unicodedata.normalize("NFKC", str(value))
    return " ".join(normalized.split())


def build_place_cache_key(
    query: str, countrycode: str | None, language: str, *, limit: int = 5
) -> PlaceCacheKey:
    normalized_countrycode = (
        str(countrycode).strip().lower() if countrycode is not None else None
    )
    if not normalized_countrycode:
        normalized_countrycode = None
    return PlaceCacheKey(
        query=normalize_place_query(query),
        countrycode=normalized_countrycode,
        language=str(language).strip() or "en",
        limit=int(limit),
    )


def place_cache_digest(key: PlaceCacheKey) -> str:
    payload = json.dumps(
        key.to_dict(), ensure_ascii=False, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def place_cache_path(
    key: PlaceCacheKey, *, cache_root: str | Path = NOMINATIM_CACHE_DIR
) -> Path:
    return Path(cache_root) / f"{place_cache_digest(key)}.json"


def _parse_utc(value: object) -> datetime:
    if not isinstance(value, str):
        raise TypeError("fetched_at_utc must be a string")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError("fetched_at_utc must include a timezone")
    return parsed.astimezone(timezone.utc)


def _format_utc(value: datetime) -> str:
    if value.tzinfo is None:
        raise ValueError("fetched_at_utc must include a timezone")
    return value.astimezone(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


def _read_bucket(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError("cache bucket must be an object")
    if payload.get("schema_version") != CACHE_SCHEMA_VERSION:
        raise ValueError("unsupported cache schema")
    if payload.get("hash_algorithm") != "sha256":
        raise ValueError("unsupported cache hash algorithm")
    entries = payload.get("entries")
    if not isinstance(entries, list) or not all(
        isinstance(item, dict) for item in entries
    ):
        raise TypeError("cache entries must be objects")
    return entries


def _entry_key(entry: dict[str, Any]) -> PlaceCacheKey:
    raw = entry.get("key")
    if not isinstance(raw, dict):
        raise TypeError("cache entry key must be an object")
    if raw.get("key_schema_version") != KEY_SCHEMA_VERSION:
        raise ValueError("unsupported cache key schema")
    query = raw.get("query")
    language = raw.get("language")
    provider = raw.get("provider")
    limit = raw.get("limit")
    countrycode = raw.get("countrycode")
    if not isinstance(query, str) or not isinstance(language, str):
        raise TypeError("invalid cache key text")
    if provider != PROVIDER:
        raise ValueError("unexpected cache provider")
    if not isinstance(limit, int) or isinstance(limit, bool) or limit <= 0:
        raise TypeError("invalid cache result limit")
    if countrycode is not None and not isinstance(countrycode, str):
        raise TypeError("invalid cache country code")
    return PlaceCacheKey(
        query=query,
        countrycode=countrycode,
        language=language,
        limit=limit,
        provider=provider,
    )


def load_place_cache(
    key: PlaceCacheKey, *, cache_root: str | Path = NOMINATIM_CACHE_DIR
) -> CachedPlaceResults | None:
    path = place_cache_path(key, cache_root=cache_root)
    if not path.is_file():
        return None
    try:
        entries = _read_bucket(path)
        for entry in entries:
            entry_key = _entry_key(entry)
            if place_cache_digest(entry_key) != path.stem:
                raise ValueError("cache entry digest does not match its bucket")
            if entry_key != key:
                continue
            original_query = entry.get("original_query")
            results = entry.get("results")
            if not isinstance(original_query, str):
                raise TypeError("invalid original query")
            if not isinstance(results, list) or not results:
                raise TypeError("cached result list must not be empty")
            if not all(isinstance(item, dict) for item in results):
                raise TypeError("cached results must be objects")
            return CachedPlaceResults(
                fetched_at_utc=_parse_utc(entry.get("fetched_at_utc")),
                original_query=original_query,
                results=tuple(dict(item) for item in results),
            )
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        logger.warning("Ignoring invalid Nominatim cache bucket: %s", path, exc_info=True)
    return None


def save_place_cache(
    key: PlaceCacheKey,
    results: tuple[dict[str, object], ...],
    *,
    original_query: str,
    fetched_at_utc: datetime,
    cache_root: str | Path = NOMINATIM_CACHE_DIR,
) -> Path | None:
    if not results:
        return None
    path = place_cache_path(key, cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    lock = FileLock(str(path) + ".lock")
    with lock:
        entries: list[dict[str, Any]] = []
        if path.exists():
            try:
                entries = _read_bucket(path)
                for existing in entries:
                    existing_key = _entry_key(existing)
                    if place_cache_digest(existing_key) != path.stem:
                        raise ValueError("cache entry digest does not match its bucket")
            except (OSError, TypeError, ValueError, json.JSONDecodeError):
                logger.warning(
                    "Refusing to overwrite invalid Nominatim cache bucket: %s",
                    path,
                    exc_info=True,
                )
                return None
        new_entry: dict[str, Any] = {
            "fetched_at_utc": _format_utc(fetched_at_utc),
            "key": key.to_dict(),
            "original_query": str(original_query),
            "results": [dict(item) for item in results],
        }
        replaced = False
        for index, existing in enumerate(entries):
            if _entry_key(existing) == key:
                entries[index] = new_entry
                replaced = True
                break
        if not replaced:
            entries.append(new_entry)
        payload = {
            "entries": entries,
            "hash_algorithm": "sha256",
            "schema_version": CACHE_SCHEMA_VERSION,
        }
        temporary_name: str | None = None
        try:
            with tempfile.NamedTemporaryFile(
                "w",
                encoding="utf-8",
                dir=path.parent,
                prefix=path.name + ".",
                suffix=".tmp",
                delete=False,
            ) as temporary:
                temporary_name = temporary.name
                json.dump(payload, temporary, ensure_ascii=False, indent=2, sort_keys=True)
                temporary.write("\n")
                temporary.flush()
                os.fsync(temporary.fileno())
            os.replace(temporary_name, path)
        finally:
            if temporary_name is not None:
                Path(temporary_name).unlink(missing_ok=True)
    return path


__all__ = [
    "CachedPlaceResults",
    "PlaceCacheKey",
    "build_place_cache_key",
    "load_place_cache",
    "normalize_place_query",
    "place_cache_digest",
    "place_cache_path",
    "save_place_cache",
]
