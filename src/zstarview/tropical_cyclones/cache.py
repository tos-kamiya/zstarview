from __future__ import annotations

import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..paths import TROPICAL_CYCLONE_CACHE_DIR
from .models import TropicalCycloneSnapshot

TROPICAL_CYCLONE_CACHE_TTL_SECONDS = 3 * 60 * 60
TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS = 90 * 60
TROPICAL_CYCLONE_CACHE_FILENAME = "active_hurricanes.json"


@dataclass(frozen=True, slots=True)
class TropicalCycloneCacheEntry:
    snapshot: TropicalCycloneSnapshot
    cached_at_utc: datetime

    def to_dict(self) -> dict[str, Any]:
        return {
            "cached_at_utc": self.cached_at_utc.astimezone(timezone.utc)
            .isoformat()
            .replace("+00:00", "Z"),
            "snapshot": self.snapshot.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> TropicalCycloneCacheEntry | None:
        cached_at_raw = data.get("cached_at_utc")
        snapshot_raw = data.get("snapshot")
        if not isinstance(snapshot_raw, dict):
            return None
        snapshot = TropicalCycloneSnapshot.from_dict(snapshot_raw)
        if snapshot is None:
            return None
        cached_at = _parse_datetime(cached_at_raw)
        if cached_at is None:
            return None
        return cls(snapshot=snapshot, cached_at_utc=cached_at)


def _parse_datetime(value: object) -> datetime | None:
    if not isinstance(value, str) or not value.strip():
        return None
    text = value.strip().replace("Z", "+00:00")
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def tropical_cyclone_cache_root(cache_root: str | os.PathLike[str] | None = None) -> Path:
    return Path(cache_root or TROPICAL_CYCLONE_CACHE_DIR).expanduser()


def tropical_cyclone_cache_path(
    cache_root: str | os.PathLike[str] | None = None,
) -> Path:
    return tropical_cyclone_cache_root(cache_root) / TROPICAL_CYCLONE_CACHE_FILENAME


def load_tropical_cyclone_cache(
    cache_root: str | os.PathLike[str] | None = None,
) -> TropicalCycloneCacheEntry | None:
    path = tropical_cyclone_cache_path(cache_root)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except Exception:
        return None
    if not isinstance(raw, dict):
        return None
    return TropicalCycloneCacheEntry.from_dict(raw)


def save_tropical_cyclone_cache(
    entry: TropicalCycloneCacheEntry,
    *,
    cache_root: str | os.PathLike[str] | None = None,
) -> Path:
    root = tropical_cyclone_cache_root(cache_root)
    root.mkdir(parents=True, exist_ok=True)
    path = tropical_cyclone_cache_path(root)
    path.write_text(
        json.dumps(entry.to_dict(), ensure_ascii=True, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return path


def is_tropical_cyclone_cache_stale(
    entry: TropicalCycloneCacheEntry,
    *,
    now_utc: datetime | None = None,
) -> bool:
    if now_utc is None:
        now_utc = datetime.now(timezone.utc)
    age_seconds = (now_utc - entry.cached_at_utc).total_seconds()
    return age_seconds >= float(TROPICAL_CYCLONE_CACHE_TTL_SECONDS)
