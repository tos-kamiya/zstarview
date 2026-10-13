from __future__ import annotations

import json
import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path
from uuid import uuid4

from ..config import cache_retention_seconds
from ..paths import CACHE_PATH

logger = logging.getLogger(__name__)


def save_place_launch_metadata(
    *,
    query: str,
    countrycode: str | None,
    language: str,
    location: object,
    observer_height_m: float,
    started_at_utc: datetime | None = None,
    cache_root: str | Path = CACHE_PATH,
) -> Path | None:
    """Persist the original --place query and resolved observer location."""
    now = _as_utc(started_at_utc or datetime.now(timezone.utc))
    try:
        latitude = float(getattr(location, "lat"))
        longitude = float(getattr(location, "lon"))
        if not (-90.0 <= latitude <= 90.0 and -180.0 <= longitude <= 180.0):
            raise ValueError("resolved coordinates are out of range")
        payload = {
            "schema_version": 1,
            "launch_id": str(uuid4()),
            "started_at_utc": now.isoformat().replace("+00:00", "Z"),
            "location": {
                "input": {
                    "option": "--place",
                    "query": str(query),
                    "countrycode": countrycode,
                    "language": str(language),
                },
                "resolved": {
                    "display_name": str(getattr(location, "display_name", "")),
                    "latitude_deg": latitude,
                    "longitude_deg": longitude,
                    "timezone": str(getattr(location, "tz", "UTC")),
                    "source": str(getattr(location, "kind", "unknown")),
                    "ground_elevation_m": float(
                        getattr(location, "ground_elevation_m", 0.0)
                    ),
                    "observer_height_m": float(observer_height_m),
                },
            },
        }
        launch_id = str(payload["launch_id"])
        directory = Path(cache_root) / "replay" / "v1" / "launches" / now.strftime(
            "%Y/%m/%d/%H/%M"
        )
        directory.mkdir(parents=True, exist_ok=True)
        destination = directory / f"{launch_id}.json"
        temporary = destination.with_name(f".{destination.name}.{uuid4().hex}.tmp")
        temporary.write_text(
            json.dumps(payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True),
            encoding="utf-8",
        )
        temporary.replace(destination)
        try:
            cleanup_place_launch_metadata(cache_root=cache_root, now_utc=now)
        except Exception:
            logger.warning("Failed to clean expired launch metadata", exc_info=True)
        return destination
    except Exception:
        logger.warning("Failed to save --place launch metadata", exc_info=True)
        return None


def cleanup_place_launch_metadata(
    *,
    cache_root: str | Path = CACHE_PATH,
    now_utc: datetime | None = None,
) -> None:
    root = Path(cache_root) / "replay" / "v1" / "launches"
    if not root.is_dir():
        return
    cutoff = _as_utc(now_utc or datetime.now(timezone.utc)) - timedelta(
        seconds=cache_retention_seconds()
    )
    for path in root.rglob("*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            started_at = _as_utc(
                datetime.fromisoformat(
                    str(payload["started_at_utc"]).replace("Z", "+00:00")
                )
            )
            expired = started_at < cutoff
        except Exception:
            expired = True
        if expired:
            try:
                path.unlink(missing_ok=True)
            except OSError:
                logger.warning("Failed to remove expired launch metadata")
    for directory in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if directory.is_dir():
            try:
                directory.rmdir()
            except OSError:
                pass


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)
