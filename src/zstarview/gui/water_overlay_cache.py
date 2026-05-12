from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from ..paths import CACHE_PATH
from ..water_overlay import WaterPolygonFootprint

WATER_OVERLAY_CACHE_ROOT_DIR = Path(CACHE_PATH) / "water_overlay"
WATER_OVERLAY_CACHE_FORMAT_VERSION = 1
WATER_OVERLAY_CACHE_RETENTION_SECONDS = 90 * 24 * 60 * 60


@dataclass(frozen=True, slots=True)
class WaterOverlayCacheSnapshot:
    footprints: tuple[WaterPolygonFootprint, ...]
    water_polygon_count: int
    fetched_at_utc: datetime | None = None


def water_overlay_cache_scope_key(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float,
) -> str:
    return "earth_{lat:+08.4f}_{lon:+09.4f}_r{radius:.2f}".format(
        lat=float(observer_lat_deg),
        lon=float(observer_lon_deg),
        radius=float(radius_km),
    )


def water_overlay_cache_path(
    scope_key: str,
    *,
    cache_root: str | Path = WATER_OVERLAY_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / f"{scope_key}.json"


def load_water_overlay_cache(
    scope_key: str,
    *,
    cache_root: str | Path = WATER_OVERLAY_CACHE_ROOT_DIR,
) -> WaterOverlayCacheSnapshot | None:
    path = water_overlay_cache_path(scope_key, cache_root=cache_root)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(payload, dict):
        return None
    try:
        footprints_payload = payload.get("footprints", [])
        footprints: list[WaterPolygonFootprint] = []
        if isinstance(footprints_payload, list):
            for item in footprints_payload:
                footprint = _deserialize_footprint(item)
                if footprint is not None:
                    footprints.append(footprint)
        return WaterOverlayCacheSnapshot(
            footprints=tuple(footprints),
            water_polygon_count=int(payload.get("water_polygon_count", 0)),
            fetched_at_utc=_parse_optional_utc(payload.get("fetched_at_utc")),
        )
    except Exception:
        return None


def save_water_overlay_cache(
    scope_key: str,
    snapshot: WaterOverlayCacheSnapshot,
    *,
    cache_root: str | Path = WATER_OVERLAY_CACHE_ROOT_DIR,
) -> Path:
    path = water_overlay_cache_path(scope_key, cache_root=cache_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cache_format_version": WATER_OVERLAY_CACHE_FORMAT_VERSION,
        "scope_key": str(scope_key),
        "fetched_at_utc": _serialize_optional_utc(snapshot.fetched_at_utc),
        "water_polygon_count": int(snapshot.water_polygon_count),
        "footprints": [_serialize_footprint(footprint) for footprint in snapshot.footprints],
    }
    path.write_text(json.dumps(payload, separators=(",", ":"), sort_keys=True), encoding="utf-8")
    return path


def _serialize_footprint(footprint: WaterPolygonFootprint) -> dict[str, object]:
    return {
        "water_id": str(footprint.water_id),
        "kind": str(footprint.kind),
        "outer_rings_lonlat": [
            [[float(lon), float(lat)] for lon, lat in ring]
            for ring in footprint.outer_rings_lonlat
        ],
        "inner_rings_lonlat": [
            [[float(lon), float(lat)] for lon, lat in ring]
            for ring in footprint.inner_rings_lonlat
        ],
        "source": str(footprint.source),
        "tags": dict(footprint.tags),
    }


def _deserialize_footprint(payload: object) -> WaterPolygonFootprint | None:
    if not isinstance(payload, dict):
        return None
    try:
        outer_rings = tuple(
            tuple((float(lon), float(lat)) for lon, lat in ring)
            for ring in payload.get("outer_rings_lonlat", [])
            if isinstance(ring, list)
        )
        inner_rings = tuple(
            tuple((float(lon), float(lat)) for lon, lat in ring)
            for ring in payload.get("inner_rings_lonlat", [])
            if isinstance(ring, list)
        )
        tags_obj = payload.get("tags", {})
        tags: dict[str, str] = {}
        if isinstance(tags_obj, dict):
            for key, value in tags_obj.items():
                if isinstance(key, str) and isinstance(value, str):
                    tags[key] = value
        return WaterPolygonFootprint(
            water_id=str(payload.get("water_id", "water")),
            kind=str(payload.get("kind", "water_polygon")),
            outer_rings_lonlat=outer_rings,
            inner_rings_lonlat=inner_rings,
            source=str(payload.get("source", "")),
            tags=tags,
        )
    except Exception:
        return None


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _serialize_optional_utc(value: datetime | None) -> str | None:
    if value is None:
        return None
    return _normalize_utc(value).isoformat()


def _parse_optional_utc(value: object) -> datetime | None:
    if not isinstance(value, str) or not value.strip():
        return None
    parsed = datetime.fromisoformat(value)
    return _normalize_utc(parsed)


def _string_or_none(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def water_overlay_cache_age_seconds(
    snapshot: WaterOverlayCacheSnapshot,
    *,
    now_utc: datetime,
) -> float | None:
    if snapshot.fetched_at_utc is None:
        return None
    age_seconds = (now_utc - snapshot.fetched_at_utc).total_seconds()
    return max(0.0, float(age_seconds))


def water_overlay_cache_is_recent(
    snapshot: WaterOverlayCacheSnapshot,
    *,
    now_utc: datetime,
    max_age_seconds: int = WATER_OVERLAY_CACHE_RETENTION_SECONDS,
) -> bool:
    age_seconds = water_overlay_cache_age_seconds(snapshot, now_utc=now_utc)
    if age_seconds is None:
        return False
    return age_seconds <= float(max_age_seconds)
