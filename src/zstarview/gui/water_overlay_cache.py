from __future__ import annotations

import json
import logging
import re
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from ..paths import CACHE_PATH
from ..water_overlay import WaterPolygonFootprint

logger = logging.getLogger(__name__)

WATER_OVERLAY_CACHE_ROOT_DIR = Path(CACHE_PATH) / "water_overlay"
WATER_OVERLAY_CACHE_FORMAT_VERSION = 2
WATER_OVERLAY_CACHE_RETENTION_SECONDS = 90 * 24 * 60 * 60
_LOCATION_CACHE_SCOPE_RE = re.compile(
    r"^earth_(?P<lat>[+-]\d+\.\d{3,4})_(?P<lon>[+-]\d+\.\d{3,4})_r(?P<radius>\d+\.\d{2})(?:_a\d+\.\d{2})?$"
)


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
    azimuth_step_deg: float | None = None,
) -> str:
    coordinate_decimals = 3 if float(radius_km) >= 128.0 else 4
    coordinate_width = coordinate_decimals + 4
    key = f"earth_{{lat:+0{coordinate_width}.{coordinate_decimals}f}}_{{lon:+0{coordinate_width + 1}.{coordinate_decimals}f}}_r{{radius:.2f}}".format(
        lat=float(observer_lat_deg),
        lon=float(observer_lon_deg),
        radius=float(radius_km),
    )
    if azimuth_step_deg is None:
        return key
    step = float(azimuth_step_deg)
    if abs(step - 2.0) < 1.0e-6:
        return key
    return f"{key}_a{step:.2f}"


def water_overlay_cache_path(
    scope_key: str,
    *,
    cache_root: str | Path = WATER_OVERLAY_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / f"{scope_key}_simplified.json"


def water_overlay_cache_legacy_path(
    scope_key: str,
    *,
    cache_root: str | Path = WATER_OVERLAY_CACHE_ROOT_DIR,
) -> Path:
    return Path(cache_root) / f"{scope_key}.json"


def load_water_overlay_cache(
    scope_key: str,
    *,
    cache_root: str | Path = WATER_OVERLAY_CACHE_ROOT_DIR,
    observer_lat_deg: float | None = None,
    observer_lon_deg: float | None = None,
) -> WaterOverlayCacheSnapshot | None:
    simplified_path = water_overlay_cache_path(scope_key, cache_root=cache_root)
    legacy_path = water_overlay_cache_legacy_path(scope_key, cache_root=cache_root)
    simplified_snapshot = _load_snapshot_from_path(simplified_path)
    legacy_snapshot = _load_snapshot_from_path(legacy_path)

    simplified_mtime_ns = _file_mtime_ns(simplified_path)
    legacy_mtime_ns = _file_mtime_ns(legacy_path)

    if simplified_snapshot is not None:
        if (
            legacy_snapshot is not None
            and legacy_mtime_ns is not None
            and simplified_mtime_ns is not None
            and legacy_mtime_ns > simplified_mtime_ns
        ):
            return _maybe_rebuild_simplified_cache(
                scope_key,
                cache_root=cache_root,
                snapshot=legacy_snapshot,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
            )
        return simplified_snapshot

    if legacy_snapshot is not None:
        if legacy_mtime_ns is not None and simplified_mtime_ns is not None:
            if legacy_mtime_ns <= simplified_mtime_ns:
                return legacy_snapshot
        return _maybe_rebuild_simplified_cache(
            scope_key,
            cache_root=cache_root,
            snapshot=legacy_snapshot,
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
        )

    return None


def load_water_overlay_cache_for_location(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    cache_root: str | Path = WATER_OVERLAY_CACHE_ROOT_DIR,
) -> WaterOverlayCacheSnapshot | None:
    """Return the newest cache for a rounded location, regardless of radius.

    The scan radius depends on observer height, so changing height can produce a
    different scope key for the same location.  This lookup is intentionally a
    fallback only; callers should try the exact scope first.
    """
    root = Path(cache_root)
    candidates: list[tuple[datetime, float, WaterOverlayCacheSnapshot]] = []
    for prefix_template in (
        "earth_{lat:+08.4f}_{lon:+09.4f}_r",
        "earth_{lat:+07.3f}_{lon:+08.3f}_r",
    ):
        prefix = prefix_template.format(
            lat=float(observer_lat_deg),
            lon=float(observer_lon_deg),
        )
        for path in root.glob(f"{prefix}*.json"):
            name = path.name
            scope_key = name.removesuffix("_simplified.json").removesuffix(".json")
            match = _LOCATION_CACHE_SCOPE_RE.fullmatch(scope_key)
            if match is None:
                continue
            snapshot = load_water_overlay_cache(
                scope_key,
                cache_root=root,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
            )
            if snapshot is None or snapshot.fetched_at_utc is None:
                continue
            candidates.append(
                (
                    _normalize_utc(snapshot.fetched_at_utc),
                    float(match.group("radius")),
                    snapshot,
                )
            )
    if not candidates:
        return None
    # Prefer the broadest valid snapshot, then the newest one.
    return max(candidates, key=lambda item: (item[1], item[0]))[2]


def save_water_overlay_cache(
    scope_key: str,
    snapshot: WaterOverlayCacheSnapshot,
    *,
    cache_root: str | Path = WATER_OVERLAY_CACHE_ROOT_DIR,
) -> Path:
    path = water_overlay_cache_path(scope_key, cache_root=cache_root)
    return save_water_overlay_cache_to_path(path, scope_key=scope_key, snapshot=snapshot)


def save_water_overlay_cache_to_path(
    path: Path,
    *,
    scope_key: str,
    snapshot: WaterOverlayCacheSnapshot,
) -> Path:
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


def _load_snapshot_from_path(path: Path) -> WaterOverlayCacheSnapshot | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(payload, dict):
        return None
    if int(payload.get("cache_format_version", 0)) != WATER_OVERLAY_CACHE_FORMAT_VERSION:
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


def _file_mtime_ns(path: Path) -> int | None:
    try:
        return path.stat().st_mtime_ns
    except FileNotFoundError:
        return None


def _maybe_rebuild_simplified_cache(
    scope_key: str,
    *,
    cache_root: str | Path,
    snapshot: WaterOverlayCacheSnapshot,
    observer_lat_deg: float | None,
    observer_lon_deg: float | None,
) -> WaterOverlayCacheSnapshot:
    if observer_lat_deg is None or observer_lon_deg is None:
        return snapshot
    from ..water_overlay import simplify_water_footprints_for_observer

    started_at = time.monotonic()
    logger.info(
        "Water overlay cache simplification started: scope=%s footprints=%d",
        scope_key,
        len(snapshot.footprints),
    )
    simplified_footprints = simplify_water_footprints_for_observer(
        snapshot.footprints,
        observer_lat_deg=float(observer_lat_deg),
        observer_lon_deg=float(observer_lon_deg),
    )
    rebuilt = WaterOverlayCacheSnapshot(
        footprints=simplified_footprints,
        water_polygon_count=len(simplified_footprints),
        fetched_at_utc=snapshot.fetched_at_utc,
    )
    save_water_overlay_cache(
        scope_key,
        rebuilt,
        cache_root=cache_root,
    )
    elapsed_s = max(0.0, time.monotonic() - started_at)
    logger.info(
        "Water overlay cache simplification finished: scope=%s input=%d output=%d removed=%d elapsed=%.3fs",
        scope_key,
        len(snapshot.footprints),
        len(rebuilt.footprints),
        max(0, len(snapshot.footprints) - len(rebuilt.footprints)),
        elapsed_s,
    )
    return rebuilt


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
