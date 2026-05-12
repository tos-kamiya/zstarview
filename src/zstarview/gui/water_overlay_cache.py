from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from ..paths import CACHE_PATH
from ..water_overlay import WaterOverlayPoint

WATER_OVERLAY_CACHE_ROOT_DIR = Path(CACHE_PATH) / "water_overlay"
WATER_OVERLAY_CACHE_FORMAT_VERSION = 1
WATER_OVERLAY_CACHE_RETENTION_SECONDS = 90 * 24 * 60 * 60


@dataclass(frozen=True, slots=True)
class WaterOverlayCacheSnapshot:
    points: tuple[WaterOverlayPoint, ...]
    water_polygon_count: int
    water_point_count: int
    fetched_at_utc: datetime | None = None


def water_overlay_cache_scope_key(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    radius_km: float,
    sample_step_m: float,
    azimuth_step_deg: float,
) -> str:
    return "earth_{lat:+08.4f}_{lon:+09.4f}_{height:+07.1f}_r{radius:.2f}_s{sample:.3f}_a{azimuth:.2f}".format(
        lat=float(observer_lat_deg),
        lon=float(observer_lon_deg),
        height=float(observer_height_m),
        radius=float(radius_km),
        sample=float(sample_step_m),
        azimuth=float(azimuth_step_deg),
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
        points_payload = payload.get("points", [])
        points: list[WaterOverlayPoint] = []
        if isinstance(points_payload, list):
            for item in points_payload:
                if not isinstance(item, dict):
                    continue
                points.append(
                    WaterOverlayPoint(
                        water_id=str(item.get("water_id", "water")),
                        alt_deg=float(item.get("alt_deg", 0.0)),
                        az_deg=float(item.get("az_deg", 0.0)),
                        distance_km=float(item.get("distance_km", 0.0)),
                        alpha_scale=float(item.get("alpha_scale", 1.0)),
                    )
                )
        return WaterOverlayCacheSnapshot(
            points=tuple(points),
            water_polygon_count=int(payload.get("water_polygon_count", 0)),
            water_point_count=int(payload.get("water_point_count", len(points))),
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
        "water_point_count": int(snapshot.water_point_count),
        "points": [_serialize_point(point) for point in snapshot.points],
    }
    path.write_text(json.dumps(payload, separators=(",", ":"), sort_keys=True), encoding="utf-8")
    return path


def _serialize_point(point: WaterOverlayPoint) -> dict[str, float | str]:
    return {
        "water_id": str(point.water_id),
        "alt_deg": float(point.alt_deg),
        "az_deg": float(point.az_deg),
        "distance_km": float(point.distance_km),
        "alpha_scale": float(point.alpha_scale),
    }


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
