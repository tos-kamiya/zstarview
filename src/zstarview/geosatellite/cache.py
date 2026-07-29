from __future__ import annotations

import datetime as dt
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

from ..paths import GEOSATELLITE_CACHE_ROOT_DIR
from .types import GeoSatelliteDownloadResult

RAW_CACHE_DIRNAME = "raw"
PROJECTION_CACHE_DIRNAME = "projection"
AVAILABLE_CACHE_DIRNAME = "available"
AVAILABLE_CACHE_MAX_AGE_SECONDS = 300.0


def geo_satellite_cache_root() -> Path:
    return Path(GEOSATELLITE_CACHE_ROOT_DIR)


def geo_satellite_raw_cache_dir() -> Path:
    return geo_satellite_cache_root() / RAW_CACHE_DIRNAME


def geo_satellite_projection_cache_dir() -> Path:
    return geo_satellite_cache_root() / PROJECTION_CACHE_DIRNAME


def geo_satellite_available_cache_dir() -> Path:
    return geo_satellite_cache_root() / AVAILABLE_CACHE_DIRNAME


def raw_png_path() -> Path:
    return geo_satellite_raw_cache_dir() / "latest.png"


def raw_metadata_path() -> Path:
    return geo_satellite_raw_cache_dir() / "latest.json"


def _legacy_raw_paths(*, kind: str) -> list[Path]:
    raw_dir = geo_satellite_raw_cache_dir()
    paths: list[Path] = []
    if not raw_dir.exists():
        return paths
    paths.extend(raw_dir.glob(f"*_{kind}.png"))
    paths.extend(raw_dir.glob(f"*_{kind}.json"))
    return paths


def purge_legacy_raw_cache(*, kind: str) -> None:
    current_paths = {raw_png_path(), raw_metadata_path()}
    for path in _legacy_raw_paths(kind=kind):
        if path in current_paths:
            continue
        try:
            path.unlink()
        except FileNotFoundError:
            continue


def purge_intermediate_cache() -> None:
    shutil.rmtree(geo_satellite_cache_root() / "intermediate", ignore_errors=True)


def projection_latest_path() -> Path:
    return geo_satellite_projection_cache_dir() / "latest.npz"


def projection_metadata_path() -> Path:
    return geo_satellite_projection_cache_dir() / "latest.json"


def write_projection_cache(
    *,
    cache_key: str,
    source_shape: tuple[int, int],
    x_src: np.ndarray,
    y_src: np.ndarray,
    valid_mask: np.ndarray,
) -> tuple[Path, Path]:
    npz_path = projection_latest_path()
    metadata_path = projection_metadata_path()
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        npz_path,
        x_src=np.asarray(x_src, dtype=np.float32),
        y_src=np.asarray(y_src, dtype=np.float32),
        valid_mask=np.asarray(valid_mask, dtype=np.uint8),
    )
    metadata_path.write_text(
        json.dumps(
            {
                "cache_key": cache_key,
                "source_shape": [int(source_shape[0]), int(source_shape[1])],
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return npz_path, metadata_path


def read_projection_cache(
    *,
    cache_key: str,
    source_shape: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    npz_path = projection_latest_path()
    metadata_path = projection_metadata_path()
    if not npz_path.exists() or not metadata_path.exists():
        return None
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if str(metadata.get("cache_key")) != cache_key:
        return None
    stored_shape = metadata.get("source_shape")
    if not isinstance(stored_shape, list) or len(stored_shape) != 2:
        return None
    if int(stored_shape[0]) != int(source_shape[0]) or int(stored_shape[1]) != int(source_shape[1]):
        return None
    with np.load(npz_path, allow_pickle=False) as data:
        x_src = np.asarray(data["x_src"], dtype=np.float32)
        y_src = np.asarray(data["y_src"], dtype=np.float32)
        valid_mask = np.asarray(data["valid_mask"], dtype=np.uint8).astype(bool, copy=False)
    return x_src, y_src, valid_mask


def write_raw_cache(result: GeoSatelliteDownloadResult) -> tuple[Path, Path]:
    image_time_utc = result.captured_at_utc or result.fetched_at_utc
    png_path = raw_png_path()
    metadata_path = raw_metadata_path()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    purge_legacy_raw_cache(kind=result.kind)
    png_path.write_bytes(result.png_bytes)
    metadata_path.write_text(
        json.dumps(
            {
                "captured_at_utc": image_time_utc.isoformat(),
                "fetched_at_utc": result.fetched_at_utc.isoformat(),
                "kind": result.kind,
                "source_url": result.source_url,
                "content_type": result.content_type,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return png_path, metadata_path


def read_raw_cache(*, image_time_utc: dt.datetime, kind: str) -> GeoSatelliteDownloadResult | None:
    png_path = raw_png_path()
    metadata_path = raw_metadata_path()
    if not png_path.exists() or not metadata_path.exists():
        return None
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if str(metadata.get("kind")) != kind:
        return None
    png_bytes = png_path.read_bytes()
    captured_at_text = metadata.get("captured_at_utc") or metadata.get("fetched_at_utc")
    normalized_expected = image_time_utc.astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
    normalized_captured = (
        dt.datetime.fromisoformat(str(captured_at_text)).astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
        if captured_at_text
        else None
    )
    if normalized_captured is None or normalized_captured != normalized_expected:
        return None
    return GeoSatelliteDownloadResult(
        fetched_at_utc=dt.datetime.fromisoformat(str(metadata["fetched_at_utc"])),
        captured_at_utc=normalized_captured,
        kind=str(metadata["kind"]),
        source_url=str(metadata["source_url"]),
        png_bytes=png_bytes,
        content_type=metadata.get("content_type"),
        cache_path=png_path,
        metadata_path=metadata_path,
    )


def compute_digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def latest_available_path(*, area: str, kind: str, size: str = "normal") -> Path:
    filename = f"{area}_{kind}_{size}.json"
    return geo_satellite_available_cache_dir() / filename


def write_latest_available_cache(
    *,
    area: str,
    kind: str,
    size: str,
    available_time_utc: dt.datetime,
    fetched_at_utc: dt.datetime,
    source_url: str,
    entry_count: int,
) -> Path:
    path = latest_available_path(area=area, kind=kind, size=size)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "area": area,
                "kind": kind,
                "size": size,
                "available_time_utc": available_time_utc.isoformat(),
                "fetched_at_utc": fetched_at_utc.isoformat(),
                "source_url": source_url,
                "entry_count": int(entry_count),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return path


def read_latest_available_cache(
    *,
    area: str,
    kind: str,
    size: str = "normal",
) -> dict[str, object] | None:
    path = latest_available_path(area=area, kind=kind, size=size)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if str(payload.get("area")) != area or str(payload.get("kind")) != kind or str(payload.get("size")) != size:
        return None
    return payload


def latest_available_cache_is_recent(
    payload: dict[str, object],
    *,
    now_utc: dt.datetime,
    max_age_seconds: float = AVAILABLE_CACHE_MAX_AGE_SECONDS,
) -> bool:
    fetched_at_text = payload.get("fetched_at_utc")
    if not fetched_at_text:
        return False
    try:
        fetched_at = dt.datetime.fromisoformat(str(fetched_at_text))
    except Exception:
        return False
    if fetched_at.tzinfo is None:
        fetched_at = fetched_at.replace(tzinfo=dt.timezone.utc)
    age_seconds = (now_utc.astimezone(dt.timezone.utc) - fetched_at.astimezone(dt.timezone.utc)).total_seconds()
    return age_seconds >= 0.0 and age_seconds <= float(max_age_seconds)
