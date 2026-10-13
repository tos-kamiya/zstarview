from __future__ import annotations

import datetime as dt
import hashlib
import json
import shutil
from pathlib import Path
from uuid import uuid4

import numpy as np

from ..config import cache_retention_seconds
from ..paths import GEOSATELLITE_CACHE_ROOT_DIR
from .types import GeoSatelliteDownloadResult

RAW_CACHE_DIRNAME = "raw"
PROJECTION_CACHE_DIRNAME = "projection"
AVAILABLE_CACHE_DIRNAME = "available"
AVAILABLE_CACHE_MAX_AGE_SECONDS = 300.0
RAW_HISTORY_DIRNAME = "history"
RAW_HISTORY_SCHEMA_VERSION = 1
RAW_HISTORY_RETENTION_SECONDS = 24 * 60 * 60  # Default setting.
RAW_HISTORY_MAX_TARGET_AGE_SECONDS = 60 * 60


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


def raw_history_root() -> Path:
    return geo_satellite_raw_cache_dir() / RAW_HISTORY_DIRNAME


def raw_history_paths(
    *, image_time_utc: dt.datetime, kind: str
) -> tuple[Path, Path]:
    image_time = _normalize_image_time(image_time_utc)
    directory = raw_history_root() / image_time.strftime("%Y/%m/%d/%H/%M")
    stem = str(kind).strip().lower()
    if stem not in {"infrared", "visible"}:
        raise ValueError("unsupported Geo-satellite image kind")
    return directory / f"{stem}.png", directory / f"{stem}.json"


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
    image_time_utc = _normalize_image_time(
        result.captured_at_utc or result.fetched_at_utc
    )
    _archive_latest_raw_cache(now_utc=result.fetched_at_utc)
    history_png_path, history_metadata_path = raw_history_paths(
        image_time_utc=image_time_utc,
        kind=result.kind,
    )
    history_png_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = _raw_metadata_payload(
        result,
        image_time_utc=image_time_utc,
        image_path=history_png_path,
    )
    _atomic_write_bytes(history_png_path, result.png_bytes)
    _atomic_write_json(history_metadata_path, metadata)

    # Keep the one-frame live cache for compatibility with existing callers.
    png_path = raw_png_path()
    metadata_path = raw_metadata_path()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    purge_legacy_raw_cache(kind=result.kind)
    _atomic_write_bytes(png_path, result.png_bytes)
    _atomic_write_json(
        metadata_path,
        _raw_metadata_payload(result, image_time_utc=image_time_utc),
    )
    return png_path, metadata_path


def read_raw_cache(*, image_time_utc: dt.datetime, kind: str) -> GeoSatelliteDownloadResult | None:
    history_png_path, history_metadata_path = raw_history_paths(
        image_time_utc=image_time_utc,
        kind=kind,
    )
    cached = _read_raw_cache_pair(
        history_png_path,
        history_metadata_path,
        kind=kind,
        expected_time_utc=image_time_utc,
    )
    if cached is not None:
        return cached
    return _read_raw_cache_pair(
        raw_png_path(),
        raw_metadata_path(),
        kind=kind,
        expected_time_utc=image_time_utc,
    )


def read_raw_cache_for_time(
    *,
    target_time_utc: dt.datetime,
    kind: str,
    max_target_age_seconds: int = RAW_HISTORY_MAX_TARGET_AGE_SECONDS,
) -> GeoSatelliteDownloadResult | None:
    target = _as_utc(target_time_utc)
    raw_history_paths(image_time_utc=target, kind=kind)
    max_age = max(0, int(max_target_age_seconds))
    retention_seconds = cache_retention_seconds()
    now = dt.datetime.now(dt.timezone.utc)
    candidates: list[tuple[dt.datetime, Path, bool]] = []
    history_root = raw_history_root()
    if history_root.is_dir():
        for metadata_path in history_root.rglob(f"{str(kind).strip().lower()}.json"):
            captured_at = _read_candidate_time(metadata_path, kind=kind)
            if captured_at is not None:
                age_seconds = (target - captured_at).total_seconds()
                if 0.0 <= age_seconds <= max_age:
                    candidates.append((captured_at, metadata_path, False))

    # Import a legacy latest-only cache lazily so upgrades can replay its one
    # retained frame before a newer live download replaces it.
    legacy_time = _read_candidate_time(raw_metadata_path(), kind=kind)
    if legacy_time is not None:
        age_seconds = (target - legacy_time).total_seconds()
        if 0.0 <= age_seconds <= max_age:
            candidates.append((legacy_time, raw_metadata_path(), True))

    for candidate_time, metadata_path, is_legacy in sorted(
        candidates, key=lambda item: item[0], reverse=True
    ):
        result = _read_raw_cache_pair(
            metadata_path.with_suffix(".png"), metadata_path, kind=kind
        )
        if result is None or result.captured_at_utc is None:
            continue
        fetched_age = (
            now - _as_utc(result.fetched_at_utc)
        ).total_seconds()
        if fetched_age < 0.0 or fetched_age > retention_seconds:
            continue
        if _normalize_image_time(result.captured_at_utc) != _normalize_image_time(
            candidate_time
        ):
            continue
        if is_legacy:
            archived_png, archived_metadata = raw_history_paths(
                image_time_utc=result.captured_at_utc,
                kind=result.kind,
            )
            if not archived_png.is_file() or not archived_metadata.is_file():
                _write_history_raw_cache(result)
        return result
    return None


def _read_candidate_time(
    metadata_path: Path, *, kind: str
) -> dt.datetime | None:
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if not isinstance(metadata, dict):
            return None
        if str(metadata.get("kind", "")).strip().lower() != str(kind).strip().lower():
            return None
        return _normalize_image_time(
            _parse_utc(metadata.get("captured_at_utc") or metadata.get("fetched_at_utc"))
        )
    except Exception:
        return None


def cleanup_geo_satellite_history(
    *,
    now_utc: dt.datetime | None = None,
    max_age_seconds: int | None = None,
) -> None:
    history_root = raw_history_root()
    if not history_root.is_dir():
        return
    retention_seconds = (
        cache_retention_seconds()
        if max_age_seconds is None
        else max(0, int(max_age_seconds))
    )
    cutoff = _as_utc(now_utc or dt.datetime.now(dt.timezone.utc)) - dt.timedelta(
        seconds=retention_seconds
    )
    for metadata_path in history_root.rglob("*.json"):
        image_path = metadata_path.with_suffix(".png")
        if not image_path.is_file():
            metadata_path.unlink(missing_ok=True)
            continue
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            fetched_at = _parse_utc(metadata.get("fetched_at_utc"))
            expired = fetched_at < cutoff
        except Exception:
            expired = True
        if expired:
            image_path.unlink(missing_ok=True)
            metadata_path.unlink(missing_ok=True)
    for image_path in history_root.rglob("*.png"):
        if image_path.with_suffix(".json").is_file():
            continue
        try:
            if image_path.stat().st_mtime < cutoff.timestamp():
                image_path.unlink(missing_ok=True)
        except OSError:
            pass
    for temporary_path in history_root.rglob("*.tmp"):
        try:
            if temporary_path.stat().st_mtime < cutoff.timestamp():
                temporary_path.unlink(missing_ok=True)
        except OSError:
            pass
    for path in sorted(
        history_root.rglob("*"), key=lambda item: len(item.parts), reverse=True
    ):
        if path.is_dir():
            try:
                path.rmdir()
            except OSError:
                pass


def _write_history_raw_cache(
    result: GeoSatelliteDownloadResult,
) -> tuple[Path, Path]:
    image_time = _normalize_image_time(
        result.captured_at_utc or result.fetched_at_utc
    )
    image_path, metadata_path = raw_history_paths(
        image_time_utc=image_time, kind=result.kind
    )
    image_path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_write_bytes(image_path, result.png_bytes)
    _atomic_write_json(
        metadata_path,
        _raw_metadata_payload(result, image_time_utc=image_time, image_path=image_path),
    )
    return image_path, metadata_path


def _archive_latest_raw_cache(*, now_utc: dt.datetime) -> None:
    latest = _read_raw_cache_pair(raw_png_path(), raw_metadata_path())
    if latest is not None and (
        _as_utc(now_utc) - _as_utc(latest.fetched_at_utc)
    ).total_seconds() <= cache_retention_seconds():
        history_png, history_metadata = raw_history_paths(
            image_time_utc=latest.captured_at_utc or latest.fetched_at_utc,
            kind=latest.kind,
        )
        if not history_png.is_file() or not history_metadata.is_file():
            _write_history_raw_cache(latest)


def _raw_metadata_payload(
    result: GeoSatelliteDownloadResult,
    *,
    image_time_utc: dt.datetime,
    image_path: Path | None = None,
) -> dict[str, object]:
    metadata: dict[str, object] = {
        "schema_version": RAW_HISTORY_SCHEMA_VERSION,
        "captured_at_utc": _isoformat_utc(image_time_utc),
        "fetched_at_utc": _isoformat_utc(result.fetched_at_utc),
        "kind": str(result.kind),
        "area": "europe",
        "size": "normal",
        "source_url": str(result.source_url),
        "content_type": result.content_type,
        "sha256": compute_digest(result.png_bytes),
    }
    if image_path is not None:
        metadata["image_path"] = image_path.resolve().relative_to(
            geo_satellite_cache_root().resolve()
        ).as_posix()
    return metadata


def _read_raw_cache_pair(
    png_path: Path,
    metadata_path: Path,
    *,
    kind: str | None = None,
    expected_time_utc: dt.datetime | None = None,
) -> GeoSatelliteDownloadResult | None:
    if not png_path.is_file() or not metadata_path.is_file():
        return None
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if not isinstance(metadata, dict):
            return None
        stored_kind = str(metadata.get("kind", "")).strip().lower()
        if kind is not None and stored_kind != str(kind).strip().lower():
            return None
        schema_version = metadata.get("schema_version")
        if schema_version is not None and int(schema_version) != RAW_HISTORY_SCHEMA_VERSION:
            return None
        captured_at = _parse_utc(
            metadata.get("captured_at_utc") or metadata.get("fetched_at_utc")
        )
        fetched_at = _parse_utc(metadata.get("fetched_at_utc"))
        if expected_time_utc is not None:
            if _normalize_image_time(captured_at) != _normalize_image_time(
                expected_time_utc
            ):
                return None
        png_bytes = png_path.read_bytes()
        expected_digest = metadata.get("sha256")
        if expected_digest is not None and compute_digest(png_bytes) != str(expected_digest):
            return None
        if "image_path" in metadata:
            root = geo_satellite_cache_root().resolve()
            recorded_path = (root / str(metadata["image_path"])).resolve()
            if root not in recorded_path.parents or recorded_path != png_path.resolve():
                return None
        return GeoSatelliteDownloadResult(
            fetched_at_utc=fetched_at,
            captured_at_utc=_normalize_image_time(captured_at),
            kind=stored_kind,  # type: ignore[arg-type]
            source_url=str(metadata.get("source_url", "cache")),
            png_bytes=png_bytes,
            content_type=(
                str(metadata["content_type"])
                if metadata.get("content_type") is not None
                else None
            ),
            cache_path=png_path,
            metadata_path=metadata_path,
        )
    except Exception:
        return None


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    temporary_path = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    temporary_path.write_bytes(data)
    temporary_path.replace(path)


def _atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    temporary_path = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    temporary_path.write_text(
        json.dumps(payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True),
        encoding="utf-8",
    )
    temporary_path.replace(path)


def _normalize_image_time(value: dt.datetime) -> dt.datetime:
    return _as_utc(value).replace(second=0, microsecond=0)


def _as_utc(value: dt.datetime) -> dt.datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=dt.timezone.utc)
    return value.astimezone(dt.timezone.utc)


def _parse_utc(value: object) -> dt.datetime:
    if not isinstance(value, str):
        raise ValueError("missing UTC timestamp")
    return _as_utc(dt.datetime.fromisoformat(value.replace("Z", "+00:00")))


def _isoformat_utc(value: dt.datetime) -> str:
    return _as_utc(value).isoformat().replace("+00:00", "Z")


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
