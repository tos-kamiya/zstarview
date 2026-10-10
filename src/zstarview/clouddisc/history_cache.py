from __future__ import annotations

import hashlib
import json
import logging
import pickle
from datetime import datetime, timedelta, timezone
from pathlib import Path
from uuid import uuid4

from .altaz_grid import CloudAltAzGrid
from .types import CloudSourceData

logger = logging.getLogger(__name__)
CLOUD_HISTORY_DIR_NAME = "source_history"
CLOUD_HISTORY_RETENTION = timedelta(hours=24)
CLOUD_HISTORY_MAX_TARGET_AGE = timedelta(minutes=30)
CLOUD_HISTORY_SCHEMA_VERSION = 1


def save_cloud_source_history(
    source: CloudSourceData,
    *,
    fetched_at_utc: datetime,
    cache_root: str | Path,
) -> Path | None:
    """Save a completed normalized source for later offline rendering."""
    grid = source.altaz_grid
    if not isinstance(grid, CloudAltAzGrid):
        logger.warning("Skipping cloud history source without a completed grid")
        return None
    root = Path(cache_root).resolve()
    fetched_at = _as_utc(fetched_at_utc)
    directory = root / CLOUD_HISTORY_DIR_NAME / fetched_at.strftime("%Y/%m/%d/%H/%M")
    record_id = uuid4().hex
    artifact_path = directory / f"{record_id}.pkl"
    metadata_path = directory / f"{record_id}.json"
    artifact_tmp = artifact_path.with_name(f".{artifact_path.name}.{uuid4().hex}.tmp")
    metadata_tmp = metadata_path.with_name(f".{metadata_path.name}.{uuid4().hex}.tmp")
    try:
        directory.mkdir(parents=True, exist_ok=True)
        data = pickle.dumps(source, protocol=pickle.HIGHEST_PROTOCOL)
        artifact_tmp.write_bytes(data)
        artifact_tmp.replace(artifact_path)
        source_key = source.source_key
        metadata = {
            "schema_version": CLOUD_HISTORY_SCHEMA_VERSION,
            "record_id": record_id,
            "fetched_at_utc": _isoformat(fetched_at),
            "source_time_utc": _isoformat(_as_utc(source.time_utc)),
            "observer_latitude_deg": float(grid.observer_lat),
            "observer_longitude_deg": float(grid.observer_lon),
            "source_key": {
                "satellite": str(getattr(source_key, "satellite", source.satellite)),
                "provider": getattr(source_key, "provider", None),
                "timeslot_utc": _isoformat(
                    _as_utc(getattr(source_key, "timeslot_utc", source.time_utc))
                ),
            },
            "source_paths": [str(path) for path in source.src_paths],
            "artifact_path": artifact_path.relative_to(root).as_posix(),
            "artifact_sha256": hashlib.sha256(data).hexdigest(),
        }
        metadata_tmp.write_text(
            json.dumps(metadata, ensure_ascii=True, separators=(",", ":"), sort_keys=True),
            encoding="utf-8",
        )
        metadata_tmp.replace(metadata_path)
        return metadata_path
    except Exception:
        logger.warning("Failed to save cloud source history", exc_info=True)
        artifact_tmp.unlink(missing_ok=True)
        metadata_tmp.unlink(missing_ok=True)
        artifact_path.unlink(missing_ok=True)
        return None


def load_cloud_source_for_time(
    *,
    observer_lat: float,
    observer_lon: float,
    target_time_utc: datetime,
    cache_root: str | Path,
    max_target_age: timedelta = CLOUD_HISTORY_MAX_TARGET_AGE,
) -> CloudSourceData | None:
    """Load the closest completed source no later than target, without network."""
    root = Path(cache_root).resolve()
    history_root = root / CLOUD_HISTORY_DIR_NAME
    target = _as_utc(target_time_utc)
    candidates: list[tuple[datetime, Path, dict[str, object]]] = []
    if not history_root.is_dir():
        return None
    for metadata_path in history_root.rglob("*.json"):
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            if not isinstance(metadata, dict):
                continue
            if int(metadata.get("schema_version", -1)) != CLOUD_HISTORY_SCHEMA_VERSION:
                continue
            fetched_at = _parse_time(metadata.get("fetched_at_utc"))
            source_time = _parse_time(metadata.get("source_time_utc"))
            if fetched_at > target or source_time > target:
                continue
            if target - source_time > max_target_age:
                continue
            if abs(float(metadata["observer_latitude_deg"]) - float(observer_lat)) > 0.00005:
                continue
            if abs(float(metadata["observer_longitude_deg"]) - float(observer_lon)) > 0.00005:
                continue
            raw_artifact = metadata.get("artifact_path")
            if not isinstance(raw_artifact, str):
                continue
            artifact_path = (root / raw_artifact).resolve()
            if root not in artifact_path.parents or not artifact_path.is_file():
                continue
            candidates.append((source_time, artifact_path, metadata))
        except Exception:
            logger.warning("Skipping invalid cloud history metadata", exc_info=True)
    for _, artifact_path, metadata in sorted(candidates, key=lambda item: item[0], reverse=True):
        try:
            data = artifact_path.read_bytes()
            if hashlib.sha256(data).hexdigest() != metadata.get("artifact_sha256"):
                continue
            source = pickle.loads(data)
            if not isinstance(source, CloudSourceData):
                continue
            grid = source.altaz_grid
            if not isinstance(grid, CloudAltAzGrid):
                continue
            if abs(float(grid.observer_lat) - float(observer_lat)) > 0.00005:
                continue
            if abs(float(grid.observer_lon) - float(observer_lon)) > 0.00005:
                continue
            if _as_utc(source.time_utc) > target:
                continue
            return source
        except Exception:
            logger.warning("Skipping invalid cloud history artifact", exc_info=True)
    return None


def cleanup_cloud_source_history(
    *,
    cache_root: str | Path,
    now_utc: datetime | None = None,
) -> None:
    root = Path(cache_root) / CLOUD_HISTORY_DIR_NAME
    if not root.is_dir():
        return
    cutoff = _as_utc(now_utc or datetime.now(timezone.utc)) - CLOUD_HISTORY_RETENTION
    referenced_artifacts: set[Path] = set()
    for metadata_path in root.rglob("*.json"):
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            artifact_path = (Path(cache_root) / str(metadata.get("artifact_path", ""))).resolve()
            if Path(cache_root).resolve() in artifact_path.parents:
                referenced_artifacts.add(artifact_path)
            fetched_at = _parse_time(metadata.get("fetched_at_utc"))
            if fetched_at >= cutoff:
                continue
            if Path(cache_root).resolve() in artifact_path.parents:
                artifact_path.unlink(missing_ok=True)
            metadata_path.unlink(missing_ok=True)
        except Exception:
            logger.warning("Failed to prune cloud history record", exc_info=True)
    for artifact_path in root.rglob("*.pkl"):
        if artifact_path.resolve() in referenced_artifacts:
            continue
        try:
            artifact_mtime = datetime.fromtimestamp(
                artifact_path.stat().st_mtime, tz=timezone.utc
            )
            if artifact_mtime < cutoff:
                artifact_path.unlink(missing_ok=True)
        except OSError:
            logger.warning("Failed to prune orphaned cloud history artifact")
    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if path.is_dir():
            try:
                path.rmdir()
            except OSError:
                pass


def _parse_time(value: object) -> datetime:
    if not isinstance(value, str):
        raise ValueError("missing UTC timestamp")
    return _as_utc(datetime.fromisoformat(value.replace("Z", "+00:00")))


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _isoformat(value: datetime) -> str:
    return _as_utc(value).isoformat().replace("+00:00", "Z")
