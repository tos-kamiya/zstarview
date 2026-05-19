#!/usr/bin/env python3
"""Download Overture building data around a point and store derived tiles."""

from __future__ import annotations

import argparse
import json
import logging
import math
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Sequence
from urllib.request import Request, urlopen

from zstarview.clouddisc.types import DownloadCancelledError
from zstarview.data import build_derived_tile_index
from zstarview.paths import CACHE_PATH
from zstarview.utils.latlon_format import (
    LAT_LON_CACHE_DECIMALS,
    format_lat_lon_cache_segment,
)

EARTH_RADIUS_KM = 6371.0088
DEFAULT_FETCH_RADIUS_KM = 2.5
DEFAULT_MIN_BUILDING_HEIGHT_M = 0.0
DEFAULT_STOREY_HEIGHT_M = 3.5
OVERTURE_CACHE_TTL_DAYS = 30
CACHE_METADATA_FILENAME = "cache_meta.json"
OVERTURE_RELEASE_CHECK_METADATA_FILENAME = "overture_buildings_release_check_meta.json"
OVERTURE_RELEASE_CHECK_MAX_AGE = timedelta(hours=24)
OVERTURE_RELEASE_CATALOG_URL = "https://stac.overturemaps.org/catalog.json"
OVERTURE_RELEASE_CATALOG_TIMEOUT_SECONDS = 10.0

logger = logging.getLogger(__name__)


def _read_json_dict(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if isinstance(loaded, dict):
        return loaded
    return None


def _metadata_path_for_derived_dir(derived_dir: Path) -> Path:
    return derived_dir.parent / CACHE_METADATA_FILENAME


def release_check_metadata_path(*, cache_root_dir: Path | None = None) -> Path:
    return Path(cache_root_dir or CACHE_PATH) / OVERTURE_RELEASE_CHECK_METADATA_FILENAME


def _parse_optional_utc(text: object) -> datetime | None:
    if not isinstance(text, str) or not text.strip():
        return None
    return _normalize_utc(datetime.fromisoformat(text.replace("Z", "+00:00")))


def read_overture_release_check_metadata(
    *,
    cache_root_dir: Path | None = None,
) -> dict[str, object] | None:
    return _read_json_dict(release_check_metadata_path(cache_root_dir=cache_root_dir))


def write_overture_release_check_metadata(
    *,
    payload: dict[str, object],
    cache_root_dir: Path | None = None,
) -> Path:
    metadata_path = release_check_metadata_path(cache_root_dir=cache_root_dir)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return metadata_path


def fetch_latest_overture_release(
    *,
    catalog_url: str = OVERTURE_RELEASE_CATALOG_URL,
    timeout_seconds: float = OVERTURE_RELEASE_CATALOG_TIMEOUT_SECONDS,
) -> str:
    request = Request(catalog_url, headers={"User-Agent": "zstarview/1.0"})
    with urlopen(request, timeout=float(timeout_seconds)) as response:  # nosec: B310
        payload = json.loads(response.read().decode("utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("Unexpected Overture catalog payload")
    latest = payload.get("latest")
    if not isinstance(latest, str) or not latest.strip():
        raise ValueError("Overture catalog payload is missing latest release")
    return latest.strip()


def resolve_overture_release_for_cache_root(
    *,
    cache_root_dir: Path | None = None,
    now_utc: datetime | None = None,
    catalog_url: str = OVERTURE_RELEASE_CATALOG_URL,
    timeout_seconds: float = OVERTURE_RELEASE_CATALOG_TIMEOUT_SECONDS,
    abort_event: threading.Event | None = None,
) -> str | None:
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    metadata = read_overture_release_check_metadata(cache_root_dir=cache_root_dir) or {}
    last_release_check_utc = _parse_optional_utc(metadata.get("last_release_check_utc"))
    last_seen_overture_release = metadata.get("last_seen_overture_release")
    if (
        last_release_check_utc is not None
        and now - last_release_check_utc < OVERTURE_RELEASE_CHECK_MAX_AGE
    ):
        if isinstance(last_seen_overture_release, str) and last_seen_overture_release.strip():
            return last_seen_overture_release.strip()
        return None

    if abort_event is not None and abort_event.is_set():
        return None

    checked_payload: dict[str, object] = dict(metadata)
    checked_payload["last_release_check_utc"] = now.isoformat()
    checked_payload.setdefault("checked_source", "stac")
    try:
        latest_release = fetch_latest_overture_release(
            catalog_url=catalog_url,
            timeout_seconds=timeout_seconds,
        )
    except Exception as exc:
        checked_payload["checked_success"] = False
        checked_payload["checked_error"] = str(exc)
        write_overture_release_check_metadata(
            payload=checked_payload,
            cache_root_dir=cache_root_dir,
        )
        logger.warning("Failed to resolve latest Overture release: %s", exc)
        return None

    checked_payload["checked_success"] = True
    checked_payload["checked_error"] = None
    checked_payload["last_seen_overture_release"] = latest_release
    write_overture_release_check_metadata(
        payload=checked_payload,
        cache_root_dir=cache_root_dir,
    )
    return latest_release


def _run_download_command(command: list[str], *, abort_event: threading.Event | None = None) -> subprocess.CompletedProcess[str]:
    if abort_event is None:
        return subprocess.run(command, check=False)

    proc = subprocess.Popen(command)
    try:
        while True:
            if abort_event.is_set():
                proc.terminate()
                try:
                    proc.wait(timeout=1.0)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                raise DownloadCancelledError("Cancelled while running overturemaps download")
            returncode = proc.poll()
            if returncode is not None:
                return subprocess.CompletedProcess(command, returncode)
            time.sleep(0.1)
    finally:
        if proc.poll() is None and abort_event is not None and abort_event.is_set():
            proc.kill()
            proc.wait()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download Overture building data around a lat/lon point and write "
            "derived urban-outline tiles into the app cache directory."
        )
    )
    parser.add_argument("--lat", type=float, required=True, help="Center latitude in degrees.")
    parser.add_argument("--lon", type=float, required=True, help="Center longitude in degrees.")
    parser.add_argument(
        "--radius-km",
        type=float,
        default=DEFAULT_FETCH_RADIUS_KM,
        help="Download radius in kilometers used to build the bbox (default: 2.5).",
    )
    parser.add_argument(
        "--derived-root-dir",
        type=Path,
        default=None,
        help="Root directory where derived datasets will be written.",
    )
    parser.add_argument(
        "--min-building-height-m",
        type=float,
        default=DEFAULT_MIN_BUILDING_HEIGHT_M,
        help="Ignore buildings lower than this height in meters.",
    )
    parser.add_argument(
        "--feature-type",
        default="building",
        choices=("building", "building_part"),
        help="Overture feature type to download (default: building).",
    )
    parser.add_argument(
        "--format",
        default="geojsonseq",
        choices=("geojson", "geojsonseq"),
        help="Intermediate download format passed to overturemaps (default: geojsonseq).",
    )
    parser.add_argument(
        "--overturemaps-bin",
        default="overturemaps",
        help="Path to the overturemaps CLI executable (default: overturemaps).",
    )
    parser.add_argument(
        "--dataset-name",
        default=None,
        help="Optional dataset directory name. Defaults to a lat/lon-based cache key.",
    )
    parser.add_argument(
        "--keep-download",
        type=Path,
        default=None,
        help="Optional path where the raw downloaded file will be copied for inspection.",
    )
    parser.add_argument(
        "--no-stac",
        action="store_true",
        help="Pass --no-stac to the overturemaps CLI.",
    )
    return parser


def bbox_from_point(lat_deg: float, lon_deg: float, radius_km: float) -> tuple[float, float, float, float]:
    if radius_km <= 0:
        raise ValueError("radius_km must be positive")
    if not -90.0 <= lat_deg <= 90.0:
        raise ValueError("lat must be between -90 and 90 degrees")
    if not -180.0 <= lon_deg <= 180.0:
        raise ValueError("lon must be between -180 and 180 degrees")

    angular_distance = radius_km / EARTH_RADIUS_KM
    lat_delta_deg = math.degrees(angular_distance)
    cos_lat = math.cos(math.radians(lat_deg))
    if abs(cos_lat) < 1e-9:
        lon_delta_deg = 180.0
    else:
        lon_delta_deg = math.degrees(angular_distance / abs(cos_lat))

    return (
        max(-180.0, lon_deg - lon_delta_deg),
        max(-90.0, lat_deg - lat_delta_deg),
        min(180.0, lon_deg + lon_delta_deg),
        min(90.0, lat_deg + lat_delta_deg),
    )


def derive_dataset_name(
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    feature_type: str,
    min_building_height_m: float = DEFAULT_MIN_BUILDING_HEIGHT_M,
) -> str:
    return (
        f"overture_{feature_type}"
        f"_lat{format_lat_lon_cache_segment(lat_deg, decimals=LAT_LON_CACHE_DECIMALS)}"
        f"_lon{format_lat_lon_cache_segment(lon_deg, decimals=LAT_LON_CACHE_DECIMALS)}"
        f"_r{radius_km:.1f}km"
        f"_h{float(min_building_height_m):.1f}m"
    ).replace("-", "m").replace(".", "p")


def output_suffix_for_format(fmt: str) -> str:
    return {
        "geojson": ".geojson",
        "geojsonseq": ".geojsonseq",
    }[fmt]


def build_download_command(
    *,
    overturemaps_bin: str,
    bbox: tuple[float, float, float, float],
    feature_type: str,
    fmt: str,
    output_path: Path,
    no_stac: bool,
) -> list[str]:
    west, south, east, north = bbox
    command = [
        overturemaps_bin,
        "download",
        f"--bbox={west:.8f},{south:.8f},{east:.8f},{north:.8f}",
        "-f",
        fmt,
        "--type",
        feature_type,
        "-o",
        str(output_path),
    ]
    if no_stac:
        command.append("--no-stac")
    return command


def derived_dataset_metadata_path(derived_dir: Path) -> Path:
    return _metadata_path_for_derived_dir(derived_dir)


def read_derived_dataset_metadata(derived_dir: Path) -> dict[str, object] | None:
    if not derived_dir.exists():
        return None
    return _read_json_dict(derived_dataset_metadata_path(derived_dir))


def read_derived_dataset_overture_release(
    derived_dir: Path,
) -> str | None:
    payload = read_derived_dataset_metadata(derived_dir)
    if payload is None:
        return None
    raw_release = payload.get("overture_release")
    if isinstance(raw_release, str) and raw_release.strip():
        return raw_release.strip()
    return None


def read_derived_dataset_fetched_at_utc(
    derived_dir: Path,
    *,
    now_utc: datetime | None = None,
    migrate_missing: bool = True,
) -> datetime | None:
    metadata_path = derived_dataset_metadata_path(derived_dir)
    payload = read_derived_dataset_metadata(derived_dir) or {}
    raw_fetched_at = payload.get("fetched_at_utc")
    if isinstance(raw_fetched_at, str):
        try:
            return _parse_utc(raw_fetched_at)
        except Exception:
            pass
    if not migrate_missing:
        return None
    fallback_mtime = _path_mtime_utc(metadata_path)
    if fallback_mtime is None and derived_dir.exists():
        fallback_mtime = _path_mtime_utc(derived_dir)
    fetched_at_utc = fallback_mtime or _normalize_utc(now_utc or datetime.now(timezone.utc))
    logger.info(
        "Urban-outline cache metadata missing; recording provisional fetched_at_utc for offline reuse: %s",
        derived_dir,
    )
    payload["fetched_at_utc"] = fetched_at_utc.isoformat()
    payload.setdefault("dataset_name", derived_dir.parent.name)
    write_derived_dataset_metadata(derived_dir, payload=payload)
    return fetched_at_utc


def is_derived_dataset_stale(
    derived_dir: Path,
    *,
    ttl_days: int = OVERTURE_CACHE_TTL_DAYS,
    now_utc: datetime | None = None,
    expected_overture_release: str | None = None,
) -> bool:
    fetched_at_utc = read_derived_dataset_fetched_at_utc(derived_dir, now_utc=now_utc)
    if fetched_at_utc is None:
        return True
    if expected_overture_release is not None:
        cached_release = read_derived_dataset_overture_release(derived_dir)
        if cached_release is not None and cached_release != expected_overture_release:
            return True
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    return (now - fetched_at_utc) > timedelta(days=max(0, int(ttl_days)))


def write_derived_dataset_metadata(
    derived_dir: Path,
    *,
    payload: dict[str, object],
) -> Path:
    metadata_path = derived_dataset_metadata_path(derived_dir)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return metadata_path


def import_overture_buildings(
    *,
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    derived_root_dir: Path,
    min_building_height_m: float,
    feature_type: str,
    fmt: str,
    overturemaps_bin: str,
    dataset_name: str | None,
    keep_download: Path | None,
    no_stac: bool,
    overture_release: str | None = None,
    skip_release_lookup: bool = False,
    now_utc: datetime | None = None,
    quiet: bool = False,
    abort_event: threading.Event | None = None,
) -> Path:
    overturemaps_path = shutil.which(overturemaps_bin)
    if overturemaps_path is None:
        raise FileNotFoundError(
            f"Could not find overturemaps CLI: {overturemaps_bin!r}. "
            "Install it separately or pass --overturemaps-bin."
        )

    bbox = bbox_from_point(lat_deg, lon_deg, radius_km)
    return import_overture_buildings_for_bbox(
        bbox=bbox,
        derived_root_dir=derived_root_dir,
        min_building_height_m=min_building_height_m,
        feature_type=feature_type,
        fmt=fmt,
        overturemaps_bin=overturemaps_path,
        dataset_name=dataset_name
        or derive_dataset_name(
            lat_deg,
            lon_deg,
            radius_km,
            feature_type,
            min_building_height_m,
        ),
        keep_download=keep_download,
        no_stac=no_stac,
        overture_release=overture_release,
        skip_release_lookup=skip_release_lookup,
        query_lat_deg=lat_deg,
        query_lon_deg=lon_deg,
        query_radius_km=radius_km,
        now_utc=now_utc,
        quiet=quiet,
        abort_event=abort_event,
    )


def import_overture_buildings_for_bbox(
    *,
    bbox: tuple[float, float, float, float],
    derived_root_dir: Path,
    min_building_height_m: float,
    feature_type: str,
    fmt: str,
    overturemaps_bin: str,
    dataset_name: str,
    keep_download: Path | None,
    no_stac: bool,
    overture_release: str | None = None,
    skip_release_lookup: bool = False,
    query_lat_deg: float | None = None,
    query_lon_deg: float | None = None,
    query_radius_km: float | None = None,
    now_utc: datetime | None = None,
    quiet: bool = False,
    abort_event: threading.Event | None = None,
) -> Path:
    overturemaps_path = shutil.which(overturemaps_bin) or overturemaps_bin
    fetched_at_utc = _normalize_utc(now_utc or datetime.now(timezone.utc))
    if skip_release_lookup:
        current_overture_release = overture_release
    else:
        current_overture_release = (
            overture_release
            or resolve_overture_release_for_cache_root(
                cache_root_dir=Path(CACHE_PATH),
                now_utc=fetched_at_utc,
                abort_event=abort_event,
            )
        )
    dataset_dir_name = dataset_name or derive_dataset_name(
        query_lat_deg or 0.0,
        query_lon_deg or 0.0,
        query_radius_km or 0.0,
        feature_type,
        min_building_height_m,
    )
    derived_dir = derived_root_dir / dataset_dir_name / "bldg"
    tile_path = derived_dir / f"{dataset_dir_name}.json"

    with tempfile.TemporaryDirectory(prefix="overture-import-") as temp_dir_str:
        temp_dir = Path(temp_dir_str)
        download_path = temp_dir / f"download{output_suffix_for_format(fmt)}"
        command = build_download_command(
            overturemaps_bin=overturemaps_path,
            bbox=bbox,
            feature_type=feature_type,
            fmt=fmt,
            output_path=download_path,
            no_stac=no_stac,
        )
        completed = _run_download_command(command, abort_event=abort_event)
        if completed.returncode != 0:
            raise RuntimeError(f"overturemaps download failed with return code {completed.returncode}")
        if keep_download is not None:
            keep_download.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(download_path, keep_download)

        payload = build_derived_tile_payload(
            download_path,
            dataset_dir_name=dataset_dir_name,
            bbox=bbox,
            min_building_height_m=min_building_height_m,
            feature_type=feature_type,
            query_lat_deg=query_lat_deg,
            query_lon_deg=query_lon_deg,
            query_radius_km=query_radius_km,
            download_format=fmt,
        )

    tile_path.parent.mkdir(parents=True, exist_ok=True)
    tile_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    build_derived_tile_index.main(["--derived-dir", str(derived_dir)], quiet=quiet)
    write_derived_dataset_metadata(
        derived_dir,
        payload={
            "bbox": {
                "west": float(bbox[0]),
                "south": float(bbox[1]),
                "east": float(bbox[2]),
                "north": float(bbox[3]),
            },
            "dataset_name": dataset_dir_name,
            "feature_type": str(feature_type),
            "fetched_at_utc": fetched_at_utc.isoformat(),
            "overture_release": current_overture_release,
            "min_building_height_m": float(min_building_height_m),
            "query_lat_deg": None if query_lat_deg is None else float(query_lat_deg),
            "query_lon_deg": None if query_lon_deg is None else float(query_lon_deg),
            "query_radius_km": None if query_radius_km is None else float(query_radius_km),
        },
    )
    return derived_dir


def _parse_utc(text: str) -> datetime:
    return _normalize_utc(datetime.fromisoformat(str(text).replace("Z", "+00:00")))


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _path_mtime_utc(path: Path) -> datetime | None:
    try:
        mtime = path.stat().st_mtime
    except FileNotFoundError:
        return None
    except OSError:
        return None
    return _normalize_utc(datetime.fromtimestamp(mtime, tz=timezone.utc))


def build_derived_tile_payload(
    download_path: Path,
    *,
    dataset_dir_name: str,
    bbox: tuple[float, float, float, float],
    min_building_height_m: float,
    feature_type: str,
    query_lat_deg: float | None,
    query_lon_deg: float | None,
    query_radius_km: float | None,
    download_format: str,
) -> dict[str, object]:
    buildings: list[dict[str, object]] = []

    for feature in iter_download_features(download_path, fmt=download_format):
        building = convert_feature_to_building(feature, min_building_height_m=min_building_height_m)
        if building is not None:
            buildings.append(building)

    west, south, east, north = bbox
    query_payload: dict[str, object] = {
        "bbox": {
            "west": west,
            "south": south,
            "east": east,
            "north": north,
        },
        "format": download_format,
    }
    if query_lat_deg is not None and query_lon_deg is not None and query_radius_km is not None:
        query_payload.update(
            {
                "lat_deg": query_lat_deg,
                "lon_deg": query_lon_deg,
                "radius_km": query_radius_km,
            }
        )

    return {
        "schema_version": 1,
        "source": {
            "format": "Overture-Buildings",
            "feature_type": feature_type,
            "query": query_payload,
        },
        "tile": {
            "id": dataset_dir_name,
            "bbox": {
                "min_lat": south,
                "min_lon": west,
                "max_lat": north,
                "max_lon": east,
            },
        },
        "filters": {
            "min_height_m": float(min_building_height_m),
            "storey_height_m": DEFAULT_STOREY_HEIGHT_M,
        },
        "buildings": buildings,
    }


def iter_download_features(download_path: Path, *, fmt: str) -> Sequence[dict[str, object]]:
    if fmt == "geojsonseq":
        features: list[dict[str, object]] = []
        with download_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                row = json.loads(line)
                if isinstance(row, dict):
                    features.append(row)
        return tuple(features)

    payload = json.loads(download_path.read_text(encoding="utf-8"))
    if isinstance(payload, dict) and isinstance(payload.get("features"), list):
        return tuple(row for row in payload["features"] if isinstance(row, dict))
    return ()


def convert_feature_to_building(
    feature: dict[str, object],
    *,
    min_building_height_m: float,
) -> dict[str, object] | None:
    geometry = feature.get("geometry")
    if not isinstance(geometry, dict):
        return None
    rings = geometry_to_rings(geometry)
    if not rings:
        return None

    properties = feature.get("properties")
    if not isinstance(properties, dict):
        properties = {}
    height_m = resolve_feature_height_m(properties)
    if height_m is None or height_m < float(min_building_height_m):
        return None

    bbox = compute_building_bbox(rings)
    building_id = feature.get("id")
    if not isinstance(building_id, str) or not building_id:
        building_id = str(properties.get("id") or f"feature-{len(rings)}")

    return {
        "id": building_id,
        "height_m": height_m,
        "min_height_m": resolve_feature_min_height_m(properties),
        "height_source": detect_height_source(properties),
        "parent_building_id": resolve_parent_building_id(properties),
        "bbox": {
            "min_lat": bbox[0],
            "min_lon": bbox[1],
            "max_lat": bbox[2],
            "max_lon": bbox[3],
        },
        "rings": [
            [[lon, lat] for lon, lat in ring]
            for ring in rings
        ],
    }


def geometry_to_rings(geometry: dict[str, object]) -> tuple[tuple[tuple[float, float], ...], ...]:
    geom_type = geometry.get("type")
    coordinates = geometry.get("coordinates")
    if geom_type == "Polygon":
        return parse_polygon_rings(coordinates)
    if geom_type == "MultiPolygon":
        if not isinstance(coordinates, list):
            return ()
        rings: list[tuple[tuple[float, float], ...]] = []
        for polygon in coordinates:
            rings.extend(parse_polygon_rings(polygon))
        return tuple(rings)
    return ()


def parse_polygon_rings(raw_polygon: object) -> tuple[tuple[tuple[float, float], ...], ...]:
    if not isinstance(raw_polygon, list):
        return ()
    rings: list[tuple[tuple[float, float], ...]] = []
    for raw_ring in raw_polygon:
        ring = parse_ring(raw_ring)
        if ring:
            rings.append(ring)
    return tuple(rings)


def parse_ring(raw_ring: object) -> tuple[tuple[float, float], ...]:
    if not isinstance(raw_ring, list):
        return ()
    points: list[tuple[float, float]] = []
    for raw_point in raw_ring:
        if not isinstance(raw_point, list) or len(raw_point) < 2:
            continue
        try:
            lon = float(raw_point[0])
            lat = float(raw_point[1])
        except (TypeError, ValueError):
            continue
        points.append((lon, lat))
    if len(points) < 4:
        return ()
    if points[0] != points[-1]:
        points.append(points[0])
    return tuple(points)


def resolve_feature_height_m(properties: dict[str, object]) -> float | None:
    for key in ("height", "roof_height"):
        value = parse_numeric(properties.get(key))
        if value is not None and value > 0:
            return value
    for key in ("num_floors", "num_floors_above_ground"):
        value = parse_numeric(properties.get(key))
        if value is not None and value > 0:
            return value * DEFAULT_STOREY_HEIGHT_M
    return None


def resolve_feature_min_height_m(properties: dict[str, object]) -> float:
    value = parse_numeric(properties.get("min_height"))
    if value is None or value <= 0.0:
        return 0.0
    return value


def detect_height_source(properties: dict[str, object]) -> str:
    if parse_numeric(properties.get("height")) is not None:
        return "height"
    if parse_numeric(properties.get("roof_height")) is not None:
        return "roof_height"
    for key in ("num_floors", "num_floors_above_ground"):
        if parse_numeric(properties.get(key)) is not None:
            return f"{key}*{DEFAULT_STOREY_HEIGHT_M}"
    return "unknown"


def resolve_parent_building_id(properties: dict[str, object]) -> str | None:
    parent_building_id = properties.get("building_id")
    if isinstance(parent_building_id, str) and parent_building_id:
        return parent_building_id
    return None


def parse_numeric(value: object) -> float | None:
    try:
        if value is None:
            return None
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def compute_building_bbox(
    rings_lonlat: Sequence[Sequence[tuple[float, float]]],
) -> tuple[float, float, float, float]:
    min_lat = float("inf")
    min_lon = float("inf")
    max_lat = float("-inf")
    max_lon = float("-inf")
    for ring in rings_lonlat:
        for lon, lat in ring:
            min_lat = min(min_lat, lat)
            min_lon = min(min_lon, lon)
            max_lat = max(max_lat, lat)
            max_lon = max(max_lon, lon)
    return (min_lat, min_lon, max_lat, max_lon)


def main(argv: Sequence[str] | None = None, *, quiet: bool = False) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    from zstarview.paths import OVERTURE_DERIVED_ROOT_DIR

    derived_root_dir = args.derived_root_dir or Path(OVERTURE_DERIVED_ROOT_DIR)
    derived_dir = import_overture_buildings(
        lat_deg=float(args.lat),
        lon_deg=float(args.lon),
        radius_km=float(args.radius_km),
        derived_root_dir=derived_root_dir,
        min_building_height_m=float(args.min_building_height_m),
        feature_type=args.feature_type,
        fmt=args.format,
        overturemaps_bin=args.overturemaps_bin,
        dataset_name=args.dataset_name,
        keep_download=args.keep_download,
        no_stac=bool(args.no_stac),
    )
    if not quiet:
        print(f"[ok] imported overture buildings -> {derived_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
