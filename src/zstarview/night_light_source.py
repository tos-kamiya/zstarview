from __future__ import annotations

import json
import hashlib
import logging
import math
import os
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import rasterio
from pyproj import Geod, Transformer

from .paths import NIGHT_LIGHTS_CACHE_DIR
from .night_lights_constants import (
    NIGHT_LIGHTS_DATASET_VERSION,
    NIGHT_LIGHTS_MANIFEST_URL,
    NIGHT_LIGHTS_TILE_NAMES,
    NIGHT_LIGHTS_TILE_WIDTH_DEG,
)
from .user_agent import build_user_agent

logger = logging.getLogger(__name__)
_GEOD = Geod(ellps="WGS84")

class NightLightsError(RuntimeError):
    """Base error for night-light cache or sampling failures."""

class NightLightsManifestError(NightLightsError):
    """Raised when the GitHub Release manifest is invalid."""

class NightLightsDownloadError(NightLightsError):
    """Raised when a night-light tile cannot be downloaded or validated."""

def _cache_root(cache_root: str | os.PathLike[str] | None = None) -> Path:
    return Path(cache_root or NIGHT_LIGHTS_CACHE_DIR).expanduser()


def _dataset_root(cache_root: str | os.PathLike[str] | None = None) -> Path:
    return _cache_root(cache_root) / NIGHT_LIGHTS_DATASET_VERSION


def _manifest_path(cache_root: str | os.PathLike[str] | None = None) -> Path:
    return _dataset_root(cache_root) / "manifest.json"


def _tile_path(
    tile_name: str,
    manifest: dict[str, object],
    cache_root: str | os.PathLike[str] | None = None,
) -> Path:
    tile = manifest["tiles"][tile_name]
    if not isinstance(tile, dict) or not isinstance(tile.get("path"), str):
        raise NightLightsManifestError(f"Night light manifest has invalid tile {tile_name}.")
    return _dataset_root(cache_root) / Path(tile["path"]).name


def _now_utc() -> datetime:
    return datetime.now(timezone.utc)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_url(url: str, *, timeout_s: float = 60.0, accept: str = "*/*") -> bytes:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": build_user_agent("night-lights"),
            "Accept": accept,
        },
    )
    with urllib.request.urlopen(request, timeout=timeout_s) as response:
        payload = response.read()
    return payload


def _validate_manifest(manifest: object) -> dict[str, object]:
    if not isinstance(manifest, dict):
        raise NightLightsManifestError("Night light manifest must be a JSON object.")
    if manifest.get("dataset_version") != NIGHT_LIGHTS_DATASET_VERSION:
        raise NightLightsManifestError("Night light manifest has an unexpected dataset version.")
    tiles = manifest.get("tiles")
    if not isinstance(tiles, dict):
        raise NightLightsManifestError("Night light manifest has no tiles object.")
    for tile_name in NIGHT_LIGHTS_TILE_NAMES:
        entry = tiles.get(tile_name)
        if not isinstance(entry, dict):
            raise NightLightsManifestError(f"Night light manifest missing tile {tile_name}.")
        path = entry.get("path")
        digest = entry.get("sha256")
        url = entry.get("url")
        if not isinstance(path, str) or Path(path).name != path:
            raise NightLightsManifestError(f"Night light manifest has invalid path for {tile_name}.")
        if not isinstance(digest, str) or len(digest) != 64:
            raise NightLightsManifestError(f"Night light manifest has invalid SHA-256 for {tile_name}.")
        if any(character not in "0123456789abcdefABCDEF" for character in digest):
            raise NightLightsManifestError(f"Night light manifest has invalid SHA-256 for {tile_name}.")
        if not isinstance(url, str) or urlparse(url).scheme != "https":
            raise NightLightsManifestError(f"Night light manifest has invalid URL for {tile_name}.")
    return manifest


def _read_or_build_manifest(
    *,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
) -> dict[str, object]:
    manifest_file = _manifest_path(cache_root)
    if manifest_file.exists():
        try:
            return _validate_manifest(json.loads(manifest_file.read_text(encoding="utf-8")))
        except Exception as exc:
            logger.warning("Night light manifest is unreadable; rebuilding: %s", exc)
        logger.info("Night light manifest is stale or incompatible; rebuilding.")
    try:
        payload = _read_url(NIGHT_LIGHTS_MANIFEST_URL, timeout_s=timeout_s, accept="application/json")
        manifest = _validate_manifest(json.loads(payload.decode("utf-8")))
    except Exception as exc:
        raise NightLightsManifestError(f"Failed to fetch night light manifest: {exc}") from exc
    manifest = dict(manifest)
    manifest["manifest_url"] = NIGHT_LIGHTS_MANIFEST_URL
    manifest["fetched_at_utc"] = _now_utc().isoformat()
    manifest_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file = manifest_file.with_suffix(".tmp")
    tmp_file.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    tmp_file.replace(manifest_file)
    return manifest


def _download_file(url: str, destination: Path, *, timeout_s: float = 300.0) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": build_user_agent("night-lights")},
    )
    digest = hashlib.sha256()
    with tempfile.NamedTemporaryFile(dir=str(destination.parent), delete=False) as tmp:
        tmp_path = Path(tmp.name)
        try:
            with urllib.request.urlopen(request, timeout=timeout_s) as response:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    tmp.write(chunk)
                    digest.update(chunk)
            tmp.flush()
            os.fsync(tmp.fileno())
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise
    tmp_path.replace(destination)
    return digest.hexdigest()


def _validate_geotiff(path: Path, tile_name: str, entry: dict[str, object]) -> None:
    try:
        with rasterio.open(path) as dataset:
            if dataset.count != 1 or dataset.width <= 0 or dataset.height <= 0:
                raise NightLightsDownloadError(f"Invalid GeoTIFF dimensions: {path}")
            if dataset.crs is None or dataset.crs.to_epsg() != 4326:
                raise NightLightsDownloadError(f"Night light tile {tile_name} is not EPSG:4326")
            expected_width = entry.get("width")
            expected_height = entry.get("height")
            if dataset.width != expected_width or dataset.height != expected_height:
                raise NightLightsDownloadError(f"Night light tile {tile_name} has unexpected dimensions")
            expected_resolution = entry.get("resolution_degrees")
            if not isinstance(expected_resolution, list) or len(expected_resolution) != 2:
                raise NightLightsManifestError(f"Night light manifest has invalid resolution for {tile_name}.")
            if not np.allclose(dataset.res, np.asarray(expected_resolution, dtype=np.float64), rtol=0, atol=1.0e-9):
                raise NightLightsDownloadError(f"Night light tile {tile_name} has unexpected resolution")
    except Exception as exc:
        raise NightLightsDownloadError(f"Invalid night light GeoTIFF: {path}") from exc


def _ensure_night_light_tiles(
    *,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    tile_names: set[str] | None = None,
) -> dict[str, Path]:
    dataset_root = _dataset_root(cache_root)
    manifest = _read_or_build_manifest(cache_root=cache_root, timeout_s=timeout_s)
    wanted_tiles = tile_names or {str(name).upper() for name in NIGHT_LIGHTS_TILE_NAMES}
    unknown_tiles = wanted_tiles.difference({str(name).upper() for name in NIGHT_LIGHTS_TILE_NAMES})
    if unknown_tiles:
        raise NightLightsManifestError(f"Unknown night light tiles: {', '.join(sorted(unknown_tiles))}")
    tiles_obj = manifest.get("tiles", {})
    if not isinstance(tiles_obj, dict):
        raise NightLightsManifestError("Night light manifest has invalid tiles.")
    dataset_root.mkdir(parents=True, exist_ok=True)
    resolved_paths: dict[str, Path] = {}
    for tile_name in NIGHT_LIGHTS_TILE_NAMES:
        tile = str(tile_name).upper()
        if tile not in wanted_tiles:
            continue
        entry = tiles_obj.get(tile)
        if not isinstance(entry, dict):
            raise NightLightsManifestError(f"Night light manifest missing tile {tile}.")
        path = _tile_path(tile, manifest, cache_root)
        expected_sha256 = str(entry["sha256"]).lower()
        if not path.exists():
            logger.info("Downloading night light tile %s...", tile)
            try:
                actual_sha256 = _download_file(str(entry["url"]), path, timeout_s=download_timeout_s)
            except Exception as exc:
                raise NightLightsDownloadError(f"Failed to download night light tile {tile}: {exc}") from exc
            if actual_sha256.lower() != expected_sha256:
                path.unlink(missing_ok=True)
                raise NightLightsDownloadError(f"SHA-256 mismatch for night light tile {tile}")
        elif _sha256_file(path).lower() != expected_sha256:
            path.unlink(missing_ok=True)
            raise NightLightsDownloadError(f"Cached SHA-256 mismatch for night light tile {tile}")
        try:
            _validate_geotiff(path, tile, entry)
        except Exception:
            path.unlink(missing_ok=True)
            raise
        resolved_paths[tile] = path
    return resolved_paths


def _tile_name_for_latlon(lat_deg: float, lon_deg: float) -> str:
    lon = float(lon_deg)
    while lon < -180.0:
        lon += 360.0
    while lon >= 180.0:
        lon -= 360.0
    lat = max(-90.0, min(90.0, float(lat_deg)))
    col = int(math.floor((lon + 180.0) / NIGHT_LIGHTS_TILE_WIDTH_DEG))
    col = max(0, min(3, col))
    row = 0 if lat >= 0.0 else 1
    return f"{'ABCD'[col]}{1 + row}"


def _required_tile_names(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    max_distance_km: float,
) -> set[str]:
    """Return tiles touched by a radial footprint around the observer."""
    azimuths = np.arange(0.0, 360.0, 5.0, dtype=np.float64)
    distances_m = np.full(azimuths.shape, max(0.0, float(max_distance_km)) * 1000.0)
    longitudes, latitudes, _ = _GEOD.fwd(
        np.full(azimuths.shape, float(observer_lon_deg), dtype=np.float64),
        np.full(azimuths.shape, float(observer_lat_deg), dtype=np.float64),
        azimuths,
        distances_m,
    )
    names = {_tile_name_for_latlon(float(observer_lat_deg), float(observer_lon_deg))}
    names.update(
        _tile_name_for_latlon(float(lat), float(lon))
        for lat, lon in zip(latitudes.tolist(), longitudes.tolist())
    )
    return names


def _open_dataset_cached(path_str: str, mtime_ns: int) -> rasterio.io.DatasetReader:
    _ = mtime_ns
    return rasterio.open(path_str)


def _sample_dataset_points(
    dataset: rasterio.io.DatasetReader,
    points_lon_lat: list[tuple[float, float]],
) -> np.ndarray:
    if not points_lon_lat:
        return np.empty(0, dtype=np.float64)
    coords = points_lon_lat
    if dataset.crs is not None and str(dataset.crs).upper() != "EPSG:4326":
        transformer = Transformer.from_crs("EPSG:4326", dataset.crs, always_xy=True)
        xs, ys = transformer.transform(
            np.asarray([pt[0] for pt in coords], dtype=np.float64),
            np.asarray([pt[1] for pt in coords], dtype=np.float64),
        )
        coords = list(zip(xs.tolist(), ys.tolist()))
    values: list[float] = []
    for sample in dataset.sample(coords, indexes=1, masked=True):
        if np.ma.is_masked(sample):
            values.append(0.0)
            continue
        value = float(np.asarray(sample, dtype=np.float64).reshape(-1)[0])
        if not math.isfinite(value):
            value = 0.0
        values.append(max(0.0, value))
    return np.asarray(values, dtype=np.float64)


def _sample_ray_night_light_samples(
    *,
    tile_paths: dict[str, Path],
    observer_lat_deg: float,
    observer_lon_deg: float,
    azimuth_deg: float,
    distances_m: np.ndarray,
) -> np.ndarray:
    if distances_m.size == 0:
        return np.zeros(0, dtype=np.float64)
    lon0 = float(observer_lon_deg)
    lat0 = float(observer_lat_deg)
    az = np.full(distances_m.shape, float(azimuth_deg), dtype=np.float64)
    lon_arr, lat_arr, _ = _GEOD.fwd(
        np.full(distances_m.shape, lon0, dtype=np.float64),
        np.full(distances_m.shape, lat0, dtype=np.float64),
        az,
        distances_m,
    )
    grouped_points: dict[str, list[tuple[float, float]]] = {}
    grouped_indices: dict[str, list[int]] = {}
    for index, (lon, lat) in enumerate(zip(lon_arr.tolist(), lat_arr.tolist())):
        tile_name = _tile_name_for_latlon(float(lat), float(lon))
        grouped_points.setdefault(tile_name, []).append((float(lon), float(lat)))
        grouped_indices.setdefault(tile_name, []).append(index)

    samples = np.zeros(distances_m.shape, dtype=np.float64)
    for tile_name, coords in grouped_points.items():
        path = tile_paths.get(tile_name)
        if path is None:
            continue
        stat = path.stat()
        dataset = _open_dataset_cached(str(path), int(stat.st_mtime_ns))
        tile_samples = _sample_dataset_points(dataset, coords)
        indices = grouped_indices[tile_name]
        samples[np.asarray(indices, dtype=np.int64)] = tile_samples
    return samples
