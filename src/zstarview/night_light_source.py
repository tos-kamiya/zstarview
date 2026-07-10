from __future__ import annotations

import json
import logging
import math
import os
import re
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Geod, Transformer

from .paths import NIGHT_LIGHTS_CACHE_DIR
from .night_lights_constants import (
    NIGHT_LIGHTS_DATASET_VERSION,
    NIGHT_LIGHTS_PAGE_URL,
    NIGHT_LIGHTS_TILE_NAMES,
    NIGHT_LIGHTS_TILE_WIDTH_DEG,
)
from .user_agent import build_user_agent

logger = logging.getLogger(__name__)
_GEOD = Geod(ellps="WGS84")
_TILE_URL_RE = re.compile(
    r'href="(?P<url>[^"]*BlackMarble_2016_(?P<tile>[A-D][12])_geo_gray\.tif)"',
    re.IGNORECASE,
)

class NightLightsError(RuntimeError):
    """Base error for night-light cache or sampling failures."""

class NightLightsManifestError(NightLightsError):
    """Raised when the NASA map page cannot be parsed."""

class NightLightsDownloadError(NightLightsError):
    """Raised when a night-light tile cannot be downloaded or validated."""

def _cache_root(cache_root: str | os.PathLike[str] | None = None) -> Path:
    return Path(cache_root or NIGHT_LIGHTS_CACHE_DIR).expanduser()


def _dataset_root(cache_root: str | os.PathLike[str] | None = None) -> Path:
    return _cache_root(cache_root) / NIGHT_LIGHTS_DATASET_VERSION


def _manifest_path(cache_root: str | os.PathLike[str] | None = None) -> Path:
    return _dataset_root(cache_root) / "manifest.json"


def _tile_filename(tile_name: str) -> str:
    return f"BlackMarble_2016_{tile_name}_geo_gray.tif"


def _tile_path(tile_name: str, cache_root: str | os.PathLike[str] | None = None) -> Path:
    return _dataset_root(cache_root) / _tile_filename(tile_name)


def _now_utc() -> datetime:
    return datetime.now(timezone.utc)


def _read_url(url: str, *, timeout_s: float = 60.0) -> str:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": build_user_agent("night-lights"),
            "Accept": "text/html,application/xhtml+xml",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout_s) as response:
        payload = response.read()
    return payload.decode("utf-8", errors="replace")


def _parse_tile_urls(page_html: str) -> dict[str, str]:
    matched: dict[str, str] = {}
    wanted = {str(name).upper() for name in NIGHT_LIGHTS_TILE_NAMES}
    for match in _TILE_URL_RE.finditer(page_html):
        tile_name = match.group("tile").upper()
        if tile_name in wanted:
            matched[tile_name] = match.group("url")
    missing = sorted(wanted.difference(matched))
    if missing:
        raise NightLightsManifestError(
            f"Night light page did not expose expected tile URLs: {', '.join(missing)}"
        )
    return matched


def _read_or_build_manifest(
    *,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
) -> dict[str, object]:
    manifest_file = _manifest_path(cache_root)
    if manifest_file.exists():
        try:
            manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
            if (
                isinstance(manifest, dict)
                and manifest.get("dataset_version") == NIGHT_LIGHTS_DATASET_VERSION
                and isinstance(manifest.get("tile_urls"), dict)
            ):
                return manifest
        except Exception as exc:
            logger.warning("Night light manifest is unreadable; rebuilding: %s", exc)
        logger.info("Night light manifest is stale or incompatible; rebuilding.")
    page_html = _read_url(NIGHT_LIGHTS_PAGE_URL, timeout_s=timeout_s)
    tile_urls = _parse_tile_urls(page_html)
    manifest = {
        "dataset_version": NIGHT_LIGHTS_DATASET_VERSION,
        "source_page_url": NIGHT_LIGHTS_PAGE_URL,
        "tile_urls": tile_urls,
        "fetched_at_utc": _now_utc().isoformat(),
    }
    manifest_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file = manifest_file.with_suffix(".tmp")
    tmp_file.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    tmp_file.replace(manifest_file)
    return manifest


def _download_file(url: str, destination: Path, *, timeout_s: float = 300.0) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": build_user_agent("night-lights")},
    )
    with tempfile.NamedTemporaryFile(dir=str(destination.parent), delete=False) as tmp:
        tmp_path = Path(tmp.name)
        try:
            with urllib.request.urlopen(request, timeout=timeout_s) as response:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    tmp.write(chunk)
            tmp.flush()
            os.fsync(tmp.fileno())
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise
    tmp_path.replace(destination)


def _validate_geotiff(path: Path) -> None:
    try:
        with rasterio.open(path) as dataset:
            if dataset.count < 1 or dataset.width <= 0 or dataset.height <= 0:
                raise NightLightsDownloadError(f"Invalid GeoTIFF dimensions: {path}")
    except Exception as exc:
        raise NightLightsDownloadError(f"Invalid night light GeoTIFF: {path}") from exc


def _ensure_night_light_tiles(
    *,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
) -> dict[str, Path]:
    dataset_root = _dataset_root(cache_root)
    manifest = _read_or_build_manifest(cache_root=cache_root, timeout_s=timeout_s)
    tile_urls_obj = manifest.get("tile_urls", {})
    if not isinstance(tile_urls_obj, dict):
        raise NightLightsManifestError("Night light manifest has invalid tile_urls.")
    dataset_root.mkdir(parents=True, exist_ok=True)
    resolved_paths: dict[str, Path] = {}
    for tile_name in NIGHT_LIGHTS_TILE_NAMES:
        tile = str(tile_name).upper()
        if tile not in tile_urls_obj:
            raise NightLightsManifestError(f"Night light manifest missing tile {tile}.")
        path = _tile_path(tile, cache_root)
        if not path.exists():
            logger.info("Downloading night light tile %s...", tile)
            try:
                _download_file(str(tile_urls_obj[tile]), path, timeout_s=download_timeout_s)
            except Exception as exc:
                raise NightLightsDownloadError(f"Failed to download night light tile {tile}: {exc}") from exc
        _validate_geotiff(path)
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

