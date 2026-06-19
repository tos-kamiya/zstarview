from __future__ import annotations

import functools
import json
import logging
import math
import os
import re
import tempfile
import time
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np
import rasterio
from pyproj import Geod, Transformer

from .paths import NIGHT_LIGHTS_CACHE_DIR
from .terrain.horizon import DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM
from .terrain.horizon import EARTH_MEAN_RADIUS_M
from .terrain.horizon import compute_apparent_altitudes
from .user_agent import build_user_agent

logger = logging.getLogger(__name__)

NIGHT_LIGHTS_DATASET_VERSION = "2016_grayscale_500m"
NIGHT_LIGHTS_PAGE_URL = "https://science.nasa.gov/earth/earth-observatory/earth-at-night/maps/"
NIGHT_LIGHTS_TILE_NAMES = ("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2")
NIGHT_LIGHTS_TILE_COUNT = len(NIGHT_LIGHTS_TILE_NAMES)
NIGHT_LIGHTS_TILE_WIDTH_DEG = 90.0
NIGHT_LIGHTS_TILE_HEIGHT_DEG = 90.0
NIGHT_LIGHTS_MAX_DISTANCE_KM = 128.0
NIGHT_LIGHTS_DISTANCE_STEP_KM = 3.0
NIGHT_LIGHTS_AZIMUTH_SMOOTHING_KERNEL_NIGHT = np.asarray([1, 14, 62, 102, 62, 14, 1], dtype=np.float64)
NIGHT_LIGHTS_BAND_CENTER_OFFSET_DEG = 1.5
NIGHT_LIGHTS_BAND_HALF_WIDTH_DEG = 1.5
NIGHT_LIGHTS_MAX_ALPHA = 0.48
NIGHT_LIGHTS_GLOW_RGB = (244, 246, 248)
NIGHT_LIGHTS_RGB = NIGHT_LIGHTS_GLOW_RGB
NIGHT_LIGHTS_SUN_BLEND_START_ALT_DEG = -6.0
NIGHT_LIGHTS_DISTANCE_BAND_EDGES_KM = DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM[1:]
NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG = 8.0
NIGHT_LIGHTS_NEIGHBORHOOD_CHUNK_SIZE = 4096
NIGHT_LIGHTS_NEIGHBORHOOD_WEIGHT_STEP_DEG = 0.5
NIGHT_LIGHTS_ALTITUDE_MIN_DEG = -90.0
NIGHT_LIGHTS_ALTITUDE_MAX_DEG = 90.0
NIGHT_LIGHTS_ALTITUDE_STEP_DEG = 1.0

_TILE_URL_RE = re.compile(
    r'href="(?P<url>[^"]*BlackMarble_2016_(?P<tile>[A-D][12])_geo_gray\.tif)"',
    re.IGNORECASE,
)
_GEOD = Geod(ellps="WGS84")


class NightLightsError(RuntimeError):
    """Base error for night-light cache or sampling failures."""


class NightLightsManifestError(NightLightsError):
    """Raised when the NASA map page cannot be parsed."""


class NightLightsDownloadError(NightLightsError):
    """Raised when a night-light tile cannot be downloaded or validated."""


@dataclass(frozen=True)
class NightLightGlowSample:
    azimuth_deg: float
    horizon_alt_deg: float
    strength: float


@dataclass(frozen=True)
class NightLightGlowProfile:
    samples: tuple[NightLightGlowSample, ...]
    sun_alt_deg: float
    band_center_offset_deg: float = NIGHT_LIGHTS_BAND_CENTER_OFFSET_DEG
    band_half_width_deg: float = NIGHT_LIGHTS_BAND_HALF_WIDTH_DEG
    altitude_bins_deg: tuple[float, ...] = ()
    alpha_grid: tuple[tuple[float, ...], ...] = ()


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


def _terrain_profile_key(
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
) -> tuple[tuple[float, float], ...]:
    if not terrain_profile_altaz:
        return ()
    return tuple(
        (round(float(alt_deg), 3), round(float(az_deg) % 360.0, 3))
        for alt_deg, az_deg in terrain_profile_altaz
    )


def _terrain_context_key(
    *,
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
    terrain_profile_distances_m: Sequence[float] | None,
    terrain_secondary_ridges_altaz_layers: Sequence[Sequence[tuple[float, float]]] | None,
    terrain_secondary_ridges_distances_m_layers: Sequence[Sequence[float]] | None,
    terrain_sample_azimuths_deg: Sequence[float] | None = None,
    terrain_sample_distances_m: Sequence[float] | np.ndarray | None = None,
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None = None,
) -> tuple[
    tuple[tuple[float, float], ...],
    tuple[float, ...],
    tuple[tuple[tuple[float, float], ...], ...],
    tuple[tuple[float, ...], ...],
    tuple[int, int, int] | None,
]:
    return (
        _terrain_profile_key(terrain_profile_altaz),
        tuple(round(float(distance_m), 3) for distance_m in terrain_profile_distances_m or ()),
        tuple(
            tuple(
                (round(float(alt_deg), 3), round(float(az_deg) % 360.0, 3))
                for alt_deg, az_deg in layer
            )
            for layer in terrain_secondary_ridges_altaz_layers or ()
        ),
        tuple(
            tuple(round(float(distance_m), 3) for distance_m in layer)
            for layer in terrain_secondary_ridges_distances_m_layers or ()
        ),
        _terrain_sample_grid_key(
            terrain_sample_azimuths_deg=terrain_sample_azimuths_deg,
            terrain_sample_distances_m=terrain_sample_distances_m,
            terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
        ),
    )


def _terrain_sample_grid_key(
    *,
    terrain_sample_azimuths_deg: Sequence[float] | None,
    terrain_sample_distances_m: Sequence[float] | np.ndarray | None,
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None,
) -> tuple[int, int, int] | None:
    if (
        terrain_sample_azimuths_deg is None
        or terrain_sample_distances_m is None
        or terrain_sample_terrain_elevation_m is None
    ):
        return None
    return (
        id(terrain_sample_azimuths_deg),
        id(terrain_sample_distances_m),
        id(terrain_sample_terrain_elevation_m),
    )


@functools.lru_cache(maxsize=8)
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


def _night_light_distance_attenuation(distances_m: np.ndarray) -> np.ndarray:
    distances_km = np.maximum(np.asarray(distances_m, dtype=np.float64) / 1000.0, 1.0)
    return 1.0 / np.square(distances_km)


def _apply_night_light_sample_floor(
    samples: np.ndarray,
    visibility_mask: np.ndarray | None,
    *,
    floor_value: float,
) -> np.ndarray:
    if samples.size == 0:
        return np.zeros(0, dtype=np.float64)
    result = np.asarray(samples, dtype=np.float64)
    floor = max(0.0, float(floor_value))
    if floor <= 0.0:
        return result if visibility_mask is None else np.where(np.asarray(visibility_mask, dtype=bool), result, 0.0)
    if visibility_mask is None:
        return result + floor
    mask = np.asarray(visibility_mask, dtype=bool)
    if mask.shape != result.shape:
        return result + floor
    return np.where(mask, result + floor, 0.0)


def _wrap_azimuth_delta_deg(left_deg: np.ndarray, right_deg: np.ndarray) -> np.ndarray:
    """Return the signed shortest azimuth delta in degrees."""
    left = np.asarray(left_deg, dtype=np.float64)
    right = np.asarray(right_deg, dtype=np.float64)
    return np.remainder(left - right + 180.0, 360.0) - 180.0


@functools.lru_cache(maxsize=8)
def _gaussian_weight_lut(sigma_deg: float, step_deg: float) -> np.ndarray:
    sigma = max(1.0e-6, float(sigma_deg))
    step = max(1.0e-6, float(step_deg))
    max_delta_deg = 180.0
    delta_bins = int(math.ceil(max_delta_deg / step))
    deltas = np.arange(delta_bins + 1, dtype=np.float64) * step
    return np.exp(-0.5 * np.square(deltas / sigma))


def _lookup_gaussian_weights(delta_deg: np.ndarray, *, sigma_deg: float, step_deg: float) -> np.ndarray:
    delta = np.abs(np.asarray(delta_deg, dtype=np.float64))
    lut = _gaussian_weight_lut(sigma_deg, step_deg)
    step = max(1.0e-6, float(step_deg))
    indices = np.clip(np.rint(delta / step).astype(np.int64), 0, lut.size - 1)
    return lut[indices]


def _accumulate_local_glow_strengths(
    *,
    source_azimuths_deg: np.ndarray,
    source_altitudes_deg: np.ndarray,
    source_strengths: np.ndarray,
    target_azimuths_deg: np.ndarray,
    target_altitudes_deg: np.ndarray,
    sigma_deg: float,
    chunk_size: int = NIGHT_LIGHTS_NEIGHBORHOOD_CHUNK_SIZE,
) -> np.ndarray:
    """Accumulate glow by summing nearby source points in az/alt space."""
    source_azimuths = np.asarray(source_azimuths_deg, dtype=np.float64).reshape(-1)
    source_altitudes = np.asarray(source_altitudes_deg, dtype=np.float64).reshape(-1)
    source_strengths_arr = np.asarray(source_strengths, dtype=np.float64).reshape(-1)
    target_azimuths = np.asarray(target_azimuths_deg, dtype=np.float64).reshape(-1)
    target_altitudes = np.asarray(target_altitudes_deg, dtype=np.float64).reshape(-1)
    if (
        source_azimuths.size == 0
        or source_altitudes.size == 0
        or source_strengths_arr.size == 0
        or target_azimuths.size == 0
        or target_altitudes.size == 0
    ):
        return np.zeros(target_azimuths.shape, dtype=np.float64)
    if not (
        source_azimuths.size == source_altitudes.size == source_strengths_arr.size
        and target_azimuths.size == target_altitudes.size
    ):
        raise ValueError("source and target arrays must have matching lengths")

    chunk = max(1, int(chunk_size))
    accumulated = np.zeros(target_azimuths.shape, dtype=np.float64)
    for start in range(0, source_strengths_arr.size, chunk):
        end = min(source_strengths_arr.size, start + chunk)
        source_az_chunk = source_azimuths[start:end][:, None]
        source_alt_chunk = source_altitudes[start:end][:, None]
        source_strength_chunk = source_strengths_arr[start:end][:, None]
        delta_az = _wrap_azimuth_delta_deg(source_az_chunk, target_azimuths[None, :])
        az_weights = _lookup_gaussian_weights(
            delta_az,
            sigma_deg=sigma_deg,
            step_deg=NIGHT_LIGHTS_NEIGHBORHOOD_WEIGHT_STEP_DEG,
        )
        delta_alt = source_alt_chunk - target_altitudes[None, :]
        alt_weights = _lookup_gaussian_weights(
            delta_alt,
            sigma_deg=sigma_deg,
            step_deg=NIGHT_LIGHTS_NEIGHBORHOOD_WEIGHT_STEP_DEG,
        )
        weights = az_weights * alt_weights
        accumulated += np.sum(source_strength_chunk * weights, axis=0)
    return accumulated


def _target_altitude_bins() -> np.ndarray:
    return np.arange(
        NIGHT_LIGHTS_ALTITUDE_MIN_DEG,
        NIGHT_LIGHTS_ALTITUDE_MAX_DEG + (0.5 * NIGHT_LIGHTS_ALTITUDE_STEP_DEG),
        NIGHT_LIGHTS_ALTITUDE_STEP_DEG,
        dtype=np.float64,
    )


def _accumulate_local_glow_field(
    *,
    source_azimuths_deg: np.ndarray,
    source_altitudes_deg: np.ndarray,
    source_strengths: np.ndarray,
    target_azimuths_deg: np.ndarray,
    target_altitudes_deg: np.ndarray,
    sigma_deg: float,
    chunk_size: int = NIGHT_LIGHTS_NEIGHBORHOOD_CHUNK_SIZE,
) -> np.ndarray:
    """Accumulate glow onto a target altitude/azimuth grid."""
    source_azimuths = np.asarray(source_azimuths_deg, dtype=np.float64).reshape(-1)
    source_altitudes = np.asarray(source_altitudes_deg, dtype=np.float64).reshape(-1)
    source_strengths_arr = np.asarray(source_strengths, dtype=np.float64).reshape(-1)
    target_azimuths = np.asarray(target_azimuths_deg, dtype=np.float64).reshape(-1)
    target_altitudes = np.asarray(target_altitudes_deg, dtype=np.float64).reshape(-1)
    if (
        source_azimuths.size == 0
        or source_altitudes.size == 0
        or source_strengths_arr.size == 0
        or target_azimuths.size == 0
        or target_altitudes.size == 0
    ):
        return np.zeros((target_altitudes.size, target_azimuths.size), dtype=np.float64)
    if not source_azimuths.size == source_altitudes.size == source_strengths_arr.size:
        raise ValueError("source arrays must have matching lengths")

    chunk = max(1, int(chunk_size))
    field = np.zeros((target_altitudes.size, target_azimuths.size), dtype=np.float64)
    for start in range(0, source_strengths_arr.size, chunk):
        end = min(source_strengths_arr.size, start + chunk)
        source_az_chunk = source_azimuths[start:end]
        source_alt_chunk = source_altitudes[start:end]
        source_strength_chunk = np.clip(source_strengths_arr[start:end], 0.0, None)
        if not np.any(source_strength_chunk > 0.0):
            continue
        delta_az = _wrap_azimuth_delta_deg(source_az_chunk[:, None], target_azimuths[None, :])
        az_weights = _lookup_gaussian_weights(
            delta_az,
            sigma_deg=sigma_deg,
            step_deg=NIGHT_LIGHTS_NEIGHBORHOOD_WEIGHT_STEP_DEG,
        )
        delta_alt = source_alt_chunk[:, None] - target_altitudes[None, :]
        alt_weights = _lookup_gaussian_weights(
            delta_alt,
            sigma_deg=sigma_deg,
            step_deg=NIGHT_LIGHTS_NEIGHBORHOOD_WEIGHT_STEP_DEG,
        )
        field += alt_weights.T @ (source_strength_chunk[:, None] * az_weights)
    return field


def _flatten_glow_source_matrix(
    values: np.ndarray,
    sample_altitudes: np.ndarray,
    sample_azimuths: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    matrix = np.asarray(values, dtype=np.float64)
    altitudes_arr = np.asarray(sample_altitudes, dtype=np.float64)
    azimuths = np.asarray(sample_azimuths, dtype=np.float64).reshape(-1, 1)
    if matrix.ndim != 2:
        raise ValueError("values must be a 2D array")
    if matrix.shape[0] != azimuths.shape[0]:
        raise ValueError("matrix shape must match sample azimuths and altitudes")
    altitudes: np.ndarray
    if altitudes_arr.ndim == 1:
        altitudes = altitudes_arr.reshape(1, -1)
        if matrix.shape[1] != altitudes.shape[1]:
            raise ValueError("matrix shape must match sample azimuths and altitudes")
        altitudes = np.repeat(altitudes, matrix.shape[0], axis=0)
    elif altitudes_arr.ndim == 2:
        altitudes = altitudes_arr
        if matrix.shape != altitudes.shape:
            raise ValueError("matrix shape must match sample azimuths and altitudes")
    else:
        raise ValueError("sample_altitudes must be 1D or 2D")
    return (
        np.repeat(azimuths, matrix.shape[1], axis=1).reshape(-1),
        altitudes.reshape(-1),
        matrix.reshape(-1),
    )


def _sample_ray_brightness_curve(
    *,
    tile_paths: dict[str, Path],
    observer_lat_deg: float,
    observer_lon_deg: float,
    azimuth_deg: float,
    distances_m: np.ndarray,
    visibility_mask: np.ndarray | None = None,
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

    samples = _apply_night_light_sample_floor(
        samples,
        visibility_mask,
        floor_value=0.0,
    )
    attenuation = _night_light_distance_attenuation(distances_m)
    return np.cumsum(samples * attenuation)


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


def _sample_ray_brightness(
    *,
    tile_paths: dict[str, Path],
    observer_lat_deg: float,
    observer_lon_deg: float,
    azimuth_deg: float,
    distances_m: np.ndarray,
) -> float:
    curve = _sample_ray_brightness_curve(
        tile_paths=tile_paths,
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        azimuth_deg=azimuth_deg,
        distances_m=distances_m,
    )
    if curve.size == 0:
        return 0.0
    return float(curve[-1])


def _surface_point_apparent_altitudes(
    distances_m: np.ndarray,
    *,
    observer_height_m: float,
    refraction_coefficient: float,
) -> np.ndarray:
    if distances_m.size == 0:
        return np.zeros(0, dtype=np.float64)
    target_elevation_m = np.zeros_like(np.asarray(distances_m, dtype=np.float64))
    return compute_apparent_altitudes(
        observer_elevation_m=max(0.0, float(observer_height_m)),
        target_elevation_m=target_elevation_m,
        surface_distance_m=np.asarray(distances_m, dtype=np.float64),
        earth_radius_m=EARTH_MEAN_RADIUS_M,
        refraction_coefficient=float(refraction_coefficient),
    )


def _terrain_sample_source_altitude_rows(
    *,
    terrain_sample_distances_m: Sequence[float] | Sequence[Sequence[float]] | np.ndarray | None,
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None,
    source_distances_m: np.ndarray,
    observer_height_m: float,
    refraction_coefficient: float,
) -> np.ndarray | None:
    if terrain_sample_distances_m is None or terrain_sample_terrain_elevation_m is None:
        return None
    terrain_distances_arr = np.asarray(terrain_sample_distances_m, dtype=np.float64)
    terrain_elevation = np.asarray(terrain_sample_terrain_elevation_m, dtype=np.float64)
    if terrain_distances_arr.size == 0 or terrain_elevation.size == 0:
        return None
    if terrain_elevation.ndim != 2:
        raise ValueError("terrain_sample_terrain_elevation_m must be a 2D array")
    if terrain_distances_arr.ndim == 1:
        terrain_distances = terrain_distances_arr.reshape(-1)
        if terrain_elevation.shape[1] != terrain_distances.size:
            raise ValueError(
                "terrain_sample_terrain_elevation_m must match terrain_sample_distances_m"
            )
    elif terrain_distances_arr.ndim == 2:
        if terrain_distances_arr.shape != terrain_elevation.shape:
            raise ValueError(
                "terrain_sample_distances_m must match terrain_sample_terrain_elevation_m"
            )
        terrain_distances = np.asarray(terrain_distances_arr[0], dtype=np.float64).reshape(-1)
        if not np.allclose(terrain_distances_arr, terrain_distances[np.newaxis, :], equal_nan=True):
            terrain_distances = None
    else:
        raise ValueError("terrain_sample_distances_m must be 1D or 2D")
    source_distances = np.asarray(source_distances_m, dtype=np.float64).reshape(-1)
    if source_distances.size == 0:
        return np.zeros((terrain_elevation.shape[0], 0), dtype=np.float64)

    if terrain_distances_arr.ndim == 2 and terrain_distances is None:
        surface_distances = terrain_distances_arr
    else:
        surface_distances = np.asarray(terrain_distances, dtype=np.float64)[np.newaxis, :]
    terrain_apparent_altitudes = compute_apparent_altitudes(
        observer_elevation_m=max(0.0, float(observer_height_m)),
        target_elevation_m=terrain_elevation,
        surface_distance_m=surface_distances,
        earth_radius_m=EARTH_MEAN_RADIUS_M,
        refraction_coefficient=float(refraction_coefficient),
    )
    if terrain_distances is not None and np.array_equal(terrain_distances, source_distances):
        return np.asarray(terrain_apparent_altitudes, dtype=np.float64)
    rows = [
        np.interp(
            source_distances,
            terrain_distances if terrain_distances is not None else np.asarray(terrain_distances_arr[row_index], dtype=np.float64).reshape(-1),
            np.asarray(row, dtype=np.float64),
            left=float(row[0]),
            right=float(row[-1]),
        )
        for row_index, row in enumerate(np.asarray(terrain_apparent_altitudes, dtype=np.float64))
    ]
    return np.asarray(rows, dtype=np.float64)


def _terrain_visibility_threshold_curve(
    *,
    azimuth_deg: float,
    distances_m: np.ndarray,
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
    terrain_profile_distances_m: Sequence[float] | None,
    terrain_secondary_ridges_altaz_layers: Sequence[Sequence[tuple[float, float]]] | None,
    terrain_secondary_ridges_distances_m_layers: Sequence[Sequence[float]] | None,
) -> np.ndarray | None:
    if distances_m.size == 0:
        return None

    az_key = round(float(azimuth_deg) % 360.0, 3)
    terrain_points: list[tuple[float, float]] = []
    if terrain_profile_altaz and terrain_profile_distances_m and len(terrain_profile_altaz) == len(terrain_profile_distances_m):
        for (alt_deg, az_deg), distance_m in zip(terrain_profile_altaz, terrain_profile_distances_m, strict=True):
            if round(float(az_deg) % 360.0, 3) != az_key:
                continue
            terrain_points.append((float(distance_m), float(alt_deg)))
    if terrain_secondary_ridges_altaz_layers and terrain_secondary_ridges_distances_m_layers:
        if len(terrain_secondary_ridges_altaz_layers) == len(terrain_secondary_ridges_distances_m_layers):
            for layer, layer_distances_m in zip(
                terrain_secondary_ridges_altaz_layers,
                terrain_secondary_ridges_distances_m_layers,
                strict=True,
            ):
                if len(layer) != len(layer_distances_m):
                    continue
                for (alt_deg, az_deg), distance_m in zip(layer, layer_distances_m, strict=True):
                    if round(float(az_deg) % 360.0, 3) != az_key:
                        continue
                    terrain_points.append((float(distance_m), float(alt_deg)))

    if not terrain_points:
        if terrain_profile_altaz:
            return None
        return np.full(distances_m.shape, -np.inf, dtype=np.float64)

    ordered_points = sorted(
        (distance_m, altitude_deg)
        for distance_m, altitude_deg in terrain_points
        if math.isfinite(float(distance_m)) and math.isfinite(float(altitude_deg))
    )
    if not ordered_points:
        return np.full(distances_m.shape, -np.inf, dtype=np.float64)

    threshold = np.full(distances_m.shape, -np.inf, dtype=np.float64)
    running_max = -np.inf
    point_index = 0
    for index, distance_m in enumerate(np.asarray(distances_m, dtype=np.float64)):
        while point_index < len(ordered_points) and float(ordered_points[point_index][0]) <= float(distance_m) + 1.0e-9:
            running_max = max(running_max, float(ordered_points[point_index][1]))
            point_index += 1
        threshold[index] = running_max
    return threshold


def _build_azimuth_grid(
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
) -> tuple[np.ndarray, np.ndarray]:
    if terrain_profile_altaz:
        az_values = np.asarray([float(az) % 360.0 for _, az in terrain_profile_altaz], dtype=np.float64)
        alt_values = np.asarray([float(alt) for alt, _ in terrain_profile_altaz], dtype=np.float64)
        order = np.argsort(az_values)
        return az_values[order], alt_values[order]
    az_values = np.linspace(0.0, 360.0, num=180, endpoint=False, dtype=np.float64)
    alt_values = np.zeros_like(az_values)
    return az_values, alt_values


def _terrain_band_target_mask(
    *,
    az_grid: np.ndarray,
    target_altitudes: np.ndarray,
    band_distance_m: float,
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
    terrain_profile_distances_m: Sequence[float] | None,
    terrain_secondary_ridges_altaz_layers: Sequence[Sequence[tuple[float, float]]] | None,
    terrain_secondary_ridges_distances_m_layers: Sequence[Sequence[float]] | None,
) -> np.ndarray:
    az_values = np.asarray(az_grid, dtype=np.float64).reshape(-1)
    target_altitudes_arr = np.asarray(target_altitudes, dtype=np.float64).reshape(-1)
    if az_values.size == 0 or target_altitudes_arr.size == 0:
        return np.zeros(0, dtype=bool)
    if az_values.size != target_altitudes_arr.size:
        raise ValueError("az_grid and target_altitudes must have the same length")
    mask = np.zeros(az_values.shape, dtype=bool)
    distances = np.asarray([float(band_distance_m)], dtype=np.float64)
    for index, az in enumerate(az_values.tolist()):
        threshold = _terrain_visibility_threshold_curve(
            azimuth_deg=float(az),
            distances_m=distances,
            terrain_profile_altaz=terrain_profile_altaz,
            terrain_profile_distances_m=terrain_profile_distances_m,
            terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
            terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
        )
        if threshold is None:
            mask[index] = True
            continue
        mask[index] = bool(target_altitudes_arr[index] > float(np.asarray(threshold, dtype=np.float64)[0]))
    return mask


def _terrain_band_target_altaz_mask(
    *,
    az_grid: np.ndarray,
    target_altitudes: np.ndarray,
    band_distance_m: float,
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
    terrain_profile_distances_m: Sequence[float] | None,
    terrain_secondary_ridges_altaz_layers: Sequence[Sequence[tuple[float, float]]] | None,
    terrain_secondary_ridges_distances_m_layers: Sequence[Sequence[float]] | None,
) -> np.ndarray:
    az_values = np.asarray(az_grid, dtype=np.float64).reshape(-1)
    alt_values = np.asarray(target_altitudes, dtype=np.float64).reshape(-1)
    if az_values.size == 0 or alt_values.size == 0:
        return np.zeros((alt_values.size, az_values.size), dtype=bool)
    mask = np.zeros((alt_values.size, az_values.size), dtype=bool)
    distances = np.asarray([float(band_distance_m)], dtype=np.float64)
    for index, az in enumerate(az_values.tolist()):
        threshold = _terrain_visibility_threshold_curve(
            azimuth_deg=float(az),
            distances_m=distances,
            terrain_profile_altaz=terrain_profile_altaz,
            terrain_profile_distances_m=terrain_profile_distances_m,
            terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
            terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
        )
        if threshold is None:
            mask[:, index] = True
            continue
        mask[:, index] = alt_values > float(np.asarray(threshold, dtype=np.float64)[0])
    return mask


def _inner_ring_mask_distance_m(distances_m: np.ndarray, sample_index: int) -> float:
    if sample_index <= 0:
        return 0.0
    distances = np.asarray(distances_m, dtype=np.float64).reshape(-1)
    index = max(0, min(int(sample_index) - 1, distances.size - 1))
    return float(distances[index])


def _circular_weighted_smooth(values: np.ndarray, kernel_weights: np.ndarray) -> np.ndarray:
    kernel = np.asarray(kernel_weights, dtype=np.float64)
    if values.size == 0 or kernel.size == 0:
        return values
    if kernel.ndim != 1 or kernel.size % 2 == 0:
        raise ValueError("kernel_weights must be a 1D array with an odd length")
    kernel = kernel / float(np.sum(kernel))
    radius = kernel.size // 2
    padded = np.pad(np.asarray(values, dtype=np.float64), radius, mode="wrap")
    smoothed = np.convolve(padded, kernel, mode="same")
    return smoothed[radius:-radius]


def _circular_smooth(values: np.ndarray) -> np.ndarray:
    return _circular_weighted_smooth(values, NIGHT_LIGHTS_AZIMUTH_SMOOTHING_KERNEL_NIGHT)


def _circular_weighted_smooth_matrix(values: np.ndarray, kernel_weights: np.ndarray) -> np.ndarray:
    kernel = np.asarray(kernel_weights, dtype=np.float64)
    if values.ndim != 2 or values.size == 0 or kernel.size == 0:
        return values
    if kernel.ndim != 1 or kernel.size % 2 == 0:
        raise ValueError("kernel_weights must be a 1D array with an odd length")
    kernel = kernel / float(np.sum(kernel))
    radius = kernel.size // 2
    padded = np.pad(np.asarray(values, dtype=np.float64), ((radius, radius), (0, 0)), mode="wrap")
    smoothed = np.zeros_like(values, dtype=np.float64)
    for offset, weight in enumerate(kernel):
        smoothed += float(weight) * padded[offset : offset + values.shape[0], :]
    return smoothed


def _smoothstep(edge0: float, edge1: float, value: float) -> float:
    lo = float(edge0)
    hi = float(edge1)
    x = float(value)
    if hi <= lo:
        return 0.0 if x <= lo else 1.0
    t = max(0.0, min(1.0, (x - lo) / (hi - lo)))
    return t * t * (3.0 - 2.0 * t)


def night_light_strength_factor(sun_alt_deg: float) -> float:
    sun_alt = float(sun_alt_deg)
    if sun_alt >= 0.0:
        return 0.0
    return 1.0 - _smoothstep(NIGHT_LIGHTS_SUN_BLEND_START_ALT_DEG, 0.0, sun_alt)


def _distance_band_ranges_km(max_distance_km: float) -> tuple[tuple[float, float], ...]:
    max_distance = float(max_distance_km)
    if max_distance <= 0.0:
        return ()
    band_edges = [
        float(edge)
        for edge in NIGHT_LIGHTS_DISTANCE_BAND_EDGES_KM
        if math.isfinite(float(edge)) and 0.0 < float(edge) <= max_distance + 1.0e-9
    ]
    if not band_edges or band_edges[-1] < max_distance - 1.0e-9:
        band_edges.append(max_distance)
    band_ranges: list[tuple[float, float]] = []
    band_start = 0.0
    for band_end in band_edges:
        if band_end <= band_start:
            continue
        band_ranges.append((band_start, band_end))
        band_start = band_end
    return tuple(band_ranges)


def _build_night_light_glow_profile_from_samples(
    *,
    az_grid: np.ndarray,
    horizon_alt_values: np.ndarray,
    distances_m: np.ndarray,
    sample_matrix: np.ndarray,
    source_altitudes: np.ndarray,
    terrain_profile_key: tuple[tuple[float, float], ...],
    terrain_profile_distances_key: tuple[float, ...],
    terrain_secondary_ridges_key: tuple[tuple[tuple[float, float], ...], ...],
    terrain_secondary_ridges_distances_key: tuple[tuple[float, ...], ...],
    max_distance_km: float,
    smooth_strengths: bool = True,
) -> NightLightGlowProfile | None:
    band_ranges_km = _distance_band_ranges_km(max_distance_km)
    if not band_ranges_km:
        return None
    source_altitudes_arr = np.asarray(source_altitudes, dtype=np.float64)
    if source_altitudes_arr.ndim == 1:
        source_altitudes_arr = np.repeat(source_altitudes_arr[np.newaxis, :], az_grid.size, axis=0)
    elif source_altitudes_arr.ndim != 2:
        raise ValueError("source_altitudes must be 1D or 2D")
    if source_altitudes_arr.shape != sample_matrix.shape:
        raise ValueError("source_altitudes must match sample_matrix")
    band_distance_indices = [
        max(
            0,
            min(
                distances_m.size - 1,
                int(np.searchsorted(distances_m, float(distance_km) * 1000.0, side="right") - 1),
            ),
        )
        for _band_min_km, distance_km in band_ranges_km
    ]
    night_sample_matrix = _apply_night_light_sample_floor(
        sample_matrix,
        None,
        floor_value=0.0,
    )
    distance_attenuation = _night_light_distance_attenuation(distances_m)
    night_weighted_matrix = night_sample_matrix * distance_attenuation[None, :]
    raw_strengths_by_band: list[np.ndarray] = []
    target_altitudes = _target_altitude_bins()
    raw_fields_by_band: list[np.ndarray] = []
    band_start_index = 0
    for distance_index in band_distance_indices:
        band_end_index = int(distance_index)
        if band_end_index < band_start_index:
            raw_strengths_by_band.append(np.zeros_like(az_grid, dtype=np.float64))
            raw_fields_by_band.append(np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64))
            band_start_index = band_end_index + 1
            continue

        band_matrix_night = night_weighted_matrix[:, band_start_index : band_end_index + 1]
        band_altitudes = source_altitudes_arr[:, band_start_index : band_end_index + 1]
        if band_matrix_night.size == 0 or band_altitudes.size == 0:
            raw_strengths_by_band.append(np.zeros_like(az_grid, dtype=np.float64))
            raw_fields_by_band.append(np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64))
            band_start_index = band_end_index + 1
            continue

        band_strengths = np.zeros_like(az_grid, dtype=np.float64)
        band_field = np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64)
        for sample_index in range(band_start_index, band_end_index + 1):
            sample_matrix_night = night_weighted_matrix[:, sample_index : sample_index + 1]
            sample_altitudes = source_altitudes_arr[:, sample_index : sample_index + 1]
            source_azimuths, source_altitudes, source_strengths = _flatten_glow_source_matrix(
                sample_matrix_night,
                sample_altitudes,
                az_grid,
            )
            sample_strengths = _accumulate_local_glow_strengths(
                source_azimuths_deg=source_azimuths,
                source_altitudes_deg=source_altitudes,
                source_strengths=source_strengths,
                target_azimuths_deg=az_grid,
                target_altitudes_deg=horizon_alt_values,
                sigma_deg=NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG,
            )
            sample_field = _accumulate_local_glow_field(
                source_azimuths_deg=source_azimuths,
                source_altitudes_deg=source_altitudes,
                source_strengths=source_strengths,
                target_azimuths_deg=az_grid,
                target_altitudes_deg=target_altitudes,
                sigma_deg=NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG,
            )
            mask_distance_m = _inner_ring_mask_distance_m(distances_m, sample_index)
            sample_mask = _terrain_band_target_mask(
                az_grid=az_grid,
                target_altitudes=horizon_alt_values,
                band_distance_m=mask_distance_m,
                terrain_profile_altaz=terrain_profile_key,
                terrain_profile_distances_m=terrain_profile_distances_key,
                terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_key,
                terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_key,
            )
            sample_field_mask = _terrain_band_target_altaz_mask(
                az_grid=az_grid,
                target_altitudes=target_altitudes,
                band_distance_m=mask_distance_m,
                terrain_profile_altaz=terrain_profile_key,
                terrain_profile_distances_m=terrain_profile_distances_key,
                terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_key,
                terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_key,
            )
            band_strengths += np.where(sample_mask, np.clip(sample_strengths, 0.0, None), 0.0)
            band_field += np.where(sample_field_mask, np.clip(sample_field, 0.0, None), 0.0)
        raw_strengths_by_band.append(band_strengths)
        raw_fields_by_band.append(band_field)
        band_start_index = band_end_index + 1

    if raw_strengths_by_band:
        full_raw_strengths = np.sum(np.stack(raw_strengths_by_band, axis=0), axis=0)
    else:
        full_raw_strengths = np.zeros_like(az_grid, dtype=np.float64)
    if raw_fields_by_band:
        full_raw_field = np.sum(np.stack(raw_fields_by_band, axis=0), axis=0)
    else:
        full_raw_field = np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64)

    scale_source = full_raw_field[full_raw_field > 0.0]
    scale = float(np.max(scale_source)) if scale_source.size else float(np.percentile(full_raw_strengths, 95))
    if not math.isfinite(scale) or scale <= 0.0:
        return None
    log_scale = float(np.log1p(scale))
    if not math.isfinite(log_scale) or log_scale <= 0.0:
        return None
    full_strengths = np.clip(np.log1p(np.clip(full_raw_strengths, 0.0, None)) / log_scale, 0.0, 1.0)
    normalized_band_strengths = [
        np.clip(np.log1p(np.clip(raw_strengths, 0.0, None)) / log_scale, 0.0, 1.0)
        for raw_strengths in raw_strengths_by_band
    ]
    full_field = np.clip(np.log1p(np.clip(full_raw_field, 0.0, None)) / log_scale, 0.0, 1.0)
    if smooth_strengths:
        full_strengths = np.clip(_circular_smooth(full_strengths), 0.0, 1.0)
        normalized_band_strengths = [
            np.clip(_circular_smooth(strengths), 0.0, 1.0)
            for strengths in normalized_band_strengths
        ]
    if not np.any(full_field > 0.0) and (
        not np.any(full_strengths > 0.0) or not any(np.any(strengths > 0.0) for strengths in normalized_band_strengths)
    ):
        return None
    return NightLightGlowProfile(
        samples=tuple(
            NightLightGlowSample(
                azimuth_deg=float(az_grid[index]) % 360.0,
                horizon_alt_deg=float(horizon_alt_values[index]),
                strength=float(full_strengths[index]),
            )
            for index in range(az_grid.size)
        ),
        sun_alt_deg=0.0,
        altitude_bins_deg=tuple(float(value) for value in target_altitudes.tolist()),
        alpha_grid=tuple(tuple(float(value) for value in row.tolist()) for row in full_field),
    )


@functools.lru_cache(maxsize=64)
def _compute_night_light_base_profile(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    terrain_refraction_coefficient: float,
    terrain_context_key: tuple[
        tuple[tuple[float, float], ...],
        tuple[float, ...],
        tuple[tuple[tuple[float, float], ...], ...],
        tuple[tuple[float, ...], ...],
        tuple[int, int, int] | None,
    ],
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    distance_step_km: float = NIGHT_LIGHTS_DISTANCE_STEP_KM,
) -> NightLightGlowProfile | None:
    started_at = time.perf_counter()
    (
        terrain_profile_key,
        terrain_profile_distances_key,
        terrain_secondary_ridges_key,
        terrain_secondary_ridges_distances_key,
        _terrain_sample_grid_key,
    ) = terrain_context_key
    az_grid, horizon_alt_values = _build_azimuth_grid(terrain_profile_key)
    if az_grid.size == 0:
        return None
    band_ranges_km = _distance_band_ranges_km(max_distance_km)
    if not band_ranges_km:
        return None
    distances_m = np.arange(
        max(500.0, float(distance_step_km) * 1000.0),
        float(max_distance_km) * 1000.0 + 0.5,
        float(distance_step_km) * 1000.0,
        dtype=np.float64,
    )
    sample_altitudes = _surface_point_apparent_altitudes(
        distances_m,
        observer_height_m=float(observer_height_m),
        refraction_coefficient=float(terrain_refraction_coefficient),
    )
    tile_paths = _ensure_night_light_tiles(
        cache_root=cache_root,
        timeout_s=timeout_s,
        download_timeout_s=download_timeout_s,
    )
    sample_matrix = np.vstack(
        [
            _sample_ray_night_light_samples(
                tile_paths=tile_paths,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                azimuth_deg=float(az),
                distances_m=distances_m,
            )
            for az in az_grid.tolist()
        ]
    ).astype(np.float64, copy=False)
    profile = _build_night_light_glow_profile_from_samples(
        az_grid=az_grid,
        horizon_alt_values=horizon_alt_values,
        distances_m=distances_m,
        sample_matrix=sample_matrix,
        source_altitudes=sample_altitudes,
        terrain_profile_key=terrain_profile_key,
        terrain_profile_distances_key=terrain_profile_distances_key,
        terrain_secondary_ridges_key=terrain_secondary_ridges_key,
        terrain_secondary_ridges_distances_key=terrain_secondary_ridges_distances_key,
        max_distance_km=max_distance_km,
        smooth_strengths=True,
    )
    if profile is not None:
        logger.info(
            "Night light alpha grid computed: az=%d alt=%d elapsed=%.3fs",
            len(profile.samples),
            len(profile.alpha_grid),
            time.perf_counter() - started_at,
        )
    return profile


def _compute_night_light_base_profile_with_terrain_samples(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    terrain_refraction_coefficient: float,
    terrain_context_key: tuple[
        tuple[tuple[float, float], ...],
        tuple[float, ...],
        tuple[tuple[tuple[float, float], ...], ...],
        tuple[tuple[float, ...], ...],
        tuple[int, int, int] | None,
    ],
    terrain_sample_distances_m: Sequence[float] | np.ndarray,
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    distance_step_km: float = NIGHT_LIGHTS_DISTANCE_STEP_KM,
) -> NightLightGlowProfile | None:
    started_at = time.perf_counter()
    (
        terrain_profile_key,
        terrain_profile_distances_key,
        terrain_secondary_ridges_key,
        terrain_secondary_ridges_distances_key,
        _terrain_sample_grid_key,
    ) = terrain_context_key
    az_grid, horizon_alt_values = _build_azimuth_grid(terrain_profile_key)
    if az_grid.size == 0:
        return None
    band_ranges_km = _distance_band_ranges_km(max_distance_km)
    if not band_ranges_km:
        return None
    distances_m = np.arange(
        max(500.0, float(distance_step_km) * 1000.0),
        float(max_distance_km) * 1000.0 + 0.5,
        float(distance_step_km) * 1000.0,
        dtype=np.float64,
    )
    source_altitude_rows = _terrain_sample_source_altitude_rows(
        terrain_sample_distances_m=terrain_sample_distances_m,
        terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
        source_distances_m=distances_m,
        observer_height_m=float(observer_height_m),
        refraction_coefficient=float(terrain_refraction_coefficient),
    )
    if source_altitude_rows is None or source_altitude_rows.shape[0] != az_grid.size:
        source_altitude_rows = np.repeat(
            _surface_point_apparent_altitudes(
                distances_m,
                observer_height_m=float(observer_height_m),
                refraction_coefficient=float(terrain_refraction_coefficient),
            )[np.newaxis, :],
            az_grid.size,
            axis=0,
        )
    tile_paths = _ensure_night_light_tiles(
        cache_root=cache_root,
        timeout_s=timeout_s,
        download_timeout_s=download_timeout_s,
    )
    sample_matrix = np.vstack(
        [
            _sample_ray_night_light_samples(
                tile_paths=tile_paths,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                azimuth_deg=float(az),
                distances_m=distances_m,
            )
            for az in az_grid.tolist()
        ]
    ).astype(np.float64, copy=False)
    profile = _build_night_light_glow_profile_from_samples(
        az_grid=az_grid,
        horizon_alt_values=horizon_alt_values,
        distances_m=distances_m,
        sample_matrix=sample_matrix,
        source_altitudes=source_altitude_rows,
        terrain_profile_key=terrain_profile_key,
        terrain_profile_distances_key=terrain_profile_distances_key,
        terrain_secondary_ridges_key=terrain_secondary_ridges_key,
        terrain_secondary_ridges_distances_key=terrain_secondary_ridges_distances_key,
        max_distance_km=max_distance_km,
        smooth_strengths=True,
    )
    if profile is not None:
        logger.info(
            "Night light alpha grid computed: az=%d alt=%d elapsed=%.3fs",
            len(profile.samples),
            len(profile.alpha_grid),
            time.perf_counter() - started_at,
        )
    return profile


def compute_night_light_glow_profile(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float = 0.0,
    sun_alt_deg: float,
    terrain_profile_altaz: Sequence[tuple[float, float]] | None = None,
    terrain_profile_distances_m: Sequence[float] | None = None,
    terrain_secondary_ridges_altaz_layers: Sequence[Sequence[tuple[float, float]]] | None = None,
    terrain_secondary_ridges_distances_m_layers: Sequence[Sequence[float]] | None = None,
    terrain_sample_distances_m: Sequence[float] | np.ndarray | None = None,
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None = None,
    terrain_refraction_coefficient: float = 0.13,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    distance_step_km: float = NIGHT_LIGHTS_DISTANCE_STEP_KM,
) -> NightLightGlowProfile | None:
    if night_light_strength_factor(sun_alt_deg) <= 0.0:
        return None
    terrain_context_key = _terrain_context_key(
        terrain_profile_altaz=terrain_profile_altaz,
        terrain_profile_distances_m=terrain_profile_distances_m,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
        terrain_sample_azimuths_deg=None,
        terrain_sample_distances_m=terrain_sample_distances_m,
        terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
    )
    compute_with_samples = _compute_night_light_base_profile_with_terrain_samples
    if terrain_sample_distances_m is not None and terrain_sample_terrain_elevation_m is not None:
        return compute_with_samples(
            observer_lat_deg=float(observer_lat_deg),
            observer_lon_deg=float(observer_lon_deg),
            observer_height_m=float(observer_height_m),
            terrain_refraction_coefficient=float(terrain_refraction_coefficient),
            terrain_context_key=terrain_context_key,
            terrain_sample_distances_m=terrain_sample_distances_m,
            terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
            cache_root=cache_root,
            timeout_s=timeout_s,
            download_timeout_s=download_timeout_s,
            max_distance_km=float(max_distance_km),
            distance_step_km=float(distance_step_km),
        )
    return _compute_night_light_base_profile(
        observer_lat_deg=float(observer_lat_deg),
        observer_lon_deg=float(observer_lon_deg),
        observer_height_m=float(observer_height_m),
        terrain_refraction_coefficient=float(terrain_refraction_coefficient),
        terrain_context_key=terrain_context_key,
        cache_root=str(cache_root) if cache_root is not None else None,
        timeout_s=timeout_s,
        download_timeout_s=download_timeout_s,
        max_distance_km=float(max_distance_km),
        distance_step_km=float(distance_step_km),
    )


def is_night_light_enabled(sun_alt_deg: float) -> bool:
    return float(sun_alt_deg) < 0.0
