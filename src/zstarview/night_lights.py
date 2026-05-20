from __future__ import annotations

import functools
import json
import logging
import math
import os
import re
import tempfile
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

logger = logging.getLogger(__name__)

NIGHT_LIGHTS_DATASET_VERSION = "2016_grayscale_500m"
NIGHT_LIGHTS_PAGE_URL = "https://science.nasa.gov/earth/earth-observatory/earth-at-night/maps/"
NIGHT_LIGHTS_TILE_NAMES = ("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2")
NIGHT_LIGHTS_TILE_COUNT = len(NIGHT_LIGHTS_TILE_NAMES)
NIGHT_LIGHTS_TILE_WIDTH_DEG = 90.0
NIGHT_LIGHTS_TILE_HEIGHT_DEG = 90.0
NIGHT_LIGHTS_MAX_DISTANCE_KM = 128.0
NIGHT_LIGHTS_DISTANCE_STEP_KM = 0.5
NIGHT_LIGHTS_ATTENUATION_SCALE_M = 22_000.0
NIGHT_LIGHTS_AZIMUTH_SMOOTHING_WIDTH = 2
NIGHT_LIGHTS_BAND_CENTER_OFFSET_DEG = 1.5
NIGHT_LIGHTS_BAND_HALF_WIDTH_DEG = 3.0
NIGHT_LIGHTS_MAX_ALPHA = 0.48
NIGHT_LIGHTS_RGB = (240, 173, 122)
NIGHT_LIGHTS_SUN_BLEND_START_ALT_DEG = -6.0

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
class NightLightDistanceBandProfile:
    min_distance_km: float
    max_distance_km: float
    samples: tuple[NightLightGlowSample, ...]


@dataclass(frozen=True)
class NightLightGlowProfile:
    samples: tuple[NightLightGlowSample, ...]
    sun_alt_deg: float
    band_center_offset_deg: float = NIGHT_LIGHTS_BAND_CENTER_OFFSET_DEG
    band_half_width_deg: float = NIGHT_LIGHTS_BAND_HALF_WIDTH_DEG
    band_profiles: tuple[NightLightDistanceBandProfile, ...] = ()


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
            "User-Agent": "zstarview-night-lights/1.0",
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
        headers={"User-Agent": "zstarview-night-lights/1.0"},
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


def _terrain_profile_cache_key(
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
) -> tuple[tuple[float, float], ...]:
    if not terrain_profile_altaz:
        return ()
    return tuple(
        (round(float(alt_deg), 3), round(float(az_deg) % 360.0, 3))
        for alt_deg, az_deg in terrain_profile_altaz
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


def _sample_ray_brightness_curve(
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

    attenuation = np.exp(-distances_m / float(NIGHT_LIGHTS_ATTENUATION_SCALE_M)) / (
        np.maximum(distances_m / 1000.0, 1.0) + 1.0
    )
    return np.cumsum(samples * attenuation)


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


def _build_azimuth_grid(
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
) -> tuple[np.ndarray, np.ndarray]:
    if terrain_profile_altaz:
        az_values = np.asarray([float(az) % 360.0 for _, az in terrain_profile_altaz], dtype=np.float64)
        order = np.argsort(az_values)
        return az_values[order], np.zeros_like(az_values[order])
    az_values = np.linspace(0.0, 360.0, num=180, endpoint=False, dtype=np.float64)
    alt_values = np.zeros_like(az_values)
    return az_values, alt_values


def _circular_smooth(values: np.ndarray) -> np.ndarray:
    width = NIGHT_LIGHTS_AZIMUTH_SMOOTHING_WIDTH
    if values.size == 0 or width <= 0:
        return values
    kernel = np.asarray([1.0, 2.0, 3.0, 2.0, 1.0], dtype=np.float64)
    kernel /= float(np.sum(kernel))
    padded = np.concatenate([values[-width:], values, values[:width]])
    smoothed = np.convolve(padded, kernel, mode="same")
    return smoothed[width:-width]


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
        for edge in DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM
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


@functools.lru_cache(maxsize=64)
def _compute_night_light_base_profile(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    terrain_profile_key: tuple[tuple[float, float], ...],
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    distance_step_km: float = NIGHT_LIGHTS_DISTANCE_STEP_KM,
) -> NightLightGlowProfile | None:
    tile_paths = _ensure_night_light_tiles(
        cache_root=cache_root,
        timeout_s=timeout_s,
        download_timeout_s=download_timeout_s,
    )
    az_grid, horizon_alt_values = _build_azimuth_grid(terrain_profile_key)
    if az_grid.size == 0:
        return None
    band_ranges_km = _distance_band_ranges_km(max_distance_km)
    if not band_ranges_km:
        return None
    distances_m = np.arange(
        max(1000.0, float(distance_step_km) * 1000.0),
        float(max_distance_km) * 1000.0 + 0.5,
        float(distance_step_km) * 1000.0,
        dtype=np.float64,
    )
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
    previous_band_distance_indices = [
        -1,
        *band_distance_indices[:-1],
    ]
    full_raw_strengths = np.zeros(az_grid.size, dtype=np.float64)
    raw_strengths_by_band = [
        np.zeros(az_grid.size, dtype=np.float64)
        for _ in band_ranges_km
    ]
    for index, az in enumerate(az_grid.tolist()):
        brightness_curve = _sample_ray_brightness_curve(
            tile_paths=tile_paths,
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
            azimuth_deg=float(az),
            distances_m=distances_m,
        )
        if brightness_curve.size == 0:
            continue
        full_raw_strengths[index] = float(brightness_curve[-1])
        for band_index, (distance_index, previous_distance_index) in enumerate(
            zip(band_distance_indices, previous_band_distance_indices, strict=True)
        ):
            band_value = float(brightness_curve[distance_index])
            if previous_distance_index >= 0:
                band_value -= float(brightness_curve[previous_distance_index])
            raw_strengths_by_band[band_index][index] = max(0.0, band_value)

    full_raw_strengths = _circular_smooth(full_raw_strengths)
    raw_strengths_by_band = [
        _circular_smooth(raw_strengths)
        for raw_strengths in raw_strengths_by_band
    ]
    scale = float(np.percentile(full_raw_strengths, 95))
    if not math.isfinite(scale) or scale <= 0.0:
        return None
    full_strengths = np.clip(np.sqrt(np.clip(full_raw_strengths / scale, 0.0, None)) * 1.35, 0.0, 1.0)
    band_strengths = [
        np.clip(np.sqrt(np.clip(raw_strengths / scale, 0.0, None)) * 1.35, 0.0, 1.0)
        for raw_strengths in raw_strengths_by_band
    ]
    if not np.any(full_strengths > 0.0) or not any(np.any(strengths > 0.0) for strengths in band_strengths):
        return None
    band_profiles = tuple(
        NightLightDistanceBandProfile(
            min_distance_km=float(band_min_km),
            max_distance_km=float(distance_km),
            samples=tuple(
                NightLightGlowSample(
                    azimuth_deg=float(az_grid[index]) % 360.0,
                    horizon_alt_deg=float(horizon_alt_values[index]),
                    strength=float(strengths[index]),
                )
                for index in range(az_grid.size)
            ),
        )
        for (band_min_km, distance_km), strengths in zip(band_ranges_km, band_strengths, strict=True)
    )
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
        band_profiles=band_profiles,
    )


def compute_night_light_glow_profile(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    sun_alt_deg: float,
    terrain_profile_altaz: Sequence[tuple[float, float]] | None = None,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    distance_step_km: float = NIGHT_LIGHTS_DISTANCE_STEP_KM,
) -> NightLightGlowProfile | None:
    if night_light_strength_factor(sun_alt_deg) <= 0.0:
        return None
    return _compute_night_light_base_profile(
        observer_lat_deg=float(observer_lat_deg),
        observer_lon_deg=float(observer_lon_deg),
        terrain_profile_key=_terrain_profile_cache_key(terrain_profile_altaz),
        cache_root=cache_root,
        timeout_s=timeout_s,
        download_timeout_s=download_timeout_s,
        max_distance_km=float(max_distance_km),
        distance_step_km=float(distance_step_km),
    )


def is_night_light_enabled(sun_alt_deg: float) -> bool:
    return float(sun_alt_deg) < 0.0
