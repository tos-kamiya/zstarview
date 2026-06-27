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
NIGHT_LIGHTS_BAND_CENTER_OFFSET_DEG = 1.5
NIGHT_LIGHTS_BAND_HALF_WIDTH_DEG = 1.5
NIGHT_LIGHTS_MAX_ALPHA = 0.48
NIGHT_LIGHTS_GLOW_RGB = (244, 246, 248)
NIGHT_LIGHTS_RGB = NIGHT_LIGHTS_GLOW_RGB
NIGHT_LIGHTS_SUN_BLEND_START_ALT_DEG = -9.0
NIGHT_LIGHTS_SUN_BLEND_END_ALT_DEG = -4.0
NIGHT_LIGHTS_DISTANCE_BAND_EDGES_KM = DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM[1:]
NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG = 6.0
NIGHT_LIGHTS_DISTANCE_SIGMA_GAMMA = 0.65
NIGHT_LIGHTS_DISTANCE_SIGMA_REFERENCE_M = NIGHT_LIGHTS_DISTANCE_STEP_KM * 1000.0
NIGHT_LIGHTS_NEIGHBORHOOD_MAX_AZ_DELTA_DEG = 25.0
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
    edge_alpha_grid: tuple[tuple[float, ...], ...] = ()


@dataclass(frozen=True)
class NightLightTerrainContext:
    terrain_profile_key: tuple[tuple[float, float], ...] = ()
    terrain_profile_distances_key: tuple[float, ...] = ()
    terrain_secondary_ridges_key: tuple[tuple[tuple[float, float], ...], ...] = ()
    terrain_secondary_ridges_distances_key: tuple[tuple[float, ...], ...] = ()
    terrain_sample_grid_key: tuple[int, int, int] | None = None
    terrain_sample_distances_m: Sequence[float] | np.ndarray | None = None
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None = None
    terrain_sample_distances_key: tuple[float, ...] = ()
    terrain_sample_terrain_elevation_key: tuple[tuple[float, ...], ...] = ()
    terrain_sample_present: bool = False

    @classmethod
    def from_inputs(
        cls,
        *,
        terrain_profile_altaz: Sequence[tuple[float, float]] | None,
        terrain_profile_distances_m: Sequence[float] | None,
        terrain_secondary_ridges_altaz_layers: Sequence[Sequence[tuple[float, float]]] | None,
        terrain_secondary_ridges_distances_m_layers: Sequence[Sequence[float]] | None,
        terrain_sample_azimuths_deg: Sequence[float] | None = None,
        terrain_sample_distances_m: Sequence[float] | np.ndarray | None = None,
        terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None = None,
    ) -> "NightLightTerrainContext":
        return cls(
            terrain_profile_key=_terrain_profile_key(terrain_profile_altaz),
            terrain_profile_distances_key=_float_sequence_key(terrain_profile_distances_m),
            terrain_secondary_ridges_key=tuple(
                tuple(
                    (round(float(alt_deg), 3), round(float(az_deg) % 360.0, 3))
                    for alt_deg, az_deg in layer
                )
                for layer in terrain_secondary_ridges_altaz_layers or ()
            ),
            terrain_secondary_ridges_distances_key=_float_sequence_layers_key(
                terrain_secondary_ridges_distances_m_layers
            ),
            terrain_sample_grid_key=_terrain_sample_grid_key(
                terrain_sample_azimuths_deg=terrain_sample_azimuths_deg,
                terrain_sample_distances_m=terrain_sample_distances_m,
                terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
            ),
            terrain_sample_distances_m=terrain_sample_distances_m,
            terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
            terrain_sample_distances_key=_float_sequence_key(terrain_sample_distances_m),
            terrain_sample_terrain_elevation_key=_terrain_sample_terrain_elevation_key(
                terrain_sample_terrain_elevation_m
            ),
            terrain_sample_present=(
                terrain_sample_distances_m is not None
                and terrain_sample_terrain_elevation_m is not None
            ),
        )

    @property
    def has_sample_grid(self) -> bool:
        return self.terrain_sample_present


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


def _float_sequence_key(values: Sequence[float] | np.ndarray | None) -> tuple[float, ...]:
    if values is None:
        return ()
    return tuple(round(float(value), 3) for value in np.asarray(values, dtype=np.float64).reshape(-1))


def _float_sequence_layers_key(
    layers: Sequence[Sequence[float]] | np.ndarray | None,
) -> tuple[tuple[float, ...], ...]:
    if layers is None:
        return ()
    layer_arrays = np.asarray(layers, dtype=np.float64)
    if layer_arrays.ndim == 1:
        return (tuple(round(float(value), 3) for value in layer_arrays.reshape(-1)),)
    if layer_arrays.ndim != 2:
        raise ValueError("terrain distance layers must be 1D or 2D")
    return tuple(
        tuple(round(float(value), 3) for value in row.tolist())
        for row in layer_arrays
    )


def _terrain_sample_distances_key(
    terrain_sample_distances_m: Sequence[float] | np.ndarray | None,
) -> tuple[float, ...]:
    if terrain_sample_distances_m is None:
        return ()
    return tuple(round(float(distance_m), 3) for distance_m in np.asarray(terrain_sample_distances_m, dtype=np.float64).reshape(-1))


def _terrain_sample_terrain_elevation_key(
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None,
) -> tuple[tuple[float, ...], ...]:
    if terrain_sample_terrain_elevation_m is None:
        return ()
    elevations = np.asarray(terrain_sample_terrain_elevation_m, dtype=np.float64)
    if elevations.ndim == 1:
        elevations = elevations.reshape(1, -1)
    if elevations.ndim != 2:
        raise ValueError("terrain_sample_terrain_elevation_m must be a 2D array")
    return tuple(
        tuple(round(float(elevation_m), 3) for elevation_m in row.tolist())
        for row in elevations
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


def _night_light_distance_boost(
    distances_m: np.ndarray,
    *,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
) -> np.ndarray:
    """Return a linear distance boost that grows toward the far edge of the band."""
    distances = np.asarray(distances_m, dtype=np.float64)
    if distances.size == 0:
        return np.zeros(0, dtype=np.float64)
    max_distance_m = max(1.0, float(max_distance_km) * 1000.0)
    ramp = np.clip(distances, 0.0, max_distance_m) / max_distance_m
    return 1.0 + ramp


def _ridge_glow_distance_gain(
    *,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    target_strength_at_max_distance: float = 255.0,
) -> float:
    max_distance_m = max(1.0, float(max_distance_km) * 1000.0)
    boost_at_max_distance = float(
        _night_light_distance_boost(
            np.asarray([max_distance_m], dtype=np.float64),
            max_distance_km=max_distance_km,
        )[0]
    )
    if not math.isfinite(boost_at_max_distance) or boost_at_max_distance <= 0.0:
        return 0.0
    return max(0.0, float(target_strength_at_max_distance)) / boost_at_max_distance


def _night_light_distance_sigma_deg(
    distance_m: float,
    *,
    sigma_deg: float = NIGHT_LIGHTS_NEIGHBORHOOD_SIGMA_DEG,
    reference_distance_m: float = NIGHT_LIGHTS_DISTANCE_SIGMA_REFERENCE_M,
    gamma: float = NIGHT_LIGHTS_DISTANCE_SIGMA_GAMMA,
) -> float:
    distance = max(1.0, float(distance_m))
    reference = max(1.0, float(reference_distance_m))
    base_sigma = max(1.0e-6, float(sigma_deg))
    exponent = -max(0.0, float(gamma))
    scaled_sigma = base_sigma * ((distance / reference) ** exponent)
    return max(1.0e-6, float(scaled_sigma))


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


def _lookup_gaussian_weights(
    delta_deg: np.ndarray,
    *,
    sigma_deg: float,
    step_deg: float,
    max_delta_deg: float | None = None,
) -> np.ndarray:
    delta = np.abs(np.asarray(delta_deg, dtype=np.float64))
    lut = _gaussian_weight_lut(sigma_deg, step_deg)
    step = max(1.0e-6, float(step_deg))
    indices = np.clip(np.rint(delta / step).astype(np.int64), 0, lut.size - 1)
    weights = lut[indices]
    if max_delta_deg is None:
        return weights
    return np.where(delta <= float(max_delta_deg), weights, 0.0)


@functools.lru_cache(maxsize=16)
def _azimuth_weight_matrix(
    source_azimuths_key: tuple[float, ...],
    target_azimuths_key: tuple[float, ...],
    sigma_deg: float,
    max_delta_deg: float | None,
) -> np.ndarray:
    source_azimuths = np.asarray(source_azimuths_key, dtype=np.float64).reshape(-1)
    target_azimuths = np.asarray(target_azimuths_key, dtype=np.float64).reshape(-1)
    if source_azimuths.size == 0 or target_azimuths.size == 0:
        return np.zeros((source_azimuths.size, target_azimuths.size), dtype=np.float64)
    delta_az = _wrap_azimuth_delta_deg(source_azimuths[:, None], target_azimuths[None, :])
    return _lookup_gaussian_weights(
        delta_az,
        sigma_deg=sigma_deg,
        step_deg=NIGHT_LIGHTS_NEIGHBORHOOD_WEIGHT_STEP_DEG,
        max_delta_deg=max_delta_deg,
    )


def _accumulate_local_glow_strengths(
    *,
    source_azimuths_deg: np.ndarray,
    source_altitudes_deg: np.ndarray,
    source_strengths: np.ndarray,
    target_azimuths_deg: np.ndarray,
    target_altitudes_deg: np.ndarray,
    sigma_deg: float,
    azimuth_weights: np.ndarray | None = None,
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
    azimuth_weight_matrix = None
    if azimuth_weights is not None:
        azimuth_weight_matrix = np.asarray(azimuth_weights, dtype=np.float64)
        if azimuth_weight_matrix.shape != (source_azimuths.size, target_azimuths.size):
            raise ValueError("azimuth_weights must match source and target azimuth lengths")
    for start in range(0, source_strengths_arr.size, chunk):
        end = min(source_strengths_arr.size, start + chunk)
        source_az_chunk = source_azimuths[start:end][:, None]
        source_strength_chunk = source_strengths_arr[start:end][:, None]
        if azimuth_weight_matrix is None:
            delta_az = _wrap_azimuth_delta_deg(source_az_chunk, target_azimuths[None, :])
            az_weights = _lookup_gaussian_weights(
                delta_az,
                sigma_deg=sigma_deg,
                step_deg=NIGHT_LIGHTS_NEIGHBORHOOD_WEIGHT_STEP_DEG,
                max_delta_deg=NIGHT_LIGHTS_NEIGHBORHOOD_MAX_AZ_DELTA_DEG,
            )
        else:
            az_weights = azimuth_weight_matrix[start:end, :]
        source_alt_chunk = source_altitudes[start:end][:, None]
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
    azimuth_weights: np.ndarray | None = None,
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
    azimuth_weight_matrix = None
    if azimuth_weights is not None:
        azimuth_weight_matrix = np.asarray(azimuth_weights, dtype=np.float64)
        if azimuth_weight_matrix.shape != (source_azimuths.size, target_azimuths.size):
            raise ValueError("azimuth_weights must match source and target azimuth lengths")
    for start in range(0, source_strengths_arr.size, chunk):
        end = min(source_strengths_arr.size, start + chunk)
        source_az_chunk = source_azimuths[start:end]
        source_alt_chunk = source_altitudes[start:end]
        source_strength_chunk = np.clip(source_strengths_arr[start:end], 0.0, None)
        if not np.any(source_strength_chunk > 0.0):
            continue
        if azimuth_weight_matrix is None:
            delta_az = _wrap_azimuth_delta_deg(source_az_chunk[:, None], target_azimuths[None, :])
            az_weights = _lookup_gaussian_weights(
                delta_az,
                sigma_deg=sigma_deg,
                step_deg=NIGHT_LIGHTS_NEIGHBORHOOD_WEIGHT_STEP_DEG,
                max_delta_deg=NIGHT_LIGHTS_NEIGHBORHOOD_MAX_AZ_DELTA_DEG,
            )
        else:
            az_weights = azimuth_weight_matrix[start:end, :]
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
    return np.cumsum(samples)


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


def _terrain_sample_edge_strength_rows(
    *,
    terrain_sample_distances_m: Sequence[float] | Sequence[Sequence[float]] | np.ndarray | None,
    terrain_sample_terrain_elevation_m: Sequence[Sequence[float]] | np.ndarray | None,
    source_distances_m: np.ndarray,
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
    terrain_rel_elevation = np.ones_like(np.asarray(terrain_elevation, dtype=np.float64), dtype=np.float64)

    if terrain_distances is not None and np.array_equal(terrain_distances, source_distances):
        return terrain_rel_elevation
    rows = [
        np.interp(
            source_distances,
            terrain_distances if terrain_distances is not None else np.asarray(terrain_distances_arr[row_index], dtype=np.float64).reshape(-1),
            np.asarray(row, dtype=np.float64),
            left=float(row[0]),
            right=float(row[-1]),
        )
        for row_index, row in enumerate(np.asarray(terrain_rel_elevation, dtype=np.float64))
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


def _terrain_visibility_threshold_grid(
    *,
    az_grid: np.ndarray,
    distances_m: np.ndarray,
    terrain_profile_altaz: Sequence[tuple[float, float]] | None,
    terrain_profile_distances_m: Sequence[float] | None,
    terrain_secondary_ridges_altaz_layers: Sequence[Sequence[tuple[float, float]]] | None,
    terrain_secondary_ridges_distances_m_layers: Sequence[Sequence[float]] | None,
) -> np.ndarray:
    az_values = np.asarray(az_grid, dtype=np.float64).reshape(-1)
    distances = np.asarray(distances_m, dtype=np.float64).reshape(-1)
    if az_values.size == 0 or distances.size == 0:
        return np.zeros((az_values.size, distances.size), dtype=np.float64)
    threshold_grid = np.full((az_values.size, distances.size), -np.inf, dtype=np.float64)
    for az_index, az in enumerate(az_values.tolist()):
        threshold = _terrain_visibility_threshold_curve(
            azimuth_deg=float(az),
            distances_m=distances,
            terrain_profile_altaz=terrain_profile_altaz,
            terrain_profile_distances_m=terrain_profile_distances_m,
            terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
            terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
        )
        if threshold is None:
            continue
        threshold_grid[az_index, :] = np.asarray(threshold, dtype=np.float64).reshape(-1)
    return threshold_grid


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
        return np.zeros(0, dtype=np.float64)
    if az_values.size != target_altitudes_arr.size:
        raise ValueError("az_grid and target_altitudes must have the same length")
    threshold_grid = _terrain_visibility_threshold_grid(
        az_grid=az_values,
        distances_m=np.asarray([float(band_distance_m)], dtype=np.float64),
        terrain_profile_altaz=terrain_profile_altaz,
        terrain_profile_distances_m=terrain_profile_distances_m,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
    )
    threshold_column = threshold_grid[:, 0] if threshold_grid.size else np.zeros_like(az_values)
    mask = np.zeros(az_values.shape, dtype=np.float64)
    fade_width_deg = max(1.0e-6, 0.5 * float(NIGHT_LIGHTS_ALTITUDE_STEP_DEG))
    finite_thresholds = np.isfinite(threshold_column)
    mask[~finite_thresholds] = np.where(threshold_column[~finite_thresholds] < 0.0, 1.0, 0.0)
    mask[finite_thresholds] = np.clip(
        (
            target_altitudes_arr[finite_thresholds]
            - (threshold_column[finite_thresholds] - fade_width_deg)
        )
        / fade_width_deg,
        0.0,
        1.0,
    )
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
        return np.zeros((alt_values.size, az_values.size), dtype=np.float64)
    mask = np.zeros((alt_values.size, az_values.size), dtype=np.float64)
    fade_width_deg = max(1.0e-6, 0.5 * float(NIGHT_LIGHTS_ALTITUDE_STEP_DEG))
    threshold_grid = _terrain_visibility_threshold_grid(
        az_grid=az_values,
        distances_m=np.asarray([float(band_distance_m)], dtype=np.float64),
        terrain_profile_altaz=terrain_profile_altaz,
        terrain_profile_distances_m=terrain_profile_distances_m,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
    )
    threshold_column = threshold_grid[:, 0] if threshold_grid.size else np.zeros_like(az_values)
    finite_thresholds = np.isfinite(threshold_column)
    mask[:, ~finite_thresholds] = np.where(threshold_column[~finite_thresholds] < 0.0, 1.0, 0.0)
    mask[:, finite_thresholds] = np.clip(
        (
            alt_values[:, None]
            - (threshold_column[finite_thresholds][None, :] - fade_width_deg)
        )
        / fade_width_deg,
        0.0,
        1.0,
    )
    return mask

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
    if sun_alt >= NIGHT_LIGHTS_SUN_BLEND_END_ALT_DEG:
        return 0.0
    return 1.0 - _smoothstep(
        NIGHT_LIGHTS_SUN_BLEND_START_ALT_DEG,
        NIGHT_LIGHTS_SUN_BLEND_END_ALT_DEG,
        sun_alt,
    )


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


def _build_night_light_glow_fields_from_samples(
    *,
    az_grid: np.ndarray,
    horizon_alt_values: np.ndarray,
    distances_m: np.ndarray,
    source_matrix: np.ndarray,
    source_altitudes: np.ndarray,
    terrain_context: NightLightTerrainContext,
    max_distance_km: float,
    smooth_strengths: bool = True,
    terrain_visibility_threshold_grid: np.ndarray | None = None,
    azimuth_weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray] | None:
    band_ranges_km = _distance_band_ranges_km(max_distance_km)
    if not band_ranges_km:
        return None
    azimuth_values = tuple(float(value) for value in np.asarray(az_grid, dtype=np.float64).tolist())
    source_altitudes_arr = np.asarray(source_altitudes, dtype=np.float64)
    if source_altitudes_arr.ndim == 1:
        source_altitudes_arr = np.repeat(source_altitudes_arr[np.newaxis, :], az_grid.size, axis=0)
    elif source_altitudes_arr.ndim != 2:
        raise ValueError("source_altitudes must be 1D or 2D")
    source_matrix_arr = np.asarray(source_matrix, dtype=np.float64)
    if source_altitudes_arr.shape != source_matrix_arr.shape:
        raise ValueError("source_altitudes must match source_matrix")
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
    night_light_source_matrix = _apply_night_light_sample_floor(source_matrix_arr, None, floor_value=0.0)
    raw_strengths_by_band: list[np.ndarray] = []
    target_altitudes = _target_altitude_bins()
    target_altitudes_arr = np.asarray(target_altitudes, dtype=np.float64)
    fade_width_deg = max(1.0e-6, 0.5 * float(NIGHT_LIGHTS_ALTITUDE_STEP_DEG))
    raw_fields_by_band: list[np.ndarray] = []
    threshold_grid = (
        _terrain_visibility_threshold_grid(
            az_grid=az_grid,
            distances_m=distances_m,
            terrain_profile_altaz=terrain_context.terrain_profile_key,
            terrain_profile_distances_m=terrain_context.terrain_profile_distances_key,
            terrain_secondary_ridges_altaz_layers=terrain_context.terrain_secondary_ridges_key,
            terrain_secondary_ridges_distances_m_layers=terrain_context.terrain_secondary_ridges_distances_key,
        )
        if terrain_visibility_threshold_grid is None
        else np.asarray(terrain_visibility_threshold_grid, dtype=np.float64)
    )
    if threshold_grid.shape != (az_grid.size, distances_m.size):
        raise ValueError("terrain_visibility_threshold_grid must match az_grid and distances_m")
    azimuth_weight_matrix = (
        np.asarray(azimuth_weights, dtype=np.float64)
        if azimuth_weights is not None
        else None
    )
    if azimuth_weight_matrix is not None and azimuth_weight_matrix.shape != (az_grid.size, az_grid.size):
        raise ValueError("azimuth_weights must be a square matrix for az_grid")
    band_start_index = 0
    for distance_index in band_distance_indices:
        band_end_index = int(distance_index)
        if band_end_index < band_start_index:
            raw_strengths_by_band.append(np.zeros_like(az_grid, dtype=np.float64))
            raw_fields_by_band.append(np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64))
            band_start_index = band_end_index + 1
            continue

        band_source_matrix = night_light_source_matrix[:, band_start_index : band_end_index + 1]
        if band_source_matrix.size == 0:
            raw_strengths_by_band.append(np.zeros_like(az_grid, dtype=np.float64))
            raw_fields_by_band.append(np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64))
            band_start_index = band_end_index + 1
            continue

        band_strengths = np.zeros_like(az_grid, dtype=np.float64)
        band_field = np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64)
        for sample_index in range(band_start_index, band_end_index + 1):
            sigma_deg = _night_light_distance_sigma_deg(float(distances_m[sample_index]))
            sample_azimuth_weights = (
                _azimuth_weight_matrix(
                    azimuth_values,
                    azimuth_values,
                    float(sigma_deg),
                    float(NIGHT_LIGHTS_NEIGHBORHOOD_MAX_AZ_DELTA_DEG),
                )
                if azimuth_weight_matrix is None
                else azimuth_weight_matrix
            )
            sample_source_column = night_light_source_matrix[:, sample_index : sample_index + 1]
            sample_altitudes = source_altitudes_arr[:, sample_index : sample_index + 1]
            source_azimuths, source_altitudes, source_strengths = _flatten_glow_source_matrix(
                sample_source_column,
                sample_altitudes,
                az_grid,
            )
            sample_strengths = _accumulate_local_glow_strengths(
                source_azimuths_deg=source_azimuths,
                source_altitudes_deg=source_altitudes,
                source_strengths=source_strengths,
                target_azimuths_deg=az_grid,
                target_altitudes_deg=horizon_alt_values,
                sigma_deg=sigma_deg,
                azimuth_weights=sample_azimuth_weights,
            )
            sample_field = _accumulate_local_glow_field(
                source_azimuths_deg=source_azimuths,
                source_altitudes_deg=source_altitudes,
                source_strengths=source_strengths,
                target_azimuths_deg=az_grid,
                target_altitudes_deg=target_altitudes,
                sigma_deg=sigma_deg,
                azimuth_weights=sample_azimuth_weights,
            )
            threshold_column = threshold_grid[:, sample_index]
            finite_thresholds = np.isfinite(threshold_column)
            sample_field_mask = np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64)
            sample_field_mask[:, ~finite_thresholds] = np.where(
                threshold_column[~finite_thresholds] < 0.0,
                1.0,
                0.0,
            )
            sample_field_mask[:, finite_thresholds] = np.clip(
                (
                    target_altitudes_arr[:, None]
                    - (threshold_column[finite_thresholds][None, :] - fade_width_deg)
                )
                / fade_width_deg,
                0.0,
                1.0,
            )
            band_strengths += np.clip(sample_strengths, 0.0, None)
            band_field += np.clip(sample_field, 0.0, None) * np.asarray(sample_field_mask, dtype=np.float64)
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
    if scale_source.size:
        scale = float(np.max(scale_source))
    else:
        positive_strengths = full_raw_strengths[full_raw_strengths > 0.0]
        scale = float(np.max(positive_strengths)) if positive_strengths.size else 0.0
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
    if (
        not np.any(full_field > 0.0)
        and not np.any(full_strengths > 0.0)
        and not any(np.any(strengths > 0.0) for strengths in normalized_band_strengths)
    ):
        return None
    return (full_strengths, full_field)


def _build_night_light_glow_profile_from_samples(
    *,
    az_grid: np.ndarray,
    horizon_alt_values: np.ndarray,
    distances_m: np.ndarray,
    source_altitudes: np.ndarray,
    terrain_context: NightLightTerrainContext,
    max_distance_km: float,
    night_light_source_matrix: np.ndarray | None = None,
    ridge_glow_source_matrix: np.ndarray | None = None,
    smooth_strengths: bool = True,
) -> NightLightGlowProfile | None:
    terrain_visibility_threshold_grid = _terrain_visibility_threshold_grid(
        az_grid=az_grid,
        distances_m=distances_m,
        terrain_profile_altaz=terrain_context.terrain_profile_key,
        terrain_profile_distances_m=terrain_context.terrain_profile_distances_key,
        terrain_secondary_ridges_altaz_layers=terrain_context.terrain_secondary_ridges_key,
        terrain_secondary_ridges_distances_m_layers=terrain_context.terrain_secondary_ridges_distances_key,
    )
    target_altitudes = _target_altitude_bins()
    if night_light_source_matrix is None:
        base_strengths = np.zeros(az_grid.size, dtype=np.float64)
        base_field = np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64)
    else:
        base_fields = _build_night_light_glow_fields_from_samples(
            az_grid=az_grid,
            horizon_alt_values=horizon_alt_values,
            distances_m=distances_m,
            source_matrix=night_light_source_matrix,
            source_altitudes=source_altitudes,
            terrain_context=terrain_context,
            max_distance_km=max_distance_km,
            smooth_strengths=smooth_strengths,
            terrain_visibility_threshold_grid=terrain_visibility_threshold_grid,
        )
        if base_fields is None:
            base_strengths = np.zeros(az_grid.size, dtype=np.float64)
            base_field = np.zeros((target_altitudes.size, az_grid.size), dtype=np.float64)
        else:
            base_strengths, base_field = base_fields

    distance_boost = _night_light_distance_boost(distances_m, max_distance_km=max_distance_km)
    ridge_glow_distance_gain = _ridge_glow_distance_gain(max_distance_km=max_distance_km)
    if ridge_glow_source_matrix is None:
        edge_field = np.zeros_like(base_field, dtype=np.float64)
    else:
        edge_fields = _build_night_light_glow_fields_from_samples(
            az_grid=az_grid,
            horizon_alt_values=horizon_alt_values,
            distances_m=distances_m,
            source_matrix=(
                np.asarray(ridge_glow_source_matrix, dtype=np.float64)
                * distance_boost[None, :]
                * ridge_glow_distance_gain
            ),
            source_altitudes=source_altitudes,
            terrain_context=terrain_context,
            max_distance_km=max_distance_km,
            smooth_strengths=smooth_strengths,
            terrain_visibility_threshold_grid=terrain_visibility_threshold_grid,
        )
        edge_field = (
            np.clip(edge_fields[1], 0.0, 1.0)
            if edge_fields is not None
            else np.zeros_like(base_field, dtype=np.float64)
        )

    if (
        not np.any(base_field > 0.0)
        and not np.any(edge_field > 0.0)
        and not np.any(base_strengths > 0.0)
    ):
        return None

    return NightLightGlowProfile(
        samples=tuple(
            NightLightGlowSample(
                azimuth_deg=float(az_grid[index]) % 360.0,
                horizon_alt_deg=float(horizon_alt_values[index]),
                strength=float(base_strengths[index]),
            )
            for index in range(az_grid.size)
        ),
        sun_alt_deg=0.0,
        altitude_bins_deg=tuple(float(value) for value in _target_altitude_bins().tolist()),
        alpha_grid=tuple(tuple(float(value) for value in row.tolist()) for row in base_field),
        edge_alpha_grid=tuple(tuple(float(value) for value in row.tolist()) for row in edge_field),
    )


@functools.lru_cache(maxsize=64)
def _compute_night_light_base_profile(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    terrain_refraction_coefficient: float,
    terrain_context: NightLightTerrainContext,
    include_night_light_tiles: bool = True,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    distance_step_km: float = NIGHT_LIGHTS_DISTANCE_STEP_KM,
) -> NightLightGlowProfile | None:
    az_grid, horizon_alt_values = _build_azimuth_grid(terrain_context.terrain_profile_key)
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
    night_light_source_matrix = None
    if include_night_light_tiles:
        tile_paths = _ensure_night_light_tiles(
            cache_root=cache_root,
            timeout_s=timeout_s,
            download_timeout_s=download_timeout_s,
        )
        night_light_source_matrix = np.vstack(
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
    return _build_night_light_glow_profile_from_samples(
        az_grid=az_grid,
        horizon_alt_values=horizon_alt_values,
        distances_m=distances_m,
        night_light_source_matrix=night_light_source_matrix,
        source_altitudes=sample_altitudes,
        terrain_context=terrain_context,
        max_distance_km=max_distance_km,
        smooth_strengths=True,
    )


def _compute_night_light_base_profile_with_terrain_samples(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    terrain_refraction_coefficient: float,
    terrain_context: NightLightTerrainContext,
    include_night_light_tiles: bool = True,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    distance_step_km: float = NIGHT_LIGHTS_DISTANCE_STEP_KM,
) -> NightLightGlowProfile | None:
    az_grid, horizon_alt_values = _build_azimuth_grid(terrain_context.terrain_profile_key)
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
        terrain_sample_distances_m=terrain_context.terrain_sample_distances_m,
        terrain_sample_terrain_elevation_m=terrain_context.terrain_sample_terrain_elevation_m,
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
    edge_strength_rows = _terrain_sample_edge_strength_rows(
        terrain_sample_distances_m=terrain_context.terrain_sample_distances_m,
        terrain_sample_terrain_elevation_m=terrain_context.terrain_sample_terrain_elevation_m,
        source_distances_m=distances_m,
    )
    if edge_strength_rows is None or edge_strength_rows.shape[0] != az_grid.size:
        edge_strength_rows = np.zeros_like(source_altitude_rows, dtype=np.float64)
    night_light_source_matrix = None
    if include_night_light_tiles:
        tile_paths = _ensure_night_light_tiles(
            cache_root=cache_root,
            timeout_s=timeout_s,
            download_timeout_s=download_timeout_s,
        )
        night_light_source_matrix = np.vstack(
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
    return _build_night_light_glow_profile_from_samples(
        az_grid=az_grid,
        horizon_alt_values=horizon_alt_values,
        distances_m=distances_m,
        night_light_source_matrix=night_light_source_matrix,
        source_altitudes=source_altitude_rows,
        ridge_glow_source_matrix=edge_strength_rows,
        terrain_context=terrain_context,
        max_distance_km=max_distance_km,
        smooth_strengths=True,
    )


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
    include_night_light_tiles: bool = True,
    cache_root: str | os.PathLike[str] | None = None,
    timeout_s: float = 60.0,
    download_timeout_s: float = 300.0,
    max_distance_km: float = NIGHT_LIGHTS_MAX_DISTANCE_KM,
    distance_step_km: float = NIGHT_LIGHTS_DISTANCE_STEP_KM,
) -> NightLightGlowProfile | None:
    if night_light_strength_factor(sun_alt_deg) <= 0.0:
        return None
    terrain_context = NightLightTerrainContext.from_inputs(
        terrain_profile_altaz=terrain_profile_altaz,
        terrain_profile_distances_m=terrain_profile_distances_m,
        terrain_secondary_ridges_altaz_layers=terrain_secondary_ridges_altaz_layers,
        terrain_secondary_ridges_distances_m_layers=terrain_secondary_ridges_distances_m_layers,
        terrain_sample_azimuths_deg=None,
        terrain_sample_distances_m=terrain_sample_distances_m,
        terrain_sample_terrain_elevation_m=terrain_sample_terrain_elevation_m,
    )
    if terrain_context.has_sample_grid:
        return _compute_night_light_base_profile_with_terrain_samples(
            observer_lat_deg=float(observer_lat_deg),
            observer_lon_deg=float(observer_lon_deg),
            observer_height_m=float(observer_height_m),
            terrain_refraction_coefficient=float(terrain_refraction_coefficient),
            terrain_context=terrain_context,
            include_night_light_tiles=bool(include_night_light_tiles),
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
        terrain_context=terrain_context,
        include_night_light_tiles=bool(include_night_light_tiles),
        cache_root=str(cache_root) if cache_root is not None else None,
        timeout_s=timeout_s,
        download_timeout_s=download_timeout_s,
        max_distance_km=float(max_distance_km),
        distance_step_km=float(distance_step_km),
    )


def is_night_light_enabled(sun_alt_deg: float) -> bool:
    return float(sun_alt_deg) < NIGHT_LIGHTS_SUN_BLEND_END_ALT_DEG
