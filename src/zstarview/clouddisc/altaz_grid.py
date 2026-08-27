"""Observer-centric (altitude, azimuth) cloud grid.

This module defines `CloudAltAzGrid`, the camera-independent intermediate
representation for cloud rendering, plus functions to build it from a
`CloudSourceData` and to read/write it from the on-disk cache.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .altaz_constants import (
    ALT_AZ_GEO_SAMPLE_EXTENT_DEG,
    ALT_AZ_GEO_SAMPLE_STEP_DEG,
    ALT_AZ_GRID_ALT_BINS,
    ALT_AZ_GRID_ALT_MAX_DEG,
    ALT_AZ_GRID_ALT_MIN_DEG,
    ALT_AZ_GRID_AZ_BINS,
    ALT_AZ_GRID_AZ_MAX_DEG,
    ALT_AZ_GRID_AZ_MIN_DEG,
    ALT_AZ_MISSING_NEIGHBORHOOD_CELLS,
)
from .altaz_projection import altaz_to_bin_indices, altaz_to_dir_ecef_array
from .projectors.az import geodetic_to_ecef
from .render.grayscale import _bt_to_weight, _suppress_low_cloud_weight
from .sampling.bt_sampler import build_bt_sampler
from .sampling.b13_b16 import (
    MAX_REDISTRIBUTION_STRENGTH,
    collocate_b13_b16,
    high_cloud_score,
)
from .sampling.estimate_bt_warm_cold import (
    estimate_bt_cold_hybrid,
    estimate_bt_warm_from_equator_band,
    estimate_bt_warm_hybrid,
)
from .types import CloudSourceData
from .workers.constants import DEFAULT_CLOUD_SHELLS_KM

logger = logging.getLogger(__name__)

ALT_AZ_CACHE_SUBDIR: str = "altaz_grids"
ALT_AZ_CACHE_META_SUFFIX: str = ".json"
ALT_AZ_CACHE_DATA_SUFFIX: str = ".npz"
ALT_AZ_GRID_ALGORITHM_B13_ONLY = "b13-only-v1"
ALT_AZ_GRID_ALGORITHM_B13_B16 = "b13-b16-v2"


@dataclass(frozen=True)
class CloudAltAzGrid:
    """Camera-independent cloud coverage as seen by an observer.

    Attributes:
        amount: (alt_bins, az_bins) float32 array, 0.0 (clear) to 1.0 (thick).
        missing_mask: (alt_bins, az_bins) uint8 array; 255 means no data.
        alt_min_deg / alt_max_deg: altitude range covered by the grid.
        az_min_deg / az_max_deg: azimuth range covered by the grid.
        observer_lat / observer_lon: observer used to build the grid.
        satellite / product / time_utc: source metadata.
        shells_km: cloud shell radii (Earth-center km) used for ingestion.
        source_key: key identifying the satellite/timeslot source.
        coverage_ratio: ratio of alt/az cells that have valid data.
        source_completeness_ratio: optional source tile coverage ratio.
        grid_resolution_deg: representative angular resolution of the grid.
        shell_amounts: optional per-shell amount fields in low-to-high order.
    """

    amount: np.ndarray
    missing_mask: np.ndarray
    alt_min_deg: float
    alt_max_deg: float
    az_min_deg: float
    az_max_deg: float
    observer_lat: float
    observer_lon: float
    satellite: str
    product: str
    time_utc: dt.datetime
    shells_km: tuple[float, ...]
    source_key: object
    coverage_ratio: float
    source_completeness_ratio: float | None = None
    grid_resolution_deg: float = 0.5
    shell_amounts: tuple[np.ndarray, ...] | None = None
    algorithm_version: str = ALT_AZ_GRID_ALGORITHM_B13_ONLY

    def __post_init__(self) -> None:
        if self.amount.shape != self.missing_mask.shape:
            raise ValueError(
                f"amount shape {self.amount.shape} != missing_mask shape {self.missing_mask.shape}"
            )
        if self.shell_amounts is not None:
            for shell_amount in self.shell_amounts:
                if shell_amount.shape != self.amount.shape:
                    raise ValueError(
                        "shell_amounts must have the same shape as amount"
                    )


def _smoothstep(edge0: float, edge1: float, x: float) -> float:
    """Return a smooth 0..1 ramp between thresholds."""
    denom = max(1e-6, float(edge1) - float(edge0))
    t = float(np.clip((float(x) - float(edge0)) / denom, 0.0, 1.0))
    return t * t * (3.0 - 2.0 * t)


def _blend_cloud_shell_weights(cloud_amount: float) -> tuple[float, ...]:
    """Return the 3-shell blend weights used by the legacy renderer."""
    low_amount = 0.25
    high_amount = 0.65
    low_weights = (0.0, 1.0, 0.0)
    default_weights = (0.20, 0.60, 0.20)
    t = _smoothstep(low_amount, high_amount, cloud_amount)
    blended = np.asarray(low_weights, dtype=np.float64) * (1.0 - t) + np.asarray(
        default_weights, dtype=np.float64
    ) * t
    return tuple(float(v) for v in blended)


def _estimate_scene_cloud_amount(bt: np.ndarray, bt_warm: float, bt_cold: float) -> float:
    """Estimate a single scene-wide cloudiness value from BT samples."""
    finite = bt[np.isfinite(bt)]
    if finite.size == 0:
        return 0.0
    weight = _bt_to_weight(finite, bt_warm, bt_cold)
    weight = _suppress_low_cloud_weight(weight)
    return float(np.mean(weight))


def _altaz_grid_centers(
    alt_bins: int,
    az_bins: int,
    *,
    alt_min_deg: float = ALT_AZ_GRID_ALT_MIN_DEG,
    alt_max_deg: float = ALT_AZ_GRID_ALT_MAX_DEG,
    az_min_deg: float = ALT_AZ_GRID_AZ_MIN_DEG,
    az_max_deg: float = ALT_AZ_GRID_AZ_MAX_DEG,
) -> tuple[np.ndarray, np.ndarray]:
    """Return 2D grids of bin centers for the observer-centric sky grid."""
    alt_edges = np.linspace(alt_min_deg, alt_max_deg, alt_bins + 1)
    az_edges = np.linspace(az_min_deg, az_max_deg, az_bins + 1)
    alt_centers = (alt_edges[:-1] + alt_edges[1:]) * 0.5
    az_centers = (az_edges[:-1] + az_edges[1:]) * 0.5
    return np.meshgrid(alt_centers, az_centers, indexing="ij")


def _intersect_altaz_rays_with_shell(
    alt_grid: np.ndarray,
    az_grid: np.ndarray,
    *,
    observer_lat: float,
    observer_lon: float,
    shell_km: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project alt/az rays onto a spherical shell and return lon/lat samples."""
    observer_pos = geodetic_to_ecef(observer_lat, observer_lon)
    ray_dirs = altaz_to_dir_ecef_array(alt_grid, az_grid, observer_lat, observer_lon)

    b_quad = 2.0 * np.sum(observer_pos * ray_dirs, axis=-1)
    c_quad = float(np.dot(observer_pos, observer_pos)) - float(shell_km * shell_km)
    discriminant = b_quad * b_quad - 4.0 * c_quad
    valid_intersection = discriminant >= 0.0

    sqrt_disc = np.sqrt(np.maximum(discriminant, 0.0))
    t1 = (-b_quad - sqrt_disc) / 2.0
    t2 = (-b_quad + sqrt_disc) / 2.0
    t = np.where(t1 > 1.0e-6, t1, np.where(t2 > 1.0e-6, t2, np.nan))
    valid = valid_intersection & np.isfinite(t)

    points = observer_pos + ray_dirs * t[..., None]
    x_int = points[..., 0]
    y_int = points[..., 1]
    z_int = points[..., 2]

    lon_grid = np.full(alt_grid.shape, np.nan, dtype=np.float32)
    lat_grid = np.full(alt_grid.shape, np.nan, dtype=np.float32)
    lon_grid[valid] = np.degrees(np.arctan2(y_int, x_int))[valid]
    hyp = np.hypot(x_int, y_int)
    lat_grid[valid] = np.degrees(np.arctan2(z_int, hyp))[valid]
    return lon_grid, lat_grid, valid


def build_altaz_grid(
    source: CloudSourceData,
    lat: float,
    lon: float,
    *,
    shells_km: Sequence[float] | None = None,
    az_bins: int = ALT_AZ_GRID_AZ_BINS,
    alt_bins: int = ALT_AZ_GRID_ALT_BINS,
    geo_sample_step_deg: float = ALT_AZ_GEO_SAMPLE_STEP_DEG,
    geo_sample_extent_deg: float = ALT_AZ_GEO_SAMPLE_EXTENT_DEG,
) -> CloudAltAzGrid:
    """Build a `CloudAltAzGrid` from satellite source data.

    Args:
        source: satellite brightness-temperature source data.
        lat: observer latitude in degrees.
        lon: observer longitude in degrees.
        shells_km: cloud shell radii from Earth's center in km. Defaults to
            the project-wide constant `CLOUD_SHELLS_KM`.
        az_bins / alt_bins: grid resolution.
        geo_sample_step_deg: geographic sampling step used to rasterize the
            source data before projecting to alt/az.
        geo_sample_extent_deg: half-width of the geographic sample box centred
            on the observer.

    Returns:
        A camera-independent `CloudAltAzGrid`.
    """
    if shells_km is None:
        shells_km = DEFAULT_CLOUD_SHELLS_KM

    sampler = source.sampler
    if sampler is None:
        sampler = build_bt_sampler(source.data_array)
    b16_source = source.auxiliary_bands.get("B16")
    b16_sampler = None
    if b16_source is not None:
        try:
            b16_sampler = b16_source.sampler or build_bt_sampler(b16_source.data_array)
        except Exception as exc:
            logger.info("B16 sampler unavailable; using B13-only grid: %s", exc)

    # 1. Build a local geographic sample grid for threshold estimation and
    #    scene cloudiness blending.
    extent = max(0.1, float(geo_sample_extent_deg))
    step = max(0.05, float(geo_sample_step_deg))
    local_lats = np.arange(lat - extent, lat + extent + step * 0.5, step, dtype=np.float64)
    local_lons = np.arange(lon - extent, lon + extent + step * 0.5, step, dtype=np.float64)
    local_lon_grid, local_lat_grid = np.meshgrid(local_lons, local_lats)

    # 2. Estimate warm/cold thresholds using the equatorial band and a small
    #    local view sample.  This mirrors the legacy renderer's approach.
    _, eq_samples = estimate_bt_warm_from_equator_band(
        source.data_array,
        lon_center_deg=lon,
        delta_lon=60.0,
        equator_lat=0.0,
        warm_p=97.0,
        half=5,
        equator_lat_half_band_deg=5.0,
    )

    # Use the central part of the geographic sample as the local view.
    central_half_extent = extent * 0.25
    central_mask = (
        (np.abs(local_lat_grid - lat) <= central_half_extent)
        & (np.abs(local_lon_grid - lon) <= central_half_extent)
    )
    central_bt = sampler(local_lon_grid[central_mask], local_lat_grid[central_mask])
    central_mask_1d = np.ones_like(central_bt, dtype=bool)

    bt_warm = estimate_bt_warm_hybrid(
        central_bt,
        central_mask_1d,
        eq_samples,
        fallback_bt_warm=310.0,
    )
    bt_cold = estimate_bt_cold_hybrid(
        central_bt,
        central_mask_1d,
        eq_samples,
        bt_warm,
    )

    # 3. Compute shell blend weights from the overall scene cloudiness.
    all_bt = sampler(local_lon_grid, local_lat_grid)
    cloud_amount = _estimate_scene_cloud_amount(all_bt, bt_warm, bt_cold)
    blend_weights = _blend_cloud_shell_weights(cloud_amount)

    logger.debug(
        "Alt/az grid thresholds: warm=%.2f cold=%.2f cloud_amount=%.3f weights=%s",
        bt_warm,
        bt_cold,
        cloud_amount,
        blend_weights,
    )

    # 4. Build the dense observer-centric alt/az grid and sample each shell.
    shell_amounts = [
        np.zeros((alt_bins, az_bins), dtype=np.float32)
        for _ in shells_km
    ]
    shell_delta_sum = np.zeros((alt_bins, az_bins), dtype=np.float64)
    shell_delta_count = np.zeros((alt_bins, az_bins), dtype=np.int32)
    sample_count = np.zeros((alt_bins, az_bins), dtype=np.int32)
    alt_grid, az_grid = _altaz_grid_centers(
        alt_bins,
        az_bins,
        alt_min_deg=ALT_AZ_GRID_ALT_MIN_DEG,
        alt_max_deg=ALT_AZ_GRID_ALT_MAX_DEG,
        az_min_deg=ALT_AZ_GRID_AZ_MIN_DEG,
        az_max_deg=ALT_AZ_GRID_AZ_MAX_DEG,
    )

    for shell_index, shell_km in enumerate(shells_km):
        lon_samples, lat_samples, valid_intersection = _intersect_altaz_rays_with_shell(
            alt_grid,
            az_grid,
            observer_lat=lat,
            observer_lon=lon,
            shell_km=float(shell_km),
        )
        if not np.any(valid_intersection):
            continue

        shell_bt = sampler(lon_samples[valid_intersection], lat_samples[valid_intersection])
        finite_bt = np.isfinite(shell_bt)
        if not np.any(finite_bt):
            continue

        shell_amount = _bt_to_weight(shell_bt, bt_warm, bt_cold)
        shell_amount = _suppress_low_cloud_weight(shell_amount)

        delta_bt = None
        if b16_sampler is not None:
            bt16 = b16_sampler(
                lon_samples[valid_intersection], lat_samples[valid_intersection]
            )
            pair = collocate_b13_b16(shell_bt, bt16)
            delta_bt = pair.delta_bt_k

        alt_idx, az_idx = altaz_to_bin_indices(
            alt_grid[valid_intersection],
            az_grid[valid_intersection],
            alt_bins=alt_bins,
            az_bins=az_bins,
        )

        # Accumulate using max to preserve thin high clouds visually.
        flat_idx = alt_idx * az_bins + az_idx
        flat_idx = flat_idx.astype(np.int64, copy=False)
        valid_amount = finite_bt & np.isfinite(shell_amount) & (shell_amount > 0.0)
        if not np.any(valid_amount):
            np.add.at(sample_count.ravel(), flat_idx.astype(np.int64, copy=False)[finite_bt], 1)
            continue

        flat_idx_finite = flat_idx.astype(np.int64, copy=False)[finite_bt]
        np.add.at(sample_count.ravel(), flat_idx_finite, 1)

        flat_idx_valid = flat_idx[valid_amount]
        amount_valid = shell_amount[valid_amount].astype(np.float32, copy=False)

        # Accumulate the maximum within this shell, while preserving shell
        # identity for per-height rendering.
        np.maximum.at(shell_amounts[shell_index].ravel(), flat_idx_valid, amount_valid)

        if delta_bt is not None:
            valid_delta = finite_bt & np.isfinite(delta_bt)
            if np.any(valid_delta):
                np.add.at(
                    shell_delta_sum.ravel(), flat_idx[valid_delta], delta_bt[valid_delta]
                )
                np.add.at(shell_delta_count.ravel(), flat_idx[valid_delta], 1)

    # Apply the B16 hint only after all shell samples are available.  The
    # middle shell provides the per-cell atmospheric hint; other shells retain
    # their B13 amount while their alpha is redistributed conservatively.
    for index, weight in enumerate(blend_weights):
        shell_amounts[index] *= float(weight)

    if b16_sampler is not None and len(shell_amounts) == 3:
        delta_grid = np.full((alt_bins, az_bins), np.nan, dtype=np.float32)
        np.divide(
            shell_delta_sum,
            shell_delta_count,
            out=delta_grid,
            where=shell_delta_count > 0,
        )
        middle_amount = shell_amounts[1]
        valid_hint = (shell_delta_count > 0) & (middle_amount > 0.0)
        scores = high_cloud_score(delta_grid)
        for alt_index, az_index in zip(*np.where(valid_hint), strict=True):
            score = float(scores[alt_index, az_index])
            signed = 2.0 * score - 1.0
            transfer = MAX_REDISTRIBUTION_STRENGTH * abs(signed)
            if transfer <= 0.0:
                continue
            layers = np.asarray(
                [shell_amount[alt_index, az_index] for shell_amount in shell_amounts],
                dtype=np.float32,
            )
            donor = 1
            recipient = 2 if signed > 0.0 else 0
            moved = min(float(layers[donor]) * transfer, float(layers[donor]))
            layers[donor] -= moved
            layers[recipient] += moved
            for index, shell_amount in enumerate(shell_amounts):
                shell_amount[alt_index, az_index] = layers[index]

    amount = np.maximum.reduce(shell_amounts)

    # 5. Build the missing mask.  A cell is missing only if it lies within the
    #    coverage neighbourhood of any valid sample but received no valid sample
    #    itself.  Cells far outside the satellite coverage are treated as clear.
    has_sample = sample_count > 0
    coverage_region = _dilate_bool(has_sample, ALT_AZ_MISSING_NEIGHBORHOOD_CELLS)
    missing_mask = (coverage_region & ~has_sample).astype(np.uint8) * 255

    covered_cells = int(np.count_nonzero(coverage_region))
    if covered_cells > 0:
        coverage_ratio = float(np.count_nonzero(has_sample)) / float(covered_cells)
    else:
        coverage_ratio = 1.0

    grid_resolution_deg = min(
        (ALT_AZ_GRID_AZ_MAX_DEG - ALT_AZ_GRID_AZ_MIN_DEG) / max(1, az_bins),
        (ALT_AZ_GRID_ALT_MAX_DEG - ALT_AZ_GRID_ALT_MIN_DEG) / max(1, alt_bins),
    )

    return CloudAltAzGrid(
        amount=amount,
        missing_mask=missing_mask,
        alt_min_deg=ALT_AZ_GRID_ALT_MIN_DEG,
        alt_max_deg=ALT_AZ_GRID_ALT_MAX_DEG,
        az_min_deg=ALT_AZ_GRID_AZ_MIN_DEG,
        az_max_deg=ALT_AZ_GRID_AZ_MAX_DEG,
        observer_lat=float(lat),
        observer_lon=float(lon),
        satellite=str(source.satellite),
        product=str(source.product),
        time_utc=source.time_utc,
        shells_km=tuple(float(v) for v in shells_km),
        source_key=source.source_key,
        coverage_ratio=coverage_ratio,
        source_completeness_ratio=source.source_completeness_ratio,
        grid_resolution_deg=float(grid_resolution_deg),
        shell_amounts=tuple(shell_amounts),
        algorithm_version=(
            ALT_AZ_GRID_ALGORITHM_B13_B16
            if b16_sampler is not None and len(shell_amounts) == 3
            else ALT_AZ_GRID_ALGORITHM_B13_ONLY
        ),
    )


def _dilate_bool(mask: np.ndarray, radius: int) -> np.ndarray:
    """Dilate a boolean mask by ``radius`` cells using a square structuring element."""
    if radius <= 0:
        return mask.astype(bool)
    # Simple square dilation via rolling maximum in 2D.
    work = mask.astype(np.uint8)
    for axis in (0, 1):
        shifted = np.stack(
            [np.roll(work, shift, axis=axis) for shift in range(-radius, radius + 1)],
            axis=-1,
        )
        work = np.any(shifted, axis=-1).astype(np.uint8)
    return work.astype(bool)


def altaz_grid_cache_key(
    observer_lat: float,
    observer_lon: float,
    *,
    satellite: str,
    product: str,
    time_utc: dt.datetime,
    source_key: object,
    shells_km: Sequence[float],
    grid_resolution_deg: float,
    algorithm_version: str = ALT_AZ_GRID_ALGORITHM_B13_ONLY,
) -> str:
    """Return a stable hash key for caching a grid on disk."""
    key_parts = [
        f"{float(observer_lat):.4f}",
        f"{float(observer_lon):.4f}",
        str(satellite),
        str(product),
        time_utc.strftime("%Y%m%d%H%M%S"),
        str(getattr(source_key, "satellite", satellite)),
        str(getattr(source_key, "timeslot_utc", time_utc)),
        ",".join(f"{float(v):.3f}" for v in shells_km),
        f"{float(grid_resolution_deg):.4f}",
        str(algorithm_version),
    ]
    raw = "|".join(key_parts)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:32]


def _grid_to_serializable_meta(grid: CloudAltAzGrid) -> dict:
    """Convert grid metadata to a JSON-serializable dict."""
    source_key = grid.source_key
    source_key_dict = {
        "satellite": getattr(source_key, "satellite", None),
        "provider": getattr(source_key, "provider", None),
        "timeslot_utc": getattr(source_key, "timeslot_utc", grid.time_utc).isoformat(),
        "sat_priority": list(getattr(source_key, "sat_priority", ("AUTO",))),
    }
    return {
        "observer_lat": grid.observer_lat,
        "observer_lon": grid.observer_lon,
        "satellite": grid.satellite,
        "product": grid.product,
        "time_utc": grid.time_utc.isoformat(),
        "shells_km": list(grid.shells_km),
        "source_key": source_key_dict,
        "coverage_ratio": grid.coverage_ratio,
        "source_completeness_ratio": grid.source_completeness_ratio,
        "grid_resolution_deg": grid.grid_resolution_deg,
        "alt_bins": grid.amount.shape[0],
        "az_bins": grid.amount.shape[1],
        "shell_amount_count": (
            len(grid.shell_amounts) if grid.shell_amounts is not None else 0
        ),
        "algorithm_version": grid.algorithm_version,
    }


def _meta_from_dict(
    meta: dict,
    amount: np.ndarray,
    missing_mask: np.ndarray,
    shell_amounts: tuple[np.ndarray, ...] | None = None,
) -> CloudAltAzGrid:
    """Reconstruct a CloudAltAzGrid from a loaded dict and arrays."""
    from .types import SourceKey

    source_key_raw = meta.get("source_key", {})
    timeslot_utc = dt.datetime.fromisoformat(source_key_raw.get("timeslot_utc", meta["time_utc"]))
    source_key = SourceKey(
        satellite=source_key_raw.get("satellite", meta["satellite"]),
        provider=source_key_raw.get("provider"),
        timeslot_utc=timeslot_utc,
        sat_priority=tuple(source_key_raw.get("sat_priority", ("AUTO",))),
    )
    time_utc = dt.datetime.fromisoformat(meta["time_utc"])
    return CloudAltAzGrid(
        amount=amount,
        missing_mask=missing_mask,
        alt_min_deg=meta.get("alt_min_deg", ALT_AZ_GRID_ALT_MIN_DEG),
        alt_max_deg=meta.get("alt_max_deg", ALT_AZ_GRID_ALT_MAX_DEG),
        az_min_deg=meta.get("az_min_deg", ALT_AZ_GRID_AZ_MIN_DEG),
        az_max_deg=meta.get("az_max_deg", ALT_AZ_GRID_AZ_MAX_DEG),
        observer_lat=meta["observer_lat"],
        observer_lon=meta["observer_lon"],
        satellite=meta["satellite"],
        product=meta["product"],
        time_utc=time_utc,
        shells_km=tuple(float(v) for v in meta["shells_km"]),
        source_key=source_key,
        coverage_ratio=meta.get("coverage_ratio", 1.0),
        source_completeness_ratio=meta.get("source_completeness_ratio"),
        grid_resolution_deg=meta.get("grid_resolution_deg", 0.5),
        shell_amounts=shell_amounts,
        algorithm_version=meta.get("algorithm_version", ALT_AZ_GRID_ALGORITHM_B13_ONLY),
    )


def save_altaz_grid(grid: CloudAltAzGrid, cache_root: Path) -> Path:
    """Persist a grid under ``cache_root / altaz_grids / <key>.npz``."""
    key = altaz_grid_cache_key(
        grid.observer_lat,
        grid.observer_lon,
        satellite=grid.satellite,
        product=grid.product,
        time_utc=grid.time_utc,
        source_key=grid.source_key,
        shells_km=grid.shells_km,
        grid_resolution_deg=grid.grid_resolution_deg,
        algorithm_version=grid.algorithm_version,
    )
    out_dir = cache_root / ALT_AZ_CACHE_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)

    data_path = out_dir / f"{key}{ALT_AZ_CACHE_DATA_SUFFIX}"
    meta_path = out_dir / f"{key}{ALT_AZ_CACHE_META_SUFFIX}"

    np.savez_compressed(
        data_path,
        amount=grid.amount,
        missing_mask=grid.missing_mask,
        **{  # type: ignore[arg-type]
            f"shell_amount_{index}": shell_amount
            for index, shell_amount in enumerate(grid.shell_amounts or ())
        },
    )
    meta_path.write_text(json.dumps(_grid_to_serializable_meta(grid), indent=2), encoding="utf-8")

    logger.info("Saved CloudAltAzGrid to %s", data_path)
    return data_path


def load_altaz_grid(cache_root: Path, key: str) -> CloudAltAzGrid | None:
    """Load a cached grid if it exists and is readable."""
    out_dir = cache_root / ALT_AZ_CACHE_SUBDIR
    data_path = out_dir / f"{key}{ALT_AZ_CACHE_DATA_SUFFIX}"
    meta_path = out_dir / f"{key}{ALT_AZ_CACHE_META_SUFFIX}"

    if not data_path.exists() or not meta_path.exists():
        return None

    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        with np.load(data_path, allow_pickle=False) as npz:
            amount = npz["amount"].astype(np.float32, copy=False)
            missing_mask = npz["missing_mask"].astype(np.uint8, copy=False)
            shell_count = int(meta.get("shell_amount_count", 0))
            shell_amounts = None
            if shell_count > 0:
                names = [f"shell_amount_{index}" for index in range(shell_count)]
                if all(name in npz.files for name in names):
                    shell_amounts = tuple(
                        npz[name].astype(np.float32, copy=False) for name in names
                    )
        return _meta_from_dict(meta, amount, missing_mask, shell_amounts)
    except Exception as e:
        logger.warning("Failed to load cached alt/az grid %s: %s", key, e)
        return None
