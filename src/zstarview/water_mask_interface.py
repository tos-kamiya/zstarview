from __future__ import annotations

import math
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from .clouddisc.types import DownloadCancelledError
from .location_resolver.place_projection import project_place_targets_to_altaz
from .terrain import WGS84_GEOD, build_ray_scan_grid
from .water_overlay import (
    DEFAULT_WATER_AZIMUTH_STEP_DEG,
    DEFAULT_WATER_SAMPLE_STEP_M,
    build_geometric_distance_samples,
    WaterOverlayPoint,
    bbox_from_point,
)


DEFAULT_WATER_TILES_ROOT_125M = Path(__file__).resolve().parent / "data" / "water_tiles_125m"
DEFAULT_WATER_TILES_ROOT_250M = Path(__file__).resolve().parent / "data" / "water_tiles_250m"
DEFAULT_WATER_TILES_ROOT_500M = Path(__file__).resolve().parent / "data" / "water_tiles_500m"
DEFAULT_WATER_TILES_ROOT = DEFAULT_WATER_TILES_ROOT_125M
DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM = 2.0
DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2 = 6.0
DEFAULT_WATER_REPRESENTATIVE_SWITCH_DISTANCE_KM = 12.0
DEFAULT_WATER_REPRESENTATIVE_SWITCH_DISTANCE_KM_2 = 18.0
DEFAULT_WATER_INTERFACE_BBOX_SCALE = 1.2
DEFAULT_WATER_INTERFACE_POINT_STRIDE = 1
DEFAULT_WATER_INTERFACE_SAMPLE_MIN_DISTANCE_M = 125.0


def _raise_if_abort_requested(abort_event: threading.Event | None) -> None:
    if abort_event is not None and abort_event.is_set():
        raise DownloadCancelledError("water surface sampling cancelled")


def _cooperative_yield(
    abort_event: threading.Event | None,
    *,
    interval: int,
    iteration_index: int,
) -> None:
    if interval <= 0 or iteration_index % interval != 0:
        return
    _raise_if_abort_requested(abort_event)
    time.sleep(0)


def _tile_paths(tile_root: Path) -> tuple[Path, ...]:
    if not tile_root.exists():
        return ()
    return tuple(sorted(path for path in tile_root.glob("tile_y*_x*") if path.suffix in {".tif", ".0", ".1"}))


def _bbox_intersection(
    a: tuple[float, float, float, float],
    b: tuple[float, float, float, float],
) -> tuple[float, float, float, float] | None:
    west = max(float(a[0]), float(b[0]))
    south = max(float(a[1]), float(b[1]))
    east = min(float(a[2]), float(b[2]))
    north = min(float(a[3]), float(b[3]))
    if west >= east or south >= north:
        return None
    return west, south, east, north


def _tile_key_from_path(tile_path: Path) -> tuple[int, int] | None:
    stem_parts = tile_path.stem.split("_")
    if len(stem_parts) != 3:
        return None
    row_part = stem_parts[1]
    col_part = stem_parts[2]
    if not row_part.startswith("y") or not col_part.startswith("x"):
        return None
    try:
        row = int(row_part[1:])
        col = int(col_part[1:])
    except ValueError:
        return None
    return row, col


def _tile_marker_value(tile_path: Path) -> bool | None:
    if tile_path.suffix == ".1":
        return True
    if tile_path.suffix == ".0":
        return False
    return None


def _collapse_tile_points_for_root(
    boundary_points: tuple[tuple[float, float], ...],
    *,
    tile_root: Path,
    tile_bounds: tuple[float, float, float, float],
) -> tuple[tuple[float, float], ...]:
    if tile_root not in {DEFAULT_WATER_TILES_ROOT_250M, DEFAULT_WATER_TILES_ROOT_500M}:
        return boundary_points
    if not boundary_points:
        return ()
    center_lon = 0.5 * (float(tile_bounds[0]) + float(tile_bounds[2]))
    center_lat = 0.5 * (float(tile_bounds[1]) + float(tile_bounds[3]))
    chosen_lonlat = min(
        boundary_points,
        key=lambda point: (
            (float(point[0]) - center_lon) ** 2 + (float(point[1]) - center_lat) ** 2,
            float(point[1]),
            float(point[0]),
        ),
    )
    return (chosen_lonlat,)


def _band_category_for_tile_root(tile_root: Path) -> str:
    if tile_root == DEFAULT_WATER_TILES_ROOT_125M:
        return "sea-125"
    if tile_root == DEFAULT_WATER_TILES_ROOT_250M:
        return "sea-250"
    if tile_root == DEFAULT_WATER_TILES_ROOT_500M:
        return "sea-500"
    return "sea"


def _water_band_specs(
    *,
    tile_root: Path | None,
    max_distance_km: float,
) -> list[tuple[Path, float, float | None, int]]:
    if tile_root is None or tile_root == DEFAULT_WATER_TILES_ROOT:
        return [
            (
                DEFAULT_WATER_TILES_ROOT_125M,
                0.0,
                min(float(max_distance_km), DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM),
                1,
            ),
            (
                DEFAULT_WATER_TILES_ROOT_250M,
                DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM,
                min(float(max_distance_km), DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2),
                1,
            ),
            (
                DEFAULT_WATER_TILES_ROOT_500M,
                DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2,
                min(float(max_distance_km), DEFAULT_WATER_REPRESENTATIVE_SWITCH_DISTANCE_KM),
                1,
            ),
            (
                DEFAULT_WATER_TILES_ROOT_500M,
                DEFAULT_WATER_REPRESENTATIVE_SWITCH_DISTANCE_KM,
                min(float(max_distance_km), DEFAULT_WATER_REPRESENTATIVE_SWITCH_DISTANCE_KM_2),
                2,
            ),
            (
                DEFAULT_WATER_TILES_ROOT_500M,
                DEFAULT_WATER_REPRESENTATIVE_SWITCH_DISTANCE_KM_2,
                float(max_distance_km),
                4,
            ),
        ]
    return [(tile_root, 0.0, float(max_distance_km), 1)]


@dataclass(frozen=True, slots=True)
class WaterSurfaceBandStats:
    band_name: str
    loaded_tile_count: int
    raw_point_count: int
    collapsed_point_count: int
    visible_point_count: int


def _tile_key_for_lonlat(lon_deg: float, lat_deg: float) -> tuple[int, int]:
    row = int(math.floor((90.0 - float(lat_deg)) / 45.0))
    col = int(math.floor((float(lon_deg) + 180.0) / 45.0))
    row = max(0, min(3, row))
    col = max(0, min(7, col))
    return row, col


def _sample_water_mask_for_lonlat_points(
    lonlat_points: list[tuple[float, float]],
    *,
    tile_root: Path | None = None,
) -> list[bool]:
    tile_root = DEFAULT_WATER_TILES_ROOT if tile_root is None else tile_root
    if not lonlat_points:
        return []

    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover - dependency mismatch
        raise RuntimeError("GeoTIFF support requires 'rasterio'.") from exc

    points_by_tile: dict[tuple[int, int], list[tuple[int, float, float]]] = {}
    for point_index, (lon_deg, lat_deg) in enumerate(lonlat_points):
        tile_key = _tile_key_for_lonlat(lon_deg, lat_deg)
        points_by_tile.setdefault(tile_key, []).append((point_index, float(lon_deg), float(lat_deg)))

    tile_paths = {
        key: path
        for path in _tile_paths(tile_root)
        if (key := _tile_key_from_path(path)) is not None
    }
    water_flags = [False] * len(lonlat_points)
    for tile_key, indexed_points in points_by_tile.items():
        tile_path = tile_paths.get(tile_key)
        if tile_path is None:
            continue
        marker_value = _tile_marker_value(tile_path)
        if marker_value is not None:
            for point_index, _lon_deg, _lat_deg in indexed_points:
                water_flags[point_index] = marker_value
            continue
        with rasterio.open(tile_path) as dataset:
            bounded_points: list[tuple[int, float, float]] = []
            for point_index, lon_deg, lat_deg in indexed_points:
                if (
                    dataset.bounds.left <= lon_deg <= dataset.bounds.right
                    and dataset.bounds.bottom <= lat_deg <= dataset.bounds.top
                ):
                    bounded_points.append((point_index, lon_deg, lat_deg))
            if not bounded_points:
                continue
            coords = [(lon_deg, lat_deg) for _point_index, lon_deg, lat_deg in bounded_points]
            try:
                samples = list(dataset.sample(coords))
            except Exception:
                continue
            for (point_index, _lon_deg, _lat_deg), sample in zip(bounded_points, samples):
                if sample.size > 0 and float(sample[0]) > 0.0:
                    water_flags[point_index] = True
    return water_flags


def _extract_lonlat_points_from_mask(
    mask: np.ndarray,
    *,
    transform: object,
    stride: int = DEFAULT_WATER_INTERFACE_POINT_STRIDE,
    representative_block_size: int = 1,
) -> tuple[tuple[float, float], ...]:
    water = np.asarray(mask, dtype=bool)
    if water.ndim != 2:
        raise ValueError("mask must be a 2D array")
    rows, cols = np.nonzero(water)
    if rows.size == 0:
        return ()
    stride = max(1, int(stride))
    if stride > 1:
        keep = ((rows + cols) % stride) == 0
        rows = rows[keep]
        cols = cols[keep]
        if rows.size == 0:
            return ()
    representative_block_size = max(1, int(representative_block_size))
    if representative_block_size > 1:
        selected_indices: dict[tuple[int, int], tuple[int, int, float]] = {}
        for row, col in zip(rows.tolist(), cols.tolist(), strict=False):
            block_row = int(row) // representative_block_size
            block_col = int(col) // representative_block_size
            center_row = block_row * representative_block_size + (representative_block_size - 1) / 2.0
            center_col = block_col * representative_block_size + (representative_block_size - 1) / 2.0
            score = (float(row) - center_row) ** 2 + (float(col) - center_col) ** 2
            key = (block_row, block_col)
            current = selected_indices.get(key)
            if current is None or score < current[2]:
                selected_indices[key] = (int(row), int(col), float(score))
        if not selected_indices:
            return ()
        ordered = [selected_indices[key] for key in sorted(selected_indices)]
        rows = np.asarray([row for row, _col, _score in ordered], dtype=np.int64)
        cols = np.asarray([col for _row, col, _score in ordered], dtype=np.int64)
    try:
        import rasterio.transform
    except ImportError as exc:  # pragma: no cover - dependency mismatch
        raise RuntimeError("GeoTIFF support requires 'rasterio'.") from exc

    lon_values, lat_values = rasterio.transform.xy(transform, rows, cols, offset="center")
    points: list[tuple[float, float]] = []
    for lon, lat in zip(lon_values, lat_values):
        try:
            points.append((float(lon), float(lat)))
        except (TypeError, ValueError):
            continue
    return tuple(points)


def _mask_candidate_count(
    mask: np.ndarray,
    *,
    stride: int = DEFAULT_WATER_INTERFACE_POINT_STRIDE,
) -> int:
    water = np.asarray(mask, dtype=bool)
    if water.ndim != 2:
        raise ValueError("mask must be a 2D array")
    rows, cols = np.nonzero(water)
    if rows.size == 0:
        return 0
    stride = max(1, int(stride))
    if stride > 1:
        keep = ((rows + cols) % stride) == 0
        rows = rows[keep]
        cols = cols[keep]
    return int(rows.size)


def _load_water_surface_interface_lonlat_points_for_root(
    *,
    center_lat_deg: float,
    center_lon_deg: float,
    radius_km: float,
    tile_root: Path,
    bbox_scale: float = DEFAULT_WATER_INTERFACE_BBOX_SCALE,
    stride: int = DEFAULT_WATER_INTERFACE_POINT_STRIDE,
    min_distance_km: float = 0.0,
    max_distance_km: float | None = None,
    representative_block_size: int = 1,
    abort_event: threading.Event | None = None,
) -> tuple[tuple[float, float], ...]:
    points, _stats = _load_water_surface_interface_lonlat_points_for_root_with_stats(
        center_lat_deg=center_lat_deg,
        center_lon_deg=center_lon_deg,
        radius_km=radius_km,
        tile_root=tile_root,
        bbox_scale=bbox_scale,
        stride=stride,
        min_distance_km=min_distance_km,
        max_distance_km=max_distance_km,
        representative_block_size=representative_block_size,
        abort_event=abort_event,
    )
    return points


def _load_water_surface_interface_lonlat_points_for_root_with_stats(
    *,
    center_lat_deg: float,
    center_lon_deg: float,
    radius_km: float,
    tile_root: Path,
    bbox_scale: float = DEFAULT_WATER_INTERFACE_BBOX_SCALE,
    stride: int = DEFAULT_WATER_INTERFACE_POINT_STRIDE,
    min_distance_km: float = 0.0,
    max_distance_km: float | None = None,
    representative_block_size: int = 1,
    abort_event: threading.Event | None = None,
) -> tuple[tuple[tuple[float, float], ...], WaterSurfaceBandStats]:
    if radius_km <= 0.0:
        raise ValueError("radius_km must be positive")
    if bbox_scale <= 0.0:
        raise ValueError("bbox_scale must be positive")
    if min_distance_km < 0.0:
        raise ValueError("min_distance_km must be non-negative")
    if max_distance_km is not None and max_distance_km <= min_distance_km:
        raise ValueError("max_distance_km must be greater than min_distance_km")
    representative_block_size = max(1, int(representative_block_size))

    request_bbox = bbox_from_point(float(center_lat_deg), float(center_lon_deg), float(radius_km))
    expanded_bbox = bbox_from_point(
        float(center_lat_deg),
        float(center_lon_deg),
        float(radius_km) * float(bbox_scale),
    )

    try:
        import rasterio
        from rasterio.windows import from_bounds
    except ImportError as exc:  # pragma: no cover - dependency mismatch
        raise RuntimeError("GeoTIFF support requires 'rasterio'.") from exc

    points: set[tuple[float, float]] = set()
    loaded_tile_count = 0
    raw_point_count = 0
    collapsed_point_count = 0
    for tile_index, tile_path in enumerate(_tile_paths(tile_root)):
        _cooperative_yield(abort_event, interval=4, iteration_index=tile_index)
        if _tile_marker_value(tile_path) is not None:
            continue
        with rasterio.open(tile_path) as dataset:
            _raise_if_abort_requested(abort_event)
            overlap = _bbox_intersection(
                expanded_bbox,
                (dataset.bounds.left, dataset.bounds.bottom, dataset.bounds.right, dataset.bounds.top),
            )
            if overlap is None:
                continue
            window = from_bounds(*overlap, transform=dataset.transform)
            try:
                mask = dataset.read(1, window=window, boundless=False)
            except Exception:
                continue
            if mask.size == 0:
                continue
            loaded_tile_count += 1
            raw_point_count += _mask_candidate_count(
                mask,
                stride=stride,
            )
            boundary_points = _extract_lonlat_points_from_mask(
                mask,
                transform=dataset.window_transform(window),
                stride=stride,
                representative_block_size=representative_block_size,
            )
            collapsed_point_count += len(boundary_points)
            for lon, lat in boundary_points:
                if request_bbox[0] <= lon <= request_bbox[2] and request_bbox[1] <= lat <= request_bbox[3]:
                    _, _, distance_m = WGS84_GEOD.inv(
                        float(center_lon_deg),
                        float(center_lat_deg),
                        float(lon),
                        float(lat),
                    )
                    distance_km = float(distance_m) / 1000.0
                    if distance_km < float(min_distance_km):
                        continue
                    if max_distance_km is not None and distance_km > float(max_distance_km):
                        continue
                    points.add((round(float(lon), 6), round(float(lat), 6)))

    band_category = _band_category_for_tile_root(tile_root)
    if band_category == "sea-125":
        band_name = "125m"
    elif band_category == "sea-250":
        band_name = "250m"
    elif band_category == "sea-500":
        band_name = "500m"
    else:
        band_name = band_category
    if representative_block_size > 1 and band_name == "500m":
        band_name = f"{band_name} x{representative_block_size}"
    visible_point_count = len(points)
    stats = WaterSurfaceBandStats(
        band_name=band_name,
        loaded_tile_count=loaded_tile_count,
        raw_point_count=raw_point_count,
        collapsed_point_count=collapsed_point_count,
        visible_point_count=visible_point_count,
    )
    return tuple(sorted(points)), stats


def _sample_water_surface_interface_ray_points_for_root_with_stats(
    *,
    center_lat_deg: float,
    center_lon_deg: float,
    observer_height_m: float,
    radius_km: float,
    tile_root: Path,
    target_ground_elevation_m_sampler: Callable[[float, float], float] | None = None,
    azimuth_step_deg: float = DEFAULT_WATER_AZIMUTH_STEP_DEG,
    sample_step_m: float = DEFAULT_WATER_SAMPLE_STEP_M,
    abort_event: threading.Event | None = None,
) -> tuple[tuple[WaterOverlayPoint, ...], WaterSurfaceBandStats]:
    if radius_km <= 0.0:
        raise ValueError("radius_km must be positive")
    if azimuth_step_deg <= 0.0:
        raise ValueError("azimuth_step_deg must be positive")
    if sample_step_m <= 0.0:
        raise ValueError("sample_step_m must be positive")

    distances_m = build_geometric_distance_samples(
        float(radius_km),
        float(sample_step_m),
    )
    distances_m = distances_m[distances_m >= float(DEFAULT_WATER_INTERFACE_SAMPLE_MIN_DISTANCE_M)]
    ray_scan = build_ray_scan_grid(
        geod=WGS84_GEOD,
        observer_latitude_deg=float(center_lat_deg),
        observer_longitude_deg=float(center_lon_deg),
        azimuth_step_deg=float(azimuth_step_deg),
        distance_samples_m=distances_m,
    )

    overlay_points: list[WaterOverlayPoint] = []
    band_category = _band_category_for_tile_root(tile_root)
    loaded_tile_count = len(_tile_paths(tile_root))
    raw_point_count = 0
    visible_point_count = 0

    for row_index, azimuth_deg in enumerate(ray_scan.azimuths_deg):
        _cooperative_yield(abort_event, interval=4, iteration_index=row_index)
        row_lonlat_points: list[tuple[float, float]] = []
        row_meta: list[tuple[int, float]] = []
        for col_index, distance_m in enumerate(ray_scan.distance_grid_m[row_index]):
            if col_index % 128 == 0:
                _raise_if_abort_requested(abort_event)
                time.sleep(0)
            lon = float(ray_scan.ray_lon_deg[row_index, col_index])
            lat = float(ray_scan.ray_lat_deg[row_index, col_index])
            row_lonlat_points.append((lon, lat))
            row_meta.append((int(col_index), float(distance_m)))
        raw_point_count += len(row_lonlat_points)
        if not row_lonlat_points:
            continue
        row_water_flags = _sample_water_mask_for_lonlat_points(
            row_lonlat_points,
            tile_root=tile_root,
        )
        water_meta: list[tuple[int, float, float, float, float]] = []
        for (col_index, distance_m), (lon_deg, lat_deg), is_water in zip(
            row_meta,
            row_lonlat_points,
            row_water_flags,
        ):
            if not is_water:
                continue
            target_height_m = 0.0
            if target_ground_elevation_m_sampler is not None:
                try:
                    target_height_m = float(
                        target_ground_elevation_m_sampler(float(lat_deg), float(lon_deg))
                    )
                except Exception:
                    target_height_m = 0.0
            water_meta.append(
                (
                    int(col_index),
                    float(distance_m),
                    float(lon_deg),
                    float(lat_deg),
                    float(target_height_m),
                )
            )
        if not water_meta:
            continue
        projections = project_place_targets_to_altaz(
            observer_latitude_deg=float(center_lat_deg),
            observer_longitude_deg=float(center_lon_deg),
            observer_height_m=float(observer_height_m),
            target_latitude_deg=[item[3] for item in water_meta],
            target_longitude_deg=[item[2] for item in water_meta],
            target_height_m=[item[4] for item in water_meta],
        )
        for (col_index, _distance_m, lon_deg, lat_deg, _target_height_m), projection in zip(
            water_meta,
            projections,
            strict=False,
        ):
            overlay_points.append(
                WaterOverlayPoint(
                    water_id="water-mask",
                    alt_deg=float(projection.alt_deg),
                    az_deg=float(projection.az_deg),
                    distance_km=float(projection.distance_km),
                    alpha_scale=1.0,
                    scan_azimuth_index=int(row_index),
                    scan_distance_index=int(col_index),
                    water_category=band_category,
                )
            )
            visible_point_count += 1

    if band_category == "sea-125":
        band_name = "125m"
    elif band_category == "sea-250":
        band_name = "250m"
    elif band_category == "sea-500":
        band_name = "500m"
    else:
        band_name = band_category

    stats = WaterSurfaceBandStats(
        band_name=band_name,
        loaded_tile_count=loaded_tile_count,
        raw_point_count=raw_point_count,
        collapsed_point_count=visible_point_count,
        visible_point_count=visible_point_count,
    )
    return tuple(overlay_points), stats


def _normalize_points_and_count(
    result: tuple[tuple[tuple[float, float], ...], int] | tuple[tuple[float, float], ...],
) -> tuple[tuple[tuple[float, float], ...], int]:
    if (
        isinstance(result, tuple)
        and len(result) == 2
        and isinstance(result[1], int)
    ):
        points, loaded_tile_count = result
        if isinstance(points, tuple):
            return points, int(loaded_tile_count)
    if isinstance(result, tuple):
        return result, 0
    return tuple(result), 0


def load_water_surface_interface_lonlat_points(
    *,
    center_lat_deg: float,
    center_lon_deg: float,
    radius_km: float,
    tile_root: Path | None = None,
    bbox_scale: float = DEFAULT_WATER_INTERFACE_BBOX_SCALE,
    stride: int = DEFAULT_WATER_INTERFACE_POINT_STRIDE,
) -> tuple[tuple[float, float], ...]:
    points: set[tuple[float, float]] = set()
    for band_root, min_distance_km, max_distance_km, representative_block_size in _water_band_specs(
        tile_root=tile_root,
        max_distance_km=float(radius_km),
    ):
        if max_distance_km is not None and max_distance_km <= min_distance_km:
            continue
        band_points = _load_water_surface_interface_lonlat_points_for_root(
            center_lat_deg=float(center_lat_deg),
            center_lon_deg=float(center_lon_deg),
            radius_km=float(max_distance_km),
            tile_root=band_root,
            bbox_scale=bbox_scale,
            stride=stride,
            min_distance_km=float(min_distance_km),
            max_distance_km=max_distance_km,
            representative_block_size=representative_block_size,
        )
        points.update(band_points)
    return tuple(sorted(points))


def sample_water_surface_interface_points_with_stats(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    max_distance_km: float,
    target_ground_elevation_m_sampler: Callable[[float, float], float] | None = None,
    tile_root: Path | None = None,
    bbox_scale: float = DEFAULT_WATER_INTERFACE_BBOX_SCALE,
    stride: int = DEFAULT_WATER_INTERFACE_POINT_STRIDE,
    abort_event: threading.Event | None = None,
) -> tuple[tuple[WaterOverlayPoint, ...], tuple[WaterSurfaceBandStats, ...]]:
    overlay_points: list[WaterOverlayPoint] = []
    band_stats: list[WaterSurfaceBandStats] = []
    for band_root, min_distance_km, max_distance_km_band, representative_block_size in _water_band_specs(
        tile_root=tile_root,
        max_distance_km=float(max_distance_km),
    ):
        if max_distance_km_band is not None and max_distance_km_band <= min_distance_km:
            continue
        band_points, stats = _sample_water_surface_interface_ray_points_for_root_with_stats(
            center_lat_deg=float(observer_lat_deg),
            center_lon_deg=float(observer_lon_deg),
            observer_height_m=float(observer_height_m),
            radius_km=float(max_distance_km_band),
            tile_root=band_root,
            target_ground_elevation_m_sampler=target_ground_elevation_m_sampler,
            abort_event=abort_event,
        )
        band_stats.append(stats)
        if not band_points:
            continue
        overlay_points.extend(band_points)
    return tuple(overlay_points), tuple(band_stats)


def sample_water_surface_interface_points(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    max_distance_km: float,
    tile_root: Path | None = None,
    bbox_scale: float = DEFAULT_WATER_INTERFACE_BBOX_SCALE,
    stride: int = DEFAULT_WATER_INTERFACE_POINT_STRIDE,
    abort_event: threading.Event | None = None,
) -> tuple[WaterOverlayPoint, ...]:
    points, _loaded_tile_counts = sample_water_surface_interface_points_with_stats(
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        observer_height_m=observer_height_m,
        max_distance_km=max_distance_km,
        tile_root=tile_root,
        bbox_scale=bbox_scale,
        stride=stride,
        abort_event=abort_event,
    )
    return points


def sample_water_surface_horizon_points(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    horizon_profile_altaz: list[tuple[float, float]],
    horizon_profile_distances_m: list[float],
    tile_root: Path | None = None,
) -> tuple[WaterOverlayPoint, ...]:
    if not horizon_profile_altaz or not horizon_profile_distances_m:
        return ()
    if len(horizon_profile_altaz) != len(horizon_profile_distances_m):
        return ()

    near_target_lonlat_points: list[tuple[float, float]] = []
    near_valid_indices: list[int] = []
    mid_target_lonlat_points: list[tuple[float, float]] = []
    mid_valid_indices: list[int] = []
    far_target_lonlat_points: list[tuple[float, float]] = []
    far_valid_indices: list[int] = []
    for index, ((alt_deg, az_deg), distance_m) in enumerate(
        zip(horizon_profile_altaz, horizon_profile_distances_m)
    ):
        if not (math.isfinite(float(alt_deg)) and math.isfinite(float(az_deg)) and math.isfinite(float(distance_m))):
            continue
        lon_deg, lat_deg, _ = WGS84_GEOD.fwd(
            float(observer_lon_deg),
            float(observer_lat_deg),
            float(az_deg),
            float(distance_m),
        )
        if float(distance_m) <= DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM * 1000.0:
            near_target_lonlat_points.append((float(lon_deg), float(lat_deg)))
            near_valid_indices.append(index)
        elif float(distance_m) <= DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2 * 1000.0:
            mid_target_lonlat_points.append((float(lon_deg), float(lat_deg)))
            mid_valid_indices.append(index)
        else:
            far_target_lonlat_points.append((float(lon_deg), float(lat_deg)))
            far_valid_indices.append(index)

    if not near_target_lonlat_points and not mid_target_lonlat_points and not far_target_lonlat_points:
        return ()

    if tile_root is None or tile_root == DEFAULT_WATER_TILES_ROOT:
        near_tile_root = DEFAULT_WATER_TILES_ROOT_125M
        mid_tile_root = DEFAULT_WATER_TILES_ROOT_250M
        far_tile_root = DEFAULT_WATER_TILES_ROOT_500M
    else:
        near_tile_root = tile_root
        mid_tile_root = tile_root
        far_tile_root = tile_root

    near_water_flags = (
        _sample_water_mask_for_lonlat_points(
            near_target_lonlat_points,
            tile_root=near_tile_root,
        )
        if near_target_lonlat_points
        else []
    )
    mid_water_flags = (
        _sample_water_mask_for_lonlat_points(
            mid_target_lonlat_points,
            tile_root=mid_tile_root,
        )
        if mid_target_lonlat_points
        else []
    )
    far_water_flags = (
        _sample_water_mask_for_lonlat_points(
            far_target_lonlat_points,
            tile_root=far_tile_root,
        )
        if far_target_lonlat_points
        else []
    )
    overlay_points: list[WaterOverlayPoint] = []
    for point_index, is_water in zip(near_valid_indices, near_water_flags):
        if not is_water:
            continue
        alt_deg, az_deg = horizon_profile_altaz[point_index]
        distance_m = horizon_profile_distances_m[point_index]
        overlay_points.append(
            WaterOverlayPoint(
                water_id="water-horizon",
                alt_deg=float(alt_deg),
                az_deg=float(az_deg),
                distance_km=float(distance_m) / 1000.0,
                water_category="sea",
            )
        )
    for point_index, is_water in zip(mid_valid_indices, mid_water_flags):
        if not is_water:
            continue
        alt_deg, az_deg = horizon_profile_altaz[point_index]
        distance_m = horizon_profile_distances_m[point_index]
        overlay_points.append(
            WaterOverlayPoint(
                water_id="water-horizon",
                alt_deg=float(alt_deg),
                az_deg=float(az_deg),
                distance_km=float(distance_m) / 1000.0,
                water_category="sea",
            )
        )
    for point_index, is_water in zip(far_valid_indices, far_water_flags):
        if not is_water:
            continue
        alt_deg, az_deg = horizon_profile_altaz[point_index]
        distance_m = horizon_profile_distances_m[point_index]
        overlay_points.append(
            WaterOverlayPoint(
                water_id="water-horizon",
                alt_deg=float(alt_deg),
                az_deg=float(az_deg),
                distance_km=float(distance_m) / 1000.0,
                water_category="sea",
            )
        )
    return tuple(overlay_points)


def sample_water_surface_horizon_layers_points(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    horizon_profile_altaz: list[tuple[float, float]],
    horizon_profile_distances_m: list[float],
    secondary_horizon_profile_altaz_layers: list[list[tuple[float, float]]] | None = None,
    secondary_horizon_profile_distances_m_layers: list[list[float]] | None = None,
    tile_root: Path | None = None,
) -> tuple[WaterOverlayPoint, ...]:
    points = list(
        sample_water_surface_horizon_points(
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
            horizon_profile_altaz=horizon_profile_altaz,
            horizon_profile_distances_m=horizon_profile_distances_m,
            tile_root=tile_root,
        )
    )
    if not secondary_horizon_profile_altaz_layers or not secondary_horizon_profile_distances_m_layers:
        return tuple(points)
    if len(secondary_horizon_profile_altaz_layers) != len(secondary_horizon_profile_distances_m_layers):
        return tuple(points)

    seen: set[tuple[float, float, float]] = {
        (round(float(point.alt_deg), 6), round(float(point.az_deg), 6), round(float(point.distance_km), 6))
        for point in points
    }
    for layer_altaz, layer_distances in zip(
        secondary_horizon_profile_altaz_layers,
        secondary_horizon_profile_distances_m_layers,
    ):
        layer_points = sample_water_surface_horizon_points(
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
            horizon_profile_altaz=list(layer_altaz),
            horizon_profile_distances_m=list(layer_distances),
            tile_root=tile_root,
        )
        for point in layer_points:
            key = (round(float(point.alt_deg), 6), round(float(point.az_deg), 6), round(float(point.distance_km), 6))
            if key in seen:
                continue
            seen.add(key)
            points.append(point)
    return tuple(points)
