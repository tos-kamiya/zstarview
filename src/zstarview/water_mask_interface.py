from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .location_resolver.place_projection import project_place_target_to_altaz
from .terrain import WGS84_GEOD
from .water_overlay import WaterOverlayPoint, bbox_from_point


DEFAULT_WATER_TILES_ROOT_125M = Path(__file__).resolve().parent / "data" / "water_tiles_125m"
DEFAULT_WATER_TILES_ROOT_250M = Path(__file__).resolve().parent / "data" / "water_tiles_250m"
DEFAULT_WATER_TILES_ROOT_500M = Path(__file__).resolve().parent / "data" / "water_tiles_500m"
DEFAULT_WATER_TILES_ROOT = DEFAULT_WATER_TILES_ROOT_125M
DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM = 2.0
DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2 = 6.0
DEFAULT_WATER_INTERFACE_BBOX_SCALE = 1.2
DEFAULT_WATER_INTERFACE_POINT_STRIDE = 1


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
) -> tuple[tuple[tuple[float, float], ...], WaterSurfaceBandStats]:
    if radius_km <= 0.0:
        raise ValueError("radius_km must be positive")
    if bbox_scale <= 0.0:
        raise ValueError("bbox_scale must be positive")
    if min_distance_km < 0.0:
        raise ValueError("min_distance_km must be non-negative")
    if max_distance_km is not None and max_distance_km <= min_distance_km:
        raise ValueError("max_distance_km must be greater than min_distance_km")

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
    for tile_path in _tile_paths(tile_root):
        if _tile_marker_value(tile_path) is not None:
            continue
        with rasterio.open(tile_path) as dataset:
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
            boundary_points = _extract_lonlat_points_from_mask(
                mask,
                transform=dataset.window_transform(window),
                stride=stride,
            )
            raw_point_count += len(boundary_points)
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
    visible_point_count = len(points)
    stats = WaterSurfaceBandStats(
        band_name=band_name,
        loaded_tile_count=loaded_tile_count,
        raw_point_count=raw_point_count,
        collapsed_point_count=collapsed_point_count,
        visible_point_count=visible_point_count,
    )
    return tuple(sorted(points)), stats


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
    tile_root = DEFAULT_WATER_TILES_ROOT if tile_root is None else tile_root
    if tile_root == DEFAULT_WATER_TILES_ROOT:
        near_tile_root = DEFAULT_WATER_TILES_ROOT_125M
        mid_tile_root = DEFAULT_WATER_TILES_ROOT_250M
        far_tile_root = DEFAULT_WATER_TILES_ROOT_500M
    else:
        near_tile_root = tile_root
        mid_tile_root = tile_root
        far_tile_root = tile_root

    points: set[tuple[float, float]] = set()
    near_radius_km = min(float(radius_km), DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM)
    near_points, _near_loaded_tile_count = _normalize_points_and_count(
        _load_water_surface_interface_lonlat_points_for_root(
        center_lat_deg=float(center_lat_deg),
        center_lon_deg=float(center_lon_deg),
        radius_km=near_radius_km,
        tile_root=near_tile_root,
        bbox_scale=bbox_scale,
        stride=stride,
        min_distance_km=0.0,
        max_distance_km=DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM,
    )
    )
    points.update(near_points)
    mid_radius_km = min(float(radius_km), DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2)
    if float(radius_km) > DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM and mid_tile_root != near_tile_root:
        mid_points, _mid_loaded_tile_count = _normalize_points_and_count(
            _load_water_surface_interface_lonlat_points_for_root(
            center_lat_deg=float(center_lat_deg),
            center_lon_deg=float(center_lon_deg),
            radius_km=mid_radius_km,
            tile_root=mid_tile_root,
            bbox_scale=bbox_scale,
            stride=stride,
            min_distance_km=DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM,
            max_distance_km=DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2,
        )
        )
        points.update(mid_points)
    if float(radius_km) > DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2 and far_tile_root != near_tile_root:
        far_points, _far_loaded_tile_count = _normalize_points_and_count(
            _load_water_surface_interface_lonlat_points_for_root(
            center_lat_deg=float(center_lat_deg),
            center_lon_deg=float(center_lon_deg),
            radius_km=float(radius_km),
            tile_root=far_tile_root,
            bbox_scale=bbox_scale,
            stride=stride,
            min_distance_km=DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2,
            max_distance_km=float(radius_km),
        )
        )
        points.update(far_points)
    return tuple(sorted(points))


def sample_water_surface_interface_points_with_stats(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_height_m: float,
    max_distance_km: float,
    tile_root: Path | None = None,
    bbox_scale: float = DEFAULT_WATER_INTERFACE_BBOX_SCALE,
    stride: int = DEFAULT_WATER_INTERFACE_POINT_STRIDE,
) -> tuple[tuple[WaterOverlayPoint, ...], tuple[WaterSurfaceBandStats, ...]]:
    overlay_points: list[WaterOverlayPoint] = []
    band_stats: list[WaterSurfaceBandStats] = []
    if tile_root is None or tile_root == DEFAULT_WATER_TILES_ROOT:
        band_roots: list[tuple[Path, float, float | None]] = [
            (
                DEFAULT_WATER_TILES_ROOT_125M,
                0.0,
                min(float(max_distance_km), DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM),
            ),
        ]
        if float(max_distance_km) > DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM:
            band_roots.append(
                (
                    DEFAULT_WATER_TILES_ROOT_250M,
                    DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM,
                    min(float(max_distance_km), DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2),
                )
            )
        if float(max_distance_km) > DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2:
            band_roots.append(
                (
                    DEFAULT_WATER_TILES_ROOT_500M,
                    DEFAULT_WATER_TILE_SWITCH_DISTANCE_KM_2,
                    float(max_distance_km),
                )
            )
    else:
        band_roots = [(tile_root, 0.0, float(max_distance_km))]

    for band_root, min_distance_km, max_distance_km_band in band_roots:
        lonlat_points, stats = _load_water_surface_interface_lonlat_points_for_root_with_stats(
            center_lat_deg=float(observer_lat_deg),
            center_lon_deg=float(observer_lon_deg),
            radius_km=float(max_distance_km_band),
            tile_root=band_root,
            bbox_scale=bbox_scale,
            stride=stride,
            min_distance_km=float(min_distance_km),
            max_distance_km=max_distance_km_band,
        )
        band_stats.append(stats)
        if not lonlat_points:
            continue
        band_category = _band_category_for_tile_root(band_root)
        for lon_deg, lat_deg in lonlat_points:
            projection = project_place_target_to_altaz(
                observer_latitude_deg=float(observer_lat_deg),
                observer_longitude_deg=float(observer_lon_deg),
                observer_height_m=float(observer_height_m),
                target_latitude_deg=float(lat_deg),
                target_longitude_deg=float(lon_deg),
                target_height_m=0.0,
            )
            overlay_points.append(
                WaterOverlayPoint(
                    water_id="water-mask",
                    alt_deg=float(projection.alt_deg),
                    az_deg=float(projection.az_deg),
                    distance_km=float(projection.distance_km),
                    water_category=band_category,
                )
            )
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
) -> tuple[WaterOverlayPoint, ...]:
    points, _loaded_tile_counts = sample_water_surface_interface_points_with_stats(
        observer_lat_deg=observer_lat_deg,
        observer_lon_deg=observer_lon_deg,
        observer_height_m=observer_height_m,
        max_distance_km=max_distance_km,
        tile_root=tile_root,
        bbox_scale=bbox_scale,
        stride=stride,
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
