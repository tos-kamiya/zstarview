#!/usr/bin/env python3
"""Experimental water-run scan demo that reuses the DEM cache.

This is a standalone validation script. It reads the same Copernicus DEM cache
directory used by zstarview, scans outward from an observer location, and marks
scan samples as water when they fall inside a supplied water polygon JSON file.

The purpose is to test the idea of drawing water as scan-aligned blue runs
instead of a separate full-screen polygon layer.
"""

# ruff: noqa: E402

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.paths import COPERNICUS_DEM_CACHE_DIR
from zstarview.terrain import (
    EARTH_MEAN_RADIUS_M,
    GeoTiffDem,
    ObserverLocation,
    WGS84_GEOD,
    build_download_bbox,
    compute_horizon_layers,
    sample_ground_elevation,
)
from zstarview.terrain.dem import collect_copernicus_tile_keys
from zstarview.water_overlay import WaterPolygonFootprint


@dataclass(frozen=True)
class PreparedWaterPolygon:
    water_id: str
    kind: str
    outer_rings_xy: tuple[tuple[tuple[float, float], ...], ...]
    inner_rings_xy: tuple[tuple[tuple[float, float], ...], ...]
    bounds_xy: tuple[float, float, float, float]
    source: str
    tags: dict[str, str]


@dataclass(frozen=True)
class WaterRunSegment:
    azimuth_deg: float
    start_distance_m: float
    end_distance_m: float
    sample_count: int
    points_xy: tuple[tuple[float, float], ...]


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan DEM samples and color water points using an existing cache.",
    )
    parser.add_argument("--lat", type=float, required=True)
    parser.add_argument("--lon", type=float, required=True)
    parser.add_argument("--water-json", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path(COPERNICUS_DEM_CACHE_DIR),
        help="Copernicus DEM cache directory to reuse.",
    )
    parser.add_argument(
        "--fetch-dem",
        action="store_true",
        help="Allow downloading missing DEM tiles into the cache directory.",
    )
    parser.add_argument("--max-distance-km", type=float, default=32.0)
    parser.add_argument("--download-margin-km", type=float, default=5.0)
    parser.add_argument(
        "--sample-step-m",
        type=float,
        dest="sample_step_m_legacy",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--sample-step-near-m",
        type=float,
        default=20.0,
        help="Ray sample spacing in the near field, in meters (default: 20.0).",
    )
    parser.add_argument(
        "--sample-step-far-m",
        type=float,
        default=100.0,
        help="Ray sample spacing beyond the near field, in meters (default: 100.0).",
    )
    parser.add_argument(
        "--sample-step-transition-km",
        type=float,
        default=2.0,
        help="Distance from the observer where sampling switches from near to far spacing, in km (default: 2.0).",
    )
    parser.add_argument("--azimuth-step-deg", type=float, default=2.0)
    parser.add_argument("--observer-height-m", type=float, default=1.7)
    parser.add_argument(
        "--water-extent-km",
        type=float,
        default=1.0,
        help="Half-width of the square water overlay window centered on the observer in km (default: 1.0).",
    )
    parser.add_argument("--dem-resampling", choices=("bilinear", "nearest"), default="bilinear")
    parser.add_argument("--canvas-size", type=int, default=1200)
    parser.add_argument("--canvas-padding", type=float, default=80.0)
    args = parser.parse_args(argv)
    if getattr(args, "sample_step_m_legacy", None) is not None:
        args.sample_step_far_m = float(args.sample_step_m_legacy)
    delattr(args, "sample_step_m_legacy")
    return args


def _load_water_footprints(path: Path) -> tuple[WaterPolygonFootprint, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    footprints: list[WaterPolygonFootprint] = []
    for item in payload.get("water_polygons", []):
        if not isinstance(item, dict):
            continue
        footprints.append(
            WaterPolygonFootprint(
                water_id=str(item.get("osm_id", "")),
                kind=str(item.get("kind", "water_polygon")),
                outer_rings_lonlat=_load_rings(item.get("outer_rings")),
                inner_rings_lonlat=_load_rings(item.get("inner_rings")),
                source=str(item.get("source", "unknown")),
                tags=_load_tags(item.get("tags")),
            )
        )
    if not footprints:
        raise ValueError(f"No water_polygons were found in {path}")
    return tuple(footprints)


def _load_rings(value: object) -> tuple[tuple[tuple[float, float], ...], ...]:
    if not isinstance(value, list):
        return ()
    rings: list[tuple[tuple[float, float], ...]] = []
    for ring in value:
        if not isinstance(ring, list):
            continue
        coords: list[tuple[float, float]] = []
        for point in ring:
            if not isinstance(point, list) or len(point) != 2:
                continue
            try:
                lon = float(point[0])
                lat = float(point[1])
            except (TypeError, ValueError):
                continue
            coords.append((lon, lat))
        if len(coords) >= 4:
            rings.append(tuple(coords))
    return tuple(rings)


def _load_tags(value: object) -> dict[str, str]:
    if not isinstance(value, dict):
        return {}
    result: dict[str, str] = {}
    for key, item in value.items():
        if isinstance(key, str) and isinstance(item, str):
            result[key] = item
    return result


def _ring_body(points: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
    ring = list(points)
    if len(ring) >= 2 and ring[0] == ring[-1]:
        ring = ring[:-1]
    if len(ring) < 3:
        return []
    return ring


def _ring_bounds(points: Iterable[tuple[float, float]]) -> tuple[float, float, float, float]:
    min_x = float("inf")
    min_y = float("inf")
    max_x = float("-inf")
    max_y = float("-inf")
    for x, y in points:
        min_x = min(min_x, float(x))
        min_y = min(min_y, float(y))
        max_x = max(max_x, float(x))
        max_y = max(max_y, float(y))
    if not math.isfinite(min_x):
        return 0.0, 0.0, 0.0, 0.0
    return min_x, min_y, max_x, max_y


def _point_in_ring_xy(point: tuple[float, float], ring: Iterable[tuple[float, float]]) -> bool:
    body = _ring_body(ring)
    if len(body) < 3:
        return False
    x, y = point
    inside = False
    for index in range(len(body)):
        x0, y0 = body[index]
        x1, y1 = body[(index + 1) % len(body)]
        if ((y0 > y) != (y1 > y)) and (y1 != y0):
            x_cross = ((x1 - x0) * (y - y0) / (y1 - y0)) + x0
            if x < x_cross:
                inside = not inside
    return inside


def _footprint_contains_point_xy(
    point: tuple[float, float],
    outer_rings_xy: tuple[tuple[tuple[float, float], ...], ...],
    inner_rings_xy: tuple[tuple[tuple[float, float], ...], ...],
) -> bool:
    for outer_ring in outer_rings_xy:
        if not _point_in_ring_xy(point, outer_ring):
            continue
        if any(_point_in_ring_xy(point, inner_ring) for inner_ring in inner_rings_xy):
            return False
        return True
    return False


def _build_polygon_list(
    footprints: Iterable[WaterPolygonFootprint],
    *,
    observer: ObserverLocation,
    water_extent_km: float,
) -> tuple[PreparedWaterPolygon, ...]:
    extent_m = float(water_extent_km) * 1000.0
    if extent_m <= 0.0:
        raise ValueError("water_extent_km must be positive")
    polygons: list[PreparedWaterPolygon] = []
    for footprint in footprints:
        outer_rings_xy: list[tuple[tuple[float, float], ...]] = []
        for ring in footprint.outer_rings_lonlat:
            projected = tuple(
                _project_local_xy(observer.longitude_deg, observer.latitude_deg, float(lon), float(lat))
                for lon, lat in ring
            )
            if len(projected) >= 4:
                outer_rings_xy.append(projected)

        inner_rings_xy: list[tuple[tuple[float, float], ...]] = []
        for ring in footprint.inner_rings_lonlat:
            projected = tuple(
                _project_local_xy(observer.longitude_deg, observer.latitude_deg, float(lon), float(lat))
                for lon, lat in ring
            )
            if len(projected) >= 4:
                inner_rings_xy.append(projected)

        if not outer_rings_xy:
            continue

        bounds_x = float("inf")
        bounds_y = float("inf")
        bounds_max_x = float("-inf")
        bounds_max_y = float("-inf")
        for ring in list(outer_rings_xy) + list(inner_rings_xy):
            min_x, min_y, max_x, max_y = _ring_bounds(ring)
            bounds_x = min(bounds_x, min_x)
            bounds_y = min(bounds_y, min_y)
            bounds_max_x = max(bounds_max_x, max_x)
            bounds_max_y = max(bounds_max_y, max_y)

        if bounds_x > extent_m or bounds_y > extent_m or bounds_max_x < -extent_m or bounds_max_y < -extent_m:
            continue

        polygons.append(
            PreparedWaterPolygon(
                water_id=str(footprint.water_id),
                kind=str(footprint.kind),
                outer_rings_xy=tuple(outer_rings_xy),
                inner_rings_xy=tuple(inner_rings_xy),
                bounds_xy=(bounds_x, bounds_y, bounds_max_x, bounds_max_y),
                source=str(footprint.source),
                tags=dict(footprint.tags),
            )
        )

    return tuple(polygons)


def _resolve_cached_dem(
    *,
    lat_deg: float,
    lon_deg: float,
    max_distance_km: float,
    margin_km: float,
    cache_dir: Path,
) -> tuple[Path, ...]:
    bbox = build_download_bbox(
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        radius_km=max_distance_km + margin_km,
    )
    tile_keys = collect_copernicus_tile_keys(bbox)
    paths = tuple(cache_dir / key for key in tile_keys)
    missing = [path for path in paths if not path.exists()]
    if missing:
        preview = "\n".join(str(path) for path in missing[:8])
        raise RuntimeError(
            "Missing Copernicus DEM cache tiles for this scan.\n"
            f"cache_dir={cache_dir}\n"
            f"missing_count={len(missing)}\n"
            f"first_missing:\n{preview}\n"
            "Run zstarview once for this area, or pass --fetch-dem to download tiles."
        )
    return paths


def _project_local_xy(
    observer_lon_deg: float,
    observer_lat_deg: float,
    lon_deg: float,
    lat_deg: float,
) -> tuple[float, float]:
    azimuth_deg, _, distance_m = WGS84_GEOD.inv(
        observer_lon_deg,
        observer_lat_deg,
        lon_deg,
        lat_deg,
    )
    azimuth_rad = math.radians(float(azimuth_deg))
    x_m = math.sin(azimuth_rad) * float(distance_m)
    y_m = math.cos(azimuth_rad) * float(distance_m)
    return x_m, y_m


def _build_distance_samples(
    *,
    max_distance_km: float,
    near_step_m: float,
    far_step_m: float,
    transition_km: float,
) -> np.ndarray:
    if max_distance_km <= 0.0:
        raise ValueError("max_distance_km must be positive")
    if near_step_m <= 0.0:
        raise ValueError("sample-step-near-m must be positive")
    if far_step_m <= 0.0:
        raise ValueError("sample-step-far-m must be positive")
    if transition_km < 0.0:
        raise ValueError("sample-step-transition-km must be non-negative")

    max_distance_m = float(max_distance_km) * 1000.0
    transition_m = min(float(transition_km) * 1000.0, max_distance_m)

    near_samples = np.arange(0.0, transition_m + float(near_step_m) * 0.5, float(near_step_m), dtype=np.float64)
    if near_samples.size == 0 or near_samples[-1] < transition_m:
        near_samples = np.append(near_samples, transition_m)

    far_start_m = transition_m + float(far_step_m)
    if far_start_m > max_distance_m:
        far_samples = np.empty(0, dtype=np.float64)
    else:
        far_samples = np.arange(far_start_m, max_distance_m + float(far_step_m) * 0.5, float(far_step_m), dtype=np.float64)

    samples = np.concatenate([near_samples, far_samples])
    if samples.size == 0 or samples[0] != 0.0:
        samples = np.insert(samples, 0, 0.0)
    if samples[-1] < max_distance_m:
        samples = np.append(samples, max_distance_m)
    return np.unique(samples)


def _sample_rays(
    *,
    observer: ObserverLocation,
    dem_grid,
    water_polygons: tuple[PreparedWaterPolygon, ...],
    max_distance_km: float,
    sample_step_near_m: float,
    sample_step_far_m: float,
    sample_step_transition_km: float,
    azimuth_step_deg: float,
    dem_resampling: str,
) -> tuple[list[WaterRunSegment], dict[str, object]]:
    azimuths = np.arange(0.0, 360.0, float(azimuth_step_deg), dtype=np.float64)
    distances_m = _build_distance_samples(
        max_distance_km=float(max_distance_km),
        near_step_m=float(sample_step_near_m),
        far_step_m=float(sample_step_far_m),
        transition_km=float(sample_step_transition_km),
    )
    az_grid_deg, distance_grid_m = np.meshgrid(azimuths, distances_m, indexing="ij")
    flat_count = az_grid_deg.size

    lon_flat, lat_flat, _ = WGS84_GEOD.fwd(
        np.full(flat_count, observer.longitude_deg, dtype=np.float64),
        np.full(flat_count, observer.latitude_deg, dtype=np.float64),
        az_grid_deg.ravel(),
        distance_grid_m.ravel(),
    )
    lon_grid_deg = np.asarray(lon_flat, dtype=np.float64).reshape(az_grid_deg.shape)
    lat_grid_deg = np.asarray(lat_flat, dtype=np.float64).reshape(az_grid_deg.shape)
    water_mask = np.zeros_like(lon_grid_deg, dtype=bool)
    for row_index in range(lon_grid_deg.shape[0]):
        for col_index in range(lon_grid_deg.shape[1]):
            x_m, y_m = _project_local_xy(
                observer.longitude_deg,
                observer.latitude_deg,
                float(lon_grid_deg[row_index, col_index]),
                float(lat_grid_deg[row_index, col_index]),
            )
            point = (x_m, y_m)
            water_mask[row_index, col_index] = any(
                _footprint_contains_point_xy(point, prepared.outer_rings_xy, prepared.inner_rings_xy)
                for prepared in water_polygons
            )

    runs: list[WaterRunSegment] = []
    for row_index, azimuth_deg in enumerate(azimuths):
        start_index: int | None = None
        for col_index, is_water in enumerate(water_mask[row_index]):
            if is_water and start_index is None:
                start_index = col_index
            elif not is_water and start_index is not None:
                runs.append(
                    _build_run_segment(
                        observer=observer,
                        azimuth_deg=float(azimuth_deg),
                        start_index=start_index,
                        end_index=col_index - 1,
                        distances_m=distances_m,
                        lon_grid_deg=lon_grid_deg[row_index],
                        lat_grid_deg=lat_grid_deg[row_index],
                    )
                )
                start_index = None
        if start_index is not None:
            runs.append(
                _build_run_segment(
                    observer=observer,
                    azimuth_deg=float(azimuth_deg),
                    start_index=start_index,
                    end_index=len(distances_m) - 1,
                    distances_m=distances_m,
                    lon_grid_deg=lon_grid_deg[row_index],
                    lat_grid_deg=lat_grid_deg[row_index],
                )
            )

    sample_total = int(water_mask.size)
    water_total = int(np.count_nonzero(water_mask))
    flat_profile = [
        {
            "azimuth_deg": float(point.azimuth_deg),
            "altitude_deg": float(point.altitude_deg),
            "distance_m": float(point.distance_m),
            "latitude_deg": float(point.latitude_deg),
            "longitude_deg": float(point.longitude_deg),
            "terrain_elevation_m": float(point.terrain_elevation_m),
        }
        for point in compute_horizon_layers(
            dem_grid=dem_grid,
            geod=WGS84_GEOD,
            observer=observer,
            azimuth_step_deg=float(azimuth_step_deg),
            distance_samples_m=distances_m,
            dem_resampling=dem_resampling,
            earth_radius_m=EARTH_MEAN_RADIUS_M,
            refraction_coefficient=0.13,
        ).main_profile
    ]
    horizon_peak = max(flat_profile, key=lambda item: float(item["altitude_deg"]))
    summary = {
        "sample_total": sample_total,
        "water_total": water_total,
        "water_ratio": float(water_total) / float(sample_total) if sample_total else 0.0,
        "horizon_peak": horizon_peak,
        "distance_max_m": float(distances_m[-1]),
        "azimuth_step_deg": float(azimuth_step_deg),
    }
    return runs, summary


def _build_run_segment(
    *,
    observer: ObserverLocation,
    azimuth_deg: float,
    start_index: int,
    end_index: int,
    distances_m: np.ndarray,
    lon_grid_deg: np.ndarray,
    lat_grid_deg: np.ndarray,
) -> WaterRunSegment:
    points_xy: list[tuple[float, float]] = []
    for sample_index in range(start_index, end_index + 1):
        x_m, y_m = _project_local_xy(
            observer.longitude_deg,
            observer.latitude_deg,
            float(lon_grid_deg[sample_index]),
            float(lat_grid_deg[sample_index]),
        )
        points_xy.append((x_m, y_m))
    return WaterRunSegment(
        azimuth_deg=float(azimuth_deg),
        start_distance_m=float(distances_m[start_index]),
        end_distance_m=float(distances_m[end_index]),
        sample_count=len(points_xy),
        points_xy=tuple(points_xy),
    )


def _collect_svg_points(
    runs: Iterable[WaterRunSegment],
    polygons: Iterable[PreparedWaterPolygon],
    *,
    max_distance_km: float,
    padding_m: float,
) -> tuple[float, float, float, float]:
    max_distance_m = float(max_distance_km) * 1000.0 + float(padding_m)
    min_x = -max_distance_m
    max_x = max_distance_m
    min_y = -max_distance_m
    max_y = max_distance_m
    for polygon in polygons:
        for ring in polygon.outer_rings_xy + polygon.inner_rings_xy:
            for x_m, y_m in ring:
                min_x = min(min_x, x_m)
                max_x = max(max_x, x_m)
                min_y = min(min_y, y_m)
                max_y = max(max_y, y_m)
    for run in runs:
        for x_m, y_m in run.points_xy:
            min_x = min(min_x, x_m)
            max_x = max(max_x, x_m)
            min_y = min(min_y, y_m)
            max_y = max(max_y, y_m)
    return min_x, max_x, min_y, max_y


def _svg_xy(
    x_m: float,
    y_m: float,
    *,
    min_x: float,
    max_y: float,
    scale: float,
    margin_px: float,
) -> tuple[float, float]:
    x_px = margin_px + (x_m - min_x) * scale
    y_px = margin_px + (max_y - y_m) * scale
    return x_px, y_px


def _svg_polyline(points: Iterable[tuple[float, float]]) -> str:
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def _polygon_to_svg_path(
    polygon: PreparedWaterPolygon,
    *,
    min_x: float,
    max_y: float,
    scale: float,
    margin_px: float,
) -> str:
    def ring_to_path(coords: Iterable[tuple[float, float]]) -> str:
        points = [
            _svg_xy(x_m, y_m, min_x=min_x, max_y=max_y, scale=scale, margin_px=margin_px)
            for x_m, y_m in coords
        ]
        if not points:
            return ""
        parts = [f"M {points[0][0]:.2f},{points[0][1]:.2f}"]
        parts.extend(f"L {x:.2f},{y:.2f}" for x, y in points[1:])
        parts.append("Z")
        return " ".join(parts)

    path_parts = [ring_to_path(ring) for ring in polygon.outer_rings_xy]
    path_parts.extend(ring_to_path(ring) for ring in polygon.inner_rings_xy)
    return " ".join(part for part in path_parts if part)


def _render_svg(
    *,
    observer: ObserverLocation,
    prepared_polygons: tuple[PreparedWaterPolygon, ...],
    runs: Iterable[WaterRunSegment],
    output_path: Path,
    max_distance_km: float,
    canvas_size: int,
    canvas_padding: float,
    summary: dict[str, object],
    azimuth_step_deg: float,
) -> None:
    min_x, max_x, min_y, max_y = _collect_svg_points(
        runs,
        prepared_polygons,
        max_distance_km=max_distance_km,
        padding_m=canvas_padding,
    )
    width_m = max_x - min_x
    height_m = max_y - min_y
    if width_m <= 0.0 or height_m <= 0.0:
        raise ValueError("SVG canvas bounds collapsed unexpectedly.")
    margin_px = 30.0
    scale = min(
        (float(canvas_size) - 2.0 * margin_px) / width_m,
        (float(canvas_size) - 2.0 * margin_px) / height_m,
    )
    max_y_px = max_y

    svg_lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{canvas_size}" height="{canvas_size}" viewBox="0 0 {canvas_size} {canvas_size}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#f8fbff"/>',
        '<g font-family="monospace" font-size="13" fill="#2c3e50">',
        '<text x="18" y="24">water-run scan demo</text>',
        f'<text x="18" y="42">observer={observer.latitude_deg:.5f},{observer.longitude_deg:.5f}</text>',
        f'<text x="18" y="60">water_runs={summary["water_run_count"]} water_samples={summary["water_total"]}/{summary["sample_total"]}</text>',
        f'<text x="18" y="78">cache={summary["cache_dir"]}</text>',
        "</g>",
    ]

    for azimuth_deg in np.arange(0.0, 360.0, float(azimuth_step_deg), dtype=np.float64):
        azimuth_rad = math.radians(float(azimuth_deg))
        x_end = math.sin(azimuth_rad) * float(max_distance_km) * 1000.0
        y_end = math.cos(azimuth_rad) * float(max_distance_km) * 1000.0
        x0, y0 = _svg_xy(0.0, 0.0, min_x=min_x, max_y=max_y_px, scale=scale, margin_px=margin_px)
        x1, y1 = _svg_xy(x_end, y_end, min_x=min_x, max_y=max_y_px, scale=scale, margin_px=margin_px)
        svg_lines.append(
            f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" stroke="#d4dde6" stroke-width="1"/>'
        )

    for polygon in prepared_polygons:
        path_data = _polygon_to_svg_path(
            polygon,
            min_x=min_x,
            max_y=max_y_px,
            scale=scale,
            margin_px=margin_px,
        )
        if path_data:
            svg_lines.append(
                f'<path d="{path_data}" fill="#64b5ff" fill-opacity="0.10" stroke="#66a8ff" stroke-width="1.2"/>'
            )

    for run in runs:
        if len(run.points_xy) < 2:
            continue
        points = [
            _svg_xy(
                x_m,
                y_m,
                min_x=min_x,
                max_y=max_y_px,
                scale=scale,
                margin_px=margin_px,
            )
            for x_m, y_m in run.points_xy
        ]
        svg_lines.append(
            f'<polyline points="{_svg_polyline(points)}" fill="none" stroke="#1b73d1" stroke-width="2.8" stroke-linecap="round" stroke-linejoin="round"/>'
        )
        start_px = _svg_xy(
            run.points_xy[0][0],
            run.points_xy[0][1],
            min_x=min_x,
            max_y=max_y_px,
            scale=scale,
            margin_px=margin_px,
        )
        end_px = _svg_xy(
            run.points_xy[-1][0],
            run.points_xy[-1][1],
            min_x=min_x,
            max_y=max_y_px,
            scale=scale,
            margin_px=margin_px,
        )
        svg_lines.append(
            f'<circle cx="{start_px[0]:.2f}" cy="{start_px[1]:.2f}" r="2.6" fill="#1b73d1"/>'
        )
        svg_lines.append(
            f'<circle cx="{end_px[0]:.2f}" cy="{end_px[1]:.2f}" r="2.6" fill="#1b73d1"/>'
        )

    observer_px = _svg_xy(0.0, 0.0, min_x=min_x, max_y=max_y_px, scale=scale, margin_px=margin_px)
    svg_lines.append(
        f'<circle cx="{observer_px[0]:.2f}" cy="{observer_px[1]:.2f}" r="5.5" fill="#e24a4a" stroke="#ffffff" stroke-width="1.4"/>'
    )
    svg_lines.append(
        f'<text x="{observer_px[0] + 8.0:.2f}" y="{observer_px[1] - 8.0:.2f}" font-family="monospace" font-size="12" fill="#8a2f2f">observer</text>'
    )
    svg_lines.append("</svg>")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(svg_lines) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    args = _parse_args(argv)
    water_footprints = _load_water_footprints(args.water_json)

    if args.fetch_dem:
        from zstarview.terrain import fetch_copernicus_dem

        download = fetch_copernicus_dem(
            observer_lat_deg=float(args.lat),
            observer_lon_deg=float(args.lon),
            max_distance_km=float(args.max_distance_km),
            margin_km=float(args.download_margin_km),
            cache_dir=args.cache_dir,
        )
        dem_paths = download.paths
        dem_source = download.source
    else:
        dem_paths = _resolve_cached_dem(
            lat_deg=float(args.lat),
            lon_deg=float(args.lon),
            max_distance_km=float(args.max_distance_km),
            margin_km=float(args.download_margin_km),
            cache_dir=args.cache_dir,
        )
        dem_source = "cache"

    dem = GeoTiffDem(dem_paths, default_elevation_m=0.0)
    try:
        bbox = build_download_bbox(
            lat_deg=float(args.lat),
            lon_deg=float(args.lon),
            radius_km=float(args.max_distance_km) + float(args.download_margin_km),
        )
        dem_grid = dem.build_grid(bbox)
        observer_ground_m = sample_ground_elevation(
            dem_grid,
            latitude_deg=float(args.lat),
            longitude_deg=float(args.lon),
            dem_resampling=str(args.dem_resampling),
        )
        observer = ObserverLocation(
            latitude_deg=float(args.lat),
            longitude_deg=float(args.lon),
            observer_ground_m=observer_ground_m,
            observer_eye_m=float(args.observer_height_m),
        )
        prepared_polygons = _build_polygon_list(
            water_footprints,
            observer=observer,
            water_extent_km=float(args.water_extent_km),
        )
        runs, scan_summary = _sample_rays(
            observer=observer,
            dem_grid=dem_grid,
            water_polygons=prepared_polygons,
            max_distance_km=float(args.max_distance_km),
            sample_step_near_m=float(args.sample_step_near_m),
            sample_step_far_m=float(args.sample_step_far_m),
            sample_step_transition_km=float(args.sample_step_transition_km),
            azimuth_step_deg=float(args.azimuth_step_deg),
            dem_resampling=str(args.dem_resampling),
        )
    finally:
        dem.close()

    summary: dict[str, object] = {
        "observer": {
            "lat": float(args.lat),
            "lon": float(args.lon),
            "ground_m": float(observer_ground_m),
            "eye_m": float(args.observer_height_m),
        },
        "cache_dir": str(args.cache_dir),
        "dem_source": dem_source,
        "water_extent_km": float(args.water_extent_km),
        "sample_step_near_m": float(args.sample_step_near_m),
        "sample_step_far_m": float(args.sample_step_far_m),
        "sample_step_transition_km": float(args.sample_step_transition_km),
        "water_polygon_count": len(prepared_polygons),
        "water_run_count": len(runs),
        "scan": scan_summary,
    }
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote_json={args.output_json}")
    if args.output_svg is not None:
        _render_svg(
            observer=observer,
            prepared_polygons=prepared_polygons,
            runs=runs,
            output_path=args.output_svg,
            max_distance_km=float(args.max_distance_km),
            canvas_size=int(args.canvas_size),
            canvas_padding=float(args.canvas_padding),
            summary={
                "water_run_count": len(runs),
                "water_total": scan_summary["water_total"],
                "sample_total": scan_summary["sample_total"],
                "cache_dir": str(args.cache_dir),
            },
            azimuth_step_deg=float(args.azimuth_step_deg),
        )
        print(f"wrote_svg={args.output_svg}")
    print(
        "summary "
        f"water_polygons={len(prepared_polygons)} "
        f"water_runs={len(runs)} "
        f"water_samples={scan_summary['water_total']}/{scan_summary['sample_total']} "
        f"water_extent_km={float(args.water_extent_km):.1f} "
        f"sample_step_near_m={float(args.sample_step_near_m):.0f} "
        f"sample_step_far_m={float(args.sample_step_far_m):.0f} "
        f"cache_dir={args.cache_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
