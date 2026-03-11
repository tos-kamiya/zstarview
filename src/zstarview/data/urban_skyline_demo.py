#!/usr/bin/env python3
"""Generate urban skyline preview images for bundled Japan tower viewpoints.

This is a standalone preprocessing sample. It reads simple building footprints
from a GeoJSON FeatureCollection, scans them from one or more observation
towers, writes azimuth-vs-altitude preview PNGs, and can export a consolidated
JSON mapping keyed by bundled tower id.

Input assumptions
-----------------
- The observer is placed at the bundled tower's ``height_m`` above the tower
  base, matching the current app's viewpoint model.
- Building input is a GeoJSON FeatureCollection with Polygon or MultiPolygon
  geometries.
- Building height is read from one of the configured height property names.
- Ground elevation is ignored in this first sample. Heights are compared in the
  same local vertical frame.

Typical usage
-------------
  uv run src/zstarview/data/urban_skyline_demo.py \
      --buildings /path/to/plateau_buildings.geojson \
      --tower "Tokyo Tower"

  uv run src/zstarview/data/urban_skyline_demo.py \
      --buildings /path/to/plateau_buildings.geojson \
      --all-japan-towers \
      --radius-km 5 \
      --output-dir src/zstarview/data/viewpoints/urban_skyline
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from PIL import Image, ImageDraw
from pyproj import CRS, Transformer

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
DATA_ROOT = SRC_ROOT / "zstarview" / "data"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.tower_viewpoints import (  # noqa: E402
    TowerViewpoint,
    load_tower_viewpoints,
    resolve_tower_viewpoint,
)


DEFAULT_HEIGHT_FIELDS = (
    "measuredHeight",
    "bldg:measuredHeight",
    "height",
    "height_m",
)
DEFAULT_CUMULATIVE_RADII_KM = (
    0.8888888888888888,
    1.3333333333333333,
    2.0,
    3.0,
    4.5,
    6.75,
    10.125,
    15.1875,
    22.78125,
)
DEFAULT_RADIUS_BAND_WIDTH_M = 90.0


@dataclass(frozen=True)
class BuildingFootprint:
    building_id: str
    height_m: float
    rings_lonlat: tuple[tuple[tuple[float, float], ...], ...]


@dataclass(frozen=True)
class SkylineSample:
    azimuth_deg: float
    altitude_deg: float


@dataclass(frozen=True)
class SkylineResult:
    tower: TowerViewpoint
    samples: tuple[SkylineSample, ...]
    buildings_considered: int
    buildings_contributing: int
    peak_altitude_deg: float
    peak_azimuth_deg: float


@dataclass(frozen=True)
class SkylineRadiusResult:
    radius_km: float
    result: SkylineResult


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate urban skyline PNG previews for bundled Japan tower "
            "viewpoints from a GeoJSON building dataset."
        )
    )
    parser.add_argument(
        "--buildings",
        type=Path,
        help="GeoJSON FeatureCollection containing building footprints.",
    )
    parser.add_argument(
        "--tower",
        action="append",
        default=[],
        help="Tower viewpoint name or wikidata:Q... identifier. May be repeated.",
    )
    parser.add_argument(
        "--all-japan-towers",
        action="store_true",
        help="Process every bundled tower whose coordinates fall in Japan.",
    )
    parser.add_argument(
        "--list-japan-towers",
        action="store_true",
        help="List bundled Japan tower names and exit.",
    )
    parser.add_argument(
        "--radius-km",
        type=float,
        default=max(DEFAULT_CUMULATIVE_RADII_KM),
        help="Maximum building search radius around each tower (default: 20.0).",
    )
    parser.add_argument(
        "--cumulative-radius-km",
        action="append",
        default=[],
        help=(
            "Skyline scan radius in km. May be repeated. Defaults to: "
            + ", ".join(str(value) for value in DEFAULT_CUMULATIVE_RADII_KM)
        ),
    )
    parser.add_argument(
        "--radius-band-width-m",
        type=float,
        default=DEFAULT_RADIUS_BAND_WIDTH_M,
        help="Radial scan band width in meters for each skyline layer (default: 90.0).",
    )
    parser.add_argument(
        "--azimuth-step",
        type=float,
        default=0.1,
        help="Azimuth sampling step in degrees (default: 0.1).",
    )
    parser.add_argument(
        "--edge-sample-step-m",
        type=float,
        default=5.0,
        help="Approximate spacing for sampling building edges in meters (default: 5.0).",
    )
    parser.add_argument(
        "--height-field",
        action="append",
        default=[],
        help=(
            "Building height property name. May be repeated. Defaults to: "
            + ", ".join(DEFAULT_HEIGHT_FIELDS)
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DATA_ROOT / "viewpoints" / "urban_skyline",
        help="Directory where PNG previews will be written.",
    )
    parser.add_argument(
        "--write-json",
        action="store_true",
        help="Also write a consolidated skyline JSON keyed by tower id.",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=None,
        help=(
            "Where to write the consolidated JSON. Defaults to "
            "<output-dir>/urban_skyline_profiles.json."
        ),
    )
    return parser


def is_japan_tower(tower: TowerViewpoint) -> bool:
    return 20.0 <= tower.latitude_deg <= 46.0 and 122.0 <= tower.longitude_deg <= 154.0


def list_japan_towers() -> tuple[TowerViewpoint, ...]:
    return tuple(tower for tower in load_tower_viewpoints() if is_japan_tower(tower))


def load_building_footprints(
    path: Path,
    *,
    height_fields: Sequence[str],
) -> tuple[BuildingFootprint, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    features = payload.get("features")
    if payload.get("type") != "FeatureCollection" or not isinstance(features, list):
        raise ValueError("Building GeoJSON must be a FeatureCollection with a features list.")

    buildings: list[BuildingFootprint] = []
    for index, feature in enumerate(features):
        if not isinstance(feature, dict):
            continue
        geometry = feature.get("geometry")
        properties = feature.get("properties")
        if not isinstance(geometry, dict) or not isinstance(properties, dict):
            continue
        rings = tuple(iter_exterior_rings_lonlat(geometry))
        if not rings:
            continue
        height_m = extract_building_height_m(properties, height_fields)
        if height_m is None or not math.isfinite(height_m) or height_m <= 0.0:
            continue
        building_id = str(
            properties.get("id")
            or properties.get("gml_id")
            or properties.get("building_id")
            or f"feature-{index}"
        )
        buildings.append(
            BuildingFootprint(
                building_id=building_id,
                height_m=float(height_m),
                rings_lonlat=rings,
            )
        )
    if not buildings:
        raise ValueError("No usable building footprints were found in the GeoJSON input.")
    return tuple(buildings)


def iter_exterior_rings_lonlat(geometry: dict[str, object]) -> Iterable[tuple[tuple[float, float], ...]]:
    geometry_type = geometry.get("type")
    coordinates = geometry.get("coordinates")
    if geometry_type == "Polygon":
        if isinstance(coordinates, list) and coordinates:
            ring = parse_ring_lonlat(coordinates[0])
            if ring:
                yield ring
        return
    if geometry_type == "MultiPolygon" and isinstance(coordinates, list):
        for polygon in coordinates:
            if isinstance(polygon, list) and polygon:
                ring = parse_ring_lonlat(polygon[0])
                if ring:
                    yield ring


def parse_ring_lonlat(raw_ring: object) -> tuple[tuple[float, float], ...]:
    if not isinstance(raw_ring, list):
        return ()
    points: list[tuple[float, float]] = []
    for raw_point in raw_ring:
        if (
            isinstance(raw_point, (list, tuple))
            and len(raw_point) >= 2
            and isinstance(raw_point[0], (int, float))
            and isinstance(raw_point[1], (int, float))
        ):
            points.append((float(raw_point[0]), float(raw_point[1])))
    if len(points) < 4:
        return ()
    if points[0] != points[-1]:
        points.append(points[0])
    return tuple(points)


def extract_building_height_m(properties: dict[str, object], height_fields: Sequence[str]) -> float | None:
    for field in height_fields:
        value = properties.get(field)
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, str):
            text = value.strip()
            if not text:
                continue
            try:
                return float(text)
            except ValueError:
                continue
    return None


def sanitize_slug(text: str) -> str:
    collapsed = re.sub(r"[^0-9A-Za-z]+", "_", text.strip())
    collapsed = collapsed.strip("_")
    return collapsed.lower() or "tower"


def make_local_transformer(tower: TowerViewpoint) -> Transformer:
    local_crs = CRS.from_proj4(
        f"+proj=aeqd +lat_0={tower.latitude_deg} +lon_0={tower.longitude_deg} "
        "+datum=WGS84 +units=m +no_defs"
    )
    return Transformer.from_crs("EPSG:4326", local_crs, always_xy=True)


def project_ring_xy(
    ring_lonlat: Sequence[tuple[float, float]],
    transformer: Transformer,
) -> np.ndarray:
    lon = np.array([point[0] for point in ring_lonlat], dtype=np.float64)
    lat = np.array([point[1] for point in ring_lonlat], dtype=np.float64)
    x, y = transformer.transform(lon, lat)
    return np.column_stack((np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)))


def bbox_min_distance_m(points_xy: np.ndarray) -> float:
    min_x = float(np.min(points_xy[:, 0]))
    max_x = float(np.max(points_xy[:, 0]))
    min_y = float(np.min(points_xy[:, 1]))
    max_y = float(np.max(points_xy[:, 1]))
    nearest_x = 0.0 if min_x <= 0.0 <= max_x else min(abs(min_x), abs(max_x))
    nearest_y = 0.0 if min_y <= 0.0 <= max_y else min(abs(min_y), abs(max_y))
    return math.hypot(nearest_x, nearest_y)


def sample_segment_points_xy(
    start_xy: np.ndarray,
    end_xy: np.ndarray,
    *,
    sample_step_m: float,
) -> np.ndarray:
    delta = end_xy - start_xy
    length = float(np.hypot(delta[0], delta[1]))
    count = max(2, int(math.ceil(length / max(sample_step_m, 0.1))) + 1)
    t = np.linspace(0.0, 1.0, num=count, dtype=np.float64)
    return start_xy[None, :] + t[:, None] * delta[None, :]


def sample_ring_points_xy(ring_xy: np.ndarray, *, sample_step_m: float) -> np.ndarray:
    segments: list[np.ndarray] = []
    for start_xy, end_xy in zip(ring_xy[:-1], ring_xy[1:]):
        segment = sample_segment_points_xy(start_xy, end_xy, sample_step_m=sample_step_m)
        if segments:
            segment = segment[1:]
        segments.append(segment)
    if not segments:
        return np.empty((0, 2), dtype=np.float64)
    return np.vstack(segments)


def iter_sampled_ring_segments_xy(
    ring_xy: np.ndarray,
    *,
    sample_step_m: float,
) -> Iterable[np.ndarray]:
    for start_xy, end_xy in zip(ring_xy[:-1], ring_xy[1:]):
        segment = sample_segment_points_xy(start_xy, end_xy, sample_step_m=sample_step_m)
        if segment.size == 0:
            continue
        yield segment


def _unwrap_azimuth_delta_deg(start_deg: float, end_deg: float) -> float:
    delta_deg = end_deg - start_deg
    if delta_deg > 180.0:
        delta_deg -= 360.0
    elif delta_deg < -180.0:
        delta_deg += 360.0
    return delta_deg


def update_altitude_bins_from_polyline(
    altitudes: np.ndarray,
    *,
    azimuth_deg: np.ndarray,
    altitude_deg: np.ndarray,
    azimuth_step_deg: float,
    closed: bool = True,
) -> bool:
    if azimuth_deg.size == 0 or altitude_deg.size == 0:
        return False
    if azimuth_deg.size != altitude_deg.size:
        raise ValueError("azimuth_deg and altitude_deg must have the same length.")

    updated = False
    bin_count = altitudes.size
    point_count = int(azimuth_deg.size)
    if point_count == 1:
        index = int(np.rint(float(azimuth_deg[0]) / azimuth_step_deg)) % bin_count
        old_value = float(altitudes[index])
        new_value = max(old_value, float(altitude_deg[0]))
        altitudes[index] = new_value
        return new_value > old_value
    if point_count < 2:
        return False
    segment_count = point_count if closed else point_count - 1
    for index in range(segment_count):
        start_az = float(azimuth_deg[index])
        end_az = float(azimuth_deg[(index + 1) % point_count])
        start_alt = float(altitude_deg[index])
        end_alt = float(altitude_deg[(index + 1) % point_count])
        delta_az = _unwrap_azimuth_delta_deg(start_az, end_az)
        steps = max(1, int(math.ceil(abs(delta_az) / azimuth_step_deg)))
        t = np.linspace(0.0, 1.0, num=steps + 1, dtype=np.float64)
        segment_az = (start_az + t * delta_az) % 360.0
        segment_alt = start_alt + t * (end_alt - start_alt)
        indices = np.rint(segment_az / azimuth_step_deg).astype(np.int64) % bin_count
        old_values = altitudes[indices].copy()
        np.maximum.at(altitudes, indices, segment_alt)
        if not updated and np.any(altitudes[indices] > old_values):
            updated = True
    return updated


def iter_true_runs(mask: np.ndarray) -> Iterable[slice]:
    if mask.ndim != 1:
        raise ValueError("mask must be 1-dimensional.")
    start: int | None = None
    for index, flag in enumerate(mask.tolist()):
        if flag:
            if start is None:
                start = index
            continue
        if start is not None:
            yield slice(start, index)
            start = None
    if start is not None:
        yield slice(start, mask.size)


def compute_band_ends_m(
    band_starts_m: np.ndarray,
    *,
    fallback_band_width_m: float,
) -> np.ndarray:
    if band_starts_m.ndim != 1 or band_starts_m.size == 0:
        raise ValueError("band_starts_m must be a non-empty 1-dimensional array.")
    if band_starts_m.size == 1:
        return band_starts_m + float(fallback_band_width_m)
    deltas_m = np.diff(band_starts_m)
    widths_m = deltas_m / 3.0
    prev_start_m = float(band_starts_m[-2])
    last_start_m = float(band_starts_m[-1])
    if prev_start_m > 0.0 and last_start_m > prev_start_m:
        ratio = last_start_m / prev_start_m
        next_start_m = last_start_m * ratio
        trailing_width_m = (next_start_m - last_start_m) / 3.0
    else:
        trailing_width_m = float(deltas_m[-1] / 3.0)
    return np.concatenate(
        (
            band_starts_m[:-1] + widths_m,
            np.array([band_starts_m[-1] + trailing_width_m], dtype=np.float64),
        )
    )


def compute_urban_skyline(
    tower: TowerViewpoint,
    buildings: Sequence[BuildingFootprint],
    *,
    radius_km: float,
    azimuth_step_deg: float,
    edge_sample_step_m: float,
) -> SkylineResult:
    if radius_km <= 0.0:
        raise ValueError("--radius-km must be positive.")
    if azimuth_step_deg <= 0.0:
        raise ValueError("--azimuth-step must be positive.")
    if edge_sample_step_m <= 0.0:
        raise ValueError("--edge-sample-step-m must be positive.")

    transformer = make_local_transformer(tower)
    azimuths = np.arange(0.0, 360.0, azimuth_step_deg, dtype=np.float64)
    altitudes = np.full(azimuths.shape, -90.0, dtype=np.float64)
    radius_m = radius_km * 1000.0
    observer_height_m = tower.observer_height_m
    buildings_considered = 0
    buildings_contributing = 0

    for building in buildings:
        projected_rings = tuple(project_ring_xy(ring, transformer) for ring in building.rings_lonlat)
        if not projected_rings:
            continue
        min_distance = min(bbox_min_distance_m(ring_xy) for ring_xy in projected_rings if ring_xy.size > 0)
        if min_distance > radius_m:
            continue
        buildings_considered += 1
        contributed = False
        for ring_xy in projected_rings:
            for sampled_points in iter_sampled_ring_segments_xy(ring_xy, sample_step_m=edge_sample_step_m):
                distances = np.hypot(sampled_points[:, 0], sampled_points[:, 1])
                valid = (distances > 0.1) & (distances <= radius_m)
                if not np.any(valid):
                    continue
                sampled_points = sampled_points[valid]
                distances = distances[valid]
                azimuth_deg = (np.degrees(np.arctan2(sampled_points[:, 0], sampled_points[:, 1])) + 360.0) % 360.0
                altitude_deg = np.degrees(np.arctan2(building.height_m - observer_height_m, distances))
                if update_altitude_bins_from_polyline(
                    altitudes,
                    azimuth_deg=azimuth_deg,
                    altitude_deg=altitude_deg,
                    azimuth_step_deg=azimuth_step_deg,
                    closed=False,
                ):
                    contributed = True
        if contributed:
            buildings_contributing += 1

    samples = tuple(
        SkylineSample(azimuth_deg=float(az), altitude_deg=float(alt))
        for az, alt in zip(azimuths, altitudes)
    )
    peak_index = int(np.argmax(altitudes))
    return SkylineResult(
        tower=tower,
        samples=samples,
        buildings_considered=buildings_considered,
        buildings_contributing=buildings_contributing,
        peak_altitude_deg=float(altitudes[peak_index]),
        peak_azimuth_deg=float(azimuths[peak_index]),
    )


def normalize_cumulative_radii_km(
    radii_km: Sequence[float] | None,
    *,
    max_radius_km: float | None = None,
) -> tuple[float, ...]:
    source = DEFAULT_CUMULATIVE_RADII_KM if not radii_km else tuple(float(value) for value in radii_km)
    normalized = tuple(sorted({float(value) for value in source if float(value) > 0.0}))
    if not normalized:
        raise ValueError("At least one positive cumulative radius is required.")
    if max_radius_km is None:
        return normalized
    limited = tuple(value for value in normalized if value <= float(max_radius_km) + 1e-9)
    if limited:
        return limited
    raise ValueError("--radius-km is smaller than the smallest cumulative radius.")


def compute_cumulative_urban_skyline(
    tower: TowerViewpoint,
    buildings: Sequence[BuildingFootprint],
    *,
    radii_km: Sequence[float],
    radius_band_width_m: float,
    azimuth_step_deg: float,
    edge_sample_step_m: float,
) -> tuple[SkylineRadiusResult, ...]:
    if azimuth_step_deg <= 0.0:
        raise ValueError("--azimuth-step must be positive.")
    if radius_band_width_m <= 0.0:
        raise ValueError("--radius-band-width-m must be positive.")
    if edge_sample_step_m <= 0.0:
        raise ValueError("--edge-sample-step-m must be positive.")

    normalized_radii_km = normalize_cumulative_radii_km(radii_km)
    band_starts_m = np.array(normalized_radii_km, dtype=np.float64) * 1000.0
    band_ends_m = compute_band_ends_m(
        band_starts_m,
        fallback_band_width_m=float(radius_band_width_m),
    )
    transformer = make_local_transformer(tower)
    azimuths = np.arange(0.0, 360.0, azimuth_step_deg, dtype=np.float64)
    altitude_layers = np.full((band_starts_m.size, azimuths.size), -90.0, dtype=np.float64)
    buildings_considered = np.zeros(band_starts_m.size, dtype=np.int64)
    buildings_contributing = np.zeros(band_starts_m.size, dtype=np.int64)
    observer_height_m = tower.observer_height_m
    max_radius_m = float(band_ends_m[-1])

    for building in buildings:
        projected_rings = tuple(project_ring_xy(ring, transformer) for ring in building.rings_lonlat)
        if not projected_rings:
            continue
        min_distance = min(bbox_min_distance_m(ring_xy) for ring_xy in projected_rings if ring_xy.size > 0)
        if min_distance > max_radius_m:
            continue
        considered = np.zeros(band_starts_m.size, dtype=bool)
        contributed = np.zeros(band_starts_m.size, dtype=bool)

        for ring_xy in projected_rings:
            for sampled_points in iter_sampled_ring_segments_xy(ring_xy, sample_step_m=edge_sample_step_m):
                distances = np.hypot(sampled_points[:, 0], sampled_points[:, 1])
                valid = (distances > 0.1) & (distances <= max_radius_m)
                if not np.any(valid):
                    continue
                sampled_points = sampled_points[valid]
                distances = distances[valid]
                azimuth_deg = (np.degrees(np.arctan2(sampled_points[:, 0], sampled_points[:, 1])) + 360.0) % 360.0
                altitude_deg = np.degrees(np.arctan2(building.height_m - observer_height_m, distances))

                for radius_index in range(band_starts_m.size):
                    mask = (distances >= band_starts_m[radius_index]) & (
                        distances <= band_ends_m[radius_index]
                    )
                    if not np.any(mask):
                        continue
                    considered[radius_index] = True
                    for run in iter_true_runs(mask):
                        if update_altitude_bins_from_polyline(
                            altitude_layers[radius_index],
                            azimuth_deg=azimuth_deg[run],
                            altitude_deg=altitude_deg[run],
                            azimuth_step_deg=azimuth_step_deg,
                            closed=False,
                        ):
                            contributed[radius_index] = True
        buildings_considered += considered.astype(np.int64)

        buildings_contributing += contributed.astype(np.int64)

    results: list[SkylineRadiusResult] = []
    for radius_index, radius_km in enumerate(normalized_radii_km):
        peak_index = int(np.argmax(altitude_layers[radius_index]))
        result = SkylineResult(
            tower=tower,
            samples=tuple(
                SkylineSample(azimuth_deg=float(az), altitude_deg=float(alt))
                for az, alt in zip(azimuths, altitude_layers[radius_index])
            ),
            buildings_considered=int(buildings_considered[radius_index]),
            buildings_contributing=int(buildings_contributing[radius_index]),
            peak_altitude_deg=float(altitude_layers[radius_index, peak_index]),
            peak_azimuth_deg=float(azimuths[peak_index]),
        )
        results.append(SkylineRadiusResult(radius_km=float(radius_km), result=result))
    return tuple(results)


def write_preview_png(path: Path, result: SkylineResult) -> None:
    width = 1440
    height = 760
    margin_left = 90
    margin_right = 30
    margin_top = 40
    margin_bottom = 90
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom

    altitudes = np.array([sample.altitude_deg for sample in result.samples], dtype=np.float64)
    min_alt = math.floor(min(-10.0, float(np.nanmin(altitudes)) - 1.0))
    max_alt = math.ceil(max(10.0, float(np.nanmax(altitudes)) + 1.0))
    if math.isclose(min_alt, max_alt):
        max_alt = min_alt + 1.0

    image = Image.new("RGB", (width, height), (249, 248, 245))
    draw = ImageDraw.Draw(image)
    plot_box = (
        margin_left,
        margin_top,
        margin_left + plot_width,
        margin_top + plot_height,
    )
    draw.rectangle(plot_box, outline=(85, 85, 85), width=1)

    for az in range(0, 361, 45):
        x = margin_left + (az / 360.0) * plot_width
        draw.line([(x, margin_top), (x, margin_top + plot_height)], fill=(220, 220, 220), width=1)
        draw.text((x + 4, margin_top + plot_height + 8), f"{az}°", fill=(60, 60, 60))

    for alt in range(int(min_alt), int(max_alt) + 1, 5):
        y = margin_top + plot_height - ((alt - min_alt) / (max_alt - min_alt)) * plot_height
        color = (165, 165, 165) if alt == 0 else (220, 220, 220)
        draw.line([(margin_left, y), (margin_left + plot_width, y)], fill=color, width=1)
        draw.text((12, y - 8), f"{alt}°", fill=(60, 60, 60))

    fill_points: list[tuple[float, float]] = [(margin_left, margin_top + plot_height)]
    polyline: list[tuple[float, float]] = []
    for sample in result.samples:
        x = margin_left + (sample.azimuth_deg / 360.0) * plot_width
        y = margin_top + plot_height - ((sample.altitude_deg - min_alt) / (max_alt - min_alt)) * plot_height
        point = (x, y)
        polyline.append(point)
        fill_points.append(point)
    fill_points.append((margin_left + plot_width, margin_top + plot_height))
    if len(fill_points) >= 3:
        draw.polygon(fill_points, fill=(80, 112, 140))
    if len(polyline) >= 2:
        draw.line(polyline, fill=(14, 43, 68), width=3)

    title = f"Urban Skyline Preview: {result.tower.name}"
    subtitle = (
        f"lat={result.tower.latitude_deg:.5f}  lon={result.tower.longitude_deg:.5f}  "
        f"observer_height={result.tower.observer_height_m:.1f} m"
    )
    summary = (
        f"buildings considered={result.buildings_considered}  "
        f"contributing={result.buildings_contributing}  "
        f"peak={result.peak_altitude_deg:.2f}° at az={result.peak_azimuth_deg:.1f}°"
    )
    draw.text((margin_left, 10), title, fill=(15, 15, 15))
    draw.text((margin_left, height - 48), subtitle, fill=(50, 50, 50))
    draw.text((margin_left, height - 24), summary, fill=(50, 50, 50))

    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def skyline_result_to_payload(result: SkylineResult) -> dict[str, object]:
    return {
        "name": result.tower.name,
        "profile": [
            {"az": round(sample.azimuth_deg, 6), "alt": round(sample.altitude_deg, 6)}
            for sample in result.samples
        ],
    }


def skyline_radius_results_to_payload(
    radius_results: Sequence[SkylineRadiusResult],
) -> dict[str, object]:
    if not radius_results:
        raise ValueError("radius_results must not be empty")
    return {
        "name": radius_results[0].result.tower.name,
        "profiles": [
            {
                "radius_km": round(radius_result.radius_km, 6),
                "profile": [
                    {"az": round(sample.azimuth_deg, 6), "alt": round(sample.altitude_deg, 6)}
                    for sample in radius_result.result.samples
                ],
            }
            for radius_result in radius_results
        ],
    }


def write_profiles_json(
    path: Path,
    results: Sequence[tuple[TowerViewpoint, Sequence[SkylineRadiusResult]]],
) -> None:
    payload = {
        tower.id: skyline_radius_results_to_payload(radius_results)
        for tower, radius_results in results
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def select_towers(args: argparse.Namespace) -> tuple[TowerViewpoint, ...]:
    if args.all_japan_towers:
        return list_japan_towers()

    selected: list[TowerViewpoint] = []
    for query in args.tower:
        tower = resolve_tower_viewpoint(query)
        if tower is None:
            raise ValueError(f"Tower viewpoint not found: {query}")
        if not is_japan_tower(tower):
            raise ValueError(f"Tower is not in Japan: {tower.name}")
        selected.append(tower)
    if not selected:
        raise ValueError("Specify --tower NAME or use --all-japan-towers.")
    return tuple(selected)


def print_result(result: SkylineResult, png_path: Path) -> None:
    print(
        f"[ok] {result.tower.name}: {png_path}  "
        f"buildings={result.buildings_considered}/{result.buildings_contributing}  "
        f"peak={result.peak_altitude_deg:.2f}deg@{result.peak_azimuth_deg:.1f}"
    )


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.list_japan_towers:
        for tower in list_japan_towers():
            print(f"{tower.name}\t{tower.latitude_deg:.6f}\t{tower.longitude_deg:.6f}\t{tower.height_m:.1f}")
        return 0

    try:
        if args.buildings is None:
            raise ValueError("Specify --buildings unless using --list-japan-towers.")
        towers = select_towers(args)
        height_fields = tuple(args.height_field) if args.height_field else DEFAULT_HEIGHT_FIELDS
        buildings = load_building_footprints(args.buildings, height_fields=height_fields)
        results: list[tuple[TowerViewpoint, Sequence[SkylineRadiusResult]]] = []
        cumulative_radii_km = normalize_cumulative_radii_km(
            [float(value) for value in args.cumulative_radius_km],
            max_radius_km=float(args.radius_km),
        )
        for tower in towers:
            radius_results = compute_cumulative_urban_skyline(
                tower,
                buildings,
                radii_km=cumulative_radii_km,
                radius_band_width_m=float(args.radius_band_width_m),
                azimuth_step_deg=float(args.azimuth_step),
                edge_sample_step_m=float(args.edge_sample_step_m),
            )
            results.append((tower, radius_results))
            result = radius_results[-1].result
            stem = sanitize_slug(str(tower.meta.get("slug") or tower.name))
            png_path = args.output_dir / f"{stem}_urban.png"
            write_preview_png(png_path, result)
            print_result(result, png_path)
        if args.write_json:
            json_output = args.json_output or (args.output_dir / "urban_skyline_profiles.json")
            write_profiles_json(json_output, results)
            print(f"[ok] skyline-json: {json_output}")
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
