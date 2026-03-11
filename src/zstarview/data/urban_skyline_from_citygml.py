#!/usr/bin/env python3
"""Build urban skyline profiles directly from PLATEAU CityGML building tiles.

This CLI scans a directory of ``udx/bldg/*.gml`` files, selects the tiles whose
envelopes intersect a radius around the observer, computes a skyline profile for
each tile, and merges the profiles with an element-wise max.
"""

from __future__ import annotations

import argparse
import os
import math
import re
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
DATA_ROOT = SRC_ROOT / "zstarview" / "data"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.tower_viewpoints import TowerViewpoint, resolve_tower_viewpoint  # noqa: E402
from zstarview.data.urban_skyline_demo import (  # noqa: E402
    BuildingFootprint,
    DEFAULT_CUMULATIVE_RADII_KM,
    DEFAULT_RADIUS_BAND_WIDTH_M,
    SkylineResult,
    SkylineRadiusResult,
    SkylineSample,
    compute_cumulative_urban_skyline,
    normalize_cumulative_radii_km,
    is_japan_tower,
    sanitize_slug,
    write_preview_png,
    write_profiles_json,
)


_ENVELOPE_RE = re.compile(r"<gml:(lowerCorner|upperCorner)>([^<]+)</gml:")
_GML_NS = {
    "bldg": "http://www.opengis.net/citygml/building/2.0",
    "gml": "http://www.opengis.net/gml",
}


@dataclass(frozen=True)
class TileEnvelope:
    path: Path
    min_lat_deg: float
    min_lon_deg: float
    max_lat_deg: float
    max_lon_deg: float


@dataclass(frozen=True)
class TileSkyline:
    envelope: TileEnvelope
    radius_results: tuple[SkylineRadiusResult, ...]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate urban skyline profiles directly from a directory of "
            "PLATEAU CityGML building tiles."
        )
    )
    parser.add_argument(
        "--citygml-dir",
        type=Path,
        required=True,
        help="Directory containing PLATEAU building GML files such as udx/bldg/.",
    )
    parser.add_argument(
        "--tower",
        required=True,
        help="Bundled tower viewpoint name or wikidata:Q... identifier.",
    )
    parser.add_argument(
        "--radius-km",
        type=float,
        default=max(DEFAULT_CUMULATIVE_RADII_KM),
        help="Maximum tile/building search radius around the tower (default: 20.0).",
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
        "--output-dir",
        type=Path,
        default=DATA_ROOT / "viewpoints" / "urban_skyline",
        help="Directory where PNG previews and JSON outputs will be written.",
    )
    parser.add_argument(
        "--write-json",
        action="store_true",
        help="Write a consolidated skyline JSON keyed by tower id.",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=None,
        help="Override consolidated JSON output path.",
    )
    parser.add_argument(
        "--print-selected-tiles",
        action="store_true",
        help="Print selected tile filenames before processing.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 1) // 2),
        help="Number of worker processes for per-tile processing (default: half of CPU cores).",
    )
    return parser


def load_tile_envelope(path: Path) -> TileEnvelope | None:
    with path.open("r", encoding="utf-8-sig", errors="ignore") as handle:
        head_lines: list[str] = []
        for _ in range(8):
            line = handle.readline()
            if not line:
                break
            head_lines.append(line)
    text = "".join(head_lines)
    matches = _ENVELOPE_RE.findall(text)
    if len(matches) < 2:
        return None
    lower = [float(value) for value in matches[0][1].split()]
    upper = [float(value) for value in matches[1][1].split()]
    return TileEnvelope(
        path=path,
        min_lat_deg=min(lower[0], upper[0]),
        min_lon_deg=min(lower[1], upper[1]),
        max_lat_deg=max(lower[0], upper[0]),
        max_lon_deg=max(lower[1], upper[1]),
    )


def approx_distance_km(
    lat0_deg: float,
    lon0_deg: float,
    lat1_deg: float,
    lon1_deg: float,
) -> float:
    dlon_km = (lon1_deg - lon0_deg) * 111.32 * math.cos(math.radians((lat0_deg + lat1_deg) * 0.5))
    dlat_km = (lat1_deg - lat0_deg) * 111.32
    return math.hypot(dlon_km, dlat_km)


def envelope_min_distance_km(
    envelope: TileEnvelope,
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
) -> float:
    nearest_lat = min(max(observer_lat_deg, envelope.min_lat_deg), envelope.max_lat_deg)
    nearest_lon = min(max(observer_lon_deg, envelope.min_lon_deg), envelope.max_lon_deg)
    return approx_distance_km(observer_lat_deg, observer_lon_deg, nearest_lat, nearest_lon)


def select_tile_envelopes(
    citygml_dir: Path,
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    radius_km: float,
) -> tuple[TileEnvelope, ...]:
    envelopes: list[TileEnvelope] = []
    for path in sorted(citygml_dir.glob("*.gml")):
        envelope = load_tile_envelope(path)
        if envelope is None:
            continue
        if envelope_min_distance_km(
            envelope,
            observer_lat_deg=observer_lat_deg,
            observer_lon_deg=observer_lon_deg,
        ) <= radius_km:
            envelopes.append(envelope)
    if not envelopes:
        raise ValueError(f"No CityGML building tiles found within {radius_km:.1f} km.")
    return tuple(envelopes)


def parse_citygml_buildings(path: Path) -> tuple[BuildingFootprint, ...]:
    root = ET.parse(path).getroot()
    buildings: list[BuildingFootprint] = []
    for index, building in enumerate(root.findall(".//bldg:Building", _GML_NS)):
        height_text = building.findtext("bldg:measuredHeight", default="", namespaces=_GML_NS)
        try:
            height_m = float(height_text)
        except ValueError:
            continue

        rings: list[tuple[tuple[float, float], ...]] = []
        for pos in building.findall(".//bldg:lod0RoofEdge//gml:posList", _GML_NS):
            ring = parse_pos_list_ring(pos.text)
            if ring:
                rings.append(ring)
        if not rings:
            continue
        building_id = building.attrib.get("{http://www.opengis.net/gml}id", f"bldg-{index}")
        buildings.append(
            BuildingFootprint(
                building_id=building_id,
                height_m=height_m,
                rings_lonlat=tuple(rings),
            )
        )
    return tuple(buildings)


def parse_pos_list_ring(text: str | None) -> tuple[tuple[float, float], ...]:
    if not text:
        return ()
    values = [float(value) for value in text.split()]
    if len(values) < 12:
        return ()
    coords: list[tuple[float, float]] = []
    for index in range(0, len(values), 3):
        lat_deg = values[index]
        lon_deg = values[index + 1]
        point = (lon_deg, lat_deg)
        if not coords or coords[-1] != point:
            coords.append(point)
    if len(coords) < 4:
        return ()
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    return tuple(coords)


def combine_tile_results(
    tower: TowerViewpoint,
    tile_results: Sequence[TileSkyline],
    *,
    radii_km: Sequence[float],
    azimuth_step_deg: float,
) -> tuple[SkylineRadiusResult, ...]:
    azimuths = np.arange(0.0, 360.0, azimuth_step_deg, dtype=np.float64)
    normalized_radii_km = normalize_cumulative_radii_km(radii_km)
    altitude_layers = np.full((len(normalized_radii_km), azimuths.size), -90.0, dtype=np.float64)
    buildings_considered = np.zeros(len(normalized_radii_km), dtype=np.int64)
    buildings_contributing = np.zeros(len(normalized_radii_km), dtype=np.int64)
    for tile_result in tile_results:
        for radius_index, radius_result in enumerate(tile_result.radius_results):
            sample_altitudes = np.full(azimuths.shape, -90.0, dtype=np.float64)
            for sample in radius_result.result.samples:
                index = int(np.rint(sample.azimuth_deg / azimuth_step_deg)) % azimuths.size
                sample_altitudes[index] = max(sample_altitudes[index], sample.altitude_deg)
            np.maximum(altitude_layers[radius_index], sample_altitudes, out=altitude_layers[radius_index])
            buildings_considered[radius_index] += radius_result.result.buildings_considered
            buildings_contributing[radius_index] += radius_result.result.buildings_contributing

    combined: list[SkylineRadiusResult] = []
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
        combined.append(SkylineRadiusResult(radius_km=float(radius_km), result=result))
    return tuple(combined)


def select_tower(query: str) -> TowerViewpoint:
    tower = resolve_tower_viewpoint(query)
    if tower is None:
        raise ValueError(f"Tower viewpoint not found: {query}")
    if not is_japan_tower(tower):
        raise ValueError(f"Tower is not in Japan: {tower.name}")
    return tower


def compute_tile_skyline(
    envelope: TileEnvelope,
    *,
    tower: TowerViewpoint,
    cumulative_radii_km: Sequence[float],
    radius_band_width_m: float,
    azimuth_step_deg: float,
    edge_sample_step_m: float,
) -> TileSkyline | None:
    buildings = parse_citygml_buildings(envelope.path)
    if not buildings:
        return None
    return TileSkyline(
        envelope=envelope,
        radius_results=compute_cumulative_urban_skyline(
            tower,
            buildings,
            radii_km=cumulative_radii_km,
            radius_band_width_m=radius_band_width_m,
            azimuth_step_deg=azimuth_step_deg,
            edge_sample_step_m=edge_sample_step_m,
        ),
    )


def compute_tile_skylines(
    envelopes: Sequence[TileEnvelope],
    *,
    tower: TowerViewpoint,
    cumulative_radii_km: Sequence[float],
    radius_band_width_m: float,
    azimuth_step_deg: float,
    edge_sample_step_m: float,
    workers: int,
) -> list[TileSkyline]:
    max_workers = max(1, min(int(workers), len(envelopes)))
    if max_workers == 1:
        tile_results: list[TileSkyline] = []
        for envelope in envelopes:
            tile_result = compute_tile_skyline(
                envelope,
                tower=tower,
                cumulative_radii_km=cumulative_radii_km,
                radius_band_width_m=radius_band_width_m,
                azimuth_step_deg=azimuth_step_deg,
                edge_sample_step_m=edge_sample_step_m,
            )
            if tile_result is not None:
                tile_results.append(tile_result)
        return tile_results

    tile_results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                compute_tile_skyline,
                envelope,
                tower=tower,
                cumulative_radii_km=tuple(cumulative_radii_km),
                radius_band_width_m=radius_band_width_m,
                azimuth_step_deg=azimuth_step_deg,
                edge_sample_step_m=edge_sample_step_m,
            )
            for envelope in envelopes
        ]
        for future in futures:
            tile_result = future.result()
            if tile_result is not None:
                tile_results.append(tile_result)
    return tile_results


def print_tile_summary(tile_result: TileSkyline) -> None:
    farthest = tile_result.radius_results[-1].result
    print(
        f"[tile] {tile_result.envelope.path.name}  "
        f"buildings={farthest.buildings_considered}/{farthest.buildings_contributing}  "
        f"peak={farthest.peak_altitude_deg:.2f}deg@{farthest.peak_azimuth_deg:.1f}"
    )


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        tower = select_tower(args.tower)
        citygml_dir = args.citygml_dir
        if not citygml_dir.is_dir():
            raise ValueError(f"CityGML directory not found: {citygml_dir}")
        cumulative_radii_km = normalize_cumulative_radii_km(
            [float(value) for value in args.cumulative_radius_km],
            max_radius_km=float(args.radius_km),
        )
        envelopes = select_tile_envelopes(
            citygml_dir,
            observer_lat_deg=tower.latitude_deg,
            observer_lon_deg=tower.longitude_deg,
            radius_km=float(args.radius_km) + (float(args.radius_band_width_m) / 1000.0),
        )
        if args.print_selected_tiles:
            for envelope in envelopes:
                print(f"[select] {envelope.path.name}")
        tile_results = compute_tile_skylines(
            envelopes,
            tower=tower,
            cumulative_radii_km=cumulative_radii_km,
            radius_band_width_m=float(args.radius_band_width_m),
            azimuth_step_deg=float(args.azimuth_step),
            edge_sample_step_m=float(args.edge_sample_step_m),
            workers=int(args.workers),
        )
        if not tile_results:
            raise ValueError("No usable buildings were found in the selected CityGML tiles.")

        radius_results = combine_tile_results(
            tower,
            tile_results,
            radii_km=cumulative_radii_km,
            azimuth_step_deg=float(args.azimuth_step),
        )
        result = radius_results[-1].result
        stem = sanitize_slug(str(tower.meta.get("slug") or tower.name))
        png_path = args.output_dir / f"{stem}_urban.png"
        write_preview_png(png_path, result)
        print(
            f"[ok] {tower.name}: {png_path}  "
            f"tiles={len(tile_results)}  "
            f"buildings={result.buildings_considered}/{result.buildings_contributing}  "
            f"peak={result.peak_altitude_deg:.2f}deg@{result.peak_azimuth_deg:.1f}"
        )
        if args.write_json:
            json_output = args.json_output or (args.output_dir / "urban_skyline_profiles.json")
            write_profiles_json(json_output, ((tower, radius_results),))
            print(f"[ok] skyline-json: {json_output}")
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
