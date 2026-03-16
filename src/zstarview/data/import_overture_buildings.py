#!/usr/bin/env python3
"""Download Overture building data around a point and store derived tiles."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

from zstarview.data import build_derived_tile_index

EARTH_RADIUS_KM = 6371.0088
DEFAULT_FETCH_RADIUS_KM = 2.5
DEFAULT_MIN_BUILDING_HEIGHT_M = 0.0
DEFAULT_STOREY_HEIGHT_M = 3.5


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download Overture building data around a lat/lon point and write "
            "derived urban-outline tiles into the app cache directory."
        )
    )
    parser.add_argument("--lat", type=float, required=True, help="Center latitude in degrees.")
    parser.add_argument("--lon", type=float, required=True, help="Center longitude in degrees.")
    parser.add_argument(
        "--radius-km",
        type=float,
        default=DEFAULT_FETCH_RADIUS_KM,
        help="Download radius in kilometers used to build the bbox (default: 2.5).",
    )
    parser.add_argument(
        "--derived-root-dir",
        type=Path,
        default=None,
        help="Root directory where derived datasets will be written.",
    )
    parser.add_argument(
        "--min-building-height-m",
        type=float,
        default=DEFAULT_MIN_BUILDING_HEIGHT_M,
        help="Ignore buildings lower than this height in meters.",
    )
    parser.add_argument(
        "--feature-type",
        default="building",
        choices=("building", "building_part"),
        help="Overture feature type to download (default: building).",
    )
    parser.add_argument(
        "--format",
        default="geojsonseq",
        choices=("geojson", "geojsonseq"),
        help="Intermediate download format passed to overturemaps (default: geojsonseq).",
    )
    parser.add_argument(
        "--overturemaps-bin",
        default="overturemaps",
        help="Path to the overturemaps CLI executable (default: overturemaps).",
    )
    parser.add_argument(
        "--dataset-name",
        default=None,
        help="Optional dataset directory name. Defaults to a lat/lon-based cache key.",
    )
    parser.add_argument(
        "--keep-download",
        type=Path,
        default=None,
        help="Optional path where the raw downloaded file will be copied for inspection.",
    )
    parser.add_argument(
        "--no-stac",
        action="store_true",
        help="Pass --no-stac to the overturemaps CLI.",
    )
    return parser


def bbox_from_point(lat_deg: float, lon_deg: float, radius_km: float) -> tuple[float, float, float, float]:
    if radius_km <= 0:
        raise ValueError("radius_km must be positive")
    if not -90.0 <= lat_deg <= 90.0:
        raise ValueError("lat must be between -90 and 90 degrees")
    if not -180.0 <= lon_deg <= 180.0:
        raise ValueError("lon must be between -180 and 180 degrees")

    angular_distance = radius_km / EARTH_RADIUS_KM
    lat_delta_deg = math.degrees(angular_distance)
    cos_lat = math.cos(math.radians(lat_deg))
    if abs(cos_lat) < 1e-9:
        lon_delta_deg = 180.0
    else:
        lon_delta_deg = math.degrees(angular_distance / abs(cos_lat))

    return (
        max(-180.0, lon_deg - lon_delta_deg),
        max(-90.0, lat_deg - lat_delta_deg),
        min(180.0, lon_deg + lon_delta_deg),
        min(90.0, lat_deg + lat_delta_deg),
    )


def derive_dataset_name(
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    feature_type: str,
    min_building_height_m: float = DEFAULT_MIN_BUILDING_HEIGHT_M,
) -> str:
    return (
        f"overture_{feature_type}"
        f"_lat{lat_deg:.4f}_lon{lon_deg:.4f}_r{radius_km:.1f}km"
        f"_h{float(min_building_height_m):.1f}m"
    ).replace("-", "m").replace(".", "p")


def output_suffix_for_format(fmt: str) -> str:
    return {
        "geojson": ".geojson",
        "geojsonseq": ".geojsonseq",
    }[fmt]


def build_download_command(
    *,
    overturemaps_bin: str,
    bbox: tuple[float, float, float, float],
    feature_type: str,
    fmt: str,
    output_path: Path,
    no_stac: bool,
) -> list[str]:
    west, south, east, north = bbox
    command = [
        overturemaps_bin,
        "download",
        f"--bbox={west:.8f},{south:.8f},{east:.8f},{north:.8f}",
        "-f",
        fmt,
        "--type",
        feature_type,
        "-o",
        str(output_path),
    ]
    if no_stac:
        command.append("--no-stac")
    return command


def import_overture_buildings(
    *,
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
    derived_root_dir: Path,
    min_building_height_m: float,
    feature_type: str,
    fmt: str,
    overturemaps_bin: str,
    dataset_name: str | None,
    keep_download: Path | None,
    no_stac: bool,
) -> Path:
    overturemaps_path = shutil.which(overturemaps_bin)
    if overturemaps_path is None:
        raise FileNotFoundError(
            f"Could not find overturemaps CLI: {overturemaps_bin!r}. "
            "Install it separately or pass --overturemaps-bin."
        )

    bbox = bbox_from_point(lat_deg, lon_deg, radius_km)
    dataset_dir_name = dataset_name or derive_dataset_name(
        lat_deg,
        lon_deg,
        radius_km,
        feature_type,
        min_building_height_m,
    )
    derived_dir = derived_root_dir / dataset_dir_name / "bldg"
    tile_path = derived_dir / f"{dataset_dir_name}.json"

    with tempfile.TemporaryDirectory(prefix="overture-import-") as temp_dir_str:
        temp_dir = Path(temp_dir_str)
        download_path = temp_dir / f"download{output_suffix_for_format(fmt)}"
        command = build_download_command(
            overturemaps_bin=overturemaps_path,
            bbox=bbox,
            feature_type=feature_type,
            fmt=fmt,
            output_path=download_path,
            no_stac=no_stac,
        )
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"overturemaps download failed with return code {completed.returncode}")
        if keep_download is not None:
            keep_download.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(download_path, keep_download)

        payload = build_derived_tile_payload(
            download_path,
            dataset_dir_name=dataset_dir_name,
            bbox=bbox,
            min_building_height_m=min_building_height_m,
            feature_type=feature_type,
            query_lat_deg=lat_deg,
            query_lon_deg=lon_deg,
            query_radius_km=radius_km,
            download_format=fmt,
        )

    tile_path.parent.mkdir(parents=True, exist_ok=True)
    tile_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    build_derived_tile_index.main(["--derived-dir", str(derived_dir)])
    return derived_dir


def build_derived_tile_payload(
    download_path: Path,
    *,
    dataset_dir_name: str,
    bbox: tuple[float, float, float, float],
    min_building_height_m: float,
    feature_type: str,
    query_lat_deg: float,
    query_lon_deg: float,
    query_radius_km: float,
    download_format: str,
) -> dict[str, object]:
    buildings: list[dict[str, object]] = []

    for feature in iter_download_features(download_path, fmt=download_format):
        building = convert_feature_to_building(feature, min_building_height_m=min_building_height_m)
        if building is not None:
            buildings.append(building)

    west, south, east, north = bbox
    return {
        "schema_version": 1,
        "source": {
            "format": "Overture-Buildings",
            "feature_type": feature_type,
            "query": {
                "lat_deg": query_lat_deg,
                "lon_deg": query_lon_deg,
                "radius_km": query_radius_km,
                "format": download_format,
            },
        },
        "tile": {
            "id": dataset_dir_name,
            "bbox": {
                "min_lat": south,
                "min_lon": west,
                "max_lat": north,
                "max_lon": east,
            },
        },
        "filters": {
            "min_height_m": float(min_building_height_m),
            "storey_height_m": DEFAULT_STOREY_HEIGHT_M,
        },
        "buildings": buildings,
    }


def iter_download_features(download_path: Path, *, fmt: str) -> Sequence[dict[str, object]]:
    if fmt == "geojsonseq":
        features: list[dict[str, object]] = []
        with download_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                row = json.loads(line)
                if isinstance(row, dict):
                    features.append(row)
        return tuple(features)

    payload = json.loads(download_path.read_text(encoding="utf-8"))
    if isinstance(payload, dict) and isinstance(payload.get("features"), list):
        return tuple(row for row in payload["features"] if isinstance(row, dict))
    return ()


def convert_feature_to_building(
    feature: dict[str, object],
    *,
    min_building_height_m: float,
) -> dict[str, object] | None:
    geometry = feature.get("geometry")
    if not isinstance(geometry, dict):
        return None
    rings = geometry_to_rings(geometry)
    if not rings:
        return None

    properties = feature.get("properties")
    if not isinstance(properties, dict):
        properties = {}
    height_m = resolve_feature_height_m(properties)
    if height_m is None or height_m < float(min_building_height_m):
        return None

    bbox = compute_building_bbox(rings)
    building_id = feature.get("id")
    if not isinstance(building_id, str) or not building_id:
        building_id = str(properties.get("id") or f"feature-{len(rings)}")

    return {
        "id": building_id,
        "height_m": height_m,
        "height_source": detect_height_source(properties),
        "parent_building_id": resolve_parent_building_id(properties),
        "bbox": {
            "min_lat": bbox[0],
            "min_lon": bbox[1],
            "max_lat": bbox[2],
            "max_lon": bbox[3],
        },
        "rings": [
            [[lon, lat] for lon, lat in ring]
            for ring in rings
        ],
    }


def geometry_to_rings(geometry: dict[str, object]) -> tuple[tuple[tuple[float, float], ...], ...]:
    geom_type = geometry.get("type")
    coordinates = geometry.get("coordinates")
    if geom_type == "Polygon":
        return parse_polygon_rings(coordinates)
    if geom_type == "MultiPolygon":
        if not isinstance(coordinates, list):
            return ()
        rings: list[tuple[tuple[float, float], ...]] = []
        for polygon in coordinates:
            rings.extend(parse_polygon_rings(polygon))
        return tuple(rings)
    return ()


def parse_polygon_rings(raw_polygon: object) -> tuple[tuple[tuple[float, float], ...], ...]:
    if not isinstance(raw_polygon, list):
        return ()
    rings: list[tuple[tuple[float, float], ...]] = []
    for raw_ring in raw_polygon:
        ring = parse_ring(raw_ring)
        if ring:
            rings.append(ring)
    return tuple(rings)


def parse_ring(raw_ring: object) -> tuple[tuple[float, float], ...]:
    if not isinstance(raw_ring, list):
        return ()
    points: list[tuple[float, float]] = []
    for raw_point in raw_ring:
        if not isinstance(raw_point, list) or len(raw_point) < 2:
            continue
        try:
            lon = float(raw_point[0])
            lat = float(raw_point[1])
        except (TypeError, ValueError):
            continue
        points.append((lon, lat))
    if len(points) < 4:
        return ()
    if points[0] != points[-1]:
        points.append(points[0])
    return tuple(points)


def resolve_feature_height_m(properties: dict[str, object]) -> float | None:
    for key in ("height", "roof_height"):
        value = parse_numeric(properties.get(key))
        if value is not None and value > 0:
            return value
    for key in ("num_floors", "num_floors_above_ground"):
        value = parse_numeric(properties.get(key))
        if value is not None and value > 0:
            return value * DEFAULT_STOREY_HEIGHT_M
    return None


def detect_height_source(properties: dict[str, object]) -> str:
    if parse_numeric(properties.get("height")) is not None:
        return "height"
    if parse_numeric(properties.get("roof_height")) is not None:
        return "roof_height"
    for key in ("num_floors", "num_floors_above_ground"):
        if parse_numeric(properties.get(key)) is not None:
            return f"{key}*{DEFAULT_STOREY_HEIGHT_M}"
    return "unknown"


def resolve_parent_building_id(properties: dict[str, object]) -> str | None:
    parent_building_id = properties.get("building_id")
    if isinstance(parent_building_id, str) and parent_building_id:
        return parent_building_id
    return None


def parse_numeric(value: object) -> float | None:
    try:
        if value is None:
            return None
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def compute_building_bbox(
    rings_lonlat: Sequence[Sequence[tuple[float, float]]],
) -> tuple[float, float, float, float]:
    min_lat = float("inf")
    min_lon = float("inf")
    max_lat = float("-inf")
    max_lon = float("-inf")
    for ring in rings_lonlat:
        for lon, lat in ring:
            min_lat = min(min_lat, lat)
            min_lon = min(min_lon, lon)
            max_lat = max(max_lat, lat)
            max_lon = max(max_lon, lon)
    return (min_lat, min_lon, max_lat, max_lon)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    from zstarview.paths import OVERTURE_DERIVED_ROOT_DIR

    derived_root_dir = args.derived_root_dir or Path(OVERTURE_DERIVED_ROOT_DIR)
    derived_dir = import_overture_buildings(
        lat_deg=float(args.lat),
        lon_deg=float(args.lon),
        radius_km=float(args.radius_km),
        derived_root_dir=derived_root_dir,
        min_building_height_m=float(args.min_building_height_m),
        feature_type=args.feature_type,
        fmt=args.format,
        overturemaps_bin=args.overturemaps_bin,
        dataset_name=args.dataset_name,
        keep_download=args.keep_download,
        no_stac=bool(args.no_stac),
    )
    print(f"[ok] imported overture buildings -> {derived_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
