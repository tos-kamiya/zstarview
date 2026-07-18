#!/usr/bin/env python3
"""Probe Overture Transportation road segments around a point.

The probe downloads a small GeoJSON segment extract with the official
``overturemaps`` CLI and reports how many road segments contain ``width_rules``.
It is intentionally independent of the zstarview runtime cache and renderer.

Example:
    uv run -p .venv/bin/python dev-samples/overture_transportation_probe.py \
        --lat 35.681236 --lon 139.767125
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any


EARTH_RADIUS_KM = 6371.0088


def bbox_from_point(
    lat_deg: float, lon_deg: float, radius_km: float
) -> tuple[float, float, float, float]:
    if radius_km <= 0.0:
        raise ValueError("radius_km must be positive")
    if not -90.0 <= lat_deg <= 90.0:
        raise ValueError("lat must be between -90 and 90 degrees")
    if not -180.0 <= lon_deg <= 180.0:
        raise ValueError("lon must be between -180 and 180 degrees")

    angular_distance = radius_km / EARTH_RADIUS_KM
    lat_delta_deg = math.degrees(angular_distance)
    cos_lat = math.cos(math.radians(lat_deg))
    lon_delta_deg = (
        180.0 if abs(cos_lat) < 1e-9 else math.degrees(angular_distance / abs(cos_lat))
    )
    return (
        max(-180.0, lon_deg - lon_delta_deg),
        max(-90.0, lat_deg - lat_delta_deg),
        min(180.0, lon_deg + lon_delta_deg),
        min(90.0, lat_deg + lat_delta_deg),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Probe Overture road segments and their width_rules."
    )
    parser.add_argument(
        "--lat", type=float, required=True, help="Center latitude in degrees."
    )
    parser.add_argument(
        "--lon", type=float, required=True, help="Center longitude in degrees."
    )
    parser.add_argument(
        "--radius-km",
        type=float,
        default=0.05,
        help="Radius used to build the download bbox (default: 0.05).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Keep the downloaded GeoJSON at this path instead of using a temporary file.",
    )
    parser.add_argument(
        "--overturemaps-bin",
        default="overturemaps",
        help="Path to the overturemaps CLI executable (default: overturemaps).",
    )
    parser.add_argument(
        "--release",
        default=None,
        help="Optional Overture release identifier passed to the CLI.",
    )
    parser.add_argument(
        "--no-stac",
        action="store_true",
        help="Pass --no-stac to the overturemaps CLI.",
    )
    return parser


def _width_values(width_rules: Any) -> list[float]:
    if not isinstance(width_rules, list):
        return []
    values: list[float] = []
    for rule in width_rules:
        if not isinstance(rule, dict):
            continue
        value = rule.get("value")
        if (
            isinstance(value, (int, float))
            and math.isfinite(float(value))
            and value > 0
        ):
            values.append(float(value))
    return values


def _features(payload: Any) -> list[dict[str, Any]]:
    if not isinstance(payload, dict):
        raise ValueError("GeoJSON root is not an object")
    features = payload.get("features")
    if not isinstance(features, list):
        raise ValueError("GeoJSON payload has no features array")
    return [feature for feature in features if isinstance(feature, dict)]


def _download(
    args: argparse.Namespace, output_path: Path, bbox: tuple[float, float, float, float]
) -> None:
    west, south, east, north = bbox
    command = [
        args.overturemaps_bin,
        "download",
        f"--bbox={west:.8f},{south:.8f},{east:.8f},{north:.8f}",
        "--type",
        "segment",
        "-f",
        "geojson",
        "-o",
        str(output_path),
    ]
    if args.release:
        command.extend(["--release", args.release])
    if args.no_stac:
        command.append("--no-stac")
    print("command=" + " ".join(command))
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Overture CLI not found: {args.overturemaps_bin}") from exc


def main() -> int:
    args = build_parser().parse_args()
    bbox = bbox_from_point(args.lat, args.lon, args.radius_km)
    print(f"center_lat={args.lat:.8f}")
    print(f"center_lon={args.lon:.8f}")
    print(f"radius_km={args.radius_km:.6f}")
    print("bbox=" + ",".join(f"{value:.8f}" for value in bbox))

    temporary_path: Path | None = None
    output_path = args.output
    if output_path is None:
        temporary_file = tempfile.NamedTemporaryFile(
            prefix="overture-transportation-probe-",
            suffix=".geojson",
            delete=False,
        )
        temporary_file.close()
        temporary_path = Path(temporary_file.name)
        output_path = temporary_path
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        _download(args, output_path, bbox)
        payload = json.loads(output_path.read_text(encoding="utf-8"))
        features = _features(payload)
        roads = [
            feature
            for feature in features
            if isinstance(feature.get("properties"), dict)
            and feature["properties"].get("subtype") == "road"
        ]
        width_counts = Counter()
        width_values: list[float] = []
        for road in roads:
            properties = road["properties"]
            rules = properties.get("width_rules")
            values = _width_values(rules)
            if values:
                width_counts["with_width_rules"] += 1
                width_values.extend(values)
            else:
                width_counts["without_width_rules"] += 1

        print(f"features={len(features)}")
        print(f"road_segments={len(roads)}")
        print(f"road_segments_with_width_rules={width_counts['with_width_rules']}")
        print(
            f"road_segments_without_width_rules={width_counts['without_width_rules']}"
        )
        if roads:
            coverage = 100.0 * width_counts["with_width_rules"] / len(roads)
            print(f"width_rule_coverage_percent={coverage:.2f}")
        else:
            print("width_rule_coverage_percent=0.00")
        if width_values:
            print(f"width_rule_value_count={len(width_values)}")
            print(f"width_rule_min_m={min(width_values):.3f}")
            print(f"width_rule_max_m={max(width_values):.3f}")
            print(
                "width_rule_values_m="
                + ",".join(f"{value:.3f}" for value in width_values[:20])
            )
        else:
            print("width_rule_value_count=0")

        classes = Counter(
            str(road["properties"].get("class", "unknown"))
            for road in roads
            if _width_values(road["properties"].get("width_rules"))
        )
        print(
            "classes_with_width_rules="
            + json.dumps(dict(sorted(classes.items())), sort_keys=True)
        )
        if args.output is not None:
            print(f"output={output_path}")
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
