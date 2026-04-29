#!/usr/bin/env python3
"""Audit how much small geometry exists in derived urban-outline data.

This script scans one or more derived building directories, loads every
building footprint it can find, and reports the footprint-size distribution.

It is intended to answer questions like:
- How many buildings are represented by very small rings?
- How small are the smallest footprints in the dataset?
- How often does the data contain footprints below a chosen size threshold?

Example:
  python scripts/audit_urban_outline_shapes.py \
    --derived-dir ~/.cache/zstarview/overture_buildings \
    --threshold-m 5 10 20 50
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from pyproj import CRS, Transformer


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.data.derived_tile_cache import parse_derived_tile_buildings  # noqa: E402
from zstarview.paths import CACHE_PATH, OVERTURE_DERIVED_ROOT_DIR, OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR  # noqa: E402
from zstarview.data.urban_outline_common import BuildingFootprint  # noqa: E402


@dataclass(frozen=True)
class ShapeRecord:
    source_file: Path
    building_id: str
    height_m: float
    ring_count: int
    point_count: int
    max_dim_m: float
    width_m: float
    height_m_bbox: float
    area_m2: float
    perimeter_m: float


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit footprint sizes in derived urban-outline building tiles."
    )
    parser.add_argument(
        "--derived-dir",
        action="append",
        dest="derived_dirs",
        default=[],
        type=Path,
        help="Derived building directory to scan. May be repeated.",
    )
    parser.add_argument(
        "--threshold-m",
        action="append",
        dest="thresholds_m",
        default=[],
        type=float,
        help="Footprint size threshold in meters. May be repeated.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of smallest shapes to print in the ranked list.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON instead of text output.",
    )
    parser.add_argument(
        "--csv-out",
        type=Path,
        help="Write a histogram CSV with diameter bins to this path.",
    )
    parser.add_argument(
        "--bin-size-m",
        type=float,
        default=1.0,
        help="Histogram bin size in meters for CSV output.",
    )
    parser.add_argument(
        "--include-skyscrapers",
        action="store_true",
        help="Deprecated: the default scan already includes skyscraper caches when present.",
    )
    return parser


def _default_scan_dirs(include_skyscrapers: bool) -> tuple[Path, ...]:
    dirs: list[Path] = []
    for candidate in (
        Path(OVERTURE_DERIVED_ROOT_DIR),
        Path(OVERTURE_SKYSCRAPER_DERIVED_ROOT_DIR),
    ):
        if candidate.exists():
            dirs.append(candidate)
    if not dirs:
        cache_root = Path(CACHE_PATH)
        for candidate in (
            cache_root / "overture_buildings",
            cache_root / "overture_skyscrapers",
        ):
            if candidate.exists():
                dirs.append(candidate)
    if not include_skyscrapers:
        # The historical flag is kept for compatibility, but the default scan now
        # includes both cache roots when present.
        return tuple(dirs)
    return tuple(dirs)


def _iter_building_files(derived_dir: Path) -> Iterable[Path]:
    if derived_dir.is_file():
        yield derived_dir
        return
    if not derived_dir.exists():
        return
    for path in sorted(derived_dir.rglob("*.json")):
        if path.name == "tile_index.json" or path.name.endswith("cache_meta.json"):
            continue
        yield path


def _select_rings(building: BuildingFootprint) -> tuple[tuple[tuple[float, float], ...], ...]:
    return tuple(ring for ring in building.rings_lonlat if len(ring) >= 4)


def _project_rings(rings: Sequence[Sequence[tuple[float, float]]]) -> tuple[np.ndarray, np.ndarray]:
    lon_points: list[float] = []
    lat_points: list[float] = []
    for ring in rings:
        for lon, lat in ring:
            lon_points.append(float(lon))
            lat_points.append(float(lat))
    if not lon_points:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)

    center_lon = float(np.mean(lon_points))
    center_lat = float(np.mean(lat_points))
    crs_local = CRS.from_proj4(
        f"+proj=aeqd +lat_0={center_lat} +lon_0={center_lon} +datum=WGS84 +units=m +no_defs"
    )
    transformer = Transformer.from_crs("EPSG:4326", crs_local, always_xy=True)
    xs, ys = transformer.transform(lon_points, lat_points)
    return np.asarray(xs, dtype=np.float64), np.asarray(ys, dtype=np.float64)


def _measure_building_shape(building: BuildingFootprint, *, source_file: Path) -> ShapeRecord | None:
    rings = _select_rings(building)
    if not rings:
        return None

    xs, ys = _project_rings(rings)
    if xs.size == 0 or ys.size == 0:
        return None

    width_m = float(np.max(xs) - np.min(xs))
    height_m_bbox = float(np.max(ys) - np.min(ys))
    max_dim_m = max(width_m, height_m_bbox)
    area_m2 = 0.0
    perimeter_m = 0.0
    point_count = 0

    for ring in rings:
        ring_xs, ring_ys = _project_rings((ring,))
        if ring_xs.size < 2 or ring_ys.size < 2:
            continue
        point_count += len(ring)
        area_m2 += 0.5 * abs(
            float(
                np.sum(
                    ring_xs[:-1] * ring_ys[1:]
                    - ring_xs[1:] * ring_ys[:-1]
                )
            )
        )
        perimeter_m += float(
            np.sum(np.hypot(ring_xs[1:] - ring_xs[:-1], ring_ys[1:] - ring_ys[:-1]))
        )

    return ShapeRecord(
        source_file=source_file,
        building_id=building.building_id,
        height_m=float(building.height_m),
        ring_count=len(rings),
        point_count=point_count,
        max_dim_m=max_dim_m,
        width_m=width_m,
        height_m_bbox=height_m_bbox,
        area_m2=area_m2,
        perimeter_m=perimeter_m,
    )


def _load_shape_records(derived_dirs: Sequence[Path]) -> tuple[ShapeRecord, ...]:
    records: list[ShapeRecord] = []
    for derived_dir in derived_dirs:
        for file_path in _iter_building_files(derived_dir):
            for building in parse_derived_tile_buildings(file_path):
                record = _measure_building_shape(building, source_file=file_path)
                if record is not None:
                    records.append(record)
    return tuple(records)


def _threshold_counts(records: Sequence[ShapeRecord], thresholds_m: Sequence[float]) -> dict[str, int]:
    counts: dict[str, int] = {}
    values = np.asarray([record.max_dim_m for record in records], dtype=np.float64)
    for threshold in thresholds_m:
        key = f"<= {threshold:g} m"
        counts[key] = int(np.sum(values <= float(threshold)))
    return counts


def _histogram_rows(records: Sequence[ShapeRecord], bin_size_m: float) -> list[dict[str, float | int]]:
    if bin_size_m <= 0.0:
        raise ValueError("bin_size_m must be positive")
    if not records:
        return []

    values = np.asarray([record.max_dim_m for record in records], dtype=np.float64)
    max_bin_index = int(np.floor(np.max(values) / bin_size_m))
    counts = np.zeros(max_bin_index + 1, dtype=np.int64)
    for value in values:
        index = int(np.floor(float(value) / bin_size_m))
        if index < 0:
            index = 0
        if index > max_bin_index:
            index = max_bin_index
        counts[index] += 1

    rows: list[dict[str, float | int]] = []
    for bin_index, count in enumerate(counts.tolist()):
        bin_start_m = float(bin_index * bin_size_m)
        rows.append(
            {
                "bin_start_m": bin_start_m,
                "bin_end_m": bin_start_m + float(bin_size_m),
                "count": int(count),
            }
        )
    return rows


def _write_csv_histogram(records: Sequence[ShapeRecord], csv_path: Path, bin_size_m: float) -> None:
    rows = _histogram_rows(records, bin_size_m)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=("bin_start_m", "bin_end_m", "count"))
        writer.writeheader()
        writer.writerows(rows)


def _format_text(records: Sequence[ShapeRecord], thresholds_m: Sequence[float], top_n: int) -> str:
    if not records:
        return "No building shapes were found."

    values = np.asarray([record.max_dim_m for record in records], dtype=np.float64)
    lines = [
        f"shapes: {len(records)}",
        f"min_max_dim_m: {float(np.min(values)):.3f}",
        f"p10_max_dim_m: {float(np.percentile(values, 10)):.3f}",
        f"p50_max_dim_m: {float(np.percentile(values, 50)):.3f}",
        f"p90_max_dim_m: {float(np.percentile(values, 90)):.3f}",
        f"p99_max_dim_m: {float(np.percentile(values, 99)):.3f}",
    ]
    if thresholds_m:
        lines.append("threshold_counts:")
        for label, count in _threshold_counts(records, thresholds_m).items():
            lines.append(f"  {label}: {count}")

    lines.append(f"smallest_{max(1, int(top_n))}:")
    for record in sorted(records, key=lambda item: (item.max_dim_m, item.height_m, item.building_id))[: max(1, int(top_n))]:
        lines.append(
            "  "
            + f"{record.max_dim_m:.3f} m | h={record.height_m:.1f} m | "
            + f"w={record.width_m:.3f} m | hbox={record.height_m_bbox:.3f} m | "
            + f"rings={record.ring_count} | points={record.point_count} | {record.building_id}"
        )
    return "\n".join(lines)


def _format_json(records: Sequence[ShapeRecord], thresholds_m: Sequence[float], top_n: int) -> str:
    values = np.asarray([record.max_dim_m for record in records], dtype=np.float64)
    payload: dict[str, object] = {
        "shape_count": len(records),
        "threshold_counts": _threshold_counts(records, thresholds_m) if thresholds_m else {},
        "summary": {},
        "smallest": [],
    }
    if records:
        payload["summary"] = {
            "min_max_dim_m": float(np.min(values)),
            "p10_max_dim_m": float(np.percentile(values, 10)),
            "p50_max_dim_m": float(np.percentile(values, 50)),
            "p90_max_dim_m": float(np.percentile(values, 90)),
            "p99_max_dim_m": float(np.percentile(values, 99)),
        }
        payload["smallest"] = [
            {
                "building_id": record.building_id,
                "source_file": str(record.source_file),
                "height_m": record.height_m,
                "max_dim_m": record.max_dim_m,
                "width_m": record.width_m,
                "height_bbox_m": record.height_m_bbox,
                "area_m2": record.area_m2,
                "perimeter_m": record.perimeter_m,
                "ring_count": record.ring_count,
                "point_count": record.point_count,
            }
            for record in sorted(records, key=lambda item: (item.max_dim_m, item.height_m, item.building_id))[: max(1, int(top_n))]
        ]
    return json.dumps(payload, ensure_ascii=False, indent=2)


def main() -> int:
    args = _build_parser().parse_args()

    derived_dirs = tuple(args.derived_dirs) or _default_scan_dirs(include_skyscrapers=bool(args.include_skyscrapers))
    if not derived_dirs:
        print("No derived directories found.", file=sys.stderr)
        return 2

    records = _load_shape_records(derived_dirs)
    thresholds_m = tuple(float(value) for value in args.thresholds_m if float(value) > 0.0)
    top_n = max(1, int(args.top_n))

    if args.csv_out is not None:
        _write_csv_histogram(records, args.csv_out, float(args.bin_size_m))
        print(f"Wrote histogram CSV to {args.csv_out}")
        return 0

    if args.json:
        print(_format_json(records, thresholds_m, top_n))
    else:
        print(f"scan_dirs: {', '.join(str(path) for path in derived_dirs)}")
        print(_format_text(records, thresholds_m, top_n))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
