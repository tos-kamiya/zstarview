#!/usr/bin/env python3
"""Compare Overture download volume for two radii around one viewpoint."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.data.import_overture_buildings import bbox_from_point, build_download_command


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download Overture building data for two radii and compare the raw "
            "download size, elapsed time, and approximate feature count."
        )
    )
    parser.add_argument("--lat", type=float, required=True, help="Center latitude in degrees.")
    parser.add_argument("--lon", type=float, required=True, help="Center longitude in degrees.")
    parser.add_argument(
        "--radii-km",
        type=float,
        nargs=2,
        default=(2.5, 7.5),
        metavar=("SMALL", "LARGE"),
        help="Two radii in km to compare (default: 2.5 7.5).",
    )
    parser.add_argument(
        "--feature-type",
        choices=("building", "building_part", "both"),
        default="both",
        help="Overture feature type(s) to compare (default: both).",
    )
    parser.add_argument(
        "--format",
        default="geojsonseq",
        choices=("geojson", "geojsonseq", "geoparquet"),
        help="Output format passed to overturemaps download (default: geojsonseq).",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="How many times to repeat each download (default: 1).",
    )
    parser.add_argument(
        "--overturemaps-bin",
        default="overturemaps",
        help="Path to the overturemaps CLI executable (default: overturemaps).",
    )
    parser.add_argument(
        "--no-stac",
        action="store_true",
        help="Pass --no-stac to the overturemaps CLI.",
    )
    parser.add_argument(
        "--keep-downloads-dir",
        type=Path,
        default=None,
        help="Optional directory where raw downloads are copied for later inspection.",
    )
    return parser


def output_suffix_for_format(fmt: str) -> str:
    return {
        "geojson": ".geojson",
        "geojsonseq": ".geojsonseq",
        "geoparquet": ".parquet",
    }[fmt]


def format_size(num_bytes: int) -> str:
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    size = float(num_bytes)
    for unit in units:
        if size < 1024.0 or unit == units[-1]:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{num_bytes} B"


def feature_types_for_mode(mode: str) -> tuple[str, ...]:
    if mode == "both":
        return ("building", "building_part")
    return (mode,)


def count_features(download_path: Path, fmt: str) -> int | None:
    if fmt == "geojsonseq":
        count = 0
        with download_path.open("r", encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    count += 1
        return count
    if fmt == "geojson":
        payload = json.loads(download_path.read_text(encoding="utf-8"))
        features = payload.get("features")
        return len(features) if isinstance(features, list) else None
    return None


def run_single_download(
    *,
    lat: float,
    lon: float,
    radius_km: float,
    feature_type: str,
    fmt: str,
    overturemaps_bin: str,
    no_stac: bool,
    keep_downloads_dir: Path | None,
    repeat_index: int,
) -> dict[str, object]:
    bbox = bbox_from_point(lat, lon, radius_km)
    with tempfile.TemporaryDirectory(prefix="overture-radius-compare-") as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        output_path = tmp_dir / f"{feature_type}_r{radius_km:g}_{repeat_index}{output_suffix_for_format(fmt)}"
        command = build_download_command(
            overturemaps_bin=overturemaps_bin,
            bbox=bbox,
            feature_type=feature_type,
            fmt=fmt,
            output_path=output_path,
            no_stac=no_stac,
        )
        started_at = time.perf_counter()
        completed = subprocess.run(command, check=False)
        elapsed_s = time.perf_counter() - started_at
        if completed.returncode != 0:
            raise RuntimeError(
                f"overturemaps download failed for feature_type={feature_type} "
                f"radius_km={radius_km} repeat={repeat_index} returncode={completed.returncode}"
            )
        size_bytes = output_path.stat().st_size if output_path.exists() else 0
        feature_count = count_features(output_path, fmt)
        kept_path = None
        if keep_downloads_dir is not None:
            keep_downloads_dir.mkdir(parents=True, exist_ok=True)
            kept_path = keep_downloads_dir / output_path.name
            shutil.copyfile(output_path, kept_path)
        return {
            "radius_km": radius_km,
            "feature_type": feature_type,
            "repeat": repeat_index,
            "bbox": bbox,
            "elapsed_s": elapsed_s,
            "size_bytes": size_bytes,
            "feature_count": feature_count,
            "kept_path": str(kept_path) if kept_path is not None else None,
        }


def summarize_results(results: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[float, str], list[dict[str, object]]] = {}
    for result in results:
        key = (float(result["radius_km"]), str(result["feature_type"]))
        grouped.setdefault(key, []).append(result)

    summaries: list[dict[str, object]] = []
    for (radius_km, feature_type), items in sorted(grouped.items()):
        size_values = [int(item["size_bytes"]) for item in items]
        elapsed_values = [float(item["elapsed_s"]) for item in items]
        feature_values = [item["feature_count"] for item in items]
        known_counts = [int(value) for value in feature_values if value is not None]
        summaries.append(
            {
                "radius_km": radius_km,
                "feature_type": feature_type,
                "runs": len(items),
                "avg_size_bytes": sum(size_values) / len(size_values),
                "avg_elapsed_s": sum(elapsed_values) / len(elapsed_values),
                "avg_feature_count": (
                    sum(known_counts) / len(known_counts) if len(known_counts) == len(items) else None
                ),
                "bbox": items[0]["bbox"],
            }
        )
    return summaries


def print_summary_table(
    *,
    summaries: list[dict[str, object]],
    radii: tuple[float, float],
    feature_types: tuple[str, ...],
) -> None:
    print("")
    print("Summary")
    print("radius_km feature_type runs size_bytes size_human elapsed_s feature_count")
    summary_by_key = {
        (float(item["radius_km"]), str(item["feature_type"])): item
        for item in summaries
    }
    for radius_km in radii:
        for feature_type in feature_types:
            summary = summary_by_key.get((float(radius_km), feature_type))
            if summary is None:
                continue
            avg_size_bytes = int(round(float(summary["avg_size_bytes"])))
            avg_elapsed_s = float(summary["avg_elapsed_s"])
            feature_count = summary["avg_feature_count"]
            feature_label = "n/a" if feature_count is None else str(int(round(float(feature_count))))
            print(
                f"{radius_km:8.2f} {feature_type:12s} {int(summary['runs']):4d} "
                f"{avg_size_bytes:10d} {format_size(avg_size_bytes):>9s} "
                f"{avg_elapsed_s:9.3f} {feature_label:>13s}"
            )

    print("")
    print("Ratios")
    print("feature_type size_ratio elapsed_ratio feature_ratio")
    small_radius, large_radius = radii
    for feature_type in feature_types:
        small = summary_by_key.get((float(small_radius), feature_type))
        large = summary_by_key.get((float(large_radius), feature_type))
        if small is None or large is None:
            continue
        size_ratio = float(large["avg_size_bytes"]) / max(float(small["avg_size_bytes"]), 1.0)
        elapsed_ratio = float(large["avg_elapsed_s"]) / max(float(small["avg_elapsed_s"]), 1e-9)
        if small["avg_feature_count"] in (None, 0) or large["avg_feature_count"] is None:
            feature_ratio = "n/a"
        else:
            ratio = float(large["avg_feature_count"]) / float(small["avg_feature_count"])
            feature_ratio = f"{ratio:.2f}x"
        print(f"{feature_type:12s} {size_ratio:10.2f}x {elapsed_ratio:12.2f}x {feature_ratio:>12s}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    overturemaps_bin = shutil.which(args.overturemaps_bin)
    if overturemaps_bin is None:
        parser.error(
            f"Could not find overturemaps CLI: {args.overturemaps_bin!r}. "
            "Install it separately or pass --overturemaps-bin."
        )

    small_radius, large_radius = (float(value) for value in args.radii_km)
    if small_radius <= 0 or large_radius <= 0:
        parser.error("--radii-km values must be positive.")
    if large_radius <= small_radius:
        parser.error("--radii-km requires SMALL < LARGE.")
    if args.repeat <= 0:
        parser.error("--repeat must be positive.")

    feature_types = feature_types_for_mode(args.feature_type)
    keep_downloads_dir = args.keep_downloads_dir.resolve() if args.keep_downloads_dir is not None else None

    print(f"center_lat={args.lat:.8f}")
    print(f"center_lon={args.lon:.8f}")
    print(f"format={args.format}")
    print(f"feature_types={','.join(feature_types)}")
    print(f"repeat={args.repeat}")
    print(f"overturemaps_bin={overturemaps_bin}")

    results: list[dict[str, object]] = []
    for radius_km in (small_radius, large_radius):
        bbox = bbox_from_point(args.lat, args.lon, radius_km)
        print("")
        print(
            f"radius_km={radius_km:.2f} "
            f"bbox={bbox[0]:.8f},{bbox[1]:.8f},{bbox[2]:.8f},{bbox[3]:.8f}"
        )
        for feature_type in feature_types:
            for repeat_index in range(1, args.repeat + 1):
                result = run_single_download(
                    lat=float(args.lat),
                    lon=float(args.lon),
                    radius_km=radius_km,
                    feature_type=feature_type,
                    fmt=str(args.format),
                    overturemaps_bin=overturemaps_bin,
                    no_stac=bool(args.no_stac),
                    keep_downloads_dir=keep_downloads_dir,
                    repeat_index=repeat_index,
                )
                results.append(result)
                feature_count = result["feature_count"]
                feature_label = "n/a" if feature_count is None else str(feature_count)
                print(
                    f"  feature_type={feature_type:12s} repeat={repeat_index:02d} "
                    f"size={int(result['size_bytes']):10d} ({format_size(int(result['size_bytes'])):>9s}) "
                    f"elapsed_s={float(result['elapsed_s']):9.3f} feature_count={feature_label}"
                )
                if result["kept_path"] is not None:
                    print(f"    kept_path={result['kept_path']}")

    summaries = summarize_results(results)
    print_summary_table(
        summaries=summaries,
        radii=(small_radius, large_radius),
        feature_types=feature_types,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
