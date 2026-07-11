#!/usr/bin/env python3
"""Build zstarview night-light tiles from an EOG VNL annual GeoTIFF.

The EOG VNL v2.2 annual download is a global EPSG:4326 GeoTIFF, commonly
distributed as a gzip-compressed file. This script keeps the source file out
of the repository and creates eight 90-degree GeoTIFF tiles plus a manifest.

Example:

    python dev-samples/build_vnl_night_lights.py \
        --source VNL_npp_2025_global_vcmslcfg_v2_c202604011200.average_masked.dat.tif.gz \
        --year 2025 \
        --output-dir build/night-lights/2025

The script can also download a source URL. If the server requires a bearer
token, pass the name of an environment variable with --bearer-token-env.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window


TILE_NAMES = ("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2")
TILE_WIDTH_DEG = 90.0
TILE_HEIGHT_DEG = 90.0
EXPECTED_RESOLUTION_DEG = 15.0 / 3600.0
CHUNK_ROWS = 1024
CHUNK_COLS = 4096


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build eight zstarview GeoTIFF tiles from an EOG VNL GeoTIFF."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--source", type=Path, help="Local .tif or .tif.gz source file.")
    source.add_argument("--url", help="URL of a .tif or .tif.gz source file.")
    parser.add_argument(
        "--bearer-token-env",
        help="Environment variable containing a bearer token for --url.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for the generated tiles and manifest.",
    )
    parser.add_argument("--year", type=int, default=2025)
    parser.add_argument(
        "--dataset-label",
        default="EOG VIIRS Nighttime Lights VNL v2.2 annual",
    )
    parser.add_argument(
        "--band",
        type=int,
        default=1,
        help="1-based source band to export (default: 1).",
    )
    parser.add_argument(
        "--source-url",
        help="Canonical source URL to record in manifest.json.",
    )
    parser.add_argument(
        "--asset-base-url",
        default=(
            "https://github.com/tos-kamiya/zstarview/releases/download/"
            "night-lights-2025/"
        ),
        help="Base URL used for tile asset URLs in manifest.json.",
    )
    parser.add_argument(
        "--conversion-commit",
        help="Git commit that produced the release assets.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_source(url: str, destination: Path, bearer_token: str | None) -> None:
    headers = {"User-Agent": "zstarview-vnl-builder/1.0"}
    if bearer_token:
        headers["Authorization"] = f"Bearer {bearer_token}"
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=300) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)


def materialize_source(source: Path, temporary_dir: Path) -> Path:
    if source.suffix.lower() != ".gz":
        return source
    materialized = temporary_dir / source.with_suffix("").name
    with gzip.open(source, "rb") as compressed, materialized.open("wb") as raw:
        shutil.copyfileobj(compressed, raw, length=1024 * 1024)
    return materialized


def materialize_url(url: str, temporary_dir: Path, bearer_token: str | None) -> Path:
    suffix = ".tif.gz" if url.lower().endswith((".gz", ".tif.gz")) else ".tif"
    downloaded = temporary_dir / f"source{suffix}"
    download_source(url, downloaded, bearer_token)
    return materialize_source(downloaded, temporary_dir)


def validate_source(dataset: rasterio.io.DatasetReader, band: int) -> None:
    if dataset.crs is None or dataset.crs.to_epsg() != 4326:
        raise ValueError(f"source CRS must be EPSG:4326, got {dataset.crs!s}")
    if not 1 <= band <= dataset.count:
        raise ValueError(f"source band must be in 1..{dataset.count}, got {band}")
    resolution_x, resolution_y = dataset.res
    tolerance = EXPECTED_RESOLUTION_DEG * 1.0e-5
    if abs(resolution_x - EXPECTED_RESOLUTION_DEG) > tolerance:
        raise ValueError(
            "source resolution must be 15 arc-seconds; "
            f"got {resolution_x:.12f} degrees"
        )
    if abs(resolution_y - EXPECTED_RESOLUTION_DEG) > tolerance:
        raise ValueError(
            "source resolution must be 15 arc-seconds; "
            f"got {resolution_y:.12f} degrees"
        )


def tile_bounds(tile_name: str) -> tuple[float, float, float, float]:
    column = "ABCD".index(tile_name[0])
    row = int(tile_name[1]) - 1
    left = -180.0 + column * TILE_WIDTH_DEG
    right = left + TILE_WIDTH_DEG
    top = 90.0 - row * TILE_HEIGHT_DEG
    bottom = top - TILE_HEIGHT_DEG
    return left, bottom, right, top


def copy_intersection(
    source: rasterio.io.DatasetReader,
    destination: rasterio.io.DatasetWriter,
    tile_bounds_value: tuple[float, float, float, float],
    band: int,
) -> None:
    tile_left, tile_bottom, tile_right, tile_top = tile_bounds_value
    source_left, source_bottom, source_right, source_top = source.bounds
    left = max(tile_left, source_left)
    bottom = max(tile_bottom, source_bottom)
    right = min(tile_right, source_right)
    top = min(tile_top, source_top)
    if left >= right or bottom >= top:
        return

    source_row = max(0, int(round((source_top - top) / source.res[1])))
    source_col = max(0, int(round((left - source_left) / source.res[0])))
    source_bottom_row = min(
        source.height,
        int(round((source_top - bottom) / source.res[1])),
    )
    source_right_col = min(
        source.width,
        int(round((right - source_left) / source.res[0])),
    )
    source_window = Window(
        source_col,
        source_row,
        max(0, source_right_col - source_col),
        max(0, source_bottom_row - source_row),
    )
    destination_col = int(round((left - tile_left) / source.res[0]))
    destination_row = int(round((tile_top - top) / source.res[1]))

    for row_offset in range(0, int(source_window.height), CHUNK_ROWS):
        rows = min(CHUNK_ROWS, int(source_window.height) - row_offset)
        for col_offset in range(0, int(source_window.width), CHUNK_COLS):
            cols = min(CHUNK_COLS, int(source_window.width) - col_offset)
            read_window = Window(
                source_window.col_off + col_offset,
                source_window.row_off + row_offset,
                cols,
                rows,
            )
            write_window = Window(
                destination_col + col_offset,
                destination_row + row_offset,
                cols,
                rows,
            )
            values = source.read(band, window=read_window, masked=True).filled(0)
            values = np.maximum(values, 0)
            destination.write(values, 1, window=write_window)


def write_tile(
    source: rasterio.io.DatasetReader,
    output_path: Path,
    tile_name: str,
    band: int,
) -> None:
    left, bottom, right, top = tile_bounds(tile_name)
    width = int(round(TILE_WIDTH_DEG / source.res[0]))
    height = int(round(TILE_HEIGHT_DEG / source.res[1]))
    transform = rasterio.transform.from_origin(left, top, source.res[0], source.res[1])
    profile = source.profile.copy()
    profile.update(
        driver="GTiff",
        count=1,
        width=width,
        height=height,
        dtype=source.dtypes[band - 1],
        crs="EPSG:4326",
        transform=transform,
        nodata=0,
        compress="deflate",
        predictor=2,
        tiled=True,
        blockxsize=256,
        blockysize=256,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(output_path, "w", **profile) as destination:
        copy_intersection(source, destination, (left, bottom, right, top), band)


def build_manifest(
    output_dir: Path,
    source_label: str,
    source_url: str | None,
    year: int,
    band: int,
    source_filename: str,
    source_sha256: str,
    asset_base_url: str,
    conversion_commit: str | None,
) -> None:
    files = {}
    for tile_name in TILE_NAMES:
        path = output_dir / f"{tile_name}.tif"
        with rasterio.open(path) as dataset:
            files[tile_name] = {
                "path": path.name,
                "url": asset_base_url.rstrip("/") + "/" + path.name,
                "sha256": sha256_file(path),
                "width": dataset.width,
                "height": dataset.height,
                "bounds": list(dataset.bounds),
                "crs": dataset.crs.to_string() if dataset.crs else None,
                "resolution_degrees": list(dataset.res),
                "dtype": dataset.dtypes[0],
                "band": 1,
            }
    manifest = {
        "dataset_version": "2025_vnl_v22_average_masked",
        "dataset": source_label,
        "year": year,
        "asset_base_url": asset_base_url,
        "source_url": source_url,
        "source_filename": source_filename,
        "source_sha256": source_sha256,
        "source_band": band,
        "conversion_script": "dev-samples/build_vnl_night_lights.py",
        "conversion_commit": conversion_commit,
        "derived_format": "GeoTIFF, eight 90-degree EPSG:4326 tiles",
        "license": "CC BY 4.0",
        "attribution": (
            "This product was made utilizing VIIRS nighttime lights data "
            "produced by the Earth Observation Group, Payne Institute for "
            "Public Policy, Colorado School of Mines."
        ),
        "conversion_note": (
            "Source GeoTIFF was split into zstarview 90-degree tiles; masked "
            "and negative values were written as zero."
        ),
        "tiles": files,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def write_notice(output_dir: Path) -> None:
    notice = """Night-light data attribution

This product was made utilizing VIIRS nighttime lights data produced by the
Earth Observation Group, Payne Institute for Public Policy, Colorado School
of Mines.

The source VNL data and this derived tile set are distributed under the
Creative Commons Attribution 4.0 International license (CC BY 4.0).

Changes made for zstarview: the source GeoTIFF was split into eight 90-degree
EPSG:4326 tiles; masked and negative values were written as zero.

Source:
https://eogdata.mines.edu/products/vnl/

License:
https://creativecommons.org/licenses/by/4.0/

EOG licensing notice:
https://eogdata.mines.edu/files/EOG_products_CC_License.pdf
"""
    (output_dir / "NOTICE-night-lights.txt").write_text(notice, encoding="ascii")


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    bearer_token = None
    if args.bearer_token_env:
        bearer_token = os.environ.get(args.bearer_token_env)
        if not bearer_token:
            raise SystemExit(f"Environment variable is empty: {args.bearer_token_env}")

    with tempfile.TemporaryDirectory(prefix="zstarview-vnl-") as temporary_name:
        temporary_dir = Path(temporary_name)
        if args.source is not None:
            source_input = args.source.resolve()
            source_path = materialize_source(source_input, temporary_dir)
            source_url = args.source_url
        else:
            source_path = materialize_url(args.url, temporary_dir, bearer_token)
            source_url = args.source_url or args.url
            source_input = temporary_dir / (
                "source.tif.gz" if args.url.lower().endswith(".gz") else "source.tif"
            )

        with rasterio.open(source_path) as source:
            validate_source(source, args.band)
            for tile_name in TILE_NAMES:
                output_path = output_dir / f"{tile_name}.tif"
                print(f"Writing {output_path}")
                write_tile(source, output_path, tile_name, args.band)

        source_filename = source_input.name
        source_sha256 = sha256_file(source_input)

    build_manifest(
        output_dir,
        args.dataset_label,
        source_url,
        args.year,
        args.band,
        source_filename,
        source_sha256,
        args.asset_base_url,
        args.conversion_commit,
    )
    write_notice(output_dir)
    print(f"Wrote manifest: {output_dir / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
