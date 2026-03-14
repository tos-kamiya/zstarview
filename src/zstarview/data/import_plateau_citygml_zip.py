#!/usr/bin/env python3
"""Import PLATEAU CityGML ZIP files into bundled derived datasets."""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence
from zipfile import ZipFile

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
DATA_ROOT = SRC_ROOT / "zstarview" / "data"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.data import build_plateau_derived_tiles  # noqa: E402
from zstarview.data import build_plateau_tile_index  # noqa: E402
from zstarview.data import urban_debug_layer_from_citygml  # noqa: E402

_CITY_DATASET_RE = re.compile(
    r"^(?P<city_code>\d{5})_(?P<city_token>.+?)(?:_city(?:_\d{4})?_citygml_|_20\d{2}_citygml_).*"
)
_CITY_SUFFIX_RE = re.compile(r"-(?:shi|ku|cho|machi|son|mura)$")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract PLATEAU CityGML ZIP files and write derived skyline datasets into the app data directory."
    )
    parser.add_argument(
        "zip_path",
        nargs="+",
        type=Path,
        help="One or more PLATEAU CityGML ZIP files.",
    )
    parser.add_argument(
        "--derived-root-dir",
        type=Path,
        default=DATA_ROOT / "plateau_derived",
        help="Root directory where derived datasets will be written.",
    )
    parser.add_argument(
        "--urban-outline-root-dir",
        type=Path,
        default=DATA_ROOT / "viewpoints" / "urban_debug_layer",
        help="Directory for optional one-off urban outline JSON exports.",
    )
    parser.add_argument(
        "--write-urban-outline-json",
        action="store_true",
        help="Also export a one-off urban outline JSON for each imported city.",
    )
    parser.add_argument(
        "--min-building-height-m",
        type=float,
        default=build_plateau_derived_tiles.DEFAULT_MIN_BUILDING_HEIGHT_M,
        help="Ignore buildings lower than this height in meters.",
    )
    parser.add_argument(
        "--radius-km",
        type=float,
        default=urban_debug_layer_from_citygml.DEFAULT_RADIUS_KM,
        help="Radius for optional urban outline export.",
    )
    parser.add_argument(
        "--edge-sample-step-m",
        type=float,
        default=urban_debug_layer_from_citygml.DEFAULT_EDGE_SAMPLE_STEP_M,
        help="Edge sampling step for optional urban outline export.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Worker count for tile conversion. Default keeps imports deterministic.",
    )
    parser.add_argument(
        "--keep-extracted-dir",
        type=Path,
        default=None,
        help="Optional directory where extracted ZIP contents will be copied for inspection.",
    )
    return parser


def detect_output_city_name(zip_path: Path) -> str:
    match = _CITY_DATASET_RE.match(zip_path.stem)
    if not match:
        raise ValueError(f"Could not detect city metadata from ZIP filename: {zip_path.name}")
    city_token = _CITY_SUFFIX_RE.sub("", match.group("city_token"))
    city_slug = city_token.replace("-", "")
    if not city_slug:
        raise ValueError(f"Could not derive output city slug from ZIP filename: {zip_path.name}")
    return f"{match.group('city_code')}_{city_slug}"


def find_citygml_building_dir(extracted_root: Path) -> Path:
    matches = sorted(path for path in extracted_root.rglob("bldg") if path.is_dir() and path.parent.name == "udx")
    if not matches:
        raise ValueError(f"No udx/bldg directory found under extracted ZIP: {extracted_root}")
    if len(matches) > 1:
        joined = ", ".join(str(path) for path in matches)
        raise ValueError(f"Multiple udx/bldg directories found under extracted ZIP: {joined}")
    return matches[0]


def build_urban_outline_output_path(urban_outline_root_dir: Path, city_dir_name: str) -> Path:
    city_slug = city_dir_name.split("_", 1)[1] if "_" in city_dir_name else city_dir_name
    return urban_outline_root_dir / f"{city_slug}_urban_outline_layer.json"


def import_zip(
    zip_path: Path,
    *,
    derived_root_dir: Path,
    urban_outline_root_dir: Path,
    write_urban_outline_json: bool,
    min_building_height_m: float,
    radius_km: float,
    edge_sample_step_m: float,
    workers: int,
    keep_extracted_dir: Path | None,
) -> tuple[Path, Path | None]:
    if not zip_path.is_file():
        raise FileNotFoundError(f"ZIP file not found: {zip_path}")

    city_dir_name = detect_output_city_name(zip_path)
    derived_dir = derived_root_dir / city_dir_name / "bldg"

    with TemporaryDirectory(prefix="plateau-import-") as temp_dir_str:
        extract_root = Path(temp_dir_str) / zip_path.stem
        extract_root.mkdir(parents=True, exist_ok=True)
        with ZipFile(zip_path) as zf:
            zf.extractall(extract_root)
        if keep_extracted_dir is not None:
            destination = keep_extracted_dir / zip_path.stem
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(extract_root, destination)
        citygml_dir = find_citygml_building_dir(extract_root)

        build_plateau_derived_tiles.main(
            [
                "--citygml-dir",
                str(citygml_dir),
                "--output-dir",
                str(derived_dir),
                "--min-building-height-m",
                str(float(min_building_height_m)),
                "--workers",
                str(max(1, int(workers))),
            ]
        )

    build_plateau_tile_index.main(
        [
            "--derived-dir",
            str(derived_dir),
        ]
    )

    outline_path: Path | None = None
    if write_urban_outline_json:
        outline_path = build_urban_outline_output_path(urban_outline_root_dir, city_dir_name)
        urban_debug_layer_from_citygml.main(
            [
                "--derived-dir",
                str(derived_dir),
                "--all-covered-towers",
                "--radius-km",
                str(float(radius_km)),
                "--min-building-height-m",
                str(float(min_building_height_m)),
                "--edge-sample-step-m",
                str(float(edge_sample_step_m)),
                "--output-json",
                str(outline_path),
            ]
        )

    return derived_dir, outline_path


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    for zip_path in args.zip_path:
        derived_dir, outline_path = import_zip(
            zip_path,
            derived_root_dir=args.derived_root_dir,
            urban_outline_root_dir=args.urban_outline_root_dir,
            write_urban_outline_json=bool(args.write_urban_outline_json),
            min_building_height_m=float(args.min_building_height_m),
            radius_km=float(args.radius_km),
            edge_sample_step_m=float(args.edge_sample_step_m),
            workers=max(1, int(args.workers)),
            keep_extracted_dir=args.keep_extracted_dir,
        )
        print(f"[ok] imported: {zip_path} -> {derived_dir}")
        if outline_path is not None:
            print(f"[ok] urban-outline-json: {outline_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
