#!/usr/bin/env python3
"""Import a PLATEAU CityGML ZIP file into bundled derived datasets."""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence
from zipfile import ZipFile

DATA_ROOT = Path(__file__).resolve().parent

from zstarview.data import build_plateau_derived_tiles  # noqa: E402
from zstarview.data import build_plateau_tile_index  # noqa: E402

_CITY_DATASET_RE = re.compile(
    r"^(?P<city_code>\d{5})_(?P<city_token>.+?)(?:_city(?:_\d{4})?_citygml_|_20\d{2}_citygml_).*"
)
_CITY_SUFFIX_RE = re.compile(r"-(?:shi|ku|cho|machi|son|mura)$")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract a PLATEAU CityGML ZIP file and write derived skyline datasets into the app data directory."
    )
    parser.add_argument(
        "zip_path",
        type=Path,
        help="PLATEAU CityGML ZIP file.",
    )
    parser.add_argument(
        "--derived-root-dir",
        type=Path,
        default=DATA_ROOT / "plateau_derived",
        help="Root directory where derived datasets will be written.",
    )
    parser.add_argument(
        "--min-building-height-m",
        type=float,
        default=build_plateau_derived_tiles.DEFAULT_MIN_BUILDING_HEIGHT_M,
        help="Ignore buildings lower than this height in meters.",
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


def import_zip(
    zip_path: Path,
    *,
    derived_root_dir: Path,
    min_building_height_m: float,
    workers: int,
    keep_extracted_dir: Path | None,
) -> Path:
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
    return derived_dir


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    derived_dir = import_zip(
        args.zip_path,
        derived_root_dir=args.derived_root_dir,
        min_building_height_m=float(args.min_building_height_m),
        workers=max(1, int(args.workers)),
        keep_extracted_dir=args.keep_extracted_dir,
    )
    print(f"[ok] imported: {args.zip_path} -> {derived_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
