"""Package experimental coastline tiles as GitHub Release assets."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import tempfile
import zipfile
from pathlib import Path


PARENT_RE = re.compile(r"tile_y(?P<y>\d+)_x(?P<x>\d+)(?P<suffix>\.geojson|\.0)$")
CHILD_RE = re.compile(
    r"tile_y(?P<y>\d+)_x(?P<x>\d+)_q(?P<q>[0-3]{2})(?P<suffix>\.geojson|\.0)$"
)
GRID_ROOT = "grid-32x16"


def _gzip_copy(source: Path, destination: Path) -> None:
    with source.open("rb") as source_handle, destination.open("wb") as destination_handle:
        with gzip.GzipFile(
            fileobj=destination_handle,
            mode="wb",
            compresslevel=9,
            mtime=0,
        ) as gzip_handle:
            shutil.copyfileobj(source_handle, gzip_handle, length=1024 * 1024)


def _zip_add(zip_handle: zipfile.ZipFile, source: Path, archive_name: str) -> None:
    # ZipFile.write streams the source instead of materializing a 250MB
    # GeoJSON file in memory. Gzip payloads are already compressed, so store
    # them without a second compression pass.
    compression = zipfile.ZIP_STORED if source.suffix == ".gz" else zipfile.ZIP_DEFLATED
    zip_handle.write(source, arcname=archive_name, compress_type=compression)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_parent(path: Path) -> tuple[int, int] | None:
    match = PARENT_RE.fullmatch(path.name)
    if match is None:
        return None
    return int(match["y"]), int(match["x"])


def _parse_child(path: Path) -> tuple[int, int, str] | None:
    match = CHILD_RE.fullmatch(path.name)
    if match is None:
        return None
    return int(match["y"]), int(match["x"]), match["q"]


def _copy_tile(source: Path, tile_dir: Path, basename: str) -> str:
    tile_dir.mkdir(parents=True, exist_ok=True)
    if source.suffix == ".0":
        destination = tile_dir / f"{basename}.0"
        destination.write_bytes(b"")
        return destination.name
    destination = tile_dir / f"{basename}.geojson.gz"
    _gzip_copy(source, destination)
    return destination.name


def _build_tree(parent_dir: Path, child_dir: Path, tree_dir: Path) -> list[dict[str, object]]:
    child_parents = {
        (parsed[0], parsed[1])
        for path in child_dir.iterdir()
        if (parsed := _parse_child(path)) is not None
    }
    tiles: list[dict[str, object]] = []
    for source in sorted(parent_dir.iterdir()):
        parsed = _parse_parent(source)
        if parsed is None:
            continue
        row, col = parsed
        tile_dir = tree_dir / f"y{row:02d}" / f"x{col:02d}"
        tile_dir.mkdir(parents=True, exist_ok=True)
        if (row, col) in child_parents:
            children: list[dict[str, str]] = []
            for child in sorted(child_dir.iterdir()):
                child_parsed = _parse_child(child)
                if child_parsed is None or child_parsed[:2] != (row, col):
                    continue
                _, _, quadrant = child_parsed
                output_name = _copy_tile(child, tile_dir / "children", f"q{quadrant}")
                children.append({"quadrant": quadrant, "file": f"children/{output_name}"})
            child_manifest = tile_dir / "children.json"
            child_manifest.write_text(
                json.dumps(
                    {
                        "schema": 1,
                        "parent": f"y{row:02d}/x{col:02d}",
                        "subdivision": "4x4",
                        "children": children,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            tiles.append(
                {
                    "y": row,
                    "x": col,
                    "layout": "4x4",
                    "children": len(children),
                }
            )
            print(f"assembled y{row:02d} x{col:02d} as 4x4", flush=True)
            continue
        output_name = _copy_tile(source, tile_dir, "tile")
        tiles.append({"y": row, "x": col, "layout": "single", "file": output_name})
        print(f"assembled y{row:02d} x{col:02d} as single", flush=True)
    return tiles


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent-dir", required=True, type=Path)
    parser.add_argument("--child-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="coastline-release-", dir=args.output_dir) as temporary:
        tree_dir = Path(temporary) / GRID_ROOT
        tree_dir.mkdir()
        tiles = _build_tree(args.parent_dir, args.child_dir, tree_dir)
        archives: list[dict[str, object]] = []
        for col in range(32):
            archive_name = f"coastline-grid-x{col:02d}-{args.version}.zip"
            archive_path = args.output_dir / archive_name
            with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
                for source in sorted(tree_dir.glob(f"y*/x{col:02d}/**/*")):
                    if source.is_file():
                        _zip_add(zip_handle, source, source.relative_to(tree_dir.parent).as_posix())
            archives.append(
                {
                    "x": col,
                    "name": archive_name,
                    "bytes": archive_path.stat().st_size,
                    "sha256": _sha256(archive_path),
                }
            )
            print(f"wrote {archive_name} ({archive_path.stat().st_size} bytes)", flush=True)
        manifest = {
            "schema": 1,
            "dataset": "coastline-vector-tiles",
            "version": args.version,
            "source": "OSM Water Polygons",
            "source_url": "https://osmdata.openstreetmap.de/data/water-polygons.html",
            "coordinate_reference_system": "EPSG:4326",
            "grid": {"columns": 32, "rows": 16, "tile_width_deg": 11.25, "tile_height_deg": 11.25},
            "encoding": "GeoJSON compressed with gzip level 9 inside ZIP assets",
            "assets": archives,
            "tiles": sorted(tiles, key=lambda item: (int(item["y"]), int(item["x"]))),
        }
        manifest_path = args.output_dir / f"coastline-manifest-{args.version}.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        checksum_path = args.output_dir / f"coastline-sha256sums-{args.version}.txt"
        assets = [*sorted(args.output_dir.glob(f"coastline-grid-x*-{args.version}.zip")), manifest_path]
        checksum_path.write_text(
            "".join(f"{_sha256(path)}  {path.name}\n" for path in assets), encoding="ascii"
        )
    print(f"Packaged {len(archives)} release assets under {args.output_dir}")


if __name__ == "__main__":
    main()
