"""Download and convert PLATEAU CityGML building data.

This module is intentionally independent from runtime source selection.  It is
the first-stage preparation pipeline used to verify that PLATEAU data can be
converted into the same derived tile shape as Overture building data.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence
from urllib.request import Request, urlopen
from zipfile import ZipFile

from tqdm import tqdm

from zstarview.__about__ import __version__
from zstarview.data import build_derived_tile_index
from zstarview.paths import CACHE_PATH
from zstarview.user_agent import build_user_agent

DEFAULT_API_URL = "https://api.plateauview.mlit.go.jp/datacatalog/citygml/{city_code}-{year}/citygml.zip"
DEFAULT_CATALOG_URL = (
    "https://api.plateauview.mlit.go.jp/datacatalog/citygml/{city_code}?types=bldg"
)
DEFAULT_STOREY_HEIGHT_M = 3.5
DEFAULT_MIN_BUILDING_HEIGHT_M = 0.0
DEFAULT_MAX_DOWNLOAD_BYTES = 10 * 1024 * 1024 * 1024
CHUNK_SIZE = 1024 * 1024
DERIVED_TILE_SCHEMA_VERSION = 1
CACHE_METADATA_SCHEMA_VERSION = 1


def build_download_url(
    city_code: str, year: str, *, api_url: str = DEFAULT_API_URL
) -> str:
    if not city_code.isdigit() or len(city_code) != 5:
        raise ValueError("city_code must be a five-digit municipality code")
    return api_url.format(city_code=city_code, year=year)


def format_binary_size(size_bytes: int) -> str:
    value = float(max(0, int(size_bytes)))
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return "0.00 B"


def _content_length(url: str) -> int | None:
    request = Request(
        url, method="HEAD", headers={"User-Agent": build_user_agent("plateau")}
    )
    try:
        with urlopen(request, timeout=30.0) as response:  # nosec: B310
            value = response.headers.get("Content-Length")
    except Exception:
        return None
    try:
        return int(value) if value is not None else None
    except ValueError:
        return None


def download_file(
    url: str,
    destination: Path,
    *,
    max_bytes: int = DEFAULT_MAX_DOWNLOAD_BYTES,
) -> int:
    request = Request(url, headers={"User-Agent": build_user_agent("plateau")})
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    total = 0
    try:
        with urlopen(request, timeout=60.0) as response:  # nosec: B310
            raw_length = response.headers.get("Content-Length")
            try:
                total_size = int(raw_length) if raw_length is not None else None
            except ValueError:
                total_size = None
            with (
                temporary.open("wb") as handle,
                tqdm(
                    total=total_size,
                    desc="Downloading CityGML",
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                ) as progress,
            ):
                while True:
                    chunk = response.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    total += len(chunk)
                    if total > max_bytes:
                        raise ValueError(f"download exceeds limit of {max_bytes} bytes")
                    handle.write(chunk)
                    progress.update(len(chunk))
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return total


def fetch_catalog(
    city_code: str, *, catalog_url: str = DEFAULT_CATALOG_URL
) -> dict[str, object]:
    url = catalog_url.format(city_code=city_code)
    request = Request(url, headers={"User-Agent": build_user_agent("plateau")})
    with urlopen(request, timeout=30.0) as response:  # nosec: B310
        payload = json.loads(response.read().decode("utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("PLATEAU catalog response is not an object")
    return payload


def catalog_file_entries(payload: dict[str, object]) -> tuple[dict[str, object], ...]:
    entries: list[dict[str, object]] = []
    cities = payload.get("cities")
    city_rows = cities if isinstance(cities, list) else [payload]
    for city in city_rows:
        if not isinstance(city, dict):
            continue
        files = city.get("files")
        if not isinstance(files, dict):
            continue
        rows = files.get("bldg")
        if not isinstance(rows, list):
            continue
        for row in rows:
            if not isinstance(row, dict):
                continue
            url = row.get("url")
            if isinstance(url, str) and url.startswith(("http://", "https://")):
                entries.append(dict(row))
    return tuple(entries)


def catalog_year(payload: dict[str, object], default: str) -> str:
    cities = payload.get("cities")
    if isinstance(cities, list) and cities and isinstance(cities[0], dict):
        value = cities[0].get("year")
        if isinstance(value, (int, str)):
            return str(value)
    return default


def catalog_archive_url(payload: dict[str, object]) -> str | None:
    cities = payload.get("cities")
    if isinstance(cities, list) and cities and isinstance(cities[0], dict):
        url = cities[0].get("url")
        if isinstance(url, str) and url.startswith(("http://", "https://")):
            return url
    return None


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _text_value(element: ET.Element, name: str) -> str | None:
    for child in element.iter():
        if _local_name(child.tag) == name and child.text:
            return child.text.strip()
    return None


def _numeric(value: str | None) -> float | None:
    try:
        parsed = float(value or "")
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) and parsed != -9999.0 else None


def _parse_ring(text: str | None) -> tuple[tuple[float, float], ...]:
    if not text:
        return ()
    try:
        values = [float(value) for value in text.split()]
    except ValueError:
        return ()
    if len(values) < 12:
        return ()
    stride = 3 if len(values) % 3 == 0 else 2
    points: list[tuple[float, float]] = []
    for index in range(0, len(values) - stride + 1, stride):
        lat, lon = values[index], values[index + 1]
        point = (lon, lat)
        if not points or points[-1] != point:
            points.append(point)
    if len(points) < 4:
        return ()
    if points[0] != points[-1]:
        points.append(points[0])
    return tuple(points)


def _building_rings(
    building: ET.Element,
) -> tuple[tuple[tuple[float, float], ...], ...]:
    for container_name in ("lod0RoofEdge", "lod0FootPrint", "lod1Solid", "lod2Solid"):
        rings = tuple(
            ring
            for container in building.iter()
            if _local_name(container.tag) == container_name
            for pos_list in container.iter()
            if _local_name(pos_list.tag) == "posList"
            for ring in (_parse_ring(pos_list.text),)
            if ring
        )
        if rings:
            return rings
    return ()


def _building_payload(
    element: ET.Element,
    index: int,
    min_building_height_m: float,
) -> dict[str, object] | None:
    measured = _numeric(_text_value(element, "measuredHeight"))
    if measured is not None:
        height_m = measured
        height_source = "measuredHeight"
    else:
        storeys = _numeric(_text_value(element, "storeysAboveGround"))
        if storeys is None or storeys <= 0 or storeys >= 9999:
            return None
        height_m = storeys * DEFAULT_STOREY_HEIGHT_M
        height_source = "storeysAboveGround*3.5"
    if height_m < min_building_height_m:
        return None
    rings = _building_rings(element)
    if not rings:
        return None
    building_id = (
        element.attrib.get("{http://www.opengis.net/gml}id") or f"bldg-{index}"
    )
    min_lat = min(point[1] for ring in rings for point in ring)
    min_lon = min(point[0] for ring in rings for point in ring)
    max_lat = max(point[1] for ring in rings for point in ring)
    max_lon = max(point[0] for ring in rings for point in ring)
    return {
        "id": building_id,
        "height_m": height_m,
        "height_source": height_source,
        "bbox": {
            "min_lat": min_lat,
            "min_lon": min_lon,
            "max_lat": max_lat,
            "max_lon": max_lon,
        },
        "rings": [[[lon, lat] for lon, lat in ring] for ring in rings],
    }


def _parse_citygml_file(
    path: Path,
    *,
    min_building_height_m: float,
) -> tuple[dict[str, float] | None, tuple[dict[str, object], ...]]:
    lower: list[float] | None = None
    upper: list[float] | None = None
    buildings: list[dict[str, object]] = []
    for index, (_, element) in enumerate(ET.iterparse(path, events=("end",))):
        name = _local_name(element.tag)
        if name in {"lowerCorner", "upperCorner"} and element.text:
            try:
                values = [float(value) for value in element.text.split()]
            except ValueError:
                values = []
            if len(values) >= 2:
                if name == "lowerCorner":
                    lower = values
                else:
                    upper = values
        elif name == "Building":
            payload = _building_payload(element, index, min_building_height_m)
            if payload is not None:
                buildings.append(payload)
            element.clear()
    tile_bbox = None
    if lower is not None and upper is not None:
        tile_bbox = {
            "min_lat": min(lower[0], upper[0]),
            "min_lon": min(lower[1], upper[1]),
            "max_lat": max(lower[0], upper[0]),
            "max_lon": max(lower[1], upper[1]),
        }
    return tile_bbox, tuple(buildings)


def parse_citygml_buildings(
    path: Path,
    *,
    min_building_height_m: float = DEFAULT_MIN_BUILDING_HEIGHT_M,
) -> tuple[dict[str, object], ...]:
    _tile_bbox_value, buildings = _parse_citygml_file(
        path, min_building_height_m=min_building_height_m
    )
    return buildings


def convert_citygml_tile(
    path: Path,
    *,
    output_path: Path,
    city_code: str,
    min_building_height_m: float,
) -> int:
    tile_bbox, buildings = _parse_citygml_file(
        path, min_building_height_m=min_building_height_m
    )
    if tile_bbox is None or not buildings:
        return 0
    payload = {
        "schema_version": DERIVED_TILE_SCHEMA_VERSION,
        "source": {
            "format": "PLATEAU-CityGML",
            "city_code": city_code,
            "path": path.name,
        },
        "tile": {"id": path.stem, "bbox": tile_bbox},
        "filters": {
            "min_height_m": float(min_building_height_m),
            "storey_height_m": DEFAULT_STOREY_HEIGHT_M,
        },
        "buildings": list(buildings),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return len(buildings)


def find_building_files(root: Path) -> tuple[Path, ...]:
    return tuple(
        sorted(
            path
            for path in root.rglob("*.gml")
            if path.is_file()
            and path.parent.name == "bldg"
            and path.parent.parent.name == "udx"
        )
    )


def extract_zip_safely(archive_path: Path, destination: Path) -> None:
    destination = destination.resolve()
    with ZipFile(archive_path) as archive:
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            try:
                target.relative_to(destination)
            except ValueError as exc:
                raise ValueError(
                    f"ZIP member escapes extraction directory: {member.filename}"
                ) from exc
        archive.extractall(destination)


def convert_extracted_citygml(
    extracted_root: Path,
    output_dir: Path,
    *,
    city_code: str,
    min_building_height_m: float,
) -> tuple[int, int]:
    files = find_building_files(extracted_root)
    written_tiles = 0
    building_count = 0
    for path in tqdm(files, desc="Converting CityGML", unit="file"):
        count = convert_citygml_tile(
            path,
            output_path=output_dir / f"{path.stem}.json",
            city_code=city_code,
            min_building_height_m=min_building_height_m,
        )
        if count:
            written_tiles += 1
            building_count += count
    if written_tiles:
        build_derived_tile_index.main(["--derived-dir", str(output_dir)])
    return written_tiles, building_count


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download or convert PLATEAU CityGML building data into derived building tiles."
    )
    parser.add_argument(
        "--city-code", required=True, help="Five-digit municipality code, e.g. 32201."
    )
    parser.add_argument(
        "--year", default="latest", help="PLATEAU preparation year or latest."
    )
    parser.add_argument(
        "--input-zip",
        type=Path,
        help="Use a local CityGML ZIP instead of downloading it.",
    )
    parser.add_argument(
        "--api-url", default=DEFAULT_API_URL, help="CityGML ZIP API URL template."
    )
    parser.add_argument(
        "--catalog-url",
        default=DEFAULT_CATALOG_URL,
        help="CityGML file catalog API URL template.",
    )
    parser.add_argument(
        "--output-root", type=Path, default=Path(CACHE_PATH) / "plateau_buildings"
    )
    parser.add_argument(
        "--min-building-height-m", type=float, default=DEFAULT_MIN_BUILDING_HEIGHT_M
    )
    parser.add_argument(
        "--max-download-bytes", type=int, default=DEFAULT_MAX_DOWNLOAD_BYTES
    )
    parser.add_argument(
        "--keep-zip",
        action="store_true",
        help="Keep the downloaded ZIP for inspection.",
    )
    parser.add_argument(
        "--yes", action="store_true", help="Confirm the download without prompting."
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace an existing prepared cache."
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    city_code = str(args.city_code)
    zip_path = args.input_zip
    temporary_dir: Path | None = None
    downloaded_total: int | None = None
    archive_url: str | None = None
    catalog_url: str | None = None
    source_spec: str | None = None
    if zip_path is None:
        catalog = fetch_catalog(city_code, catalog_url=str(args.catalog_url))
        entries = catalog_file_entries(catalog)
        if not entries:
            raise ValueError(
                f"No CityGML building files found for city code {city_code}"
            )
        dataset_year = catalog_year(catalog, str(args.year))
        archive_url = catalog_archive_url(catalog)
        if archive_url is None:
            raise ValueError(
                f"PLATEAU catalog has no CityGML archive URL for city code {city_code}"
            )
        download_url = archive_url
        catalog_url = str(args.catalog_url).format(city_code=city_code)
        cities = catalog.get("cities")
        if isinstance(cities, list) and cities and isinstance(cities[0], dict):
            spec = cities[0].get("spec")
            source_spec = str(spec) if spec is not None else None
        known_sizes = [
            int(row["fileSize"])
            for row in entries
            if isinstance(row.get("fileSize"), (int, float))
            and int(row["fileSize"]) >= 0
        ]
        if len(known_sizes) == len(entries):
            size_text = f"{format_binary_size(sum(known_sizes))} (GML files)"
        else:
            archive_size = (
                _content_length(archive_url) if archive_url is not None else None
            )
            size_text = (
                f"{format_binary_size(archive_size)} (CityGML ZIP)"
                if archive_size is not None
                else "unknown"
            )
        print(f"PLATEAU catalog: {args.catalog_url.format(city_code=city_code)}")
        print(f"CityGML files: {len(entries)}")
        print(f"Estimated download size: {size_text}")
        if not args.yes:
            answer = input("Continue with PLATEAU download? [y/N] ").strip().lower()
            if answer not in {"y", "yes"}:
                print("Download cancelled.")
                return 0
        temporary_dir = Path(tempfile.mkdtemp(prefix="plateau-download-"))
        zip_path = temporary_dir / "citygml.zip"
        downloaded_total = download_file(
            download_url,
            zip_path,
            max_bytes=max(1, int(args.max_download_bytes)),
        )
        print(f"Downloaded: {downloaded_total} bytes")
    elif not zip_path.is_file():
        raise FileNotFoundError(f"Input ZIP not found: {zip_path}")

    dataset_year = str(args.year) if args.input_zip is not None else dataset_year
    output_root = Path(args.output_root)
    dataset_dir = output_root / f"{city_code}_{dataset_year}"
    existing_metadata_path = dataset_dir / "cache_meta.json"
    if dataset_dir.exists() and not args.overwrite:
        if existing_metadata_path.is_file():
            try:
                existing_metadata = json.loads(
                    existing_metadata_path.read_text(encoding="utf-8")
                )
            except (OSError, ValueError):
                existing_metadata = {}
            if existing_metadata.get("status") == "complete":
                raise FileExistsError(
                    f"Prepared PLATEAU cache already exists: {dataset_dir} (use --overwrite to replace)"
                )
        incomplete_dir = dataset_dir.with_name(
            f"{dataset_dir.name}.incomplete-{datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')}"
        )
        dataset_dir.rename(incomplete_dir)
        print(f"Moved incomplete cache to: {incomplete_dir}")
    output_root.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{city_code}-{dataset_year}-", dir=output_root)
    )
    derived_dir = staging_dir / "bldg"
    try:
        with tempfile.TemporaryDirectory(prefix="plateau-extract-") as extracted_text:
            extracted_root = Path(extracted_text)
            extract_zip_safely(zip_path, extracted_root)
            written_tiles, building_count = convert_extracted_citygml(
                extracted_root,
                derived_dir,
                city_code=city_code,
                min_building_height_m=float(args.min_building_height_m),
            )
        if written_tiles == 0:
            raise ValueError("No PLATEAU building tiles were generated")
        metadata = {
            "metadata_schema_version": CACHE_METADATA_SCHEMA_VERSION,
            "derived_tile_schema_version": DERIVED_TILE_SCHEMA_VERSION,
            "source": "PLATEAU-CityGML",
            "source_spec": source_spec,
            "city_code": city_code,
            "year": dataset_year,
            "fetched_at_utc": datetime.now(timezone.utc).isoformat(),
            "converter": "zstarview-plateau-buildings",
            "converter_version": __version__,
            "tile_count": written_tiles,
            "building_count": building_count,
            "min_building_height_m": float(args.min_building_height_m),
            "catalog_url": catalog_url,
            "archive_url": archive_url,
            "archive_size_bytes": downloaded_total,
            "status": "complete",
        }
        (staging_dir / "cache_meta.json").write_text(
            json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        if temporary_dir is not None and args.keep_zip and zip_path.exists():
            shutil.copyfile(zip_path, staging_dir / "source-citygml.zip")
        if dataset_dir.exists():
            if args.overwrite:
                backup_dir = dataset_dir.with_name(
                    f"{dataset_dir.name}.replaced-{datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')}"
                )
                dataset_dir.rename(backup_dir)
                print(f"Moved replaced cache to: {backup_dir}")
        staging_dir.rename(dataset_dir)
        print(f"Prepared: {dataset_dir / 'bldg'}")
        print(f"Tiles: {written_tiles}; buildings: {building_count}")
        return 0
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir, ignore_errors=True)
        if temporary_dir is not None:
            shutil.rmtree(temporary_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
