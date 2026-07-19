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
from urllib.error import HTTPError
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
DERIVED_TILE_SCHEMA_VERSION = 3
CACHE_METADATA_SCHEMA_VERSION = 1
MAX_CITY_CODE_RANGE_SIZE = 1000


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


def parse_city_codes(value: str) -> tuple[str, ...]:
    """Expand one or more five-digit municipality codes."""

    codes: list[str] = []
    seen: set[str] = set()
    for item in str(value).split(","):
        token = item.strip()
        if not token:
            raise ValueError("city-code contains an empty item")
        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2:
                raise ValueError(f"invalid city-code range: {token}")
            start_text, end_text = (part.strip() for part in parts)
            _validate_city_code(start_text)
            _validate_city_code(end_text)
            start = int(start_text)
            end = int(end_text)
            if start > end:
                raise ValueError(f"city-code range is reversed: {token}")
            if end - start + 1 > MAX_CITY_CODE_RANGE_SIZE:
                raise ValueError(
                    f"city-code range is too large (maximum: {MAX_CITY_CODE_RANGE_SIZE})"
                )
            expanded = (f"{number:05d}" for number in range(start, end + 1))
        else:
            _validate_city_code(token)
            expanded = (token,)
        for code in expanded:
            if code not in seen:
                seen.add(code)
                codes.append(code)
    if not codes:
        raise ValueError("city-code must not be empty")
    return tuple(codes)


def _validate_city_code(value: str) -> None:
    if not value.isdigit() or len(value) != 5:
        raise ValueError(f"city-code must be a five-digit code: {value}")


def _catalog_size_text(
    catalog: dict[str, object], entries: tuple[dict[str, object], ...]
) -> tuple[int | None, str]:
    archive_url = catalog_archive_url(catalog)
    archive_size = _content_length(archive_url) if archive_url is not None else None
    if archive_size is not None:
        return archive_size, "CityGML ZIP"
    file_size = catalog_file_size_bytes(entries)
    if file_size is not None:
        return file_size, "GML files"
    return None, "unknown"


def _complete_cache_metadata(
    output_root: Path, city_code: str
) -> dict[str, object] | None:
    if not output_root.is_dir():
        return None
    for dataset_dir in output_root.glob(f"{city_code}_*"):
        metadata_path = dataset_dir / "cache_meta.json"
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if (
            isinstance(metadata, dict)
            and metadata.get("city_code") == city_code
            and metadata.get("status") == "complete"
        ):
            return metadata
    return None


def catalog_registration_year(payload: dict[str, object], default: str) -> str:
    cities = payload.get("cities")
    if isinstance(cities, list) and cities and isinstance(cities[0], dict):
        value = cities[0].get("registrationYear")
        if isinstance(value, (int, str)):
            return str(value)
    return default


def catalog_file_size_bytes(
    entries: tuple[dict[str, object], ...],
) -> int | None:
    sizes = [
        int(row["fileSize"])
        for row in entries
        if isinstance(row.get("fileSize"), (int, float)) and int(row["fileSize"]) >= 0
    ]
    return sum(sizes) if len(sizes) == len(entries) else None


def _catalog_source_metadata(
    catalog: dict[str, object],
    entries: tuple[dict[str, object], ...],
    default_year: str,
) -> dict[str, object]:
    cities = catalog.get("cities")
    spec: str | None = None
    if isinstance(cities, list) and cities and isinstance(cities[0], dict):
        value = cities[0].get("spec")
        if value is not None:
            spec = str(value)
    return {
        "preparation_year": catalog_year(catalog, default_year),
        "registration_year": catalog_registration_year(catalog, ""),
        "source_spec": spec,
        "source_file_size_bytes": catalog_file_size_bytes(entries),
        "source_file_count": len(entries),
    }


def _cache_matches_catalog(
    metadata: dict[str, object],
    catalog: dict[str, object],
    entries: tuple[dict[str, object], ...],
    default_year: str,
) -> bool:
    if metadata.get("metadata_schema_version") != CACHE_METADATA_SCHEMA_VERSION:
        return False
    if metadata.get("derived_tile_schema_version") != DERIVED_TILE_SCHEMA_VERSION:
        return False
    expected = _catalog_source_metadata(catalog, entries, default_year)
    return all(metadata.get(key) == value for key, value in expected.items())


def _cache_is_current(
    output_root: Path,
    city_code: str,
    catalog: dict[str, object],
    entries: tuple[dict[str, object], ...],
    default_year: str,
) -> bool:
    metadata = _complete_cache_metadata(output_root, city_code)
    return metadata is not None and _cache_matches_catalog(
        metadata, catalog, entries, default_year
    )


def _preflight_city_codes(
    args: argparse.Namespace, city_codes: tuple[str, ...]
) -> tuple[str, ...]:
    total_size = 0
    known_total = True
    available_codes: list[str] = []
    print("PLATEAU batch download estimate:")
    for city_code in city_codes:
        try:
            catalog = fetch_catalog(city_code, catalog_url=str(args.catalog_url))
        except HTTPError as exc:
            if exc.code != 404:
                raise
            print(
                f"Skipping city code {city_code}: PLATEAU catalog not found (HTTP 404)."
            )
            continue
        entries = catalog_file_entries(catalog)
        if not entries:
            print(f"Skipping city code {city_code}: no building files in catalog.")
            continue
        if not args.overwrite and _cache_is_current(
            Path(args.output_root), city_code, catalog, entries, str(args.year)
        ):
            print(f"Skipping city code {city_code}: cache is up to date.")
            continue
        size_bytes, size_kind = _catalog_size_text(catalog, entries)
        if size_bytes is None:
            known_total = False
            size_text = "unknown"
        else:
            total_size += size_bytes
            size_text = format_binary_size(size_bytes)
        print(f"  {city_code}: {len(entries)} files, {size_text} ({size_kind})")
        available_codes.append(city_code)
    if not available_codes:
        raise ValueError(
            "No PLATEAU building catalogs found for the requested city codes"
        )
    total_text = format_binary_size(total_size) if known_total else "unknown"
    print(f"Total estimated download size: {total_text}")
    if not args.yes:
        answer = input("Continue with PLATEAU batch download? [y/N] ").strip().lower()
        if answer not in {"y", "yes"}:
            print("Download cancelled.")
            return ()
    args.yes = True
    return tuple(available_codes)


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


def _parse_surface_ring(
    text: str | None,
) -> tuple[tuple[float, float, float], ...]:
    if not text:
        return ()
    try:
        values = [float(value) for value in text.split()]
    except ValueError:
        return ()
    if len(values) < 12 or len(values) % 3 != 0:
        return ()
    points: list[tuple[float, float, float]] = []
    for index in range(0, len(values), 3):
        lat, lon, elevation = values[index : index + 3]
        point = (lon, lat, elevation)
        if not points or points[-1] != point:
            points.append(point)
    if len(points) < 4:
        return ()
    if points[0] != points[-1]:
        points.append(points[0])
    return tuple(points)


def _roof_surfaces(
    building: ET.Element,
) -> tuple[tuple[tuple[float, float, float], ...], ...]:
    surfaces: list[tuple[tuple[float, float, float], ...]] = []
    for roof_surface in building.iter():
        if _local_name(roof_surface.tag) != "RoofSurface":
            continue
        for pos_list in roof_surface.iter():
            if _local_name(pos_list.tag) != "posList":
                continue
            ring = _parse_surface_ring(pos_list.text)
            if ring:
                surfaces.append(ring)
    return tuple(surfaces)


def _building_rings(
    building: ET.Element,
) -> tuple[int, tuple[tuple[tuple[float, float], ...], ...]]:
    for lod, container_names in (
        (1, {"lod1Solid"}),
        (0, {"lod0RoofEdge", "lod0FootPrint"}),
    ):
        rings = tuple(
            ring
            for container in building.iter()
            if _local_name(container.tag) in container_names
            for pos_list in container.iter()
            if _local_name(pos_list.tag) == "posList"
            for ring in (_parse_ring(pos_list.text),)
            if ring
        )
        if rings:
            if lod == 1:
                rings = (_largest_ring(rings),)
            return lod, rings
    return 0, ()


def _largest_ring(
    rings: tuple[tuple[tuple[float, float], ...], ...],
) -> tuple[tuple[float, float], ...]:
    def area(ring: tuple[tuple[float, float], ...]) -> float:
        return abs(
            sum(
                lon0 * lat1 - lon1 * lat0
                for (lon0, lat0), (lon1, lat1) in zip(ring, ring[1:])
            )
            * 0.5
        )

    return max(rings, key=area)


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
    roof_surfaces = _roof_surfaces(element)
    geometry_lod, rings = _building_rings(element)
    if roof_surfaces:
        geometry_lod = 2
    if not rings and roof_surfaces:
        rings = (
            tuple((lon, lat) for lon, lat, _elevation in roof_surfaces[0]),
        )
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
        "geometry_lod": geometry_lod,
        "roof_surfaces": [
            [[lon, lat, elevation] for lon, lat, elevation in surface]
            for surface in roof_surfaces
        ],
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


def summarize_geometry_lods(output_dir: Path) -> dict[str, int]:
    counts = {
        "lod0": 0,
        "lod1": 0,
        "lod2": 0,
        "lod3": 0,
        "lod4": 0,
        "detailed_surfaces": 0,
    }
    for path in output_dir.glob("*.json"):
        if path.name == "tile_index.json":
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        buildings = payload.get("buildings") if isinstance(payload, dict) else None
        if not isinstance(buildings, list):
            continue
        for building in buildings:
            if not isinstance(building, dict):
                continue
            try:
                lod = int(building.get("geometry_lod", 0))
            except (TypeError, ValueError):
                lod = 0
            key = f"lod{lod}" if 0 <= lod <= 4 else "lod0"
            counts[key] += 1
            raw_surfaces = building.get("roof_surfaces")
            if isinstance(raw_surfaces, list):
                counts["detailed_surfaces"] += len(raw_surfaces)
    return counts


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download or convert PLATEAU CityGML building data into derived building tiles."
    )
    parser.add_argument(
        "--city-code",
        required=True,
        help=(
            "Five-digit municipality code, range, or comma-separated codes, "
            "e.g. 32201, 13100-13122."
        ),
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


def _prepare_city_code(args: argparse.Namespace, city_code: str) -> int:
    zip_path = args.input_zip
    temporary_dir: Path | None = None
    downloaded_total: int | None = None
    archive_url: str | None = None
    catalog_url: str | None = None
    source_spec: str | None = None
    source_metadata: dict[str, object] = {}
    catalog_checked = False
    if zip_path is None:
        catalog = fetch_catalog(city_code, catalog_url=str(args.catalog_url))
        catalog_checked = True
        entries = catalog_file_entries(catalog)
        if not entries:
            raise ValueError(
                f"No CityGML building files found for city code {city_code}"
            )
        source_metadata = _catalog_source_metadata(catalog, entries, str(args.year))
        if not args.overwrite and _cache_matches_catalog(
            _complete_cache_metadata(Path(args.output_root), city_code) or {},
            catalog,
            entries,
            str(args.year),
        ):
            print(f"Skipping city code {city_code}: cache is up to date.")
            return 0
        dataset_year = catalog_year(catalog, str(args.year))
        archive_url = catalog_archive_url(catalog)
        if archive_url is None:
            raise ValueError(
                f"PLATEAU catalog has no CityGML archive URL for city code {city_code}"
            )
        download_url = archive_url
        catalog_url = str(args.catalog_url).format(city_code=city_code)
        source_spec_value = source_metadata.get("source_spec")
        source_spec = source_spec_value if isinstance(source_spec_value, str) else None
        size_bytes, size_kind = _catalog_size_text(catalog, entries)
        size_text = (
            f"{format_binary_size(size_bytes)} ({size_kind})"
            if size_bytes is not None
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
                if not catalog_checked:
                    raise FileExistsError(
                        f"Prepared PLATEAU cache already exists: {dataset_dir} (use --overwrite to replace)"
                    )
                outdated_dir = dataset_dir.with_name(
                    f"{dataset_dir.name}.outdated-{datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')}"
                )
                dataset_dir.rename(outdated_dir)
                print(f"Moved outdated cache to: {outdated_dir}")
            else:
                incomplete_dir = dataset_dir.with_name(
                    f"{dataset_dir.name}.incomplete-{datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')}"
                )
                dataset_dir.rename(incomplete_dir)
                print(f"Moved incomplete cache to: {incomplete_dir}")
        elif dataset_dir.exists():
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
        geometry_lod_counts = summarize_geometry_lods(derived_dir)
        max_geometry_lod = max(
            (lod for lod in range(4, -1, -1) if geometry_lod_counts[f"lod{lod}"]),
            default=0,
        )
        metadata = {
            "metadata_schema_version": CACHE_METADATA_SCHEMA_VERSION,
            "derived_tile_schema_version": DERIVED_TILE_SCHEMA_VERSION,
            "source": "PLATEAU-CityGML",
            "source_spec": source_spec,
            "city_code": city_code,
            "year": dataset_year,
            "preparation_year": source_metadata.get("preparation_year", dataset_year),
            "registration_year": source_metadata.get("registration_year", ""),
            "source_file_size_bytes": source_metadata.get("source_file_size_bytes"),
            "source_file_count": source_metadata.get("source_file_count"),
            "geometry_mode": (
                "lod2-roof-surfaces"
                if max_geometry_lod >= 2
                else "lod1-footprint"
                if max_geometry_lod >= 1
                else "lod0-footprint"
            ),
            "max_geometry_lod": max_geometry_lod,
            "lod0_building_count": geometry_lod_counts["lod0"],
            "lod1_building_count": geometry_lod_counts["lod1"],
            "lod2_building_count": geometry_lod_counts["lod2"],
            "detailed_surface_count": geometry_lod_counts["detailed_surfaces"],
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


def main(argv: Sequence[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    city_codes = parse_city_codes(str(args.city_code))
    if args.input_zip is not None and len(city_codes) != 1:
        raise ValueError("--input-zip can only be used with one city code")
    if args.input_zip is None and len(city_codes) > 1:
        try:
            city_codes = _preflight_city_codes(args, city_codes)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        if not city_codes:
            return 0
    for index, city_code in enumerate(city_codes, start=1):
        if len(city_codes) > 1:
            print(
                f"Preparing PLATEAU city code {city_code} ({index}/{len(city_codes)})"
            )
        try:
            result = _prepare_city_code(args, city_code)
        except HTTPError as exc:
            if exc.code != 404:
                raise
            if len(city_codes) == 1:
                print(
                    f"PLATEAU catalog request failed for city code {city_code}: "
                    "HTTP 404. PLATEAU building data may not be available "
                    "for this municipality."
                )
                return 1
            print(
                f"Skipping city code {city_code}: PLATEAU catalog not found "
                "(HTTP 404); PLATEAU building data may not be available "
                "for this municipality."
            )
            continue
        if result != 0:
            return result
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
