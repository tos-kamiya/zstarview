"""Download and convert PLATEAU CityGML building data.

This module is intentionally independent from runtime source selection.  It is
the first-stage preparation pipeline used to verify that PLATEAU data can be
converted into the same derived tile shape as Overture building data.
"""

from __future__ import annotations

import argparse
import errno
import hashlib
import json
import math
import shutil
import sys
import tempfile
import uuid
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from zipfile import BadZipFile, ZipFile

from tqdm import tqdm

from zstarview.__about__ import __version__
from zstarview.data import build_derived_tile_index
from zstarview.paths import CACHE_PATH
from zstarview.user_agent import build_user_agent

if TYPE_CHECKING:
    from zstarview.data.plateau_building_cache import PlateauCacheEntry

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
MIN_DISK_SPACE_BUFFER_BYTES = 256 * 1024 * 1024


class InsufficientPlateauDiskSpaceError(RuntimeError):
    """Raised before a PLATEAU download when its working space is insufficient."""


@dataclass(frozen=True)
class SourceZipInfo:
    """Identity of the CityGML ZIP a cache was converted from."""

    size_bytes: int
    sha256: str


def hash_source_zip(path: Path) -> SourceZipInfo:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(CHUNK_SIZE)
            if not chunk:
                break
            size += len(chunk)
            digest.update(chunk)
    return SourceZipInfo(size_bytes=size, sha256=digest.hexdigest())


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
) -> dict[str, object]:
    cities = catalog.get("cities")
    spec: str | None = None
    if isinstance(cities, list) and cities and isinstance(cities[0], dict):
        value = cities[0].get("spec")
        if value is not None:
            spec = str(value)
    return {
        "preparation_year": catalog_year(catalog, ""),
        "registration_year": catalog_registration_year(catalog, ""),
        "source_spec": spec,
        "source_file_size_bytes": catalog_file_size_bytes(entries),
        "source_file_count": len(entries),
    }


def _cache_matches_catalog(
    metadata: dict[str, object],
    catalog: dict[str, object],
    entries: tuple[dict[str, object], ...],
) -> bool:
    if metadata.get("metadata_schema_version") != CACHE_METADATA_SCHEMA_VERSION:
        return False
    if metadata.get("derived_tile_schema_version") != DERIVED_TILE_SCHEMA_VERSION:
        return False
    expected = _catalog_source_metadata(catalog, entries)
    return all(metadata.get(key) == value for key, value in expected.items())


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


def _required_working_space_bytes(
    archive_size_bytes: int | None,
    entries: tuple[dict[str, object], ...],
) -> int | None:
    if archive_size_bytes is None:
        return None
    source_size_bytes = catalog_file_size_bytes(entries)
    if source_size_bytes is None:
        return None
    working_bytes = int(archive_size_bytes) + int(source_size_bytes)
    return working_bytes + max(MIN_DISK_SPACE_BUFFER_BYTES, working_bytes // 10)


def _existing_space_path(path: Path) -> Path:
    candidate = path.expanduser()
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    return candidate


def _check_working_space_before_download(
    *,
    temp_base_dir: Path,
    archive_size_bytes: int | None,
    entries: tuple[dict[str, object], ...],
) -> None:
    required_bytes = _required_working_space_bytes(archive_size_bytes, entries)
    if required_bytes is None:
        return
    usage = shutil.disk_usage(_existing_space_path(temp_base_dir))
    if usage.free >= required_bytes:
        return
    raise InsufficientPlateauDiskSpaceError(
        "PLATEAU preparation needs about "
        f"{format_binary_size(required_bytes)} in the temporary filesystem, "
        f"but only {format_binary_size(usage.free)} is available."
    )


def download_file(
    url: str,
    destination: Path,
    *,
    max_bytes: int = DEFAULT_MAX_DOWNLOAD_BYTES,
) -> SourceZipInfo:
    """Download ``url`` to ``destination`` and return its size and SHA-256."""

    request = Request(url, headers={"User-Agent": build_user_agent("plateau")})
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    total = 0
    digest = hashlib.sha256()
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
                    digest.update(chunk)
                    progress.update(len(chunk))
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return SourceZipInfo(size_bytes=total, sha256=digest.hexdigest())


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


_SUBCOMMANDS = ("prepare", "list", "remove")
_CATALOG_FIELDS = (
    "year",
    "preparation_year",
    "registration_year",
    "source_spec",
    "source_file_size_bytes",
    "source_file_count",
    "catalog_url",
    "archive_url",
    "fetched_at_utc",
)


def default_output_root() -> Path:
    return Path(CACHE_PATH) / "plateau_buildings"


def build_arg_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--output-root",
        type=Path,
        default=default_output_root(),
        help="PLATEAU cache root directory (default: %(default)s).",
    )
    parser = argparse.ArgumentParser(
        prog="zstarview-download-plateau-buildings",
        description=(
            "Prepare, list, or remove PLATEAU building caches. Each municipality "
            "keeps one cache generation in <output-root>/<city-code>/."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    prepare = subparsers.add_parser(
        "prepare",
        parents=[common],
        help="Download (or convert a local ZIP) and replace municipality caches.",
        description=(
            "Prepare PLATEAU building caches. The PLATEAU catalog is checked first; "
            "a cache that matches the catalog is not downloaded again. A new cache "
            "is built in a staging directory and replaces the existing one only "
            "after it has been verified; older layouts and leftovers of the same "
            "municipality are then removed."
        ),
    )
    prepare.add_argument(
        "city_codes",
        nargs="+",
        metavar="CITY_CODE",
        help=(
            "Five-digit municipality code, range, or comma-separated codes, "
            "e.g. 32201, 13101-13123."
        ),
    )
    prepare.add_argument(
        "--input-zip",
        type=Path,
        help=(
            "Convert a local CityGML ZIP instead of downloading it "
            "(one city code only)."
        ),
    )
    prepare.add_argument(
        "--keep-zip",
        action="store_true",
        help="Keep the source CityGML ZIP as <cache>/source-citygml.zip.",
    )
    prepare.add_argument(
        "--temp-dir",
        type=Path,
        help=(
            "Base directory for the downloaded ZIP and extracted CityGML "
            "temporary files (default: system temporary directory)."
        ),
    )
    prepare.add_argument(
        "--max-download-bytes",
        type=int,
        default=DEFAULT_MAX_DOWNLOAD_BYTES,
        help="Abort a download larger than this many bytes.",
    )
    prepare.add_argument(
        "--yes", action="store_true", help="Do not ask for confirmation."
    )

    list_parser = subparsers.add_parser(
        "list",
        parents=[common],
        help="List caches, including unused directories that can be removed.",
        description=(
            "List PLATEAU caches without network access. Each line shows the "
            "city code, data year, status (current, legacy, removable), and path. "
            "Removable directories are never read by zstarview."
        ),
    )
    list_parser.add_argument(
        "--jsonl",
        action="store_true",
        help="Output cache metadata as JSON Lines.",
    )

    remove = subparsers.add_parser(
        "remove",
        parents=[common],
        help="Remove municipality caches or unused cache directories.",
    )
    remove.add_argument(
        "city_codes",
        nargs="*",
        metavar="CITY_CODE",
        help="Municipality codes to remove (all codes with --unused if omitted).",
    )
    remove.add_argument(
        "--unused",
        action="store_true",
        help="Remove only directories that zstarview does not read.",
    )
    remove.add_argument(
        "--yes", action="store_true", help="Do not ask for confirmation."
    )
    return parser


def _confirm(prompt: str, *, assume_yes: bool) -> bool:
    if assume_yes:
        return True
    answer = input(f"{prompt} [y/N] ").strip().lower()
    return answer in {"y", "yes"}


def _cache_entry(output_root: Path, city_code: str) -> PlateauCacheEntry | None:
    from zstarview.data.plateau_building_cache import classify_plateau_caches

    for entry in classify_plateau_caches(output_root):
        if entry.city_code == city_code:
            return entry
    return None


def _selected_metadata(
    output_root: Path, city_code: str
) -> tuple[Path | None, bool, dict[str, object] | None]:
    """Return the readable cache for ``city_code`` with its metadata."""

    from zstarview.data.plateau_building_cache import read_plateau_cache_metadata

    entry = _cache_entry(output_root, city_code)
    if entry is None or entry.selected is None:
        return None, False, None
    return entry.selected, entry.legacy, read_plateau_cache_metadata(entry.selected)


def _remove_tree(path: Path) -> bool:
    try:
        shutil.rmtree(path)
    except OSError as exc:
        print(f"Warning: could not remove {path}: {exc}", file=sys.stderr)
        return False
    print(f"Removed: {path}")
    return True


def _leftover_dirs(
    output_root: Path, city_code: str, keep: Path | None = None
) -> tuple[Path, ...]:
    from zstarview.data.plateau_building_cache import plateau_city_dirs

    current = output_root / city_code
    return tuple(
        path
        for path in plateau_city_dirs(output_root, city_code)
        if path != current and path != keep
    )


def _install_cache(output_root: Path, city_code: str, source_dir: Path) -> Path:
    """Move ``source_dir`` to ``<root>/<code>/`` and remove other directories.

    The existing ``<code>/`` is first renamed to a dot-prefixed name, which is
    never read, so an interruption between the two renames leaves either the
    old or the new cache readable (or only the legacy layout, which is read as
    a fallback).  Leftovers are removed on the next preparation.
    """

    target = output_root / city_code
    retired: Path | None = None
    if target.exists():
        retired = output_root / f".{city_code}-retired-{uuid.uuid4().hex}"
        target.rename(retired)
    try:
        source_dir.rename(target)
    except OSError:
        if retired is not None and not target.exists():
            retired.rename(target)
        raise
    print(f"Prepared: {target / 'bldg'}")
    for path in _leftover_dirs(output_root, city_code):
        _remove_tree(path)
    return target


@dataclass
class _PreparePlan:
    city_code: str
    action: str  # "skip", "migrate", "download", or "convert"
    current: Path | None
    current_legacy: bool
    removals: tuple[Path, ...]
    needs_confirmation: bool
    description: str
    catalog: dict[str, object] | None = None
    entries: tuple[dict[str, object], ...] = ()
    size_bytes: int | None = None
    size_kind: str = "unknown"
    source_metadata: dict[str, object] = field(default_factory=dict)
    zip_path: Path | None = None
    zip_info: SourceZipInfo | None = None


def _plan_catalog(args: argparse.Namespace, city_code: str) -> _PreparePlan:
    output_root = Path(args.output_root)
    catalog = fetch_catalog(city_code)
    entries = catalog_file_entries(catalog)
    if not entries:
        raise ValueError(f"No CityGML building files found for city code {city_code}")
    current, legacy, metadata = _selected_metadata(output_root, city_code)
    source_metadata = _catalog_source_metadata(catalog, entries)
    if metadata is not None and _cache_matches_catalog(metadata, catalog, entries):
        removals = _leftover_dirs(output_root, city_code, keep=current)
        if legacy:
            return _PreparePlan(
                city_code=city_code,
                action="migrate",
                current=current,
                current_legacy=True,
                removals=removals,
                needs_confirmation=True,
                description=(
                    "cache is up to date; move legacy cache to "
                    f"{output_root / city_code} without downloading"
                ),
            )
        return _PreparePlan(
            city_code=city_code,
            action="skip",
            current=current,
            current_legacy=False,
            removals=removals,
            needs_confirmation=bool(removals),
            description="cache is up to date; nothing to download",
        )
    size_bytes, size_kind = _catalog_size_text(catalog, entries)
    size_text = (
        f"{format_binary_size(size_bytes)} ({size_kind})"
        if size_bytes is not None
        else "unknown size"
    )
    year = source_metadata.get("preparation_year") or "unknown year"
    return _PreparePlan(
        city_code=city_code,
        action="download",
        current=current,
        current_legacy=legacy,
        removals=_leftover_dirs(output_root, city_code),
        needs_confirmation=True,
        description=(
            f"download {len(entries)} building files of {year}, {size_text}"
            + ("; replace current cache" if current is not None else "")
        ),
        catalog=catalog,
        entries=entries,
        size_bytes=size_bytes,
        size_kind=size_kind,
        source_metadata=source_metadata,
    )


def _plan_zip(args: argparse.Namespace, city_code: str) -> _PreparePlan:
    output_root = Path(args.output_root)
    zip_path = Path(args.input_zip)
    if not zip_path.is_file():
        raise FileNotFoundError(f"Input ZIP not found: {zip_path}")
    print(f"Computing SHA-256 of {zip_path} ...")
    zip_info = hash_source_zip(zip_path)
    current, legacy, metadata = _selected_metadata(output_root, city_code)
    removals = _leftover_dirs(output_root, city_code)
    same_source = (
        metadata is not None
        and not legacy
        and metadata.get("source_zip_sha256") == zip_info.sha256
    )
    if same_source and metadata is not None:
        source_metadata = {
            key: metadata.get(key) for key in _CATALOG_FIELDS if key in metadata
        }
        source_metadata["archive_size_bytes"] = metadata.get("archive_size_bytes")
        description = f"reconvert {zip_path} (same source ZIP as current cache)"
    else:
        source_metadata = {}
        description = f"convert {zip_path}" + (
            "; replace current cache built from a different or unknown source"
            if current is not None
            else ""
        )
    return _PreparePlan(
        city_code=city_code,
        action="convert",
        current=current,
        current_legacy=legacy,
        removals=removals,
        needs_confirmation=not same_source or bool(removals),
        description=description,
        source_metadata=source_metadata,
        zip_path=zip_path,
        zip_info=zip_info,
    )


def _print_plans(plans: Sequence[_PreparePlan]) -> None:
    print("PLATEAU preparation plan:")
    total_size = 0
    known_total = True
    downloads = 0
    for plan in plans:
        print(f"  {plan.city_code}: {plan.description}")
        if plan.current is not None:
            layout = " (legacy layout)" if plan.current_legacy else ""
            print(f"    current cache: {plan.current}{layout}")
        if plan.removals:
            print("    remove after success:")
            for path in plan.removals:
                print(f"      {path}")
        if plan.action == "download":
            downloads += 1
            if plan.size_bytes is None:
                known_total = False
            else:
                total_size += plan.size_bytes
    if downloads:
        total_text = format_binary_size(total_size) if known_total else "unknown"
        print(f"Total estimated download size: {total_text}")


def _build_cache(
    args: argparse.Namespace,
    plan: _PreparePlan,
    *,
    zip_path: Path,
    zip_info: SourceZipInfo,
    downloaded: bool,
) -> None:
    from zstarview.data.plateau_building_cache import is_valid_plateau_cache

    city_code = plan.city_code
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir.expanduser() if args.temp_dir is not None else None
    staging_dir = Path(tempfile.mkdtemp(prefix=f".{city_code}-", dir=output_root))
    derived_dir = staging_dir / "bldg"
    try:
        with tempfile.TemporaryDirectory(
            prefix="plateau-extract-",
            dir=str(temp_dir) if temp_dir is not None else None,
        ) as extracted_text:
            extracted_root = Path(extracted_text)
            extract_zip_safely(zip_path, extracted_root)
            written_tiles, building_count = convert_extracted_citygml(
                extracted_root,
                derived_dir,
                city_code=city_code,
                min_building_height_m=DEFAULT_MIN_BUILDING_HEIGHT_M,
            )
        if written_tiles == 0:
            raise ValueError("No PLATEAU building tiles were generated")
        geometry_lod_counts = summarize_geometry_lods(derived_dir)
        max_geometry_lod = max(
            (lod for lod in range(4, -1, -1) if geometry_lod_counts[f"lod{lod}"]),
            default=0,
        )
        source = plan.source_metadata
        now = datetime.now(timezone.utc).isoformat()
        preparation_year = source.get("preparation_year") or None
        metadata = {
            "metadata_schema_version": CACHE_METADATA_SCHEMA_VERSION,
            "derived_tile_schema_version": DERIVED_TILE_SCHEMA_VERSION,
            "source": "PLATEAU-CityGML",
            "source_spec": source.get("source_spec"),
            "city_code": city_code,
            "year": source.get("year", preparation_year),
            "preparation_year": source.get("preparation_year"),
            "registration_year": source.get("registration_year"),
            "source_file_size_bytes": source.get("source_file_size_bytes"),
            "source_file_count": source.get("source_file_count"),
            "source_zip_sha256": zip_info.sha256,
            "source_zip_size_bytes": zip_info.size_bytes,
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
            "fetched_at_utc": now if downloaded else source.get("fetched_at_utc"),
            "prepared_at_utc": now,
            "converter": "zstarview-plateau-buildings",
            "converter_version": __version__,
            "tile_count": written_tiles,
            "building_count": building_count,
            "min_building_height_m": DEFAULT_MIN_BUILDING_HEIGHT_M,
            "catalog_url": source.get("catalog_url"),
            "archive_url": source.get("archive_url"),
            "archive_size_bytes": (
                zip_info.size_bytes if downloaded else source.get("archive_size_bytes")
            ),
            "status": "complete",
        }
        (staging_dir / "cache_meta.json").write_text(
            json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        if args.keep_zip:
            shutil.copyfile(zip_path, staging_dir / "source-citygml.zip")
        if not is_valid_plateau_cache(staging_dir):
            raise ValueError("Prepared PLATEAU cache failed verification")
        _install_cache(output_root, city_code, staging_dir)
        print(f"Tiles: {written_tiles}; buildings: {building_count}")
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir, ignore_errors=True)


def _execute_plan(args: argparse.Namespace, plan: _PreparePlan) -> None:
    output_root = Path(args.output_root)
    if plan.action == "skip":
        print(f"Skipping city code {plan.city_code}: cache is up to date.")
        for path in _leftover_dirs(output_root, plan.city_code):
            _remove_tree(path)
        return
    if plan.action == "migrate":
        assert plan.current is not None
        _install_cache(output_root, plan.city_code, plan.current)
        return
    if plan.action == "convert":
        assert plan.zip_path is not None and plan.zip_info is not None
        _build_cache(
            args,
            plan,
            zip_path=plan.zip_path,
            zip_info=plan.zip_info,
            downloaded=False,
        )
        return
    assert plan.catalog is not None
    archive_url = catalog_archive_url(plan.catalog)
    if archive_url is None:
        raise ValueError(
            "PLATEAU catalog has no CityGML archive URL for city code "
            f"{plan.city_code}"
        )
    temp_dir = args.temp_dir.expanduser() if args.temp_dir is not None else None
    if temp_dir is not None:
        temp_dir.mkdir(parents=True, exist_ok=True)
    _check_working_space_before_download(
        temp_base_dir=(temp_dir or Path(tempfile.gettempdir())),
        archive_size_bytes=(
            plan.size_bytes if plan.size_kind == "CityGML ZIP" else None
        ),
        entries=plan.entries,
    )
    temporary_dir = Path(
        tempfile.mkdtemp(
            prefix="plateau-download-",
            dir=str(temp_dir) if temp_dir is not None else None,
        )
    )
    try:
        zip_path = temporary_dir / "citygml.zip"
        zip_info = download_file(
            archive_url,
            zip_path,
            max_bytes=max(1, int(args.max_download_bytes)),
        )
        print(f"Downloaded: {zip_info.size_bytes} bytes")
        plan.source_metadata = dict(plan.source_metadata)
        plan.source_metadata["catalog_url"] = DEFAULT_CATALOG_URL.format(
            city_code=plan.city_code
        )
        plan.source_metadata["archive_url"] = archive_url
        plan.source_metadata["year"] = (
            plan.source_metadata.get("preparation_year") or None
        )
        _build_cache(
            args, plan, zip_path=zip_path, zip_info=zip_info, downloaded=True
        )
    finally:
        shutil.rmtree(temporary_dir, ignore_errors=True)


def _missing_catalog_message(city_code: str) -> str:
    return (
        f"PLATEAU catalog not found for city code {city_code} (HTTP 404); "
        "PLATEAU building data may not be available for this municipality."
    )


def _prepare_cli(args: argparse.Namespace) -> int:
    try:
        city_codes = parse_city_codes(",".join(args.city_codes))
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    if args.input_zip is not None and len(city_codes) != 1:
        print("Error: --input-zip can only be used with one city code", file=sys.stderr)
        return 2
    plans: list[_PreparePlan] = []
    for city_code in city_codes:
        try:
            if args.input_zip is not None:
                plans.append(_plan_zip(args, city_code))
            else:
                plans.append(_plan_catalog(args, city_code))
        except HTTPError as exc:
            if exc.code != 404:
                raise
            if len(city_codes) == 1:
                print(_missing_catalog_message(city_code), file=sys.stderr)
                return 1
            print(f"Skipping city code {city_code}: " + _missing_catalog_message(city_code))
        except (ValueError, FileNotFoundError) as exc:
            if len(city_codes) == 1:
                print(f"Error: {exc}", file=sys.stderr)
                return 1
            print(f"Skipping city code {city_code}: {exc}")
    if not plans:
        print(
            "Error: No PLATEAU building catalogs found for the requested city codes",
            file=sys.stderr,
        )
        return 1
    _print_plans(plans)
    if any(plan.needs_confirmation for plan in plans):
        if not _confirm("Continue?", assume_yes=bool(args.yes)):
            print("Cancelled.")
            return 0
    for index, plan in enumerate(plans, start=1):
        if len(plans) > 1:
            print(
                f"Preparing PLATEAU city code {plan.city_code} ({index}/{len(plans)})"
            )
        try:
            _execute_plan(args, plan)
        except (ValueError, BadZipFile, URLError) as exc:
            print(
                f"PLATEAU preparation failed for city code {plan.city_code}: {exc}. "
                "The existing cache was kept.",
                file=sys.stderr,
            )
            return 1
        except OSError as exc:
            if exc.errno != errno.ENOSPC:
                raise
            print(
                "PLATEAU preparation failed: no space left on the temporary "
                "filesystem. The existing cache was kept. Try again with "
                "--temp-dir, for example:",
                file=sys.stderr,
            )
            print(
                "  zstarview-download-plateau-buildings prepare "
                f'{" ".join(args.city_codes)} --temp-dir "$HOME/zstarview-tmp"',
                file=sys.stderr,
            )
            return 1
        except InsufficientPlateauDiskSpaceError as exc:
            print(f"PLATEAU preparation stopped before download: {exc}", file=sys.stderr)
            print(
                "Use --temp-dir to select a filesystem with more free space, "
                "or free disk space and try again.",
                file=sys.stderr,
            )
            return 1
    return 0


def _list_plateau_caches(output_root: Path) -> list[dict[str, object]]:
    """Describe every municipality directory under ``output_root``."""

    from zstarview.data.plateau_building_cache import (
        classify_plateau_caches,
        read_plateau_cache_metadata,
    )

    rows: list[dict[str, object]] = []
    for entry in classify_plateau_caches(output_root):
        if entry.selected is not None:
            row = dict(read_plateau_cache_metadata(entry.selected) or {})
            row["city_code"] = entry.city_code
            row["path"] = str(entry.selected)
            row["status"] = "legacy" if entry.legacy else "current"
            rows.append(row)
        for path in entry.removable:
            rows.append(
                {"city_code": entry.city_code, "path": str(path), "status": "removable"}
            )
    return rows


def _list_cli(args: argparse.Namespace) -> int:
    rows = _list_plateau_caches(Path(args.output_root))
    if args.jsonl:
        for row in rows:
            print(json.dumps(row, ensure_ascii=True, sort_keys=True))
        return 0
    for row in rows:
        year = row.get("preparation_year") or row.get("year") or "-"
        print(f"{row['city_code']} {year} {row['status']} {row['path']}")
    return 0


def _remove_cli(args: argparse.Namespace) -> int:
    from zstarview.data.plateau_building_cache import classify_plateau_caches

    if not args.city_codes and not args.unused:
        print(
            "Error: specify city codes, or use --unused to remove unused "
            "directories of all municipalities",
            file=sys.stderr,
        )
        return 2
    try:
        city_codes = (
            set(parse_city_codes(",".join(args.city_codes)))
            if args.city_codes
            else None
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    targets: list[Path] = []
    for entry in classify_plateau_caches(Path(args.output_root)):
        if city_codes is not None and entry.city_code not in city_codes:
            continue
        if not args.unused and entry.selected is not None:
            targets.append(entry.selected)
        targets.extend(entry.removable)
    if not targets:
        print("Nothing to remove.")
        return 0
    print("PLATEAU cache directories to remove:")
    for path in targets:
        print(f"  {path}")
    if not _confirm("Remove these directories?", assume_yes=bool(args.yes)):
        print("Cancelled.")
        return 0
    failed = [path for path in targets if not _remove_tree(path)]
    return 1 if failed else 0


def main(argv: Sequence[str] | None = None) -> int:
    argv_list = list(sys.argv[1:] if argv is None else argv)
    if argv_list and argv_list[0].startswith("--") and argv_list[0] != "--help":
        print(
            "Error: this command uses subcommands since zstarview 2.2.0: "
            "'prepare CITY_CODE...', 'list', 'remove'. "
            "--city-code, --list, --year, and --overwrite were removed; "
            "see --help.",
            file=sys.stderr,
        )
        return 2
    parser = build_arg_parser()
    args = parser.parse_args(argv_list)
    if args.command == "prepare":
        return _prepare_cli(args)
    if args.command == "list":
        return _list_cli(args)
    if args.command == "remove":
        return _remove_cli(args)
    parser.print_help(sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
