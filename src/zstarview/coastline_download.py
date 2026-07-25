"""Download and install coastline vector Release assets."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath
from typing import Callable

from tqdm import tqdm

from .paths import CACHE_PATH
from .user_agent import build_user_agent

COASTLINE_DATASET_VERSION = "20260725"
COASTLINE_SCHEMA_VERSION = 1
COASTLINE_RELEASE_TAG = "coastline-data-20260725"
COASTLINE_RELEASE_BASE_URL = (
    "https://github.com/tos-kamiya/zstarview/releases/download/"
    f"{COASTLINE_RELEASE_TAG}"
)
COASTLINE_CACHE_DIR = Path(CACHE_PATH) / "coastline"
GRID_COLS = 32
GRID_ROWS = 16
TILE_WIDTH_DEG = 360.0 / GRID_COLS
_SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
_ASSET_RE = re.compile(r"^coastline-grid-x(?P<x>\d{2})-20260725\.zip$")


class CoastlineDownloadError(RuntimeError):
    """Raised when coastline Release data cannot be downloaded or installed."""


def format_binary_size(size_bytes: int) -> str:
    value = float(max(0, int(size_bytes)))
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return "0.00 B"


def dataset_root(cache_dir: str | os.PathLike[str] | None = None) -> Path:
    return (
        Path(cache_dir or COASTLINE_CACHE_DIR).expanduser()
        / "osm-water-polygons"
        / COASTLINE_DATASET_VERSION
        / "schema-1"
    )


def select_columns(
    *, lon_min: float | None = None, lon_max: float | None = None, all_columns: bool = False
) -> tuple[int, ...]:
    if all_columns:
        if lon_min is not None or lon_max is not None:
            raise ValueError("--all cannot be combined with longitude bounds")
        return tuple(range(GRID_COLS))
    if lon_min is None or lon_max is None:
        raise ValueError("both longitude bounds are required")
    if not -180.0 <= lon_min <= 180.0 or not -180.0 <= lon_max <= 180.0:
        raise ValueError("longitude must be between -180 and 180 degrees")
    if lon_min > lon_max:
        raise ValueError("lon-min must not be greater than lon-max")
    start = max(0, min(GRID_COLS - 1, int((lon_min + 180.0) // TILE_WIDTH_DEG)))
    end_lon = min(lon_max, 180.0 - 1.0e-12)
    end = max(0, min(GRID_COLS - 1, int((end_lon + 180.0) // TILE_WIDTH_DEG)))
    return tuple(range(start, end + 1))


def _read_url(url: str, *, timeout_s: float, accept: str = "*/*") -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": build_user_agent("coastline-download"), "Accept": accept},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_s) as response:
            return response.read()
    except Exception as exc:
        raise CoastlineDownloadError(f"failed to download {url}: {exc}") from exc


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def validate_manifest(payload: object) -> dict[str, object]:
    if not isinstance(payload, dict):
        raise CoastlineDownloadError("coastline manifest must be a JSON object")
    if payload.get("schema") != COASTLINE_SCHEMA_VERSION:
        raise CoastlineDownloadError("unsupported coastline manifest schema")
    if payload.get("dataset") != "coastline-vector-tiles":
        raise CoastlineDownloadError("unexpected coastline dataset")
    if payload.get("version") != COASTLINE_DATASET_VERSION:
        raise CoastlineDownloadError("unexpected coastline dataset version")
    grid = payload.get("grid")
    if not isinstance(grid, dict) or grid.get("columns") != GRID_COLS or grid.get("rows") != GRID_ROWS:
        raise CoastlineDownloadError("invalid coastline grid metadata")
    assets = payload.get("assets")
    if not isinstance(assets, list) or len(assets) != GRID_COLS:
        raise CoastlineDownloadError("coastline manifest must contain 32 assets")
    seen: set[int] = set()
    for asset in assets:
        if not isinstance(asset, dict):
            raise CoastlineDownloadError("invalid coastline asset entry")
        column = asset.get("x")
        name = asset.get("name")
        size = asset.get("bytes")
        digest = asset.get("sha256")
        if not isinstance(column, int) or not 0 <= column < GRID_COLS or column in seen:
            raise CoastlineDownloadError("invalid coastline asset column")
        match = _ASSET_RE.fullmatch(name) if isinstance(name, str) else None
        if match is None:
            raise CoastlineDownloadError("invalid coastline asset name")
        if int(match["x"]) != column:
            raise CoastlineDownloadError("coastline asset name and column disagree")
        if not isinstance(size, int) or size <= 0:
            raise CoastlineDownloadError("invalid coastline asset size")
        if not isinstance(digest, str) or _SHA256_RE.fullmatch(digest) is None:
            raise CoastlineDownloadError("invalid coastline asset SHA-256")
        seen.add(column)
    return payload


def _safe_member_name(name: str, column: int) -> PurePosixPath:
    if "\\" in name:
        raise CoastlineDownloadError("ZIP entry contains a backslash")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts:
        raise CoastlineDownloadError("ZIP entry has an unsafe path")
    parts = path.parts
    if len(parts) < 4 or parts[0] != "grid-32x16" or not re.fullmatch(r"y\d{2}", parts[1]):
        raise CoastlineDownloadError("ZIP entry is outside the coastline grid")
    if parts[2] != f"x{column:02d}":
        raise CoastlineDownloadError("ZIP entry is in the wrong longitude column")
    return path


def _extract_checked(archive: Path, destination: Path, column: int) -> None:
    seen: set[str] = set()
    with zipfile.ZipFile(archive) as zip_handle:
        for info in zip_handle.infolist():
            path = _safe_member_name(info.filename, column)
            if info.is_dir():
                continue
            if path.as_posix() in seen:
                raise CoastlineDownloadError("ZIP contains duplicate entries")
            seen.add(path.as_posix())
            mode = (info.external_attr >> 16) & 0o170000
            if mode == 0o120000:
                raise CoastlineDownloadError("ZIP contains a symbolic link")
            target = destination.joinpath(*path.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with zip_handle.open(info) as source, target.open("wb") as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)


def _replace_directory(source: Path, target: Path) -> None:
    """Replace a cache column, including when the old directory is non-empty."""
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        os.replace(source, target)
        return
    backup_root = Path(tempfile.mkdtemp(prefix=f".{target.name}.backup-", dir=target.parent))
    backup = backup_root / target.name
    try:
        os.replace(target, backup)
        try:
            os.replace(source, target)
        except Exception:
            os.replace(backup, target)
            raise
    finally:
        shutil.rmtree(backup_root, ignore_errors=True)


def _download_asset(
    *,
    url: str,
    expected_size: int,
    expected_sha256: str,
    temporary_dir: Path,
    timeout_s: float,
    description: str,
    read_url: Callable[..., bytes] | None = None,
) -> Path:
    if read_url is not None:
        payload = read_url(url, timeout_s=timeout_s)
        actual_sha256 = _sha256_bytes(payload)
        if len(payload) != expected_size or actual_sha256.lower() != expected_sha256.lower():
            raise CoastlineDownloadError("coastline asset size or SHA-256 mismatch")
        path = temporary_dir / "asset.zip"
        path.write_bytes(payload)
        return path
    request = urllib.request.Request(
        url,
        headers={"User-Agent": build_user_agent("coastline-download")},
    )
    path = temporary_dir / "asset.zip"
    digest = hashlib.sha256()
    size = 0
    try:
        with urllib.request.urlopen(request, timeout=timeout_s) as response, path.open("wb") as output, tqdm(
            total=expected_size,
            desc=description,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            ascii=True,
        ) as progress:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
                digest.update(chunk)
                size += len(chunk)
                progress.update(len(chunk))
    except Exception as exc:
        raise CoastlineDownloadError(f"failed to download coastline asset: {exc}") from exc
    if size != expected_size or digest.hexdigest().lower() != expected_sha256.lower():
        path.unlink(missing_ok=True)
        raise CoastlineDownloadError("coastline asset size or SHA-256 mismatch")
    return path


def download_coastline_data(
    *,
    lon_min: float | None = None,
    lon_max: float | None = None,
    all_columns: bool = False,
    cache_dir: str | os.PathLike[str] | None = None,
    base_url: str = COASTLINE_RELEASE_BASE_URL,
    timeout_s: float = 60.0,
    download_timeout_s: float = 600.0,
    read_url: Callable[..., bytes] | None = None,
    status_callback: Callable[[str], None] | None = None,
) -> tuple[int, ...]:
    columns = select_columns(lon_min=lon_min, lon_max=lon_max, all_columns=all_columns)
    manifest_url = f"{base_url.rstrip('/')}/coastline-manifest-{COASTLINE_DATASET_VERSION}.json"
    if read_url is None:
        manifest_payload = _read_url(manifest_url, timeout_s=timeout_s, accept="application/json")
    else:
        manifest_payload = read_url(manifest_url, timeout_s=timeout_s)
    try:
        manifest = validate_manifest(json.loads(manifest_payload.decode("utf-8")))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CoastlineDownloadError(f"invalid coastline manifest: {exc}") from exc
    assets = {int(asset["x"]): asset for asset in manifest["assets"] if isinstance(asset, dict)}
    if status_callback is not None:
        estimated_size = sum(int(assets[column]["bytes"]) for column in columns)
        status_callback(
            "Coastline columns: " + ", ".join(f"x{column:02d}" for column in columns)
        )
        status_callback(f"Estimated download size: {format_binary_size(estimated_size)}")
    root = dataset_root(cache_dir)
    root.mkdir(parents=True, exist_ok=True)
    (root / "grid-32x16").mkdir(exist_ok=True)
    # A failed multi-column update must not leave an old READY marker that
    # makes the partially updated cache look complete to the reader.
    (root / "READY").unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix="coastline-download-", dir=root) as temporary:
        temporary_dir = Path(temporary)
        for column in columns:
            asset = assets[column]
            name = str(asset["name"])
            url = f"{base_url.rstrip('/')}/{name}"
            archive = _download_asset(
                url=url,
                expected_size=int(asset["bytes"]),
                expected_sha256=str(asset["sha256"]),
                temporary_dir=temporary_dir,
                timeout_s=download_timeout_s,
                description=f"Downloading coastline x{column:02d}",
                read_url=read_url,
            )
            if status_callback is not None:
                status_callback(f"Verified coastline x{column:02d}; extracting")
            extracted = temporary_dir / f"x{column:02d}"
            extracted.mkdir()
            _extract_checked(archive, extracted, column)
            source_grid = extracted / "grid-32x16"
            for row_dir in sorted(source_grid.glob("y*/")):
                source_column = row_dir / f"x{column:02d}"
                if not source_column.is_dir():
                    continue
                target_row = root / "grid-32x16" / row_dir.name
                target_row.mkdir(parents=True, exist_ok=True)
                _replace_directory(source_column, target_row / source_column.name)
            if status_callback is not None:
                status_callback(f"Installed coastline x{column:02d}")
    local_manifest = dict(manifest)
    local_manifest["manifest_url"] = manifest_url
    (root / "manifest.json").write_text(
        json.dumps(local_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (root / "READY").write_text(
        "coastline-data-20260725\n" + " ".join(f"x{column:02d}" for column in columns) + "\n",
        encoding="ascii",
    )
    if status_callback is not None:
        status_callback(f"Coastline cache ready: {root}")
    return columns
