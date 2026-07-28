"""Download and prepare an AKARI far-infrared display asset."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
from pathlib import Path
from typing import Callable, Mapping, Protocol, cast

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from tqdm import tqdm

from .paths import CACHE_PATH
from .user_agent import build_user_agent

AKARI_DATASET = "akari-far-infrared-all-sky"
AKARI_RELEASE = "1"
AKARI_SCHEMA = 1
AKARI_SOURCE_BASE_URL = "https://lambda.gsfc.nasa.gov/data/foregrounds/akari/images"
DEFAULT_BANDS = ("90", "140", "160")
DEFAULT_WIDTH = 2048
DEFAULT_HEIGHT = 1024
DEFAULT_ZERO_RUN_MAX_WIDTH = 4
DEFAULT_ZERO_RUN_VALUE_FRACTION = 0.05
AKARI_BAND_FILENAMES = {
    "90": "akari_mollweide_WideS_1_4096.fits",
    "140": "akari_mollweide_WideL_1_4096.fits",
    "160": "akari_mollweide_160_1_4096.fits",
}
AKARI_CACHE_DIR = Path(CACHE_PATH) / "molecular-cloud"


class MolecularCloudDownloadError(RuntimeError):
    """Raised when an AKARI source cannot be downloaded or processed."""


class _HTTPResponse(Protocol):
    headers: Mapping[str, str]

    def __enter__(self) -> "_HTTPResponse": ...

    def __exit__(self, *args: object) -> None: ...

    def read(self, size: int = -1) -> bytes: ...


def format_binary_size(size_bytes: int) -> str:
    value = float(max(0, int(size_bytes)))
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return "0.00 B"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _download_file(
    url: str,
    destination: Path,
    *,
    timeout_s: float,
    description: str,
    urlopen: Callable[..., object] | None = None,
) -> int:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": build_user_agent("molecular-cloud-download")},
    )
    opener = urlopen or urllib.request.urlopen
    try:
        response = cast(_HTTPResponse, opener(request, timeout=timeout_s))
        with response, destination.open("wb") as output:
            header = response.headers.get("Content-Length")
            total = int(header) if header and header.isdigit() else None
            size = 0
            with tqdm(
                total=total,
                desc=description,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                ascii=True,
            ) as progress:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
                    size += len(chunk)
                    progress.update(len(chunk))
    except Exception as exc:
        destination.unlink(missing_ok=True)
        raise MolecularCloudDownloadError(f"failed to download {url}: {exc}") from exc
    if size <= 0:
        destination.unlink(missing_ok=True)
        raise MolecularCloudDownloadError(f"downloaded file is empty: {url}")
    return size


def _load_image(path: Path) -> tuple[np.ndarray, WCS]:
    try:
        with fits.open(path, memmap=True) as hdul:
            image_hdu = next(
                (
                    hdu
                    for hdu in hdul
                    if getattr(hdu, "data", None) is not None and np.asarray(hdu.data).ndim >= 2
                ),
                None,
            )
            if image_hdu is None:
                raise MolecularCloudDownloadError(f"no image HDU found: {path.name}")
            data = np.asarray(image_hdu.data, dtype=np.float32)
            while data.ndim > 2:
                data = data[0]
            header = image_hdu.header.copy()
        wcs = WCS(header).celestial
    except MolecularCloudDownloadError:
        raise
    except Exception as exc:
        raise MolecularCloudDownloadError(f"failed to read FITS {path.name}: {exc}") from exc
    if data.ndim != 2 or data.size == 0:
        raise MolecularCloudDownloadError(f"FITS image is not two-dimensional: {path.name}")
    ctypes = tuple(str(value).upper() for value in wcs.wcs.ctype)
    if not any("GLON" in value for value in ctypes) or not any("GLAT" in value for value in ctypes):
        raise MolecularCloudDownloadError(
            f"FITS image is not a Galactic-coordinate map: {path.name} ({ctypes})"
        )
    return data, wcs


def _sample_to_galactic_grid(
    data: np.ndarray,
    wcs: WCS,
    *,
    width: int,
    height: int,
) -> np.ndarray:
    lon = np.linspace(-180.0, 180.0, width, endpoint=False, dtype=np.float64)
    lat = np.linspace(90.0, -90.0, height, dtype=np.float64)
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    try:
        pixel_x, pixel_y = wcs.world_to_pixel_values(lon_grid, lat_grid)
    except Exception as exc:
        raise MolecularCloudDownloadError(f"failed to project Galactic grid: {exc}") from exc
    valid = (
        np.isfinite(pixel_x)
        & np.isfinite(pixel_y)
        & (pixel_x >= 0)
        & (pixel_y >= 0)
        & (pixel_x < data.shape[1])
        & (pixel_y < data.shape[0])
    )
    x = np.clip(np.rint(np.nan_to_num(pixel_x, nan=0.0)).astype(np.int64), 0, data.shape[1] - 1)
    y = np.clip(np.rint(np.nan_to_num(pixel_y, nan=0.0)).astype(np.int64), 0, data.shape[0] - 1)
    sampled = np.zeros((height, width), dtype=np.float32)
    sampled[valid] = data[y[valid], x[valid]]
    sampled[~np.isfinite(sampled)] = 0.0
    return sampled


def _repair_short_zero_runs(
    data: np.ndarray,
    *,
    max_width: int,
    value_threshold: float = 0.0,
) -> np.ndarray:
    """Interpolate short low-value runs, treating longitude as cyclic."""
    if max_width < 0:
        raise ValueError("max_width must be non-negative")
    if value_threshold < 0.0:
        raise ValueError("value_threshold must be non-negative")
    if max_width == 0:
        return data.copy()

    repaired = np.asarray(data, dtype=np.float32).copy()
    height, width = repaired.shape
    for row_index in range(height):
        row = repaired[row_index]
        zero = row <= value_threshold
        if not np.any(zero) or np.all(zero):
            continue
        starts = np.flatnonzero(zero & ~np.roll(zero, 1))
        ends = np.flatnonzero(zero & ~np.roll(zero, -1))
        for start, end in zip(starts.tolist(), ends.tolist()):
            run_width = (end - start) % width + 1
            if run_width > max_width:
                continue
            left = row[(start - 1) % width]
            right = row[(end + 1) % width]
            if left <= value_threshold or right <= value_threshold:
                continue
            positions = (start + np.arange(run_width)) % width
            fraction = (np.arange(run_width, dtype=np.float32) + 1.0) / (run_width + 1.0)
            row[positions] = left + (right - left) * fraction
    return repaired


def _normalize_display(data: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    finite = data[np.isfinite(data)]
    positive = finite[finite > 0]
    if positive.size == 0:
        return np.zeros(data.shape, dtype=np.uint16), {"low": 0.0, "high": 0.0}
    low = float(np.percentile(positive, 1.0))
    high = float(np.percentile(positive, 99.5))
    if not high > low:
        high = low + 1.0
    clipped = np.clip(data, low, high)
    normalized = np.clip((clipped - low) / (high - low), 0.0, 1.0)
    normalized[data <= 0] = 0.0
    # arcsinh preserves broad, faint structure without allowing bright clouds to dominate.
    stretched = np.arcsinh(8.0 * normalized) / np.arcsinh(8.0)
    return np.rint(stretched * 65535.0).astype(np.uint16), {"low": low, "high": high}


def _cache_root(cache_dir: str | os.PathLike[str] | None) -> Path:
    return (
        Path(cache_dir).expanduser()
        if cache_dir is not None
        else AKARI_CACHE_DIR
    ) / AKARI_DATASET / f"release-{AKARI_RELEASE}" / f"schema-{AKARI_SCHEMA}"


def prepare_akari_data(
    *,
    bands: tuple[str, ...] = DEFAULT_BANDS,
    cache_dir: str | os.PathLike[str] | None = None,
    source_base_url: str = AKARI_SOURCE_BASE_URL,
    width: int = DEFAULT_WIDTH,
    height: int = DEFAULT_HEIGHT,
    zero_run_max_width: int = DEFAULT_ZERO_RUN_MAX_WIDTH,
    zero_run_value_fraction: float = DEFAULT_ZERO_RUN_VALUE_FRACTION,
    timeout_s: float = 600.0,
    urlopen: Callable[..., object] | None = None,
) -> Path:
    """Download AKARI maps and write a compact Galactic-coordinate display cache."""
    if not bands:
        raise ValueError("at least one AKARI band is required")
    if any(band not in AKARI_BAND_FILENAMES for band in bands):
        raise ValueError("supported AKARI bands are 90, 140, and 160")
    if len(set(bands)) != len(bands):
        raise ValueError("AKARI bands must be unique")
    if width < 2 or height < 2:
        raise ValueError("output dimensions must be at least 2x2")
    if zero_run_max_width < 0:
        raise ValueError("zero_run_max_width must be non-negative")
    if not 0.0 <= zero_run_value_fraction <= 1.0:
        raise ValueError("zero_run_value_fraction must be between 0 and 1")

    root = _cache_root(cache_dir)
    root.parent.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(tempfile.mkdtemp(prefix=".akari-build-", dir=root.parent))
    try:
        arrays: list[np.ndarray] = []
        normalization: dict[str, Mapping[str, float]] = {}
        sources: list[dict[str, object]] = []
        for band in tqdm(bands, desc="Processing AKARI bands", unit="band", ascii=True):
            filename = AKARI_BAND_FILENAMES[band]
            url = f"{source_base_url.rstrip('/')}/{filename}"
            existing_source = root / filename
            if existing_source.is_file():
                source_path = existing_source
                size = source_path.stat().st_size
            else:
                source_path = temporary_root / filename
                size = _download_file(
                    url,
                    source_path,
                    timeout_s=timeout_s,
                    description=f"AKARI {band} um",
                    urlopen=urlopen,
                )
            digest = _sha256(source_path)
            data, wcs = _load_image(source_path)
            sampled = _sample_to_galactic_grid(data, wcs, width=width, height=height)
            normalized, parameters = _normalize_display(sampled)
            normalized = _repair_short_zero_runs(
                normalized,
                max_width=zero_run_max_width,
                value_threshold=65535.0 * zero_run_value_fraction,
            )
            normalized = np.clip(np.rint(normalized), 0, 65535).astype(np.uint16)
            arrays.append(normalized)
            normalization[band] = parameters
            sources.append(
                {
                    "band_um": int(band),
                    "url": url,
                    "filename": filename,
                    "bytes": size,
                    "sha256": digest,
                }
            )

        output_name = "akari-galactic-display.npz"
        output_path = temporary_root / output_name
        np.savez_compressed(
            output_path,
            bands=np.asarray([int(band) for band in bands], dtype=np.uint16),
            data=np.stack(arrays, axis=0),
            width=np.asarray(width, dtype=np.int32),
            height=np.asarray(height, dtype=np.int32),
            longitude_min=np.asarray(-180.0, dtype=np.float64),
            longitude_max=np.asarray(180.0, dtype=np.float64),
            latitude_min=np.asarray(-90.0, dtype=np.float64),
            latitude_max=np.asarray(90.0, dtype=np.float64),
        )
        manifest = {
            "schema": AKARI_SCHEMA,
            "dataset": AKARI_DATASET,
            "release": AKARI_RELEASE,
            "coordinate_system": "galactic_plate_carree",
            "projection": "equirectangular",
            "dimensions": {"width": width, "height": height},
            "bands_um": [int(band) for band in bands],
            "channel_order": list(bands),
            "encoding": "uint16_normalized_arcsinh",
            "zero_run_repair": {
                "max_width_pixels": zero_run_max_width,
                "value_fraction": zero_run_value_fraction,
                "method": "horizontal_linear_interpolation",
            },
            "normalization": normalization,
            "output": {"name": output_name, "bytes": output_path.stat().st_size, "sha256": _sha256(output_path)},
            "sources": sources,
            "attribution": {
                "source": "AKARI Far-Infrared All-Sky Survey Maps",
                "archive": "DARTS/ISAS/JAXA; mirrored by NASA LAMBDA",
                "source_url": "https://lambda.gsfc.nasa.gov/product/foreground/fg_akari_info.html",
                "note": "Display-oriented derived asset; not suitable for scientific measurement.",
            },
        }
        manifest_path = temporary_root / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        if root.exists():
            for child in temporary_root.iterdir():
                os.replace(child, root / child.name)
            temporary_root.rmdir()
        else:
            os.replace(temporary_root, root)
        return root
    except Exception:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
