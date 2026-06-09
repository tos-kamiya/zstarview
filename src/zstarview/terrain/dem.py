from __future__ import annotations

import json
import logging
import math
import socket
import threading
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import numpy as np
from pyproj import CRS, Geod, Transformer

from ..clouddisc.types import DownloadCancelledError
from ..user_agent import build_user_agent

WGS84_GEOD = Geod(ellps="WGS84")
COPERNICUS_DEM_BUCKET = "copernicus-dem-90m"
COPERNICUS_DEM_REGION = "eu-central-1"
DEM_CACHE_TTL_DAYS = 90
_DEM_CHUNK_SIZE = 1024 * 1024

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class DownloadedDem:
    paths: tuple[Path, ...]
    bbox_crs84: tuple[float, float, float, float]
    tile_keys: tuple[str, ...]
    source: str


@dataclass
class DemGrid:
    array: np.ndarray
    transform: object
    crs: object
    default_elevation_m: float | None
    _to_grid: Transformer

    def sample_lonlat(
        self,
        lon_deg: np.ndarray,
        lat_deg: np.ndarray,
        *,
        method: str = "bilinear",
    ) -> np.ndarray:
        lon_deg = np.asarray(lon_deg, dtype=np.float64)
        lat_deg = np.asarray(lat_deg, dtype=np.float64)
        sample_shape = lon_deg.shape
        if lon_deg.size == 1:
            x, y = self._to_grid.transform(float(lon_deg.reshape(())), float(lat_deg.reshape(())))
            x = np.asarray([x], dtype=np.float64)
            y = np.asarray([y], dtype=np.float64)
        else:
            x, y = self._to_grid.transform(lon_deg.ravel(), lat_deg.ravel())
        x = np.asarray(x, dtype=np.float64).reshape(sample_shape)
        y = np.asarray(y, dtype=np.float64).reshape(sample_shape)
        result = np.full(lon_deg.shape, np.nan, dtype=np.float64)
        inv_transform = ~self.transform
        cols_f, rows_f = inv_transform * (x.ravel(), y.ravel())
        cols_f = np.asarray(cols_f, dtype=np.float64).reshape(sample_shape) - 0.5
        rows_f = np.asarray(rows_f, dtype=np.float64).reshape(sample_shape) - 0.5

        inside = (
            np.isfinite(x)
            & np.isfinite(y)
            & (rows_f >= 0.0)
            & (rows_f <= self.array.shape[0] - 1)
            & (cols_f >= 0.0)
            & (cols_f <= self.array.shape[1] - 1)
        )

        if method == "nearest":
            rows = np.rint(rows_f).astype(np.int64)
            cols = np.rint(cols_f).astype(np.int64)
            if np.any(inside):
                result[inside] = self.array[rows[inside], cols[inside]]
        elif method == "bilinear":
            bilinear_mask = (
                inside
                & (rows_f < self.array.shape[0] - 1)
                & (cols_f < self.array.shape[1] - 1)
            )
            if np.any(bilinear_mask):
                row0 = np.floor(rows_f[bilinear_mask]).astype(np.int64)
                col0 = np.floor(cols_f[bilinear_mask]).astype(np.int64)
                row1 = row0 + 1
                col1 = col0 + 1
                dy = rows_f[bilinear_mask] - row0
                dx = cols_f[bilinear_mask] - col0
                v00 = self.array[row0, col0]
                v01 = self.array[row0, col1]
                v10 = self.array[row1, col0]
                v11 = self.array[row1, col1]
                result[bilinear_mask] = (
                    v00 * (1.0 - dx) * (1.0 - dy)
                    + v01 * dx * (1.0 - dy)
                    + v10 * (1.0 - dx) * dy
                    + v11 * dx * dy
                )

            edge_mask = inside & ~bilinear_mask
            if np.any(edge_mask):
                rows = np.rint(rows_f[edge_mask]).astype(np.int64)
                cols = np.rint(cols_f[edge_mask]).astype(np.int64)
                result[edge_mask] = self.array[rows, cols]
        else:
            raise ValueError(f"Unsupported DEM resampling method: {method}")

        if self.default_elevation_m is not None:
            result = np.where(np.isnan(result), float(self.default_elevation_m), result)
        return result


class GeoTiffDem:
    """Load a DEM mosaic once, then sample it from NumPy."""

    def __init__(
        self,
        dataset_paths: Sequence[Path],
        *,
        band: int = 1,
        default_elevation_m: float | None = None,
    ) -> None:
        try:
            import rasterio
        except ImportError as exc:  # pragma: no cover - depends on local setup
            raise RuntimeError(
                "GeoTIFF support requires 'rasterio'. Install it in the current environment and retry."
            ) from exc

        self._datasets = [rasterio.open(path) for path in dataset_paths]
        if not self._datasets:
            raise ValueError("At least one DEM GeoTIFF is required.")
        crs_values = {str(ds.crs) for ds in self._datasets}
        if len(crs_values) != 1:
            raise ValueError("All DEM GeoTIFFs must share the same CRS.")
        self._band = band
        self._default_elevation_m = default_elevation_m
        self._crs = self._datasets[0].crs

    def close(self) -> None:
        for dataset in self._datasets:
            dataset.close()

    def build_grid(self, bbox_wgs84: tuple[float, float, float, float]) -> DemGrid:
        import rasterio.merge
        import rasterio.warp

        bounds = rasterio.warp.transform_bounds(
            "EPSG:4326",
            self._crs,
            *bbox_wgs84,
            densify_pts=21,
        )
        merged, transform = rasterio.merge.merge(
            self._datasets,
            bounds=bounds,
            indexes=self._band,
            nodata=np.nan,
            dtype="float32",
        )
        array = merged[0].astype(np.float64, copy=False)
        return DemGrid(
            array=array,
            transform=transform,
            crs=self._crs,
            default_elevation_m=self._default_elevation_m,
            _to_grid=Transformer.from_crs(CRS.from_epsg(4326), self._crs, always_xy=True),
        )


def build_download_bbox(
    *,
    lat_deg: float,
    lon_deg: float,
    radius_km: float,
) -> tuple[float, float, float, float]:
    radius_m = radius_km * 1000.0
    lon_n, lat_n, _ = WGS84_GEOD.fwd(lon_deg, lat_deg, 0.0, radius_m)
    lon_s, lat_s, _ = WGS84_GEOD.fwd(lon_deg, lat_deg, 180.0, radius_m)
    lon_e, lat_e, _ = WGS84_GEOD.fwd(lon_deg, lat_deg, 90.0, radius_m)
    lon_w, lat_w, _ = WGS84_GEOD.fwd(lon_deg, lat_deg, 270.0, radius_m)
    min_lon = min(lon_w, lon_s, lon_n, lon_e)
    max_lon = max(lon_w, lon_s, lon_n, lon_e)
    min_lat = min(lat_s, lat_w, lat_e, lat_n)
    max_lat = max(lat_s, lat_w, lat_e, lat_n)
    return (min_lon, min_lat, max_lon, max_lat)


def collect_copernicus_tile_keys(
    bbox: tuple[float, float, float, float]
) -> list[str]:
    min_lon, min_lat, max_lon, max_lat = bbox
    min_lon_tile = math.floor(min_lon)
    max_lon_tile = math.ceil(max_lon) - 1
    min_lat_tile = math.floor(min_lat)
    max_lat_tile = math.ceil(max_lat) - 1
    keys: list[str] = []
    for lat_tile in range(min_lat_tile, max_lat_tile + 1):
        for lon_tile in range(min_lon_tile, max_lon_tile + 1):
            keys.append(copernicus_tile_key(lat_tile, lon_tile))
    return keys


def copernicus_tile_key(lat_tile: int, lon_tile: int) -> str:
    lat_prefix = "N" if lat_tile >= 0 else "S"
    lon_prefix = "E" if lon_tile >= 0 else "W"
    lat_label = f"{lat_prefix}{abs(lat_tile):02d}_00"
    lon_label = f"{lon_prefix}{abs(lon_tile):03d}_00"
    stem = f"Copernicus_DSM_COG_30_{lat_label}_{lon_label}_DEM"
    return f"{stem}/{stem}.tif"


def dem_tile_metadata_path(tile_path: Path) -> Path:
    return tile_path.with_suffix(tile_path.suffix + ".meta.json")


def read_dem_tile_fetched_at_utc(
    tile_path: Path,
    *,
    now_utc: datetime | None = None,
    migrate_missing: bool = True,
) -> datetime | None:
    if not tile_path.exists():
        return None
    if not _is_valid_dem_tile(tile_path):
        logger.warning("Discarding invalid DEM cache tile: %s", tile_path)
        _discard_dem_tile(tile_path)
        return None
    metadata_path = dem_tile_metadata_path(tile_path)
    payload: dict[str, object] = {}
    if metadata_path.exists():
        try:
            loaded = json.loads(metadata_path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                payload = loaded
            raw_fetched_at = payload.get("fetched_at_utc")
            if isinstance(raw_fetched_at, str):
                return _parse_utc(raw_fetched_at)
        except Exception:
            payload = {}
    if not migrate_missing:
        return None
    fetched_at_utc = _normalize_utc(now_utc or datetime.now(timezone.utc))
    logger.info(
        "DEM cache metadata missing; recording provisional fetched_at_utc for offline reuse: %s",
        tile_path,
    )
    payload["fetched_at_utc"] = fetched_at_utc.isoformat()
    _write_dem_tile_metadata(tile_path, payload)
    return fetched_at_utc


def is_dem_tile_stale(
    tile_path: Path,
    *,
    ttl_days: int = DEM_CACHE_TTL_DAYS,
    now_utc: datetime | None = None,
) -> bool:
    fetched_at_utc = read_dem_tile_fetched_at_utc(tile_path, now_utc=now_utc)
    if fetched_at_utc is None:
        return True
    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    return (now - fetched_at_utc) > timedelta(days=max(0, int(ttl_days)))


def _write_dem_tile_metadata(tile_path: Path, payload: dict[str, object]) -> None:
    metadata_path = dem_tile_metadata_path(tile_path)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")


def _discard_dem_tile(tile_path: Path) -> None:
    tile_path.unlink(missing_ok=True)
    dem_tile_metadata_path(tile_path).unlink(missing_ok=True)


def _is_valid_dem_tile(tile_path: Path) -> bool:
    if not tile_path.exists():
        return False
    try:
        import rasterio
    except ImportError:
        return True

    try:
        with rasterio.open(tile_path) as dataset:
            return dataset.count > 0
    except Exception:
        return False


def _parse_utc(text: str) -> datetime:
    return _normalize_utc(datetime.fromisoformat(str(text).replace("Z", "+00:00")))


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _copernicus_dem_url(key: str) -> str:
    quoted_key = quote(key, safe="/-_.~")
    return f"https://{COPERNICUS_DEM_BUCKET}.s3.{COPERNICUS_DEM_REGION}.amazonaws.com/{quoted_key}"


def _request_timeout(timeout_s: float | None) -> float | None:
    if timeout_s is None:
        return None
    timeout = float(timeout_s)
    return timeout if timeout > 0.0 else 0.0


def _download_dem_tile(
    *,
    key: str,
    dst: Path,
    abort_event: threading.Event | None,
    timeout_s: float | None = 120.0,
) -> None:
    req = Request(
        _copernicus_dem_url(key),
        method="GET",
        headers={"User-Agent": build_user_agent("copernicus-dem")},
    )
    with urlopen(req, timeout=_request_timeout(timeout_s)) as resp, dst.open("wb") as handle:
        while True:
            if abort_event is not None and abort_event.is_set():
                raise KeyboardInterrupt()
            chunk = resp.read(_DEM_CHUNK_SIZE)
            if not chunk:
                break
            handle.write(chunk)


def fetch_copernicus_dem(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    max_distance_km: float,
    margin_km: float,
    cache_dir: Path,
    ttl_days: int = DEM_CACHE_TTL_DAYS,
    now_utc: datetime | None = None,
    abort_event: threading.Event | None = None,
) -> DownloadedDem:
    if margin_km < 0.0:
        raise ValueError("margin_km must be non-negative.")

    now = _normalize_utc(now_utc or datetime.now(timezone.utc))
    radius_km = max_distance_km + margin_km
    bbox = build_download_bbox(
        lat_deg=observer_lat_deg,
        lon_deg=observer_lon_deg,
        radius_km=radius_km,
    )
    tile_keys = collect_copernicus_tile_keys(bbox)
    cache_dir.mkdir(parents=True, exist_ok=True)
    downloaded_paths: list[Path] = []
    downloaded_any = False
    stale_fallback_used = False

    for key in tile_keys:
        if abort_event is not None and abort_event.is_set():
            raise DownloadCancelledError("Cancelled while downloading Copernicus DEM tiles")
        dst = cache_dir / key
        if dst.exists():
            if not _is_valid_dem_tile(dst):
                logger.warning("Discarding invalid DEM cache tile: %s", dst)
                _discard_dem_tile(dst)
            elif not is_dem_tile_stale(dst, ttl_days=ttl_days, now_utc=now):
                downloaded_paths.append(dst)
                continue
            logger.info(
                "Refreshing expired DEM cache tile: %s (ttl=%s days)",
                dst,
                int(ttl_days),
            )
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = dst.with_suffix(dst.suffix + ".tmp")
        try:
            _download_dem_tile(key=key, dst=tmp_path, abort_event=abort_event)
            if not _is_valid_dem_tile(tmp_path):
                logger.warning(
                    "Discarding invalid downloaded DEM tile: s3://%s/%s",
                    COPERNICUS_DEM_BUCKET,
                    key,
                )
                continue
            tmp_path.replace(dst)
            _write_dem_tile_metadata(
                dst,
                {
                    "bucket": COPERNICUS_DEM_BUCKET,
                    "fetched_at_utc": now.isoformat(),
                    "key": key,
                },
            )
            downloaded_paths.append(dst)
            downloaded_any = True
        except (KeyboardInterrupt, DownloadCancelledError) as exc:
            raise DownloadCancelledError("Cancelled while downloading Copernicus DEM tiles") from exc
        except HTTPError as exc:
            if exc.code == 404:
                if dst.exists():
                    logger.warning("DEM refresh failed; using stale cached tile: %s", dst)
                    downloaded_paths.append(dst)
                    stale_fallback_used = True
                continue
            if dst.exists():
                logger.warning("DEM refresh failed; using stale cached tile: %s", dst)
                downloaded_paths.append(dst)
                stale_fallback_used = True
                continue
            raise RuntimeError(
                f"Failed to download s3://{COPERNICUS_DEM_BUCKET}/{key}: {exc}"
            ) from exc
        except URLError as exc:
            reason = getattr(exc, "reason", None)
            if isinstance(reason, (TimeoutError, socket.timeout)):
                if dst.exists():
                    logger.warning("DEM refresh failed; using stale cached tile: %s", dst)
                    downloaded_paths.append(dst)
                    stale_fallback_used = True
                    continue
                raise RuntimeError(
                    f"Failed to download s3://{COPERNICUS_DEM_BUCKET}/{key}: {exc}"
                ) from exc
            if dst.exists():
                logger.warning("DEM refresh failed; using stale cached tile: %s", dst)
                downloaded_paths.append(dst)
                stale_fallback_used = True
                continue
            raise RuntimeError(
                f"Failed to download s3://{COPERNICUS_DEM_BUCKET}/{key}: {exc}"
            ) from exc
        except Exception as exc:
            message = str(exc)
            if "404" in message or "Not Found" in message or "NoSuchKey" in message:
                if dst.exists():
                    logger.warning("DEM refresh failed; using stale cached tile: %s", dst)
                    downloaded_paths.append(dst)
                    stale_fallback_used = True
            else:
                if dst.exists():
                    logger.warning("DEM refresh failed; using stale cached tile: %s", dst)
                    downloaded_paths.append(dst)
                    stale_fallback_used = True
                    continue
                raise RuntimeError(
                    f"Failed to download s3://{COPERNICUS_DEM_BUCKET}/{key}: {exc}"
                ) from exc
        finally:
            tmp_path.unlink(missing_ok=True)

    if not downloaded_paths:
        raise RuntimeError(
            "No Copernicus DEM tiles were downloaded for the requested area."
        )

    return DownloadedDem(
        paths=tuple(downloaded_paths),
        bbox_crs84=bbox,
        tile_keys=tuple(tile_keys),
        source="download" if downloaded_any else ("cache-stale" if stale_fallback_used else "cache"),
    )


def sample_ground_elevation(
    dem_grid: DemGrid,
    *,
    latitude_deg: float,
    longitude_deg: float,
    dem_resampling: str = "bilinear",
) -> float:
    ground_m = float(
        dem_grid.sample_lonlat(
            np.array([longitude_deg], dtype=np.float64),
            np.array([latitude_deg], dtype=np.float64),
            method=dem_resampling,
        )[0]
    )
    if not math.isfinite(ground_m):
        raise ValueError("Observer location is outside the DEM or falls on nodata.")
    return ground_m
