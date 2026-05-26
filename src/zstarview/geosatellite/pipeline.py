from __future__ import annotations

import datetime as dt
import io
import logging
from dataclasses import replace
from pathlib import Path

import numpy as np
from PIL import Image

from ..paths import GEOSATELLITE_GRAY_COMMON_MASK_FILE
from .cache import (
    AVAILABLE_CACHE_MAX_AGE_SECONDS,
    compute_digest,
    latest_available_cache_is_recent,
    read_latest_available_cache,
    read_raw_cache,
    purge_intermediate_cache,
    write_latest_available_cache,
    write_raw_cache,
)
from .client import download_latest_image, fetch_latest_available_image_time
from .mask import DEFAULT_GRAY_SPREAD, DEFAULT_WHITE_BRIGHTNESS, fill_masked_regions, load_default_mask
from .projection import (
    DEFAULT_CLOUD_HEIGHT_KM,
    DEFAULT_GRID_NPZ,
    DEFAULT_RENDER_RADIUS_PX,
    EUROPE_MAX_LAT_DEG,
    EUROPE_MAX_LON_DEG,
    EUROPE_MIN_LAT_DEG,
    EUROPE_MIN_LON_DEG,
    project_gray_image_to_disc,
)
from .proxy import DEFAULT_CONTRAST_HIGH, DEFAULT_CONTRAST_LOW, DEFAULT_LOGO_MASK_FRAC, build_cloud_proxy
from .types import GeoSatelliteDownloadResult, GeoSatelliteIntermediateResult, GeoSatelliteKind, GeoSatellitePipelineResult

logger = logging.getLogger(__name__)


def is_within_europe_band(lat: float, lon: float) -> bool:
    return EUROPE_MIN_LAT_DEG <= float(lat) <= EUROPE_MAX_LAT_DEG and EUROPE_MIN_LON_DEG <= float(lon) <= EUROPE_MAX_LON_DEG


def _load_png_image(png_bytes: bytes) -> Image.Image:
    with Image.open(io.BytesIO(png_bytes)) as image:
        return image.convert("RGB")


def _load_mask_image(mask_path: Path | None) -> Image.Image:
    if mask_path is None:
        return load_default_mask()
    with Image.open(mask_path) as image:
        return image.convert("L")


def _mask_digest(mask_image: Image.Image) -> str:
    return compute_digest(np.asarray(mask_image, dtype=np.uint8).tobytes())


def build_proxy_image(
    image: Image.Image,
    *,
    mask_logo: bool = True,
    mask_frac: tuple[float, float] = DEFAULT_LOGO_MASK_FRAC,
    contrast_low: float = DEFAULT_CONTRAST_LOW,
    contrast_high: float = DEFAULT_CONTRAST_HIGH,
) -> Image.Image:
    return build_cloud_proxy(
        image,
        mask_logo=mask_logo,
        mask_frac=mask_frac,
        contrast_low=contrast_low,
        contrast_high=contrast_high,
    )


def build_inpainted_image(
    proxy_image: Image.Image,
    mask_image: Image.Image,
    *,
    mask_threshold: float = 127.0,
    gray_spread: float = DEFAULT_GRAY_SPREAD,
    white_brightness: float = DEFAULT_WHITE_BRIGHTNESS,
) -> Image.Image:
    # gray_spread and white_brightness are kept as explicit knobs for compatibility
    # with the exploratory dev-sample workflow, even though the mask image is now the
    # main driver for the hole filling stage.
    del gray_spread, white_brightness
    proxy_array = np.asarray(proxy_image.convert("L"), dtype=np.float32)
    mask_array = np.asarray(mask_image.convert("L"), dtype=np.float32) >= float(mask_threshold)
    inpainted = fill_masked_regions(proxy_array, mask_array)
    return Image.fromarray(np.rint(inpainted).astype(np.uint8), mode="L")


def build_download_result_from_cache(
    *,
    image_time_utc: dt.datetime,
    kind: GeoSatelliteKind,
) -> GeoSatelliteDownloadResult | None:
    return read_raw_cache(image_time_utc=image_time_utc, kind=kind)


def _resolve_latest_available_image_time(
    *,
    kind: GeoSatelliteKind,
    area: str = "europe",
    size: str = "normal",
    max_age_seconds: float = AVAILABLE_CACHE_MAX_AGE_SECONDS,
) -> dt.datetime:
    now_utc = dt.datetime.now(dt.timezone.utc)
    cached = read_latest_available_cache(area=area, kind=kind, size=size)
    if cached is not None and latest_available_cache_is_recent(
        cached,
        now_utc=now_utc,
        max_age_seconds=max_age_seconds,
    ):
        cached_text = cached.get("available_time_utc") or cached.get("fetched_at_utc")
        if isinstance(cached_text, str) and cached_text:
            try:
                return dt.datetime.fromisoformat(cached_text).astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
            except ValueError:
                pass
    try:
        available_time_utc, entry_count, source_url = fetch_latest_available_image_time(
            kind=kind,
            area=area,
            size=size,
        )
    except Exception:
        if cached is not None:
            cached_text = cached.get("available_time_utc") or cached.get("fetched_at_utc")
            if isinstance(cached_text, str) and cached_text:
                logger.warning(
                    "Geo-sat available list fetch failed; using cached latest slot: %s %s",
                    kind,
                    cached_text,
                )
                return dt.datetime.fromisoformat(cached_text).astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
        raise
    write_latest_available_cache(
        area=area,
        kind=kind,
        size=size,
        available_time_utc=available_time_utc,
        fetched_at_utc=now_utc,
        source_url=source_url,
        entry_count=entry_count,
    )
    return available_time_utc


def run_geo_satellite_pipeline(
    *,
    observer_lat: float,
    observer_lon: float,
    alt: float,
    az: float,
    fov_deg: float,
    kind: GeoSatelliteKind = "infrared",
    cloud_height_km: float = DEFAULT_CLOUD_HEIGHT_KM,
    radius_px: int = DEFAULT_RENDER_RADIUS_PX,
    raw_png: bytes | None = None,
    download_time_utc: dt.datetime | None = None,
    mask_path: Path | None = Path(GEOSATELLITE_GRAY_COMMON_MASK_FILE),
    grid_npz: Path = DEFAULT_GRID_NPZ,
) -> GeoSatellitePipelineResult:
    """Run the experimental Geo-satellite workflow end to end."""
    purge_intermediate_cache()
    if raw_png is None:
        if download_time_utc is None:
            download_time_utc = _resolve_latest_available_image_time(kind=kind)
        cache_slot = download_time_utc.astimezone(dt.timezone.utc).strftime("%Y%m%dT%H%MZ")
        cached = build_download_result_from_cache(image_time_utc=download_time_utc, kind=kind)
        if cached is None:
            logger.info("Geo-sat raw cache miss: %s %s", kind, cache_slot)
            download = download_latest_image(kind=kind, requested_time_utc=download_time_utc)
            download_path, metadata_path = write_raw_cache(download)
            download = replace(download, cache_path=download_path, metadata_path=metadata_path)
        else:
            logger.info("Geo-sat raw cache hit: %s %s", kind, cache_slot)
            download = cached
    else:
        if download_time_utc is None:
            download_time_utc = dt.datetime.now(dt.timezone.utc)
        download = GeoSatelliteDownloadResult(
            fetched_at_utc=download_time_utc,
            captured_at_utc=download_time_utc,
            kind=kind,
            source_url="memory",
            png_bytes=raw_png,
            content_type="image/png",
        )

    raw_digest = compute_digest(download.png_bytes)
    mask_image = _load_mask_image(mask_path)
    mask_digest = _mask_digest(mask_image)
    source_image = _load_png_image(download.png_bytes)
    proxy_image = build_proxy_image(source_image)
    inpainted_image = build_inpainted_image(proxy_image, mask_image)
    proxy_array = np.asarray(proxy_image, dtype=np.uint8)
    inpainted_array = np.asarray(inpainted_image, dtype=np.uint8)
    manifest = {
        "raw_digest": raw_digest,
        "mask_digest": mask_digest,
        "kind": kind,
        "mask_path": str(mask_path) if mask_path is not None else None,
    }

    disc_gray = project_gray_image_to_disc(
        inpainted_array,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        alt=alt,
        az=az,
        fov_deg=fov_deg,
        grid_npz=grid_npz,
        cloud_height_km=cloud_height_km,
        radius_px=radius_px,
    )

    intermediate = GeoSatelliteIntermediateResult(
        raw_digest=raw_digest,
        proxy_gray=proxy_array,
        inpainted_gray=inpainted_array,
        manifest=manifest,
        proxy_path=None,
        inpainted_path=None,
        mask_path=mask_path,
    )
    return GeoSatellitePipelineResult(download=download, intermediate=intermediate, disc_gray=disc_gray)
