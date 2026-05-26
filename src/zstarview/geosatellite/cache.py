from __future__ import annotations

import datetime as dt
import hashlib
import json
from pathlib import Path

from PIL import Image

from ..paths import GEOSATELLITE_CACHE_ROOT_DIR
from .types import GeoSatelliteDownloadResult

RAW_CACHE_DIRNAME = "raw"
INTERMEDIATE_CACHE_DIRNAME = "intermediate"
AVAILABLE_CACHE_DIRNAME = "available"
AVAILABLE_CACHE_MAX_AGE_SECONDS = 300.0


def geo_satellite_cache_root() -> Path:
    return Path(GEOSATELLITE_CACHE_ROOT_DIR)


def geo_satellite_raw_cache_dir() -> Path:
    return geo_satellite_cache_root() / RAW_CACHE_DIRNAME


def geo_satellite_intermediate_cache_dir() -> Path:
    return geo_satellite_cache_root() / INTERMEDIATE_CACHE_DIRNAME


def geo_satellite_available_cache_dir() -> Path:
    return geo_satellite_cache_root() / AVAILABLE_CACHE_DIRNAME


def _format_time_slot(timeslot_utc: dt.datetime) -> str:
    slot = timeslot_utc.astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
    return slot.strftime("%Y%m%dT%H%MZ")


def raw_cache_stem(*, image_time_utc: dt.datetime, kind: str) -> str:
    return f"{_format_time_slot(image_time_utc)}_{kind}"


def raw_png_path(*, image_time_utc: dt.datetime, kind: str) -> Path:
    return geo_satellite_raw_cache_dir() / f"{raw_cache_stem(image_time_utc=image_time_utc, kind=kind)}.png"


def raw_metadata_path(*, image_time_utc: dt.datetime, kind: str) -> Path:
    return geo_satellite_raw_cache_dir() / f"{raw_cache_stem(image_time_utc=image_time_utc, kind=kind)}.json"


def write_raw_cache(result: GeoSatelliteDownloadResult) -> tuple[Path, Path]:
    image_time_utc = result.captured_at_utc or result.fetched_at_utc
    png_path = raw_png_path(image_time_utc=image_time_utc, kind=result.kind)
    metadata_path = raw_metadata_path(image_time_utc=image_time_utc, kind=result.kind)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.write_bytes(result.png_bytes)
    metadata_path.write_text(
        json.dumps(
            {
                "captured_at_utc": image_time_utc.isoformat(),
                "fetched_at_utc": result.fetched_at_utc.isoformat(),
                "kind": result.kind,
                "source_url": result.source_url,
                "content_type": result.content_type,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return png_path, metadata_path


def read_raw_cache(*, image_time_utc: dt.datetime, kind: str) -> GeoSatelliteDownloadResult | None:
    png_path = raw_png_path(image_time_utc=image_time_utc, kind=kind)
    metadata_path = raw_metadata_path(image_time_utc=image_time_utc, kind=kind)
    if not png_path.exists() or not metadata_path.exists():
        return None
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    png_bytes = png_path.read_bytes()
    captured_at_text = metadata.get("captured_at_utc") or metadata.get("fetched_at_utc")
    return GeoSatelliteDownloadResult(
        fetched_at_utc=dt.datetime.fromisoformat(str(metadata["fetched_at_utc"])),
        captured_at_utc=dt.datetime.fromisoformat(str(captured_at_text)) if captured_at_text else None,
        kind=str(metadata["kind"]),
        source_url=str(metadata["source_url"]),
        png_bytes=png_bytes,
        content_type=metadata.get("content_type"),
        cache_path=png_path,
        metadata_path=metadata_path,
    )


def compute_digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def intermediate_cache_dir(raw_digest: str, *, mask_digest: str) -> Path:
    return geo_satellite_intermediate_cache_dir() / raw_digest / mask_digest


def intermediate_proxy_path(raw_digest: str, *, mask_digest: str) -> Path:
    return intermediate_cache_dir(raw_digest, mask_digest=mask_digest) / "proxy.png"


def intermediate_inpainted_path(raw_digest: str, *, mask_digest: str) -> Path:
    return intermediate_cache_dir(raw_digest, mask_digest=mask_digest) / "inpainted.png"


def intermediate_manifest_path(raw_digest: str, *, mask_digest: str) -> Path:
    return intermediate_cache_dir(raw_digest, mask_digest=mask_digest) / "manifest.json"


def latest_available_path(*, area: str, kind: str, size: str = "normal") -> Path:
    filename = f"{area}_{kind}_{size}.json"
    return geo_satellite_available_cache_dir() / filename


def write_latest_available_cache(
    *,
    area: str,
    kind: str,
    size: str,
    available_time_utc: dt.datetime,
    fetched_at_utc: dt.datetime,
    source_url: str,
    entry_count: int,
) -> Path:
    path = latest_available_path(area=area, kind=kind, size=size)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "area": area,
                "kind": kind,
                "size": size,
                "available_time_utc": available_time_utc.isoformat(),
                "fetched_at_utc": fetched_at_utc.isoformat(),
                "source_url": source_url,
                "entry_count": int(entry_count),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return path


def read_latest_available_cache(
    *,
    area: str,
    kind: str,
    size: str = "normal",
) -> dict[str, object] | None:
    path = latest_available_path(area=area, kind=kind, size=size)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if str(payload.get("area")) != area or str(payload.get("kind")) != kind or str(payload.get("size")) != size:
        return None
    return payload


def latest_available_cache_is_recent(
    payload: dict[str, object],
    *,
    now_utc: dt.datetime,
    max_age_seconds: float = AVAILABLE_CACHE_MAX_AGE_SECONDS,
) -> bool:
    fetched_at_text = payload.get("fetched_at_utc")
    if not fetched_at_text:
        return False
    try:
        fetched_at = dt.datetime.fromisoformat(str(fetched_at_text))
    except Exception:
        return False
    if fetched_at.tzinfo is None:
        fetched_at = fetched_at.replace(tzinfo=dt.timezone.utc)
    age_seconds = (now_utc.astimezone(dt.timezone.utc) - fetched_at.astimezone(dt.timezone.utc)).total_seconds()
    return age_seconds >= 0.0 and age_seconds <= float(max_age_seconds)


def write_intermediate_cache(
    *,
    raw_digest: str,
    mask_digest: str,
    proxy_image: Image.Image,
    inpainted_image: Image.Image,
    manifest: dict[str, object],
) -> tuple[Path, Path, Path]:
    proxy_path = intermediate_proxy_path(raw_digest, mask_digest=mask_digest)
    inpainted_path = intermediate_inpainted_path(raw_digest, mask_digest=mask_digest)
    manifest_path = intermediate_manifest_path(raw_digest, mask_digest=mask_digest)
    proxy_path.parent.mkdir(parents=True, exist_ok=True)
    proxy_image.save(proxy_path)
    inpainted_image.save(inpainted_path)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    return proxy_path, inpainted_path, manifest_path


def read_intermediate_cache(
    raw_digest: str,
    *,
    mask_digest: str,
) -> tuple[Image.Image, Image.Image, dict[str, object]] | None:
    proxy_path = intermediate_proxy_path(raw_digest, mask_digest=mask_digest)
    inpainted_path = intermediate_inpainted_path(raw_digest, mask_digest=mask_digest)
    manifest_path = intermediate_manifest_path(raw_digest, mask_digest=mask_digest)
    if not proxy_path.exists() or not inpainted_path.exists() or not manifest_path.exists():
        return None
    with Image.open(proxy_path) as proxy_image:
        proxy = proxy_image.convert("L")
    with Image.open(inpainted_path) as inpainted_image:
        inpainted = inpainted_image.convert("L")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    return proxy, inpainted, manifest
