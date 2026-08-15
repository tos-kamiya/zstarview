"""HelioViewer SDO/AIA 193 image retrieval for the Sun hover overlay."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable

from PySide6.QtGui import QImage

from .paths import CACHE_PATH
from .user_agent import build_user_agent

_CLOSEST_IMAGE_URL = "https://api.helioviewer.org/v2/getClosestImage/"
_SCREENSHOT_URL = "https://api.helioviewer.org/v2/takeScreenshot/"
_SOURCE_ID = 11
_IMAGE_SIZE = 1024
_IMAGE_SCALE_ARCSEC = 2.4
_CACHE_SUBDIR = "solar-hover"
_CACHE_VERSION = "v1"
_DEFAULT_TIMEOUT_S = 30.0


@dataclass(frozen=True, slots=True)
class SolarHoverImage:
    """A decoded AIA 193 image and the observation metadata it represents."""

    image: QImage
    time_utc: datetime
    source_radius_px: float
    image_id: int
    stale: bool = False


def normalize_solar_hover_time(value: datetime) -> datetime:
    """Floor a UTC datetime to the ten-minute cache bucket."""
    value = value.astimezone(timezone.utc).replace(second=0, microsecond=0)
    return value.replace(minute=(value.minute // 10) * 10)


def closest_image_url(value: datetime) -> str:
    """Build the metadata URL for the requested UTC datetime."""
    value = value.astimezone(timezone.utc)
    query = urllib.parse.urlencode(
        {
            "date": value.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "sourceId": _SOURCE_ID,
        }
    )
    return f"{_CLOSEST_IMAGE_URL}?{query}"


def screenshot_url(value: datetime) -> str:
    """Build a centered, north-up AIA 193 screenshot URL."""
    value = value.astimezone(timezone.utc)
    query = urllib.parse.urlencode(
        {
            "date": value.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "imageScale": _IMAGE_SCALE_ARCSEC,
            "layers": "[SDO,AIA,AIA,193,1,100]",
            "eventLabels": "false",
            "x0": 0,
            "y0": 0,
            "width": _IMAGE_SIZE,
            "height": _IMAGE_SIZE,
            "display": "true",
            "watermark": "false",
        }
    )
    return f"{_SCREENSHOT_URL}?{query}"


def _cache_root(cache_root: str | Path | None) -> Path:
    return Path(cache_root or CACHE_PATH) / _CACHE_SUBDIR / _CACHE_VERSION


def _request_bytes(
    url: str,
    *,
    timeout_s: float,
    opener: Callable[..., Any],
) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": build_user_agent("solar-hover")},
    )
    with opener(request, timeout=timeout_s) as response:
        return response.read()


def _parse_closest_image(payload: bytes) -> dict[str, object]:
    data = json.loads(payload.decode("utf-8"))
    image_id = int(data["id"])
    observed = datetime.strptime(str(data["date"]), "%Y-%m-%d %H:%M:%S").replace(
        tzinfo=timezone.utc
    )
    scale = float(data["scale"])
    radius = float(data["rsun"])
    if image_id <= 0 or scale <= 0.0 or radius <= 0.0:
        raise ValueError("HelioViewer metadata contains invalid image geometry")
    return {
        "image_id": image_id,
        "time_utc": observed.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_radius_px": radius * scale / _IMAGE_SCALE_ARCSEC,
    }


def _decode_image(payload: bytes) -> QImage:
    image = QImage()
    if not image.loadFromData(payload):
        raise ValueError("HelioViewer solar image could not be decoded")
    if image.width() != _IMAGE_SIZE or image.height() != _IMAGE_SIZE:
        raise ValueError("HelioViewer solar image is not 1024x1024")
    return image.convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)


def _load_cached_entry(metadata_path: Path, *, stale: bool = False) -> SolarHoverImage:
    data = json.loads(metadata_path.read_text(encoding="ascii"))
    image_id = int(data["image_id"])
    observed = datetime.strptime(
        str(data["time_utc"]), "%Y-%m-%dT%H:%M:%SZ"
    ).replace(tzinfo=timezone.utc)
    radius = float(data["source_radius_px"])
    if image_id <= 0 or radius <= 0.0:
        raise ValueError("Cached solar metadata contains invalid geometry")
    image = _decode_image((metadata_path.parent / f"{image_id}.png").read_bytes())
    return SolarHoverImage(image, observed, radius, image_id, stale=stale)


def _latest_cached_entry(root: Path) -> SolarHoverImage | None:
    candidates = sorted(root.glob("request-*.json"), reverse=True)
    for path in candidates:
        try:
            return _load_cached_entry(path, stale=True)
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            continue
    return None


def fetch_solar_hover_image(
    target_time: datetime,
    *,
    cache_root: str | Path | None = None,
    timeout_s: float = _DEFAULT_TIMEOUT_S,
    opener: Callable[..., Any] = urllib.request.urlopen,
) -> SolarHoverImage:
    """Fetch and decode the centered AIA 193 image nearest ``target_time``."""
    key = normalize_solar_hover_time(target_time)
    root = _cache_root(cache_root)
    root.mkdir(parents=True, exist_ok=True)
    metadata_path = root / f"request-{key.strftime('%Y%m%dT%H%MZ')}.json"
    try:
        return _load_cached_entry(metadata_path)
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        pass

    try:
        metadata_payload = _request_bytes(
            closest_image_url(target_time), timeout_s=timeout_s, opener=opener
        )
        metadata = _parse_closest_image(metadata_payload)
        observed = datetime.strptime(
            str(metadata["time_utc"]), "%Y-%m-%dT%H:%M:%SZ"
        ).replace(tzinfo=timezone.utc)
        image_payload = _request_bytes(
            screenshot_url(observed), timeout_s=timeout_s, opener=opener
        )
        image = _decode_image(image_payload)
    except Exception:
        stale = _latest_cached_entry(root)
        if stale is not None:
            return stale
        raise

    image_id = int(metadata["image_id"])
    image_path = root / f"{image_id}.png"
    image_path.write_bytes(image_payload)
    metadata_path.write_text(
        json.dumps(metadata, sort_keys=True, separators=(",", ":")),
        encoding="ascii",
    )
    return SolarHoverImage(
        image=image,
        time_utc=observed,
        source_radius_px=float(str(metadata["source_radius_px"])),
        image_id=image_id,
    )


def solar_hover_expiry(value: datetime) -> datetime:
    """Return the next ten-minute cache boundary for ``value``."""
    return normalize_solar_hover_time(value) + timedelta(minutes=10)


__all__ = [
    "SolarHoverImage",
    "closest_image_url",
    "fetch_solar_hover_image",
    "normalize_solar_hover_time",
    "screenshot_url",
    "solar_hover_expiry",
]
