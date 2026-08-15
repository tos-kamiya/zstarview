"""NASA SVS Dial-A-Moon image retrieval for the Moon hover overlay."""

from __future__ import annotations

import json
import logging
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable

from PySide6.QtGui import QImage

from .paths import CACHE_PATH

logger = logging.getLogger(__name__)

_DIAL_A_MOON_URL = "https://svs.gsfc.nasa.gov/api/dialamoon/{}"
_CACHE_SUBDIR = "moon-hover"
_CACHE_VERSION = "v1"
_IMAGE_SIZE = 730
_DEFAULT_TIMEOUT_S = 15.0


@dataclass(frozen=True, slots=True)
class MoonHoverImage:
    """A decoded NASA Moon image and the timestamp it represents."""

    image: QImage
    time_utc: datetime
    phase_percent: float | None = None
    age_days: float | None = None
    diameter_arcsec: float | None = None
    posangle_deg: float | None = None


def normalize_dialamoon_time(value: datetime) -> datetime:
    """Round a UTC datetime to the nearest whole hour."""
    value = value.astimezone(timezone.utc).replace(second=0, microsecond=0)
    if value.minute >= 30:
        value += timedelta(hours=1)
    return value.replace(minute=0)


def dialamoon_url(value: datetime) -> str:
    """Build the public Dial-A-Moon API URL for a UTC datetime."""
    normalized = normalize_dialamoon_time(value)
    return _DIAL_A_MOON_URL.format(normalized.strftime("%Y-%m-%dT%H:%M"))


def _cache_root(cache_root: str | Path | None) -> Path:
    return Path(cache_root or CACHE_PATH) / _CACHE_SUBDIR / _CACHE_VERSION


def _cache_path(root: Path, key: str, suffix: str) -> Path:
    return root / f"{key}.{suffix}"


def _request_bytes(
    url: str,
    *,
    timeout_s: float,
    opener: Callable[..., Any] = urllib.request.urlopen,
) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "zstarview/1.0 Moon hover"},
    )
    with opener(request, timeout=timeout_s) as response:
        return response.read()


def _parse_response(payload: bytes) -> tuple[str, datetime, dict[str, Any]]:
    data = json.loads(payload.decode("utf-8"))
    image_data = data.get("image")
    if not isinstance(image_data, dict) or not isinstance(image_data.get("url"), str):
        raise ValueError("Dial-A-Moon response has no image URL")
    raw_time = data.get("time")
    if not isinstance(raw_time, str):
        raise ValueError("Dial-A-Moon response has no image time")
    image_time = datetime.strptime(raw_time, "%Y-%m-%dT%H:%M").replace(
        tzinfo=timezone.utc
    )
    return str(image_data["url"]), image_time, data


def _decode_image(image_bytes: bytes) -> QImage:
    image = QImage()
    if not image.loadFromData(image_bytes):
        raise ValueError("NASA Moon image could not be decoded")
    if image.width() != _IMAGE_SIZE or image.height() != _IMAGE_SIZE:
        raise ValueError("NASA Moon image is not 730x730")
    return image.convertToFormat(QImage.Format.Format_ARGB32_Premultiplied)


def fetch_moon_hover_image(
    target_time: datetime,
    *,
    cache_root: str | Path | None = None,
    timeout_s: float = _DEFAULT_TIMEOUT_S,
    opener: Callable[..., Any] = urllib.request.urlopen,
) -> MoonHoverImage:
    """Fetch and decode the NASA image for ``target_time``.

    The function is intentionally synchronous and is expected to run in the
    shared GUI worker pool. It only writes complete response/image files after
    successful parsing and decoding.
    """
    normalized = normalize_dialamoon_time(target_time)
    request_key = normalized.strftime("%Y%m%dT%H00Z")
    root = _cache_root(cache_root)
    root.mkdir(parents=True, exist_ok=True)
    metadata_path = _cache_path(root, request_key, "json")

    try:
        metadata_bytes = metadata_path.read_bytes()
        image_url, image_time, data = _parse_response(metadata_bytes)
    except (OSError, ValueError, json.JSONDecodeError):
        response_bytes = _request_bytes(
            dialamoon_url(normalized), timeout_s=timeout_s, opener=opener
        )
        image_url, image_time, data = _parse_response(response_bytes)
        metadata_path.write_bytes(response_bytes)

    image_key = Path(image_url).name.rsplit(".", 1)[0]
    image_path = _cache_path(root, image_key, "jpg")
    try:
        image_bytes = image_path.read_bytes()
    except OSError:
        image_bytes = _request_bytes(image_url, timeout_s=timeout_s, opener=opener)
        image = _decode_image(image_bytes)
        image_path.write_bytes(image_bytes)
    else:
        image = _decode_image(image_bytes)

    return MoonHoverImage(
        image=image,
        time_utc=image_time,
        phase_percent=_optional_float(data.get("phase")),
        age_days=_optional_float(data.get("age")),
        diameter_arcsec=_optional_float(data.get("diameter")),
        posangle_deg=_optional_float(data.get("posangle")),
    )


def _optional_float(value: object) -> float | None:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


__all__ = [
    "MoonHoverImage",
    "dialamoon_url",
    "fetch_moon_hover_image",
    "normalize_dialamoon_time",
]
