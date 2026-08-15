"""SDO/HMI Continuum image retrieval for the Sun hover overlay."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np
from PySide6.QtGui import QImage

from .paths import CACHE_PATH
from .render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba
from .user_agent import build_user_agent

_CLOSEST_IMAGE_URL = "https://api.helioviewer.org/v2/getClosestImage/"
_SCREENSHOT_URL = "https://api.helioviewer.org/v2/takeScreenshot/"
_SOURCE_ID = 18
_IMAGE_SIZE = 1024
_IMAGE_SCALE_ARCSEC = 2.4
_HMI_LAYERS = "[SDO,HMI,HMI,continuum,1,100]"
_CACHE_SUBDIR = "solar-hover"
_CACHE_VERSION = "v3"
_DEFAULT_TIMEOUT_S = 30.0
_BLACK_ALPHA_THRESHOLD = 3
_LIMB_FEATHER_WIDTH_PX = 4.0


@dataclass(frozen=True, slots=True)
class SolarHoverImage:
    """A decoded HMI Continuum image and its cache metadata."""

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
    """Build the HMI Continuum metadata URL for a requested UTC datetime."""
    value = value.astimezone(timezone.utc)
    query = urllib.parse.urlencode(
        {
            "date": value.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "sourceId": _SOURCE_ID,
        }
    )
    return f"{_CLOSEST_IMAGE_URL}?{query}"


def screenshot_url(value: datetime) -> str:
    """Build a centered, north-up HMI Continuum screenshot URL."""
    value = value.astimezone(timezone.utc)
    query = urllib.parse.urlencode(
        {
            "date": value.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "imageScale": _IMAGE_SCALE_ARCSEC,
            "layers": _HMI_LAYERS,
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
        raise ValueError("HelioViewer HMI metadata contains invalid image geometry")
    return {
        "image_id": image_id,
        "time_utc": observed.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_radius_px": radius * scale / _IMAGE_SCALE_ARCSEC,
    }


def apply_solar_black_to_alpha(
    image: QImage,
    source_radius_px: float,
    *,
    feather_width_px: float = _LIMB_FEATHER_WIDTH_PX,
    black_threshold: int = _BLACK_ALPHA_THRESHOLD,
) -> QImage:
    """Make the black background transparent outside the solar disc.

    The disc remains opaque. Across the first ``feather_width_px`` beyond the
    limb, alpha transitions smoothly to a max-channel black-to-alpha mask.
    """
    radius = float(source_radius_px)
    feather = float(feather_width_px)
    threshold = int(black_threshold)
    if radius <= 0.0 or feather < 0.0 or not 0 <= threshold <= 255:
        raise ValueError("Invalid solar black-to-alpha parameters")

    rgba = qimage_to_np_rgba(image)
    height, width, _ = rgba.shape
    yy, xx = np.ogrid[:height, :width]
    center_x = (width - 1) / 2.0
    center_y = (height - 1) / 2.0
    distance = np.hypot(xx - center_x, yy - center_y)
    if feather <= 0.0:
        outside_weight = (distance > radius).astype(np.float32)
    else:
        outside_weight = np.clip((distance - radius) / feather, 0.0, 1.0)
        outside_weight = outside_weight * outside_weight * (
            3.0 - 2.0 * outside_weight
        )

    rgb = rgba[..., :3].astype(np.float32)
    max_channel = np.max(rgb, axis=2)
    visible = max_channel > float(threshold)
    outside_alpha = np.where(visible, max_channel, 0.0)
    outside_premultiplied = np.where(visible[..., None], rgb, 0.0)

    weight = outside_weight[..., None]
    alpha = 255.0 * (1.0 - outside_weight) + outside_alpha * outside_weight
    premultiplied = rgb * (1.0 - weight) + outside_premultiplied * weight
    straight_rgb = np.zeros_like(rgb)
    np.divide(
        premultiplied * 255.0,
        alpha[..., None],
        out=straight_rgb,
        where=alpha[..., None] > 0.0,
    )
    output = np.empty_like(rgba)
    output[..., :3] = np.clip(np.rint(straight_rgb), 0.0, 255.0).astype(np.uint8)
    output[..., 3] = np.clip(np.rint(alpha), 0.0, 255.0).astype(np.uint8)
    return np_rgba_to_qimage(output).convertToFormat(
        QImage.Format.Format_ARGB32_Premultiplied
    )


def _decode_image(payload: bytes, source_radius_px: float) -> QImage:
    image = QImage()
    if not image.loadFromData(payload):
        raise ValueError("HelioViewer HMI solar image could not be decoded")
    if image.width() != _IMAGE_SIZE or image.height() != _IMAGE_SIZE:
        raise ValueError("HelioViewer HMI solar image is not 1024x1024")
    return apply_solar_black_to_alpha(image, source_radius_px)


def _load_cached_entry(metadata_path: Path, *, stale: bool = False) -> SolarHoverImage:
    data = json.loads(metadata_path.read_text(encoding="ascii"))
    image_id = int(data["image_id"])
    observed = datetime.strptime(
        str(data["time_utc"]), "%Y-%m-%dT%H:%M:%SZ"
    ).replace(tzinfo=timezone.utc)
    radius = float(data["source_radius_px"])
    if image_id <= 0 or radius <= 0.0:
        raise ValueError("Cached solar metadata contains invalid geometry")
    image = _decode_image(
        (metadata_path.parent / f"{image_id}.png").read_bytes(),
        radius,
    )
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
    """Fetch and decode the HMI Continuum image nearest ``target_time``."""
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
        image = _decode_image(image_payload, float(str(metadata["source_radius_px"])))
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
    "apply_solar_black_to_alpha",
    "closest_image_url",
    "fetch_solar_hover_image",
    "normalize_solar_hover_time",
    "screenshot_url",
    "solar_hover_expiry",
]
