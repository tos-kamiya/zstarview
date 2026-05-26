from __future__ import annotations

import datetime as dt
import io
import logging
import urllib.error
import urllib.request

from PIL import Image

from .types import GeoSatelliteDownloadResult, GeoSatelliteKind

logger = logging.getLogger(__name__)

DEFAULT_GEO_SATELLITE_URL = "https://api.met.no/weatherapi/geosatellite/1.4/?area=europe&type=infrared"
DEFAULT_USER_AGENT = "zstarview/1.0 github.com/tos-kamiya/zstarview"
DEFAULT_TIMEOUT_SECONDS = 30.0


def build_geo_satellite_url(*, kind: GeoSatelliteKind = "infrared") -> str:
    return f"https://api.met.no/weatherapi/geosatellite/1.4/?area=europe&type={kind}"


def _verify_png(data: bytes) -> None:
    with Image.open(io.BytesIO(data)) as image:
        image.load()
        if image.mode not in {"RGB", "RGBA", "L"}:
            image.convert("RGBA")


def download_latest_image(
    *,
    kind: GeoSatelliteKind = "infrared",
    user_agent: str = DEFAULT_USER_AGENT,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
) -> GeoSatelliteDownloadResult:
    """Download the latest Geo-satellite image and validate it as PNG-like data."""
    source_url = build_geo_satellite_url(kind=kind)
    request = urllib.request.Request(source_url, headers={"User-Agent": user_agent})
    fetched_at_utc = dt.datetime.now(dt.timezone.utc)
    try:
        with urllib.request.urlopen(request, timeout=float(timeout_seconds)) as response:
            status = getattr(response, "status", None)
            if status != 200:
                raise RuntimeError(f"Geo-satellite download failed with HTTP status {status}")
            content_type = response.headers.get_content_type() if hasattr(response, "headers") else None
            payload = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Geo-satellite download failed with HTTP status {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Geo-satellite download failed: {exc.reason}") from exc

    try:
        _verify_png(payload)
    except Exception as exc:
        raise RuntimeError("Geo-satellite download did not contain a readable image") from exc

    logger.info("Downloaded Geo-satellite %s frame from %s", kind, source_url)
    return GeoSatelliteDownloadResult(
        fetched_at_utc=fetched_at_utc,
        kind=kind,
        source_url=source_url,
        png_bytes=payload,
        content_type=content_type,
    )

