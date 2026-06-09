from __future__ import annotations

import datetime as dt
import io
import logging
from urllib.parse import urlencode, urlparse, parse_qs
import urllib.error
import urllib.request

from PIL import Image

from ..user_agent import build_user_agent
from .types import GeoSatelliteDownloadResult, GeoSatelliteKind

logger = logging.getLogger(__name__)

DEFAULT_USER_AGENT = build_user_agent("geosatellite")
DEFAULT_TIMEOUT_SECONDS = 30.0
DEFAULT_GEO_SATELLITE_AREA = "europe"
DEFAULT_GEO_SATELLITE_SIZE = "normal"


def _format_api_time(time_utc: dt.datetime) -> str:
    return time_utc.astimezone(dt.timezone.utc).replace(second=0, microsecond=0).strftime("%Y-%m-%dT%H:%M:00Z")


def build_geo_satellite_url(
    *,
    kind: GeoSatelliteKind = "infrared",
    area: str = DEFAULT_GEO_SATELLITE_AREA,
    time_utc: dt.datetime | None = None,
) -> str:
    params = {"area": area, "type": kind}
    if time_utc is not None:
        params["time"] = _format_api_time(time_utc)
    return f"https://api.met.no/weatherapi/geosatellite/1.4/?{urlencode(params)}"


def build_geo_satellite_available_url(
    *,
    kind: GeoSatelliteKind = "infrared",
    area: str = DEFAULT_GEO_SATELLITE_AREA,
    size: str = DEFAULT_GEO_SATELLITE_SIZE,
) -> str:
    params = {"area": area, "type": kind, "size": size}
    return f"https://api.met.no/weatherapi/geosatellite/1.4/available.json?{urlencode(params)}"


def _verify_png(data: bytes) -> None:
    with Image.open(io.BytesIO(data)) as image:
        image.load()
        if image.mode not in {"RGB", "RGBA", "L"}:
            image.convert("RGBA")


def download_latest_image(
    *,
    kind: GeoSatelliteKind = "infrared",
    requested_time_utc: dt.datetime | None = None,
    user_agent: str = DEFAULT_USER_AGENT,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
) -> GeoSatelliteDownloadResult:
    """Download the latest Geo-satellite image and validate it as PNG-like data."""
    source_url = build_geo_satellite_url(kind=kind, time_utc=requested_time_utc)
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
        captured_at_utc=requested_time_utc.astimezone(dt.timezone.utc).replace(second=0, microsecond=0) if requested_time_utc is not None else None,
        kind=kind,
        source_url=source_url,
        png_bytes=payload,
        content_type=content_type,
    )


def _parse_available_time_from_item(item: object) -> dt.datetime | None:
    if isinstance(item, dict):
        params = item.get("params")
        if isinstance(params, dict):
            time_text = params.get("time")
            if isinstance(time_text, str) and time_text:
                try:
                    return dt.datetime.fromisoformat(time_text.replace("Z", "+00:00"))
                except ValueError:
                    pass
        uri = item.get("uri")
        if isinstance(uri, str) and uri:
            try:
                query = parse_qs(urlparse(uri).query)
            except Exception:
                query = {}
            time_values = query.get("time")
            if time_values:
                try:
                    return dt.datetime.fromisoformat(time_values[0].replace("Z", "+00:00"))
                except ValueError:
                    return None
    return None


def fetch_latest_available_image_time(
    *,
    kind: GeoSatelliteKind = "infrared",
    area: str = DEFAULT_GEO_SATELLITE_AREA,
    size: str = DEFAULT_GEO_SATELLITE_SIZE,
    user_agent: str = DEFAULT_USER_AGENT,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
) -> tuple[dt.datetime, int, str]:
    source_url = build_geo_satellite_available_url(kind=kind, area=area, size=size)
    request = urllib.request.Request(source_url, headers={"User-Agent": user_agent})
    try:
        with urllib.request.urlopen(request, timeout=float(timeout_seconds)) as response:
            status = getattr(response, "status", None)
            if status != 200:
                raise RuntimeError(f"Geo-satellite available list failed with HTTP status {status}")
            payload = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Geo-satellite available list failed with HTTP status {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Geo-satellite available list failed: {exc.reason}") from exc

    try:
        import json

        items = json.loads(payload.decode("utf-8"))
    except Exception as exc:
        raise RuntimeError("Geo-satellite available list did not contain readable JSON") from exc

    times: list[dt.datetime] = []
    if isinstance(items, list):
        for item in items:
            item_time = _parse_available_time_from_item(item)
            if item_time is not None:
                times.append(item_time)
    elif isinstance(items, dict):
        for value in items.values():
            if isinstance(value, list):
                for item in value:
                    item_time = _parse_available_time_from_item(item)
                    if item_time is not None:
                        times.append(item_time)
    if not times:
        raise RuntimeError("Geo-satellite available list did not contain any image times")
    latest_time = max(times).astimezone(dt.timezone.utc).replace(second=0, microsecond=0)
    logger.info("Geo-satellite available list latest %s frame is %s", kind, _format_api_time(latest_time))
    return latest_time, len(times), source_url
