"""Encode and write export-image output, including SIXEL and metadata."""

from __future__ import annotations

import json
import time

try:
    import termios
except ImportError:  # pragma: no cover - non-Unix fallback
    termios = None  # type: ignore[assignment]

from pathlib import Path

from PySide6.QtCore import QBuffer, QByteArray, QIODevice
from PySide6.QtGui import QImage, QImageWriter

from ..location_resolver import ResolvedLocation
from ..paths import OBSERVER_MAX_ALT_DEG, OBSERVER_MIN_ALT_DEG
from ..render import background as render_background
from ..render.pipeline import RenderStyle
from ..search.models import SearchJumpTarget
from ..types import CelestialData, ViewerData
from .export_image_support import (
    EXPORT_IMAGE_METADATA_SCHEMA,
    EXPORT_IMAGE_METADATA_TEXT_KEY,
    host,
    logger,
)


def _resolved_location_source(location: ResolvedLocation) -> str:
    persistence_value = location.persistence_value
    if isinstance(persistence_value, dict):
        resolver = persistence_value.get("resolver")
        if isinstance(resolver, str) and resolver.strip():
            return resolver.strip()
    if location.kind:
        return str(location.kind)
    return "unknown"

def _resolved_location_metadata(location: ResolvedLocation) -> dict[str, object]:
    return {
        "display_name": location.display_name,
        "lat_deg": float(location.lat),
        "lon_deg": float(location.lon),
        "timezone_name": location.tz,
        "source": _resolved_location_source(location),
        "kind": location.kind,
        "observer_height_m": float(location.observer_height_m),
        "ground_elevation_m": float(location.ground_elevation_m),
        "location_height_label": location.location_height_label,
        "location_height_m": float(location.location_height_m),
        "height_add_m": float(location.height_add_m),
        "cc": location.cc,
    }

def _search_target_metadata(target: SearchJumpTarget) -> dict[str, object]:
    payload: dict[str, object] = {
        "label": target.label,
        "kind": target.kind,
        "object_key": target.object_key,
        "command": target.command,
    }
    if target.subtitle:
        payload["subtitle"] = target.subtitle
    if target.alt_deg is not None:
        payload["alt_deg"] = float(target.alt_deg)
    if target.az_deg is not None:
        payload["az_deg"] = float(target.az_deg)
    if target.latitude_deg is not None:
        payload["latitude_deg"] = float(target.latitude_deg)
    if target.longitude_deg is not None:
        payload["longitude_deg"] = float(target.longitude_deg)
    return payload

def _build_export_image_metadata_payload(
    *,
    app_version: str,
    viewer_data: ViewerData,
    celestial_data: CelestialData,
    style: RenderStyle,
    place_query: str | None,
    place_location: ResolvedLocation | None,
    search_overlay_target: SearchJumpTarget | None,
    cloud_coverage_ratio: float | None,
    urban_outline_source: str | None = None,
    urban_outline_count: int | None = None,
) -> dict[str, object]:
    hud_lines = [
        " ".join(str(line).split())
        for line in render_background.format_overlay_info_lines(
            celestial_data,
            viewer_data,
            float(style.vmag_limit),
            include_vmag_limit=True,
        )
    ]
    payload: dict[str, object] = {
        "schema": EXPORT_IMAGE_METADATA_SCHEMA,
        "version": str(app_version),
        "hud": {
            "lines": hud_lines,
            "view": {
                "city_name": viewer_data.city_name,
                "lat_deg": float(viewer_data.lat_deg),
                "lon_deg": float(viewer_data.lon_deg),
                "view_center_alt_deg": float(viewer_data.view_center[0]),
                "view_center_az_deg": float(viewer_data.view_center[1]),
                "vmag_limit": float(style.vmag_limit),
            },
        },
        "place": {},
        "extra": {},
    }
    if place_query is not None and place_location is not None:
        payload["place"] = {
            "query": place_query,
            "resolved": _resolved_location_metadata(place_location),
        }
    if search_overlay_target is not None:
        payload["extra"]["search_target"] = _search_target_metadata(
            search_overlay_target
        )
    if cloud_coverage_ratio is not None:
        payload["extra"]["cloud_coverage_ratio"] = float(cloud_coverage_ratio)
    if urban_outline_source is not None or urban_outline_count is not None:
        payload["extra"]["urban_outline"] = {
            "source": urban_outline_source,
            "outline_count": urban_outline_count,
        }
    return payload

def _require_img2sixel_binary() -> str:
    executable = host().shutil.which("img2sixel")
    if executable:
        return executable
    logger.error("--sixel was requested, but 'img2sixel' was not found in PATH.")
    raise SystemExit(1)

def _response_indicates_sixel_support(response: bytes) -> bool:
    if not response.startswith(b"\x1b[") or not response.endswith(b"c"):
        return False
    for token in response[2:-1].split(b";"):
        token = token.lstrip(b"?")
        if token == b"4":
            return True
    return False

def _require_sixel_terminal_support(timeout_seconds: float = 0.25) -> None:
    if host().os.environ.get("TERM", "").startswith("yaft"):
        return
    if host().termios is None:
        logger.error(
            "--sixel was requested, but terminal control is unavailable on this platform."
        )
        raise SystemExit(1)
    tty_fd: int | None = None
    old_attrs = None
    response = bytearray()
    try:
        tty_fd = host().os.open("/dev/tty", host().os.O_RDWR | host().os.O_NOCTTY)
        old_attrs = host().termios.tcgetattr(tty_fd)
        new_attrs = host().termios.tcgetattr(tty_fd)
        new_attrs[3] &= ~(host().termios.ECHO | host().termios.ICANON)
        new_attrs[6][host().termios.VMIN] = 0
        new_attrs[6][host().termios.VTIME] = 0
        host().termios.tcsetattr(tty_fd, host().termios.TCSANOW, new_attrs)
        host().os.write(tty_fd, b"\x1b[c")
        deadline = time.monotonic() + timeout_seconds
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            ready, _, _ = host().select.select([tty_fd], [], [], remaining)
            if not ready:
                break
            chunk = host().os.read(tty_fd, 1024)
            if not chunk:
                break
            response.extend(chunk)
            if b"c" in chunk:
                break
    except OSError:
        logger.error(
            "--sixel was requested, but the terminal does not report SIXEL graphics support."
        )
        raise SystemExit(1)
    finally:
        if tty_fd is not None and old_attrs is not None:
            try:
                host().termios.tcsetattr(tty_fd, host().termios.TCSANOW, old_attrs)
            except OSError:
                pass
        if tty_fd is not None:
            try:
                host().os.close(tty_fd)
            except OSError:
                pass

    if not host()._response_indicates_sixel_support(bytes(response)):
        logger.error(
            "--sixel was requested, but the terminal does not report SIXEL graphics support."
        )
        raise SystemExit(1)

def _encode_image_as_png_bytes(
    image: QImage,
    *,
    metadata_payload: dict[str, object] | None = None,
) -> bytes:
    ba = QByteArray()
    buf = QBuffer(ba)
    if not buf.open(QIODevice.OpenModeFlag.WriteOnly):
        logger.error("Failed to open in-memory buffer for PNG encoding.")
        raise SystemExit(1)
    try:
        writer = QImageWriter(buf, b"PNG")
        if metadata_payload is not None:
            writer.setText(
                EXPORT_IMAGE_METADATA_TEXT_KEY,
                json.dumps(
                    metadata_payload,
                    ensure_ascii=False,
                    separators=(",", ":"),
                    sort_keys=True,
                ),
            )
        if not writer.write(image):
            logger.error("Failed to encode image as PNG for SIXEL output.")
            logger.error("%s", writer.errorString())
            raise SystemExit(1)
    finally:
        buf.close()
    return bytes(ba.data())

def _write_sixel_to_stdout(
    image: QImage,
    *,
    img2sixel_bin: str,
    metadata_payload: dict[str, object] | None = None,
) -> bool:
    png_bytes = host()._encode_image_as_png_bytes(
        image, metadata_payload=metadata_payload
    )
    try:
        proc = host().subprocess.run(
            [img2sixel_bin, "-"],
            input=png_bytes,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        logger.warning("SIXEL output failed: %s", exc)
        return False
    if proc.returncode != 0:
        stderr_text = proc.stderr.decode("utf-8", errors="replace").strip()
        if stderr_text:
            logger.warning("SIXEL output failed: %s", stderr_text)
        else:
            logger.warning("SIXEL output failed with exit status %s.", proc.returncode)
        return False
    host().sys.stdout.buffer.write(proc.stdout)
    host().sys.stdout.buffer.flush()
    return True

def _write_png_to_stdout(
    image: QImage,
    *,
    metadata_payload: dict[str, object] | None = None,
) -> bool:
    png_bytes = host()._encode_image_as_png_bytes(
        image, metadata_payload=metadata_payload
    )
    try:
        host().sys.stdout.buffer.write(png_bytes)
        host().sys.stdout.buffer.flush()
    except OSError as exc:
        logger.error("Failed to write PNG image to stdout: %s", exc)
        return False
    return True

def _write_png_to_path(
    image: QImage,
    output_path: Path,
    *,
    metadata_payload: dict[str, object] | None = None,
) -> bool:
    png_bytes = host()._encode_image_as_png_bytes(
        image, metadata_payload=metadata_payload
    )
    try:
        output_path.write_bytes(png_bytes)
    except OSError as exc:
        logger.error("Failed to save image: %s", exc)
        return False
    return True

def _write_export_overlay_summary_to_stderr(
    *,
    viewer_data: ViewerData,
    celestial_data: CelestialData,
    vmag_limit: float,
    cloud_coverage_ratio: float | None = None,
    search_overlay_target: SearchJumpTarget | None = None,
) -> None:
    lines = render_background.format_overlay_info_lines(
        celestial_data,
        viewer_data,
        vmag_limit,
        include_vmag_limit=True,
    )
    if cloud_coverage_ratio is None:
        lines.append("Cloud coverage n/a")
    else:
        lines.append(f"Cloud coverage {float(cloud_coverage_ratio) * 100.0:.1f}%")
    if search_overlay_target is not None:
        lines.append(host()._format_search_target_line(search_overlay_target))
    host().sys.stderr.write("\n".join(lines) + "\n")
    host().sys.stderr.flush()

def _format_search_target_line(target: SearchJumpTarget) -> str:
    parts = [
        f"Search target label={target.label}",
        f"id={target.object_key or target.command or target.label}",
        f"kind={target.kind}",
    ]
    if target.alt_deg is not None:
        parts.append(f"alt={float(target.alt_deg):.1f} deg")
    if target.az_deg is not None:
        parts.append(f"az={float(target.az_deg):.1f} deg")
    return " | ".join(parts)

def _format_search_candidate_line(target: SearchJumpTarget) -> str:
    parts = [
        f"label={target.label}",
        f"id={target.object_key or target.command or target.label}",
        f"kind={target.kind}",
    ]
    return " | ".join(parts)

def _format_search_failure_message(query: str, candidate_count: int) -> str:
    if candidate_count <= 0:
        return f"No search candidates found for {query!r}"
    return f"Search query {query!r} matched {candidate_count} candidates"

def _clamp_view_center_altitude(alt_deg: float) -> float:
    return max(
        float(OBSERVER_MIN_ALT_DEG), min(float(OBSERVER_MAX_ALT_DEG), float(alt_deg))
    )

def _search_view_center_for_target(
    *,
    base_view_center: tuple[float, float],
    target_alt_deg: float,
    target_az_deg: float,
    fixed_alt: bool,
    fixed_az: bool,
) -> tuple[float, float]:
    view_center_alt = (
        float(base_view_center[0])
        if fixed_alt
        else host()._clamp_view_center_altitude(target_alt_deg)
    )
    view_center_az = (
        float(base_view_center[1]) if fixed_az else float(target_az_deg) % 360.0
    )
    return view_center_alt, view_center_az
