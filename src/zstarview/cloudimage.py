# -*- coding: utf-8 -*-
"""Standalone cloud-image CLI backed by CloudDisc."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

from .clouddisc import CloudDisc, CloudDiscConfig
from .clouddisc.altaz_render import render_altaz_grid_circles
from .logging_utils import setup_root_logger
from .paths import CACHE_PATH, CLOUD_SHELLS_KM
from .render.qt_image import np_rgba_to_qimage

logger = logging.getLogger(__name__)

DEFAULT_OUTPUT_PATH = Path("latest_cloud.png")
DEFAULT_RENDER_RADIUS_PX = 255
DEFAULT_ALT_MIN_DEG = 0.0
MAX_FOV_DEG = 135.0


def parse_observer_spec(value: str) -> tuple[float, float]:
    """Parse an observer location in the form '@lat,lon'."""
    text = (value or "").strip()
    if not text.startswith("@"):
        raise argparse.ArgumentTypeError("Observer must be given as '@lat,lon'.")
    payload = text[1:].strip()
    parts = [part.strip() for part in payload.split(",")]
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise argparse.ArgumentTypeError("Observer must be given as '@lat,lon'.")
    try:
        lat = float(parts[0])
        lon = float(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Observer must be given as '@lat,lon'."
        ) from exc
    if not (-90.0 <= lat <= 90.0):
        raise argparse.ArgumentTypeError("Latitude must be between -90 and 90 degrees.")
    if not (-180.0 <= lon <= 180.0):
        raise argparse.ArgumentTypeError(
            "Longitude must be between -180 and 180 degrees."
        )
    return float(lat), float(lon)


def parse_fov_deg(value: str) -> float:
    """Parse the center-to-edge field of view in degrees."""
    try:
        fov = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid FOV: {value!r}") from exc
    if not (0.0 < fov <= MAX_FOV_DEG):
        raise argparse.ArgumentTypeError(
            f"FOV must be greater than 0 and at most {MAX_FOV_DEG:.0f} degrees."
        )
    return float(fov)


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for the cloud-image generator."""
    parser = argparse.ArgumentParser(
        prog="zstarview.cloudimage",
        description="Generate a standalone cloud image from GOES/Himawari data.",
    )
    parser.add_argument(
        "observer",
        type=parse_observer_spec,
        metavar="@LAT,LON",
        help="Observer position as '@lat,lon'.",
    )
    parser.add_argument(
        "--alt",
        type=float,
        required=True,
        help="View center altitude in degrees.",
    )
    parser.add_argument(
        "--az",
        type=float,
        required=True,
        help="View center azimuth in degrees.",
    )
    parser.add_argument(
        "--fov",
        type=parse_fov_deg,
        required=True,
        help="Field of view in degrees from the center to the image edge.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"Output PNG path (default: {DEFAULT_OUTPUT_PATH}).",
    )
    return parser


def _compose_on_black_background(cloud_rgba: np.ndarray) -> np.ndarray:
    """Flatten a cloud RGBA image onto an opaque black background."""
    arr = np.asarray(cloud_rgba)
    if arr.ndim != 3 or arr.shape[2] != 4:
        raise ValueError("cloud_rgba must have shape (H, W, 4)")
    if arr.dtype != np.uint8:
        raise ValueError("cloud_rgba must have dtype uint8")
    alpha = arr[..., 3:4].astype(np.uint16, copy=False)
    rgb = arr[..., :3].astype(np.uint16, copy=False)
    out = np.empty_like(arr)
    out[..., :3] = ((rgb * alpha + 127) // 255).astype(np.uint8, copy=False)
    out[..., 3] = 255
    return out


def render_cloud_image(
    *,
    observer_lat: float,
    observer_lon: float,
    alt: float,
    az: float,
    fov_deg: float,
    radius_px: int = DEFAULT_RENDER_RADIUS_PX,
) -> np.ndarray:
    """Fetch and render a cloud image as an opaque RGBA NumPy array."""
    logger.info(
        "Fetching cloud image for lat=%.4f lon=%.4f alt=%.2f az=%.2f fov=%.2f",
        observer_lat,
        observer_lon,
        alt,
        az,
        fov_deg,
    )
    clouddisc = CloudDisc(
        CloudDiscConfig(
            cache_dir=CACHE_PATH,
            sat_priority=("AUTO",),
            bt_warm_k=310.0,
            bt_cold_k=190.0,
            alt_min_deg=DEFAULT_ALT_MIN_DEG,
            search_back_minutes=120,
        )
    )
    source = clouddisc.fetch_source(lat=observer_lat, lon=observer_lon)
    source_satellite = getattr(source, "satellite", "unknown")
    source_time = getattr(source, "time_utc", None)
    source_time_text = (
        source_time.isoformat() if hasattr(source_time, "isoformat") else "unknown"
    )
    logger.info(
        "Rendering cloud disc from satellite=%s time=%s",
        source_satellite,
        source_time_text,
    )
    logger.info("Building alt/az cloud grid...")
    source.altaz_grid = clouddisc.build_altaz_grid_from_source(
        source=source,
        lat=observer_lat,
        lon=observer_lon,
        cloud_shells_km=CLOUD_SHELLS_KM,
    )
    logger.info("Alt/az cloud grid ready.")
    cloud_rgba = render_altaz_grid_circles(
        source.altaz_grid,
        width=radius_px * 2 + 1,
        height=radius_px * 2 + 1,
        center_alt_deg=alt,
        center_az_deg=az,
        edge_fov_deg=fov_deg,
        mask_fov_deg=fov_deg,
    )
    return _compose_on_black_background(cloud_rgba)


def save_png(image_rgba: np.ndarray, output_path: Path) -> None:
    """Save an RGBA image buffer as PNG."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not np_rgba_to_qimage(image_rgba).save(str(output_path), "PNG"):
        raise OSError(f"Failed to save image: {output_path}")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the cloud-image CLI."""
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    setup_root_logger()
    observer_lat, observer_lon = args.observer
    try:
        image_rgba = render_cloud_image(
            observer_lat=observer_lat,
            observer_lon=observer_lon,
            alt=float(args.alt),
            az=float(args.az),
            fov_deg=float(args.fov),
        )
        save_png(image_rgba, Path(args.output).expanduser())
        logger.info("Saved image: %s", Path(args.output).expanduser())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
