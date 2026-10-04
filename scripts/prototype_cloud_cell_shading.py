#!/usr/bin/env python3
"""Render sky and stars with experimental native satellite cloud voxels."""

from __future__ import annotations

# The source tree is added to sys.path below for direct script execution.
# ruff: noqa: E402

import argparse
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from PySide6.QtGui import QImage

from zstarview.astro import load_ephemeris
from zstarview.clouddisc import CloudDisc, CloudDiscConfig, CloudDiscError
from zstarview.clouddisc.altaz_grid import CloudAltAzGrid, build_altaz_grid
from zstarview.clouddisc.workers.cloud_source import build_cloud_source_fetch_request
from zstarview.clouddisc.workers.cloud_source_worker import (
    run_cloud_source_worker_process,
)
from zstarview.cli.export_image_support import EXPORT_IMAGE_METADATA_TEXT_KEY
from zstarview.location_resolver import LocationResolveError, resolve_launch_location
from zstarview.paths import CACHE_PATH, CLOUD_SHELLS_KM
from zstarview.render.cloud_shading import (
    AMBIENT_FLOOR,
    SUN_EXTINCTION,
    VIEW_EXTINCTION,
    shade_cloud_cells,
)
from zstarview.clouddisc.altaz_projection import altaz_to_bin_indices
from zstarview.render.cloud_voxels import (
    CLOUD_WHITENESS,
    ENVIRONMENT_LIGHT_FRACTION,
    SUNLIGHT_LEVELS,
    shade_native_voxels,
)
from zstarview.render.geometry import get_screen_geometry
from zstarview.render.ground_mask import inverse_project_disc
from zstarview.render.qt_image import np_rgba_to_qimage, qimage_to_np_rgba


def _parse_datetime(text: str | None) -> datetime:
    if text is None:
        return datetime.now(timezone.utc)
    value = text.strip()
    if value.endswith("Z"):
        value = value[:-1] + "+00:00"
    parsed = datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _parse_size(text: str) -> tuple[int, int]:
    try:
        width_text, height_text = text.lower().split("x", 1)
        size = int(width_text), int(height_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("size must be WIDTHxHEIGHT") from exc
    if min(size) < 64 or max(size) > 8192:
        raise argparse.ArgumentTypeError("image dimensions must be between 64 and 8192")
    return size


def _sun_altaz(when_utc: datetime, lat: float, lon: float) -> tuple[float, float]:
    import skyfield.api
    from skyfield.api import Topos

    ephemeris = load_ephemeris()
    time = skyfield.api.load.timescale().from_datetime(when_utc)
    observer = ephemeris["earth"] + Topos(latitude_degrees=lat, longitude_degrees=lon)
    apparent = observer.at(time).observe(ephemeris["sun"]).apparent()
    altitude, azimuth, _ = apparent.altaz()
    return float(altitude.degrees), float(azimuth.degrees)


def _base_sky_image(
    *,
    location: str,
    when_utc: datetime,
    alt: float,
    az: float,
    fov: float,
    content_fov: float,
    size: tuple[int, int],
    output: Path,
    sky_opacity: float = 1.0,
) -> None:
    args = [
        sys.executable,
        "-c",
        "from zstarview.cli.export_image import main; main()",
        location,
        "--datetime",
        when_utc.strftime("%Y-%m-%d %H:%M:%S UTC"),
        "--view-center-alt",
        str(alt),
        "--view-center-az",
        str(az),
        "--edge-fov-deg",
        str(fov),
        "--content-fov-deg",
        str(content_fov),
        "--image-size",
        f"{size[0]},{size[1]}",
        "--output",
        str(output),
        "--sky-opacity",
        str(sky_opacity),
        "--show-dso-initial",
        "false",
        "--show-asterisms-initial",
        "false",
        "--show-guidelines-initial",
        "false",
        "--sky-disc-altaz-rings",
        "off",
        "--sky-disc-altaz-rings-hover",
        "off",
        "--diffuse-sky-opacity",
        "0",
        "--cloud-opacity",
        "0",
        "--precipitation-opacity",
        "0",
        "--tropical-cyclone-opacity",
        "0",
        "--aircraft-opacity",
        "0",
        "--satellite-opacity",
        "0",
        "--meteor-trails-opacity",
        "0",
        "--terrain-horizon-opacity",
        "0",
        "--earth-guide-opacity",
        "0",
        "--water-surface-opacity",
        "0",
        "--road-light-opacity",
        "0",
        "--night-light-opacity",
        "0",
        "--ridge-glow-opacity",
        "0",
        "--urban-outline-opacity",
        "0",
    ]
    try:
        child_env = os.environ.copy()
        child_env["QT_QPA_PLATFORM"] = "offscreen"
        subprocess.run(args, cwd=REPO_ROOT, check=True, env=child_env)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"base sky image export failed with exit status {exc.returncode}"
        ) from exc


def _read_png_metadata(image: QImage) -> dict[str, object]:
    raw = image.text(EXPORT_IMAGE_METADATA_TEXT_KEY)
    if not raw:
        return {}
    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        return {}
    return value if isinstance(value, dict) else {}


def _draw_grid_edges(
    rgba: np.ndarray,
    grid: CloudAltAzGrid,
    altitudes: np.ndarray,
    azimuths: np.ndarray,
    inside: np.ndarray,
) -> None:
    """Overlay visible cell boundaries across the projected cloud-grid area."""
    height, width = inside.shape
    ai, zi = altaz_to_bin_indices(
        altitudes,
        azimuths,
        alt_bins=grid.amount.shape[0],
        az_bins=grid.amount.shape[1],
        alt_min_deg=grid.alt_min_deg,
        alt_max_deg=grid.alt_max_deg,
        az_min_deg=grid.az_min_deg,
        az_max_deg=grid.az_max_deg,
    )
    valid = (
        (altitudes >= grid.alt_min_deg)
        & (altitudes <= grid.alt_max_deg)
        & (azimuths >= grid.az_min_deg)
        & (azimuths <= grid.az_max_deg)
    )
    cells_alt = np.full((height, width), -1, dtype=np.int32)
    cells_az = np.full((height, width), -1, dtype=np.int32)
    flat_inside = np.flatnonzero(inside)
    cells_alt.flat[flat_inside[valid]] = ai[valid]
    cells_az.flat[flat_inside[valid]] = zi[valid]
    edges = np.zeros((height, width), dtype=bool)
    horizontal = (
        (cells_alt[:, 1:] >= 0)
        & (cells_alt[:, :-1] >= 0)
        & (
            (cells_alt[:, 1:] != cells_alt[:, :-1])
            | (cells_az[:, 1:] != cells_az[:, :-1])
        )
    )
    vertical = (
        (cells_alt[1:, :] >= 0)
        & (cells_alt[:-1, :] >= 0)
        & (
            (cells_alt[1:, :] != cells_alt[:-1, :])
            | (cells_az[1:, :] != cells_az[:-1, :])
        )
    )
    edges[:, 1:] |= horizontal
    edges[1:, :] |= vertical
    pixels = rgba[..., :3].astype(np.float32) / 255.0
    luminance = np.mean(pixels, axis=2)
    dark_line = np.array([0.10, 0.13, 0.17], dtype=np.float32)
    light_line = np.array([0.16, 0.62, 0.82], dtype=np.float32)
    line_colors = np.where((luminance > 0.35)[..., None], dark_line, light_line)
    alpha = 0.72
    pixels[edges] = pixels[edges] * (1.0 - alpha) + line_colors[edges] * alpha
    rgba[..., :3] = np.rint(np.clip(pixels, 0.0, 1.0) * 255.0).astype(np.uint8)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "location",
        help="City name or coordinates accepted by zstarview, such as @35.68,139.69",
    )
    parser.add_argument(
        "--view-alt", type=float, default=45.0, help="View center altitude in degrees"
    )
    parser.add_argument(
        "--view-az", type=float, default=180.0, help="View center azimuth in degrees"
    )
    parser.add_argument(
        "--fov", type=float, default=90.0, help="Edge field of view in degrees"
    )
    parser.add_argument(
        "--content-fov",
        type=float,
        default=115.0,
        help="Projection overscan field of view in degrees",
    )
    parser.add_argument(
        "--cloud-model",
        choices=("voxel", "altaz"),
        default="voxel",
        help="Native satellite voxels (default) or legacy angular cells",
    )
    parser.add_argument(
        "--alt-bins",
        type=int,
        default=18,
        help="Number of altitude cells (default: 18; 5-degree steps)",
    )
    parser.add_argument(
        "--az-bins",
        type=int,
        default=72,
        help="Number of azimuth cells (default: 72; 5-degree steps)",
    )
    parser.add_argument(
        "--size", type=_parse_size, default=(1200, 800), metavar="WIDTHxHEIGHT"
    )
    parser.add_argument(
        "--datetime", help="UTC ISO datetime; defaults to the current time"
    )
    parser.add_argument(
        "--sky-opacity", type=float, default=1.0, help="Sky color opacity from 0 to 1"
    )
    parser.add_argument(
        "--opacity", type=float, default=0.85, help="Cloud opacity from 0 to 1"
    )
    parser.add_argument(
        "--show-grid",
        action="store_true",
        help="Emphasize voxel edges, or overlay angular cell boundaries in altaz mode",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output PNG path")
    parser.add_argument(
        "--timeout",
        type=float,
        default=900.0,
        help="Cloud data worker timeout in seconds",
    )
    args = parser.parse_args(argv)
    if not 0.0 <= args.sky_opacity <= 1.0:
        parser.error("sky opacity must be between 0 and 1")
    if not 0.0 <= args.opacity <= 1.0:
        parser.error("opacity must be between 0 and 1")
    if args.alt_bins < 1 or args.az_bins < 1:
        parser.error("altitude and azimuth bin counts must be positive")

    when_utc = _parse_datetime(args.datetime).replace(microsecond=0)
    try:
        location = resolve_launch_location(args.location)
    except LocationResolveError as exc:
        parser.error(f"could not resolve location {args.location!r}: {exc}")
    resolved_token = f"@{location.lat:.8f},{location.lon:.8f}"
    sun_alt, sun_az = _sun_altaz(when_utc, location.lat, location.lon)
    clouddisc = CloudDisc(
        CloudDiscConfig(
            cache_dir=CACHE_PATH,
            sat_priority=("AUTO",),
            bt_warm_k=310.0,
            bt_cold_k=190.0,
            alt_min_deg=0.0,
            search_back_minutes=120,
        )
    )
    request = build_cloud_source_fetch_request(
        lat=location.lat,
        lon=location.lon,
        when_utc=when_utc,
        cloud_shells_km=CLOUD_SHELLS_KM,
    )
    try:
        source = run_cloud_source_worker_process(
            clouddisc,
            request,
            request_id=os.getpid(),
            timeout_s=max(1.0, float(args.timeout)),
            skip_altaz_grid=True,
        )
    except CloudDiscError as exc:
        parser.error(f"cloud data is unavailable for this location: {exc}")
    grid = None
    if args.cloud_model == "altaz":
        grid = build_altaz_grid(
            source,
            location.lat,
            location.lon,
            shells_km=CLOUD_SHELLS_KM,
            alt_bins=args.alt_bins,
            az_bins=args.az_bins,
            preserve_physical_shells=True,
        )
        if not isinstance(grid, CloudAltAzGrid):
            parser.error("could not build the prototype cloud grid")
    voxel_info = None

    output_path = args.output.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    width, height = args.size
    with tempfile.TemporaryDirectory(prefix="zstarview-cloud-prototype-") as temp_dir:
        base_path = Path(temp_dir) / "sky.png"
        try:
            _base_sky_image(
                location=resolved_token,
                when_utc=when_utc,
                alt=args.view_alt,
                az=args.view_az,
                fov=args.fov,
                content_fov=args.content_fov,
                size=args.size,
                output=base_path,
                sky_opacity=args.sky_opacity,
            )
        except RuntimeError as exc:
            parser.error(str(exc))
        base_image = QImage(str(base_path))
        if base_image.isNull():
            parser.error(f"failed to read base sky image: {base_path}")
        rgba = qimage_to_np_rgba(base_image)
        geometry = get_screen_geometry(
            width,
            height,
            args.view_alt,
            edge_fov_deg=args.fov,
            content_fov_deg=args.content_fov,
        )
        altitudes, azimuths, inside = inverse_project_disc(
            width,
            height,
            geometry,
            (args.view_alt, args.view_az),
            edge_fov_deg=args.fov,
            content_fov_deg=args.content_fov,
        )
        flat_indices = np.flatnonzero(inside)
        rgb = rgba[..., :3].reshape((-1, 3)).astype(np.float32) / 255.0
        if args.cloud_model == "voxel":
            rgb[flat_indices], voxel_info = shade_native_voxels(
                source,
                location.lat,
                location.lon,
                altitudes,
                azimuths,
                sun_alt,
                sun_az,
                rgb[flat_indices],
                opacity=args.opacity,
                show_grid=args.show_grid,
            )
        else:
            rgb[flat_indices] = shade_cloud_cells(
                grid,
                rgb[flat_indices],
                altitudes,
                azimuths,
                sun_alt,
                sun_az,
                opacity=args.opacity,
            )
        rgba[..., :3] = np.rint(
            np.clip(rgb.reshape((height, width, 3)), 0.0, 1.0) * 255.0
        ).astype(np.uint8)
        if args.show_grid and grid is not None:
            _draw_grid_edges(rgba, grid, altitudes, azimuths, inside)
        output_image = np_rgba_to_qimage(rgba)
        if not output_image.save(str(output_path), "PNG"):
            parser.error(f"failed to write output image: {output_path}")

    metadata = _read_png_metadata(base_image)
    payload: dict[str, object] = {
        "prototype": "cloud-native-voxel-v4" if voxel_info else "cloud-cell-shading-v1",
        "render_time_utc": when_utc.isoformat().replace("+00:00", "Z"),
        "cloud_observation_time_utc": source.time_utc.isoformat().replace(
            "+00:00", "Z"
        ),
        "location": {
            "name": location.display_name,
            "lat_deg": location.lat,
            "lon_deg": location.lon,
        },
        "view": {
            "alt_deg": args.view_alt,
            "az_deg": args.view_az,
            "edge_fov_deg": args.fov,
            "size": [width, height],
        },
        "sun": {"alt_deg": sun_alt, "az_deg": sun_az},
        "cloud_source": {
            "satellite": source.satellite,
            "product": source.product,
            "coverage_ratio": voxel_info["coverage_ratio"]
            if voxel_info
            else grid.coverage_ratio,
        },
        "grid": voxel_info
        or {
            "altitude_bins": int(grid.amount.shape[0]),
            "azimuth_bins": int(grid.amount.shape[1]),
            "shell_altitudes_km": [
                round(value - 6371.0, 3) for value in grid.shells_km
            ],
        },
        "grid_overlay": bool(args.show_grid),
        "model": {
            "sun_extinction": SUN_EXTINCTION,
            "view_extinction": VIEW_EXTINCTION,
            "ambient_floor": AMBIENT_FLOOR,
            "opacity": args.opacity,
            "sky_opacity": args.sky_opacity,
            "voxel_color": "sunlight-environment-whiteness-table"
            if voxel_info
            else None,
            "sunlight_levels": SUNLIGHT_LEVELS.tolist() if voxel_info else None,
            "environment_light_fraction": ENVIRONMENT_LIGHT_FRACTION
            if voxel_info
            else None,
            "cloud_whiteness": CLOUD_WHITENESS.tolist() if voxel_info else None,
            "sunlight_baseline": "max(0, sin(sun_altitude))" if voxel_info else None,
        },
        "base_image_metadata": metadata,
    }
    sidecar = output_path.with_suffix(".json")
    sidecar.write_text(
        json.dumps(payload, ensure_ascii=True, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Saved image: {output_path}")
    print(f"Saved metadata: {sidecar}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
