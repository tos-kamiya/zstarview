#!/usr/bin/env python3
"""Render a small spherical-Earth single-scattering sky experiment.

This is deliberately independent of the application renderer.  Each pixel is
an altitude/azimuth direction.  Sunlight is attenuated on its way to an
atmospheric sample, scattered toward the observer, attenuated again on the
view path, and integrated up to an atmosphere shell 500 km above the surface.
The direct solar beam is never added to the result.

Examples:
  uv run -p .venv/bin/python scripts/atmospheric_sky_experiment.py
  uv run -p .venv/bin/python scripts/atmospheric_sky_experiment.py \
      --sun-altitudes 60 15 0 -6 --outfile /tmp/sky.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from PIL import Image, ImageDraw
import numpy as np


EARTH_RADIUS_KM = 6371.0
ATMOSPHERE_TOP_KM = 500.0
SCALE_HEIGHT_KM = 8.0

# Relative Rayleigh coefficients for R, G, B.  They are intentionally
# explicit experiment inputs, not display-tuning constants.
SCATTERING_RGB = np.array([0.25, 0.58, 1.0], dtype=np.float32)
EXTINCTION_RGB = SCATTERING_RGB.copy()
SUN_RADIANCE_RGB = np.ones(3, dtype=np.float32)
OPTICAL_DEPTH_SCALE = 0.018
DEFAULT_SUN_ALTITUDES = tuple(float(value) / 2.0 for value in range(12, -13, -1))


def _ray_shell_distance(origin: np.ndarray, direction: np.ndarray, radius: float) -> np.ndarray:
    """Return the forward intersection distance with a sphere."""
    b = np.sum(origin * direction, axis=-1)
    c = np.sum(origin * origin, axis=-1) - radius * radius
    discriminant = np.maximum(0.0, b * b - c)
    return -b + np.sqrt(discriminant)


def _sun_visibility(point: np.ndarray, sun_direction: np.ndarray) -> np.ndarray:
    """Whether the ray from each point toward the Sun clears the Earth."""
    b = np.sum(point * sun_direction[None, :], axis=-1)
    closest_distance = np.linalg.norm(point, axis=-1)
    outward = b >= 0.0
    projected_distance = np.sqrt(np.maximum(0.0, closest_distance**2 - b**2))
    return outward | (projected_distance >= EARTH_RADIUS_KM)


def _density(position: np.ndarray) -> np.ndarray:
    height = np.maximum(0.0, np.linalg.norm(position, axis=-1) - EARTH_RADIUS_KM)
    return np.exp(-height / SCALE_HEIGHT_KM).astype(np.float32)


def _sun_optical_depth(point: np.ndarray, sun_direction: np.ndarray, steps: int) -> np.ndarray:
    """Numerically integrate density from points toward the Sun."""
    distances = _ray_shell_distance(point, sun_direction[None, :], EARTH_RADIUS_KM + ATMOSPHERE_TOP_KM)
    u = (np.arange(steps, dtype=np.float32) + 0.5) / steps
    samples = point[:, None, :] + distances[:, None, None] * u[None, :, None] * sun_direction[None, None, :]
    column = np.mean(_density(samples.reshape(-1, 3)).reshape(point.shape[0], steps), axis=1)
    # Return normalized column density; the RGB extinction scale is applied
    # once, together with the view-path column density, by the caller.
    return (column * distances).astype(np.float32)


def _render_panel(
    width: int,
    height: int,
    sun_alt_deg: float,
    sun_az_deg: float,
    view_steps: int,
    sun_steps: int,
) -> np.ndarray:
    """Render one altitude/azimuth sky map as uint8 RGB."""
    azimuth = np.linspace(0.0, 360.0, width, endpoint=False, dtype=np.float32)
    altitude = np.linspace(90.0, 0.0, height, dtype=np.float32)
    az_grid, alt_grid = np.meshgrid(azimuth, altitude)

    alt = np.radians(alt_grid.reshape(-1))
    az = np.radians(az_grid.reshape(-1))
    view_direction = np.stack(
        [np.cos(alt) * np.sin(az), np.cos(alt) * np.cos(az), np.sin(alt)], axis=1
    ).astype(np.float32)
    observer = np.array([0.0, 0.0, EARTH_RADIUS_KM], dtype=np.float32)
    sun_alt = math.radians(sun_alt_deg)
    sun_az = math.radians(sun_az_deg)
    sun_direction = np.array(
        [math.cos(sun_alt) * math.sin(sun_az), math.cos(sun_alt) * math.cos(sun_az), math.sin(sun_alt)],
        dtype=np.float32,
    )

    max_distance = _ray_shell_distance(
        np.broadcast_to(observer, view_direction.shape), view_direction, EARTH_RADIUS_KM + ATMOSPHERE_TOP_KM
    )
    sample_fraction = (np.arange(view_steps, dtype=np.float32) + 0.5) / view_steps
    points = observer[None, None, :] + max_distance[:, None, None] * sample_fraction[None, :, None] * view_direction[:, None, :]
    densities = _density(points.reshape(-1, 3)).reshape(-1, view_steps)
    ds = max_distance / view_steps

    # Accumulate view-path optical depth and in-scattered radiance in chunks
    # to keep this exploratory script small in memory.
    rgb = np.zeros((view_direction.shape[0], 3), dtype=np.float32)
    for start in range(0, view_direction.shape[0], 4096):
        end = min(view_direction.shape[0], start + 4096)
        chunk_points = points[start:end].reshape(-1, 3)
        chunk_density = densities[start:end]
        chunk_ds = ds[start:end]
        view_tau = np.cumsum(chunk_density * chunk_ds[:, None], axis=1) - 0.5 * chunk_density * chunk_ds[:, None]
        sun_tau = _sun_optical_depth(chunk_points, sun_direction, sun_steps).reshape(end - start, view_steps)
        visible = _sun_visibility(chunk_points, sun_direction).reshape(end - start, view_steps)

        # Rayleigh phase function. Incoming propagation is -sun_direction;
        # outgoing propagation is view_direction toward the observer.
        directions = view_direction[start:end]
        cos_theta = np.sum((-sun_direction)[None, None, :] * directions[:, None, :], axis=-1)
        phase = 3.0 / (16.0 * math.pi) * (1.0 + cos_theta**2)
        transmittance = np.exp(-OPTICAL_DEPTH_SCALE * (sun_tau[:, :, None] * EXTINCTION_RGB[None, None, :]))
        transmittance *= np.exp(-OPTICAL_DEPTH_SCALE * view_tau[:, :, None] * EXTINCTION_RGB[None, None, :])
        contribution = (
            chunk_density[:, :, None]
            * SCATTERING_RGB[None, None, :]
            * SUN_RADIANCE_RGB[None, None, :]
            * phase[:, :, None]
            * transmittance
            * visible[:, :, None]
            * chunk_ds[:, None, None]
        )
        rgb[start:end] = np.sum(contribution, axis=1)

    # Exposure is only a display conversion. It does not alter the sky model.
    rgb = np.clip(1.0 - np.exp(-rgb * 2.8), 0.0, 1.0)
    return np.round(rgb.reshape(height, width, 3) * 255.0).astype(np.uint8)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outfile", type=Path, default=Path("atmospheric_sky_experiment.png"))
    parser.add_argument(
        "--sun-altitudes",
        type=float,
        nargs="+",
        default=DEFAULT_SUN_ALTITUDES,
        help="Solar altitudes for panels (default: 6 to -6 degrees in 0.5-degree steps).",
    )
    parser.add_argument("--sun-azimuth", type=float, default=180.0)
    parser.add_argument("--width", type=int, default=240)
    parser.add_argument("--height", type=int, default=100)
    parser.add_argument("--view-steps", type=int, default=64)
    parser.add_argument("--sun-steps", type=int, default=24)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    if args.width < 1 or args.height < 1 or args.view_steps < 1 or args.sun_steps < 1:
        raise SystemExit("width, height, view-steps, and sun-steps must be positive")

    panels = [
        _render_panel(args.width, args.height, altitude, args.sun_azimuth, args.view_steps, args.sun_steps)
        for altitude in args.sun_altitudes
    ]
    margin = 24
    canvas = Image.new("RGB", (args.width * len(panels), args.height + margin), (8, 8, 12))
    for index, (altitude, panel) in enumerate(zip(args.sun_altitudes, panels)):
        image = Image.fromarray(panel, mode="RGB")
        canvas.paste(image, (index * args.width, margin))
        ImageDraw.Draw(canvas).text((index * args.width + 4, 4), f"Sun alt {altitude:g} deg", fill=(235, 235, 235))
    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(args.outfile)
    print(f"wrote {args.outfile}")


if __name__ == "__main__":
    main()
