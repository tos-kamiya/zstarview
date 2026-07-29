#!/usr/bin/env python3
"""Render the atmospheric sky experiment with an additional Mie component.

The original ``atmospheric_sky_experiment.py`` is intentionally left
unchanged.  This companion script adds a low-altitude aerosol layer and a
forward-peaked Henyey-Greenstein phase function to the Rayleigh model.

Example:
  uv run -p .venv/bin/python scripts/atmospheric_sky_mie_experiment.py \
      --outfile /tmp/atmospheric_sky_mie.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from atmospheric_sky_experiment import (
    DEFAULT_SUN_ALTITUDES,
    EARTH_RADIUS_KM,
    OPTICAL_DEPTH_SCALE,
    SCALE_HEIGHT_KM,
    SUN_RADIANCE_RGB,
    _density,
    _ray_shell_distance,
    _sun_visibility,
)


AEROSOL_SCALE_HEIGHT_KM = 1.4
AEROSOL_DENSITY_RATIO = 0.045
AEROSOL_SCATTERING_RGB = np.array([0.86, 0.90, 0.96], dtype=np.float32)
AEROSOL_EXTINCTION_RGB = AEROSOL_SCATTERING_RGB.copy()
MIE_ANISOTROPY = 0.76


def _aerosol_density(position: np.ndarray) -> np.ndarray:
    height = np.maximum(0.0, np.linalg.norm(position, axis=-1) - EARTH_RADIUS_KM)
    return (AEROSOL_DENSITY_RATIO * np.exp(-height / AEROSOL_SCALE_HEIGHT_KM)).astype(np.float32)


def _column_densities(
    point: np.ndarray,
    sun_direction: np.ndarray,
    steps: int,
    atmosphere_top_km: float,
) -> tuple[np.ndarray, np.ndarray]:
    distances = _ray_shell_distance(
        point,
        sun_direction[None, :],
        EARTH_RADIUS_KM + atmosphere_top_km,
    )
    fraction = (np.arange(steps, dtype=np.float32) + 0.5) / steps
    samples = point[:, None, :] + distances[:, None, None] * fraction[None, :, None] * sun_direction[None, None, :]
    flat = samples.reshape(-1, 3)
    shape = (point.shape[0], steps)
    return (
        (np.mean(_density(flat).reshape(shape), axis=1) * distances).astype(np.float32),
        (np.mean(_aerosol_density(flat).reshape(shape), axis=1) * distances).astype(np.float32),
    )


def _henyey_greenstein(cos_theta: np.ndarray, anisotropy: float) -> np.ndarray:
    g = float(anisotropy)
    denominator = np.power(1.0 + g * g - 2.0 * g * cos_theta, 1.5)
    return ((1.0 - g * g) / (4.0 * math.pi * denominator)).astype(np.float32)


def _render_panel(
    width: int,
    height: int,
    sun_alt_deg: float,
    sun_az_deg: float,
    view_steps: int,
    sun_steps: int,
    atmosphere_top_km: float,
) -> np.ndarray:
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
        np.broadcast_to(observer, view_direction.shape),
        view_direction,
        EARTH_RADIUS_KM + atmosphere_top_km,
    )
    fraction = (np.arange(view_steps, dtype=np.float32) + 0.5) / view_steps
    points = observer[None, None, :] + max_distance[:, None, None] * fraction[None, :, None] * view_direction[:, None, :]
    shape = (view_direction.shape[0], view_steps)
    rayleigh_density = _density(points.reshape(-1, 3)).reshape(shape)
    aerosol_density = _aerosol_density(points.reshape(-1, 3)).reshape(shape)
    ds = max_distance / view_steps
    rgb = np.zeros((view_direction.shape[0], 3), dtype=np.float32)

    for start in range(0, view_direction.shape[0], 4096):
        end = min(view_direction.shape[0], start + 4096)
        chunk_points = points[start:end].reshape(-1, 3)
        chunk_rayleigh = rayleigh_density[start:end]
        chunk_aerosol = aerosol_density[start:end]
        chunk_ds = ds[start:end]
        rayleigh_view_column = np.cumsum(chunk_rayleigh * chunk_ds[:, None], axis=1) - 0.5 * chunk_rayleigh * chunk_ds[:, None]
        aerosol_view_column = np.cumsum(chunk_aerosol * chunk_ds[:, None], axis=1) - 0.5 * chunk_aerosol * chunk_ds[:, None]
        rayleigh_sun_column, aerosol_sun_column = _column_densities(
            chunk_points,
            sun_direction,
            sun_steps,
            atmosphere_top_km,
        )
        rayleigh_sun_column = rayleigh_sun_column.reshape(end - start, view_steps)
        aerosol_sun_column = aerosol_sun_column.reshape(end - start, view_steps)
        visible = _sun_visibility(chunk_points, sun_direction).reshape(end - start, view_steps)

        direction = view_direction[start:end]
        cos_phase = np.sum(sun_direction[None, None, :] * direction[:, None, :], axis=-1)
        rayleigh_phase = 3.0 / (16.0 * math.pi) * (1.0 + cos_phase**2)
        mie_phase = _henyey_greenstein(cos_phase, MIE_ANISOTROPY)
        rayleigh_transmission = np.exp(
            -OPTICAL_DEPTH_SCALE
            * (rayleigh_sun_column[:, :, None] + rayleigh_view_column[:, :, None])
            * np.array([0.25, 0.58, 1.0], dtype=np.float32)[None, None, :]
        )
        aerosol_transmission = np.exp(
            -OPTICAL_DEPTH_SCALE
            * (aerosol_sun_column[:, :, None] + aerosol_view_column[:, :, None])
            * AEROSOL_EXTINCTION_RGB[None, None, :]
        )
        transmission = rayleigh_transmission * aerosol_transmission
        contribution = (
            chunk_rayleigh[:, :, None] * np.array([0.25, 0.58, 1.0], dtype=np.float32)[None, None, :] * rayleigh_phase[:, :, None]
            + chunk_aerosol[:, :, None] * AEROSOL_SCATTERING_RGB[None, None, :] * mie_phase[:, :, None]
        )
        rgb[start:end] = np.sum(
            contribution * SUN_RADIANCE_RGB[None, None, :] * transmission * visible[:, :, None] * chunk_ds[:, None, None],
            axis=1,
        )

    rgb = np.clip(1.0 - np.exp(-rgb * 2.8), 0.0, 1.0)
    return np.round(rgb.reshape(height, width, 3) * 255.0).astype(np.uint8)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outfile", type=Path, default=Path("atmospheric_sky_mie_experiment.png"))
    parser.add_argument("--sun-altitudes", type=float, nargs="+", default=DEFAULT_SUN_ALTITUDES)
    parser.add_argument("--sun-azimuth", type=float, default=180.0)
    parser.add_argument("--width", type=int, default=240)
    parser.add_argument("--height", type=int, default=100)
    parser.add_argument("--view-steps", type=int, default=64)
    parser.add_argument("--sun-steps", type=int, default=24)
    parser.add_argument(
        "--atmosphere-top-km",
        type=float,
        default=100.0,
        help="Atmosphere shell height in km (default: 100).",
    )
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    if min(args.width, args.height, args.view_steps, args.sun_steps, args.atmosphere_top_km) < 1:
        raise SystemExit("width, height, view-steps, sun-steps, and atmosphere-top-km must be positive")
    panels = [
        _render_panel(
            args.width,
            args.height,
            altitude,
            args.sun_azimuth,
            args.view_steps,
            args.sun_steps,
            args.atmosphere_top_km,
        )
        for altitude in args.sun_altitudes
    ]
    margin = 24
    canvas = Image.new("RGB", (args.width * len(panels), args.height + margin), (8, 8, 12))
    draw = ImageDraw.Draw(canvas)
    for index, (altitude, panel) in enumerate(zip(args.sun_altitudes, panels)):
        canvas.paste(Image.fromarray(panel, mode="RGB"), (index * args.width, margin))
        draw.text((index * args.width + 4, 4), f"Sun alt {altitude:g} deg", fill=(235, 235, 235))
    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(args.outfile)
    print(f"wrote {args.outfile}")


if __name__ == "__main__":
    main()
