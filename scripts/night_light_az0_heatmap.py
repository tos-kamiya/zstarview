#!/usr/bin/env python3
"""Render a night-light az=0 heatmap with signed altitude on the x-axis.

The plot shows distance samples on the y-axis and altitude on the x-axis from
-90 to +90 degrees, so it is easier to see whether the accumulation leaks above
or below the horizon.

Example:
  python scripts/night_light_az0_heatmap.py \
    --observer-lat 35.6581 \
    --observer-lon 139.7456 \
    --observer-eye-m 334.0 \
    --azimuth-deg 0.0 \
    --outfile night-light-az0-heatmap.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from zstarview.night_lights import _sample_ray_brightness_curve  # noqa: E402
from zstarview.paths import COPERNICUS_DEM_CACHE_DIR, NIGHT_LIGHTS_CACHE_DIR  # noqa: E402
from zstarview.terrain import (  # noqa: E402
    EARTH_MEAN_RADIUS_M,
    GeoTiffDem,
    ObserverLocation,
    WGS84_GEOD,
    build_ray_scan_grid,
    build_download_bbox,
    compute_apparent_altitudes,
)
from zstarview.terrain.dem import collect_copernicus_tile_keys, sample_ground_elevation  # noqa: E402


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observer-lat", type=float, required=True)
    parser.add_argument("--observer-lon", type=float, required=True)
    parser.add_argument("--observer-eye-m", type=float, default=0.0)
    parser.add_argument("--azimuth-deg", type=float, default=0.0)
    parser.add_argument("--outfile", default="night-light-az0-heatmap.png")
    parser.add_argument("--max-distance-km", type=float, default=128.0)
    parser.add_argument("--distance-step-km", type=float, default=1.0)
    parser.add_argument("--alt-step-deg", type=float, default=0.25)
    return parser


def _load_dem(observer_lat_deg: float, observer_lon_deg: float) -> tuple[GeoTiffDem, object]:
    dem_cache_root = Path(COPERNICUS_DEM_CACHE_DIR)
    bbox = build_download_bbox(
        lat_deg=float(observer_lat_deg),
        lon_deg=float(observer_lon_deg),
        radius_km=138.0,
    )
    dem_tile_paths = [
        dem_cache_root / key
        for key in collect_copernicus_tile_keys(bbox)
        if (dem_cache_root / key).exists()
    ]
    if not dem_tile_paths:
        raise FileNotFoundError(
            "No Copernicus DEM tiles found in cache. Download terrain data first."
        )
    dem = GeoTiffDem(dem_tile_paths, default_elevation_m=0.0)
    dem_grid = dem.build_grid(bbox)
    return dem, dem_grid


def _load_night_light_tile_paths() -> dict[str, Path]:
    night_light_root = Path(NIGHT_LIGHTS_CACHE_DIR) / "2016_grayscale_500m"
    tile_paths: dict[str, Path] = {}
    for path in sorted(night_light_root.glob("BlackMarble_2016_*_geo_gray.tif")):
        tile_name = path.stem.replace("BlackMarble_2016_", "").replace("_geo_gray", "")
        tile_paths[tile_name] = path
    if not tile_paths:
        raise FileNotFoundError(
            "No night-light GeoTIFF tiles found in cache. Download night-light data first."
        )
    return tile_paths


def _build_heatmap(
    *,
    observer_lat_deg: float,
    observer_lon_deg: float,
    observer_eye_m: float,
    azimuth_deg: float,
    max_distance_km: float,
    distance_step_km: float,
    alt_step_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dem, dem_grid = _load_dem(observer_lat_deg, observer_lon_deg)
    try:
        ground_m = sample_ground_elevation(
            dem_grid,
            latitude_deg=float(observer_lat_deg),
            longitude_deg=float(observer_lon_deg),
            dem_resampling="bilinear",
        )
        observer = ObserverLocation(
            latitude_deg=float(observer_lat_deg),
            longitude_deg=float(observer_lon_deg),
            observer_ground_m=float(ground_m),
            observer_eye_m=float(observer_eye_m),
        )
        distances_m = np.arange(
            float(distance_step_km) * 1000.0,
            float(max_distance_km) * 1000.0 + 0.5,
            float(distance_step_km) * 1000.0,
            dtype=np.float64,
        )
        ray_scan = build_ray_scan_grid(
            geod=WGS84_GEOD,
            observer_latitude_deg=float(observer_lat_deg),
            observer_longitude_deg=float(observer_lon_deg),
            azimuth_step_deg=360.0,
            distance_samples_m=distances_m,
        )
        terrain_m = dem_grid.sample_lonlat(
            ray_scan.ray_lon_deg,
            ray_scan.ray_lat_deg,
            method="bilinear",
        )
        apparent_alt_deg = compute_apparent_altitudes(
            observer_elevation_m=observer.observer_elevation_m,
            target_elevation_m=terrain_m,
            surface_distance_m=ray_scan.distance_grid_m,
            earth_radius_m=EARTH_MEAN_RADIUS_M,
            refraction_coefficient=0.13,
        )[0]
    finally:
        dem.close()

    tile_paths = _load_night_light_tile_paths()
    brightness_curve = _sample_ray_brightness_curve(
        tile_paths=tile_paths,
        observer_lat_deg=float(observer_lat_deg),
        observer_lon_deg=float(observer_lon_deg),
        azimuth_deg=float(azimuth_deg),
        distances_m=distances_m,
    )

    step_brightness = np.clip(
        np.diff(np.concatenate([[0.0], np.asarray(brightness_curve, dtype=np.float64)])),
        0.0,
        None,
    )

    alt_bins = np.arange(-90.0, 90.0 + alt_step_deg, alt_step_deg, dtype=np.float64)
    heat = np.zeros((distances_m.size, alt_bins.size), dtype=np.float32)
    nonzero = step_brightness[step_brightness > 0.0]
    scale = float(np.percentile(nonzero, 95)) if nonzero.size else 1.0
    scale = max(scale, 1.0e-12)
    sigma_deg = 0.85
    for i, (floor_alt, amp_raw) in enumerate(zip(apparent_alt_deg, step_brightness, strict=True)):
        amp = float(np.log1p(float(amp_raw) / scale * 12.0) / np.log1p(12.0))
        amp = max(0.0, min(1.0, amp))
        if amp <= 0.0:
            continue
        profile = amp * np.exp(-0.5 * np.square((alt_bins - float(floor_alt)) / sigma_deg))
        heat[i, :] = np.maximum(heat[i, :], profile.astype(np.float32))

    return distances_m / 1000.0, alt_bins, heat, apparent_alt_deg


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()

    distances_km, alt_bins, heat, apparent_alt_deg = _build_heatmap(
        observer_lat_deg=float(args.observer_lat),
        observer_lon_deg=float(args.observer_lon),
        observer_eye_m=float(args.observer_eye_m),
        azimuth_deg=float(args.azimuth_deg),
        max_distance_km=float(args.max_distance_km),
        distance_step_km=float(args.distance_step_km),
        alt_step_deg=float(args.alt_step_deg),
    )

    display = np.power(np.clip(heat, 0.0, 1.0), 0.72)
    fig, ax = plt.subplots(figsize=(13.0, 7.0), dpi=160)
    mesh = ax.imshow(
        display,
        origin="lower",
        aspect="auto",
        extent=[float(alt_bins[0]), float(alt_bins[-1]), float(distances_km[0]), float(distances_km[-1])],
        cmap="magma",
        vmin=0.0,
        vmax=1.0,
    )
    fig.colorbar(mesh, ax=ax, label="brightness")

    ax.axvline(0.0, color="#7fd3ff", linewidth=1.2, alpha=0.9)
    ax.plot(apparent_alt_deg, distances_km, color="#43ffd0", linewidth=1.4)
    ax.set_xlabel("apparent altitude (deg)")
    ax.set_ylabel("distance (km)")
    ax.set_title("az={:.1f} brightness map: distance vs signed altitude".format(float(args.azimuth_deg)))
    ax.set_xlim(float(alt_bins[0]), float(alt_bins[-1]))
    ax.set_ylim(float(distances_km[0]), float(distances_km[-1]))
    ax.grid(True, color="white", alpha=0.12, linewidth=0.7)
    fig.tight_layout()
    fig.savefig(args.outfile, facecolor="black")
    plt.close(fig)
    print(f"wrote {args.outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
