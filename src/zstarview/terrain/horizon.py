from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np


EARTH_MEAN_RADIUS_M = 6_371_008.8


@dataclass(frozen=True)
class ObserverLocation:
    latitude_deg: float
    longitude_deg: float
    observer_ground_m: float
    observer_eye_m: float

    @property
    def observer_elevation_m(self) -> float:
        return self.observer_ground_m + self.observer_eye_m


@dataclass(frozen=True)
class HorizonProfilePoint:
    azimuth_deg: float
    altitude_deg: float
    distance_m: float
    latitude_deg: float
    longitude_deg: float
    terrain_elevation_m: float


def build_distance_samples(max_distance_km: float, sample_step_m: float) -> np.ndarray:
    max_distance_m = max_distance_km * 1000.0
    if max_distance_m <= 0.0:
        raise ValueError("max_distance_km must be positive.")
    if sample_step_m <= 0.0:
        raise ValueError("sample_step_m must be positive.")
    count = max(1, int(math.ceil(max_distance_m / sample_step_m)))
    return np.linspace(sample_step_m, max_distance_m, num=count, dtype=np.float64)


def compute_apparent_altitudes(
    *,
    observer_elevation_m: float,
    target_elevation_m: np.ndarray,
    surface_distance_m: np.ndarray,
    earth_radius_m: float,
    refraction_coefficient: float,
) -> np.ndarray:
    if earth_radius_m <= 0.0:
        raise ValueError("earth_radius_m must be positive.")
    if refraction_coefficient >= 1.0:
        raise ValueError("refraction_coefficient must be < 1.0.")

    effective_radius_m = earth_radius_m / max(1.0 - refraction_coefficient, 1e-6)
    theta = np.asarray(surface_distance_m, dtype=np.float64) / effective_radius_m
    horizontal_m = effective_radius_m * np.sin(theta)
    vertical_m = (
        np.asarray(target_elevation_m, dtype=np.float64)
        - float(observer_elevation_m)
        + effective_radius_m * (np.cos(theta) - 1.0)
    )
    return np.degrees(np.arctan2(vertical_m, horizontal_m))


def compute_horizon_profile(
    *,
    dem_grid,
    geod,
    observer: ObserverLocation,
    azimuth_step_deg: float,
    distance_samples_m: np.ndarray,
    dem_resampling: str,
    earth_radius_m: float,
    refraction_coefficient: float,
) -> list[HorizonProfilePoint]:
    if azimuth_step_deg <= 0.0:
        raise ValueError("azimuth_step_deg must be positive.")

    azimuths = np.arange(0.0, 360.0, azimuth_step_deg, dtype=np.float64)
    azimuth_grid_deg, distance_grid_m = np.meshgrid(
        azimuths,
        distance_samples_m,
        indexing="ij",
    )
    flat_count = azimuth_grid_deg.size
    ray_lon_flat, ray_lat_flat, _ = geod.fwd(
        np.full(flat_count, observer.longitude_deg, dtype=np.float64),
        np.full(flat_count, observer.latitude_deg, dtype=np.float64),
        azimuth_grid_deg.ravel(),
        distance_grid_m.ravel(),
    )
    ray_lon_deg = np.asarray(ray_lon_flat, dtype=np.float64).reshape(azimuth_grid_deg.shape)
    ray_lat_deg = np.asarray(ray_lat_flat, dtype=np.float64).reshape(azimuth_grid_deg.shape)
    terrain_m = dem_grid.sample_lonlat(
        ray_lon_deg,
        ray_lat_deg,
        method=dem_resampling,
    )
    altitude_deg = compute_apparent_altitudes(
        observer_elevation_m=observer.observer_elevation_m,
        target_elevation_m=terrain_m,
        surface_distance_m=distance_grid_m,
        earth_radius_m=earth_radius_m,
        refraction_coefficient=refraction_coefficient,
    )
    valid = np.isfinite(terrain_m)
    if not np.any(valid):
        raise ValueError("The DEM did not provide any valid samples for the scan.")

    altitude_masked = np.where(valid, altitude_deg, -np.inf)
    valid_rows = np.any(valid, axis=1)
    peak_indices = np.argmax(altitude_masked, axis=1)

    points: list[HorizonProfilePoint] = []
    for row_index, azimuth_deg in enumerate(azimuths):
        if not valid_rows[row_index]:
            continue
        peak_index = int(peak_indices[row_index])
        points.append(
            HorizonProfilePoint(
                azimuth_deg=float(azimuth_deg),
                altitude_deg=float(altitude_deg[row_index, peak_index]),
                distance_m=float(distance_grid_m[row_index, peak_index]),
                latitude_deg=float(ray_lat_deg[row_index, peak_index]),
                longitude_deg=float(ray_lon_deg[row_index, peak_index]),
                terrain_elevation_m=float(terrain_m[row_index, peak_index]),
            )
        )

    if not points:
        raise ValueError("The DEM did not provide any valid samples for the scan.")
    return points


def reduce_profile_to_altaz(points: Sequence[HorizonProfilePoint]) -> list[tuple[float, float]]:
    return [(float(point.altitude_deg), float(point.azimuth_deg)) for point in points]
