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


@dataclass(frozen=True)
class HorizonLayerSet:
    main_profile: list[HorizonProfilePoint]
    secondary_layers: list[list[HorizonProfilePoint]]


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


def _select_secondary_peak_indices(
    altitude_row_deg: np.ndarray,
    distance_row_m: np.ndarray,
    *,
    main_peak_index: int,
    min_prominence_deg: float,
    min_drop_deg: float,
    min_separation_m: float,
    prominence_window: int = 8,
    max_secondary_peaks: int = 2,
) -> list[int]:
    values = np.asarray(altitude_row_deg, dtype=np.float64)
    distances = np.asarray(distance_row_m, dtype=np.float64)
    if values.size < 3:
        return []

    candidate_indices: list[int] = []
    for index in range(1, values.size - 1):
        center = float(values[index])
        if not (center >= float(values[index - 1]) and center > float(values[index + 1])):
            continue

        left_start = max(0, index - int(prominence_window))
        right_end = min(values.size, index + int(prominence_window) + 1)
        left_slice = values[left_start:index]
        right_slice = values[index + 1:right_end]
        if left_slice.size == 0 or right_slice.size == 0:
            continue

        left_min = float(np.min(left_slice))
        right_min = float(np.min(right_slice))
        prominence_deg = center - max(left_min, right_min)
        if prominence_deg < float(min_prominence_deg):
            continue
        if (center - right_min) < float(min_drop_deg):
            continue
        candidate_indices.append(index)

    if not candidate_indices:
        return []

    selected: list[int] = []
    for index in sorted(
        candidate_indices,
        key=lambda idx: (
            -float(values[idx]),
            -float(values[idx] - float(np.min(values[max(0, idx - prominence_window):idx]))),
            float(distances[idx]),
        ),
    ):
        if index == int(main_peak_index):
            continue
        if any(abs(float(distances[index]) - float(distances[existing])) < float(min_separation_m) for existing in selected):
            continue
        selected.append(index)
        if len(selected) >= int(max_secondary_peaks):
            break

    selected.sort(key=lambda idx: float(distances[idx]))
    return selected


def _prune_secondary_peak_indices_by_visibility(
    altitude_row_deg: np.ndarray,
    secondary_peak_indices: Sequence[int],
) -> list[int]:
    values = np.asarray(altitude_row_deg, dtype=np.float64)
    visible_indices: list[int] = []
    best_altitude_deg = -np.inf
    for index in sorted(int(value) for value in secondary_peak_indices):
        if index < 0 or index >= values.size:
            continue
        altitude_deg = float(values[index])
        if altitude_deg <= best_altitude_deg:
            continue
        visible_indices.append(index)
        best_altitude_deg = altitude_deg
    return visible_indices


def _should_break_secondary_ridge(
    previous_distance_m: float,
    current_distance_m: float,
    *,
    max_distance_jump_ratio: float,
) -> bool:
    prev = float(previous_distance_m)
    cur = float(current_distance_m)
    if not (math.isfinite(prev) and math.isfinite(cur)):
        return True
    if prev <= 0.0 or cur <= 0.0:
        return True
    jump_m = abs(cur - prev)
    smaller = min(prev, cur)
    if smaller > 0.0 and (jump_m / smaller) > float(max_distance_jump_ratio):
        return True
    return False


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
    layers = compute_horizon_layers(
        dem_grid=dem_grid,
        geod=geod,
        observer=observer,
        azimuth_step_deg=azimuth_step_deg,
        distance_samples_m=distance_samples_m,
        dem_resampling=dem_resampling,
        earth_radius_m=earth_radius_m,
        refraction_coefficient=refraction_coefficient,
    )
    return layers.main_profile


def compute_horizon_layers(
    *,
    dem_grid,
    geod,
    observer: ObserverLocation,
    azimuth_step_deg: float,
    distance_samples_m: np.ndarray,
    dem_resampling: str,
    earth_radius_m: float,
    refraction_coefficient: float,
    max_secondary_peaks: int = 2,
    min_prominence_deg: float = 0.08,
    min_drop_deg: float = 0.05,
    min_separation_m: float = 1_200.0,
    max_distance_jump_ratio: float = 0.25,
) -> HorizonLayerSet:
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

    active_secondary_layers: list[list[HorizonProfilePoint] | None] = [None for _ in range(max(0, int(max_secondary_peaks)))]
    secondary_layers: list[list[HorizonProfilePoint]] = []
    for row_index, azimuth_deg in enumerate(azimuths):
        if not valid_rows[row_index]:
            continue
        secondary_indices = _select_secondary_peak_indices(
            altitude_deg[row_index],
            distance_grid_m[row_index],
            main_peak_index=int(peak_indices[row_index]),
            min_prominence_deg=float(min_prominence_deg),
            min_drop_deg=float(min_drop_deg),
            min_separation_m=float(min_separation_m),
            max_secondary_peaks=max(0, int(max_secondary_peaks)),
        )
        secondary_indices = _prune_secondary_peak_indices_by_visibility(
            altitude_deg[row_index],
            secondary_indices,
        )
        for layer_index, peak_index in enumerate(secondary_indices):
            if layer_index >= len(active_secondary_layers):
                break
            point = HorizonProfilePoint(
                azimuth_deg=float(azimuth_deg),
                altitude_deg=float(altitude_deg[row_index, peak_index]),
                distance_m=float(distance_grid_m[row_index, peak_index]),
                latitude_deg=float(ray_lat_deg[row_index, peak_index]),
                longitude_deg=float(ray_lon_deg[row_index, peak_index]),
                terrain_elevation_m=float(terrain_m[row_index, peak_index]),
            )
            current_fragment = active_secondary_layers[layer_index]
            if (
                current_fragment is None
                or _should_break_secondary_ridge(
                    current_fragment[-1].distance_m,
                    point.distance_m,
                    max_distance_jump_ratio=float(max_distance_jump_ratio),
                )
            ):
                if current_fragment is not None and len(current_fragment) >= 2:
                    secondary_layers.append(current_fragment)
                current_fragment = [point]
                active_secondary_layers[layer_index] = current_fragment
            else:
                current_fragment.append(point)

    for fragment in active_secondary_layers:
        if fragment is not None and len(fragment) >= 2:
            secondary_layers.append(fragment)
    return HorizonLayerSet(main_profile=points, secondary_layers=secondary_layers)


def reduce_profile_to_altaz(points: Sequence[HorizonProfilePoint]) -> list[tuple[float, float]]:
    return [(float(point.altitude_deg), float(point.azimuth_deg)) for point in points]
