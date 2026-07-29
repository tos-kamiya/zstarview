from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

EARTH_MEAN_RADIUS_M = 6_371_008.8
DEFAULT_TERRAIN_DISTANCE_SAMPLE_STEP_M = 60.0


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
    sample_azimuths_deg: np.ndarray | None = None
    sample_distances_m: np.ndarray | None = None
    sample_terrain_elevation_m: np.ndarray | None = None


@dataclass(frozen=True)
class RayScanGrid:
    azimuths_deg: np.ndarray
    distance_grid_m: np.ndarray
    ray_lon_deg: np.ndarray
    ray_lat_deg: np.ndarray


DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0)


def build_distance_samples(max_distance_km: float, sample_step_m: float) -> np.ndarray:
    max_distance_m = max_distance_km * 1000.0
    if max_distance_m <= 0.0:
        raise ValueError("max_distance_km must be positive.")
    if sample_step_m <= 0.0:
        raise ValueError("sample_step_m must be positive.")
    count = max(1, int(math.ceil(max_distance_m / sample_step_m)))
    return np.linspace(sample_step_m, max_distance_m, num=count, dtype=np.float64)


def build_ray_scan_grid(
    *,
    geod,
    observer_latitude_deg: float,
    observer_longitude_deg: float,
    azimuth_step_deg: float,
    distance_samples_m: np.ndarray,
) -> RayScanGrid:
    if azimuth_step_deg <= 0.0:
        raise ValueError("azimuth_step_deg must be positive.")

    azimuths_deg = np.arange(0.0, 360.0, azimuth_step_deg, dtype=np.float64)
    azimuth_grid_deg, distance_grid_m = np.meshgrid(
        azimuths_deg,
        np.asarray(distance_samples_m, dtype=np.float64),
        indexing="ij",
    )
    flat_count = azimuth_grid_deg.size
    ray_lon_flat, ray_lat_flat, _ = geod.fwd(
        np.full(flat_count, float(observer_longitude_deg), dtype=np.float64),
        np.full(flat_count, float(observer_latitude_deg), dtype=np.float64),
        azimuth_grid_deg.ravel(),
        distance_grid_m.ravel(),
    )
    ray_lon_deg = np.asarray(ray_lon_flat, dtype=np.float64).reshape(azimuth_grid_deg.shape)
    ray_lat_deg = np.asarray(ray_lat_flat, dtype=np.float64).reshape(azimuth_grid_deg.shape)
    return RayScanGrid(
        azimuths_deg=azimuths_deg,
        distance_grid_m=distance_grid_m,
        ray_lon_deg=ray_lon_deg,
        ray_lat_deg=ray_lat_deg,
    )


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


def _select_distance_band_peak_index(
    altitude_row_deg: np.ndarray,
    distance_row_m: np.ndarray,
    *,
    band_min_m: float,
    band_max_m: float,
    include_upper_bound: bool,
) -> int | None:
    values = np.asarray(altitude_row_deg, dtype=np.float64)
    distances = np.asarray(distance_row_m, dtype=np.float64)
    if values.size == 0 or values.size != distances.size:
        return None
    lower = float(band_min_m)
    upper = float(band_max_m)
    if not (math.isfinite(lower) and math.isfinite(upper)) or upper <= lower:
        return None

    lower_mask = distances >= lower
    upper_mask = distances <= upper if include_upper_bound else distances < upper
    mask = lower_mask & upper_mask & np.isfinite(values) & np.isfinite(distances)
    if not np.any(mask):
        return None

    masked_values = np.where(mask, values, -np.inf)
    peak_index = int(np.argmax(masked_values))
    if not math.isfinite(float(masked_values[peak_index])):
        return None
    return peak_index


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


def _prune_secondary_peak_indices_by_main_profile(
    altitude_row_deg: np.ndarray,
    distance_row_m: np.ndarray,
    secondary_peak_indices: Sequence[int],
    *,
    main_peak_distance_m: float,
) -> list[int]:
    values = np.asarray(altitude_row_deg, dtype=np.float64)
    distances = np.asarray(distance_row_m, dtype=np.float64)
    filtered_indices: list[int] = []
    main_distance_m = float(main_peak_distance_m)
    for index in sorted(int(value) for value in secondary_peak_indices):
        if index < 0 or index >= values.size:
            continue
        distance_m = float(distances[index])
        if distance_m >= main_distance_m - 1.0e-9:
            continue
        filtered_indices.append(index)
    return filtered_indices


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
    distance_band_edges_km: Sequence[float] = DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM,
) -> HorizonLayerSet:
    ray_scan = build_ray_scan_grid(
        geod=geod,
        observer_latitude_deg=observer.latitude_deg,
        observer_longitude_deg=observer.longitude_deg,
        azimuth_step_deg=azimuth_step_deg,
        distance_samples_m=distance_samples_m,
    )
    terrain_m = dem_grid.sample_lonlat(
        ray_scan.ray_lon_deg,
        ray_scan.ray_lat_deg,
        method=dem_resampling,
    )
    altitude_deg = compute_apparent_altitudes(
        observer_elevation_m=observer.observer_elevation_m,
        target_elevation_m=terrain_m,
        surface_distance_m=ray_scan.distance_grid_m,
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
    for row_index, azimuth_deg in enumerate(ray_scan.azimuths_deg):
        if not valid_rows[row_index]:
            continue
        peak_index = int(peak_indices[row_index])
        points.append(
            HorizonProfilePoint(
                azimuth_deg=float(azimuth_deg),
                altitude_deg=float(altitude_deg[row_index, peak_index]),
                distance_m=float(ray_scan.distance_grid_m[row_index, peak_index]),
                latitude_deg=float(ray_scan.ray_lat_deg[row_index, peak_index]),
                longitude_deg=float(ray_scan.ray_lon_deg[row_index, peak_index]),
                terrain_elevation_m=float(terrain_m[row_index, peak_index]),
            )
        )

    if not points:
        raise ValueError("The DEM did not provide any valid samples for the scan.")

    band_edges_km = sorted(
        {
            float(edge)
            for edge in distance_band_edges_km
            if math.isfinite(float(edge)) and float(edge) > 0.0
        }
    )
    band_limits_m: list[tuple[float, float]] = []
    band_start_m = 0.0
    for band_edge_km in band_edges_km:
        band_end_m = float(band_edge_km) * 1000.0
        if band_end_m <= band_start_m:
            continue
        band_limits_m.append((band_start_m, band_end_m))
        band_start_m = band_end_m

    band_layers: list[list[HorizonProfilePoint]] = [[] for _ in band_limits_m]
    for row_index, azimuth_deg in enumerate(ray_scan.azimuths_deg):
        if not valid_rows[row_index]:
            continue
        for band_index, (band_min_m, band_max_m) in enumerate(band_limits_m):
            peak_index = _select_distance_band_peak_index(
                altitude_deg[row_index],
                ray_scan.distance_grid_m[row_index],
                band_min_m=band_min_m,
                band_max_m=band_max_m,
                include_upper_bound=band_index == len(band_limits_m) - 1,
            )
            if peak_index is None:
                continue
            band_layers[band_index].append(
                HorizonProfilePoint(
                    azimuth_deg=float(azimuth_deg),
                    altitude_deg=float(altitude_deg[row_index, peak_index]),
                    distance_m=float(ray_scan.distance_grid_m[row_index, peak_index]),
                    latitude_deg=float(ray_scan.ray_lat_deg[row_index, peak_index]),
                    longitude_deg=float(ray_scan.ray_lon_deg[row_index, peak_index]),
                    terrain_elevation_m=float(terrain_m[row_index, peak_index]),
                )
            )

    secondary_layers = [layer for layer in band_layers if len(layer) >= 2]
    return HorizonLayerSet(
        main_profile=points,
        secondary_layers=secondary_layers,
        sample_azimuths_deg=np.asarray(ray_scan.azimuths_deg, dtype=np.float64),
        sample_distances_m=np.asarray(ray_scan.distance_grid_m, dtype=np.float64),
        sample_terrain_elevation_m=np.asarray(terrain_m, dtype=np.float64),
    )


def compute_flat_ground_horizon_layers(
    *,
    geod,
    observer: ObserverLocation,
    azimuth_step_deg: float,
    distance_samples_m: np.ndarray,
    earth_radius_m: float,
    refraction_coefficient: float,
    distance_band_edges_km: Sequence[float] = DEFAULT_TERRAIN_DISTANCE_BAND_EDGES_KM,
) -> HorizonLayerSet:
    class _FlatDemGrid:
        def sample_lonlat(
            self,
            lon_deg: np.ndarray,
            lat_deg: np.ndarray,
            *,
            method: str = "bilinear",
        ) -> np.ndarray:
            return np.zeros_like(lon_deg, dtype=np.float64)

    return compute_horizon_layers(
        dem_grid=_FlatDemGrid(),
        geod=geod,
        observer=observer,
        azimuth_step_deg=azimuth_step_deg,
        distance_samples_m=distance_samples_m,
        dem_resampling="bilinear",
        earth_radius_m=earth_radius_m,
        refraction_coefficient=refraction_coefficient,
        distance_band_edges_km=distance_band_edges_km,
    )


def compute_sea_level_horizon_layers(
    *,
    geod,
    observer: ObserverLocation,
    azimuth_step_deg: float,
    earth_radius_m: float,
    refraction_coefficient: float,
) -> HorizonLayerSet:
    if azimuth_step_deg <= 0.0:
        raise ValueError("azimuth_step_deg must be positive.")
    if earth_radius_m <= 0.0:
        raise ValueError("earth_radius_m must be positive.")
    if refraction_coefficient >= 1.0:
        raise ValueError("refraction_coefficient must be < 1.0.")

    effective_radius_m = earth_radius_m / max(1.0 - refraction_coefficient, 1e-6)
    observer_elevation_m = max(0.0, float(observer.observer_elevation_m))
    horizon_distance_m = math.sqrt(max(0.0, observer_elevation_m * (2.0 * effective_radius_m + observer_elevation_m)))
    azimuths_deg = np.arange(0.0, 360.0, float(azimuth_step_deg), dtype=np.float64)

    points: list[HorizonProfilePoint] = []
    for azimuth_deg in azimuths_deg:
        lon_deg, lat_deg, _ = geod.fwd(
            float(observer.longitude_deg),
            float(observer.latitude_deg),
            float(azimuth_deg),
            float(horizon_distance_m),
        )
        points.append(
            HorizonProfilePoint(
                azimuth_deg=float(azimuth_deg),
                altitude_deg=0.0,
                distance_m=float(horizon_distance_m),
                latitude_deg=float(lat_deg),
                longitude_deg=float(lon_deg),
                terrain_elevation_m=0.0,
            )
        )
    return HorizonLayerSet(main_profile=points, secondary_layers=[])


def reduce_profile_to_altaz(points: Sequence[HorizonProfilePoint]) -> list[tuple[float, float]]:
    return [(float(point.altitude_deg), float(point.azimuth_deg)) for point in points]
