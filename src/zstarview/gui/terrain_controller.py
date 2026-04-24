# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
import threading
from pathlib import Path
from typing import Optional

from PySide6.QtCore import QObject, Signal

from ..terrain import (
    EARTH_MEAN_RADIUS_M,
    COPERNICUS_DEM_BUCKET,
    GeoTiffDem,
    ObserverLocation,
    WGS84_GEOD,
    build_distance_samples,
    build_download_bbox,
    compute_horizon_profile,
    fetch_copernicus_dem,
    reduce_profile_to_altaz,
    sample_ground_elevation,
)

logger = logging.getLogger(__name__)


class TerrainHorizonController(QObject):
    terrain_started = Signal(object)
    terrain_ready = Signal(object)
    terrain_failed = Signal(object)

    def __init__(
        self,
        *,
        cache_dir: Path,
        observer_eye_m: float = 1.7,
        max_distance_km: float = 120.0,
        download_margin_km: float = 10.0,
        sample_step_m: float = 90.0,
        azimuth_step_deg: float = 1.0,
        dem_resampling: str = "bilinear",
        earth_radius_m: float = EARTH_MEAN_RADIUS_M,
        refraction_coefficient: float = 0.13,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._cache_dir = Path(cache_dir)
        self._observer_eye_m = float(observer_eye_m)
        self._max_distance_km = float(max_distance_km)
        self._download_margin_km = float(download_margin_km)
        self._sample_step_m = float(sample_step_m)
        self._azimuth_step_deg = float(azimuth_step_deg)
        self._dem_resampling = str(dem_resampling)
        self._earth_radius_m = float(earth_radius_m)
        self._refraction_coefficient = float(refraction_coefficient)
        self._running = False
        self._stopping = False
        self._failed_this_session = False
        self._completed_for_location: Optional[tuple[float, float, float]] = None
        self._lock = threading.Lock()

    def shutdown(self) -> None:
        with self._lock:
            self._stopping = True

    def update(
        self,
        *,
        lat: float,
        lon: float,
        observer_height_m: float | None = None,
        reason: str = "manual",
    ) -> bool:
        eye_height_m = self._observer_eye_m if observer_height_m is None else float(observer_height_m)
        location = (float(lat), float(lon), eye_height_m)
        with self._lock:
            if self._stopping or self._running or self._failed_this_session:
                return False
            if self._completed_for_location == location:
                return False
            self._running = True

        self.terrain_started.emit({"banner": "Terrain horizon: loading DEM..."})
        worker = threading.Thread(
            target=self._run_update,
            kwargs={"lat": location[0], "lon": location[1], "observer_height_m": eye_height_m, "reason": reason},
            daemon=True,
        )
        worker.start()
        return True

    def _run_update(self, *, lat: float, lon: float, observer_height_m: float, reason: str) -> None:
        try:
            if reason == "initial":
                logger.info("Fetching initial terrain horizon data...")
            else:
                logger.info("Fetching terrain horizon data...")

            try:
                download = fetch_copernicus_dem(
                    observer_lat_deg=lat,
                    observer_lon_deg=lon,
                    max_distance_km=self._max_distance_km,
                    margin_km=self._download_margin_km,
                    cache_dir=self._cache_dir,
                )
            except RuntimeError as exc:
                if str(exc) != "No Copernicus DEM tiles were downloaded for the requested area.":
                    raise
                logger.info(
                    "No Copernicus DEM tiles available for observer area; using sea-level horizon only."
                )
                with self._lock:
                    if not self._stopping:
                        self._completed_for_location = (float(lat), float(lon), float(observer_height_m))
                    should_emit = not self._stopping
                if should_emit:
                    self.terrain_ready.emit(
                        {
                            "profile_altaz": [],
                            "source": f"{COPERNICUS_DEM_BUCKET}:ocean",
                        }
                    )
                return
            dem = GeoTiffDem(download.paths, default_elevation_m=0.0)
            try:
                bbox = build_download_bbox(
                    lat_deg=lat,
                    lon_deg=lon,
                    radius_km=self._max_distance_km + self._download_margin_km,
                )
                dem_grid = dem.build_grid(bbox)
                ground_m = sample_ground_elevation(
                    dem_grid,
                    latitude_deg=lat,
                    longitude_deg=lon,
                    dem_resampling=self._dem_resampling,
                )
                observer = ObserverLocation(
                    latitude_deg=lat,
                    longitude_deg=lon,
                    observer_ground_m=ground_m,
                    observer_eye_m=observer_height_m,
                )
                points = compute_horizon_profile(
                    dem_grid=dem_grid,
                    geod=WGS84_GEOD,
                    observer=observer,
                    azimuth_step_deg=self._azimuth_step_deg,
                    distance_samples_m=build_distance_samples(
                        self._max_distance_km,
                        self._sample_step_m,
                    ),
                    dem_resampling=self._dem_resampling,
                    earth_radius_m=self._earth_radius_m,
                    refraction_coefficient=self._refraction_coefficient,
                )
            finally:
                dem.close()

            profile_altaz = reduce_profile_to_altaz(points)
            with self._lock:
                if not self._stopping:
                    self._completed_for_location = (float(lat), float(lon), float(observer_height_m))

            with self._lock:
                should_emit = not self._stopping
            if should_emit:
                self.terrain_ready.emit(
                    {
                        "profile_altaz": profile_altaz,
                        "source": f"{COPERNICUS_DEM_BUCKET}:{download.source}",
                    }
                )
        except Exception as exc:
            logger.warning("Terrain horizon update failed: %s", exc, exc_info=True)
            with self._lock:
                self._failed_this_session = True
                should_emit = not self._stopping
            if should_emit:
                self.terrain_failed.emit({"banner": "Terrain horizon: unavailable"})
        finally:
            with self._lock:
                self._running = False
