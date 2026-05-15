from __future__ import annotations

import logging
import threading
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

from PySide6.QtCore import QObject, Signal

from ..types import ViewerData
from ..clouddisc.types import DownloadCancelledError
from ..water_overlay import (
    DEFAULT_WATER_RADIUS_KM,
    DEFAULT_WATER_SAMPLE_STEP_M,
    DEFAULT_WATER_AZIMUTH_STEP_DEG,
    resolve_water_scan_radius_km,
)
from ..water_mask_interface import sample_water_surface_interface_points
from ..terrain import GeoTiffDem
from ..terrain import build_download_bbox
from ..terrain import fetch_copernicus_dem
from ..terrain import sample_ground_elevation
from ..paths import CACHE_PATH
from .water_overlay_cache import (
    water_overlay_cache_scope_key,
)

logger = logging.getLogger(__name__)


@dataclass(slots=True)
class _WaterOverlayScopeCache:
    footprints: tuple
    fetched_at_utc: datetime | None
    sea_points: tuple | None = None
    dem_points: tuple | None = None
    dem_ground_m: float | None = None


class WaterOverlayController(QObject):
    water_started = Signal(object)
    water_ready = Signal(object)
    water_failed = Signal(object)

    def __init__(
        self,
        *,
        radius_km: float = DEFAULT_WATER_RADIUS_KM,
        sample_step_m: float = DEFAULT_WATER_SAMPLE_STEP_M,
        azimuth_step_deg: float = DEFAULT_WATER_AZIMUTH_STEP_DEG,
        endpoint: str | None = None,
        user_agent: str | None = None,
        timeout_s: float = 60.0,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._radius_km = float(radius_km)
        self._sample_step_m = float(sample_step_m)
        self._azimuth_step_deg = float(azimuth_step_deg)
        self._endpoint = endpoint
        self._user_agent = user_agent
        self._timeout_s = float(timeout_s)
        self._dem_cache_dir = Path(CACHE_PATH) / "copernicus-dem"
        self._running = False
        self._stopping = False
        self._active_key: Optional[tuple[float, float, float, float, bool]] = None
        self._completed_key: Optional[tuple[float, float, float, float, bool]] = None
        self._failed_key: Optional[tuple[float, float, float, float, bool]] = None
        self._pending_request: tuple[ViewerData, float, bool, str] | None = None
        self._scope_cache: dict[str, _WaterOverlayScopeCache] = {}
        self._active_workers: set[threading.Thread] = set()
        self._download_abort_event = threading.Event()
        self._lock = threading.Lock()

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        with self._lock:
            self._stopping = True
            self._download_abort_event.set()
        self._wait_for_workers(wait_timeout_s)

    def update(
        self,
        *,
        viewer_data: ViewerData,
        observer_ground_m: float,
        use_dem_ground: bool,
        reason: str = "manual",
    ) -> bool:
        observer_absolute_height_m = float(viewer_data.observer_height_m) + float(observer_ground_m)
        scan_radius_km = resolve_water_scan_radius_km(
            observer_absolute_height_m,
            minimum_distance_km=self._radius_km,
        )
        key = (
            float(viewer_data.lat_deg),
            float(viewer_data.lon_deg),
            float(viewer_data.observer_height_m),
            float(observer_ground_m),
            bool(use_dem_ground),
        )
        scope_key = water_overlay_cache_scope_key(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            radius_km=scan_radius_km,
        )
        with self._lock:
            in_memory_scope = self._scope_cache.get(scope_key)
        cached_scope = in_memory_scope
        if cached_scope is not None:
            cached_variant = self._select_cached_variant(
                cached_scope,
                use_dem_ground=bool(use_dem_ground),
                observer_ground_m=float(observer_ground_m),
            )
            if cached_variant is not None:
                self._emit_variant(
                    cached_variant["points"],
                    mode=cached_variant["mode"],
                    sea_points=cached_scope.sea_points,
                    dem_points=cached_scope.dem_points,
                    water_polygon_count=len(cached_scope.footprints),
                    source="Water: cache",
                )
                with self._lock:
                    self._active_key = key
                    self._completed_key = key
                    self._failed_key = None
                    self._pending_request = None
                return False
        with self._lock:
            if self._stopping or self._running:
                if self._running and self._completed_key != key and self._failed_key != key:
                    self._pending_request = (
                        viewer_data,
                        float(observer_ground_m),
                        bool(use_dem_ground),
                        reason,
                    )
                return False
            if self._completed_key == key:
                return False
            if self._failed_key == key:
                return False
            self._running = True
            self._active_key = key

        self.water_started.emit({"banner": "Water: loading..."})
        try:
            self._spawn_worker(
                target=self._run_update,
                kwargs={
                    "lat_deg": float(viewer_data.lat_deg),
                    "lon_deg": float(viewer_data.lon_deg),
                    "observer_height_m": float(viewer_data.observer_height_m),
                    "observer_ground_m": float(observer_ground_m),
                    "use_dem_ground": bool(use_dem_ground),
                    "reason": reason,
                    "key": key,
                    "scope_key": scope_key,
                    "scan_radius_km": scan_radius_km,
                    "cached_scope": cached_scope,
                },
                label="water",
            )
        except Exception:
            with self._lock:
                self._running = False
            raise
        return True

    def _spawn_worker(
        self,
        *,
        target: Callable[..., None],
        kwargs: dict[str, object],
        label: str,
    ) -> None:
        def runner() -> None:
            try:
                target(**kwargs)
            finally:
                self._unregister_worker(threading.current_thread())

        worker = threading.Thread(target=runner, name=f"WaterOverlayController-{label}", daemon=True)
        with self._lock:
            if self._stopping:
                return
            self._active_workers.add(worker)
        try:
            worker.start()
        except Exception:
            with self._lock:
                self._active_workers.discard(worker)
            raise

    def _unregister_worker(self, worker: threading.Thread) -> None:
        with self._lock:
            self._active_workers.discard(worker)

    def _wait_for_workers(self, wait_timeout_s: float | None) -> None:
        deadline = None if wait_timeout_s is None else time.monotonic() + max(0.0, float(wait_timeout_s))
        while True:
            with self._lock:
                workers = tuple(self._active_workers)
            if not workers:
                return
            if deadline is None:
                for worker in workers:
                    worker.join()
                continue
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                logger.warning(
                    "Timed out waiting for %d water-overlay worker thread(s) to finish during shutdown",
                    len(workers),
                )
                return
            for worker in workers:
                worker.join(timeout=remaining)

    def _run_update(
        self,
        *,
        lat_deg: float,
        lon_deg: float,
        observer_height_m: float,
        observer_ground_m: float,
        use_dem_ground: bool,
        reason: str,
        key: tuple[float, float, float, float, bool],
        scope_key: str,
        scan_radius_km: float,
        cached_scope: _WaterOverlayScopeCache | None,
    ) -> None:
        try:
            if reason == "initial":
                logger.info("Fetching initial water surface mask...")
            else:
                logger.info("Fetching water surface mask...")

            scope_cache = self._ensure_scope_cache(
                scope_key=scope_key,
                lat_deg=float(lat_deg),
                lon_deg=float(lon_deg),
                scan_radius_km=scan_radius_km,
                cached_scope=cached_scope,
                now_utc=datetime.now(timezone.utc),
            )
            active_points, sea_points, dem_points = self._build_requested_variants(
                scope_cache,
                observer_lat_deg=float(lat_deg),
                observer_lon_deg=float(lon_deg),
                observer_height_m=float(observer_height_m),
                observer_ground_m=float(observer_ground_m),
                use_dem_ground=bool(use_dem_ground),
                scan_radius_km=scan_radius_km,
            )
            if active_points:
                nearest_distance_km = min(float(point.distance_km) for point in active_points)
                logger.info(
                    "Water mask points: %d visible, nearest sea point %.3f km",
                    len(active_points),
                    nearest_distance_km,
                )
            mode = "sea"
            self._store_scope_cache(
                scope_key,
                scope_cache,
                sea_points=sea_points,
                dem_points=dem_points,
                dem_ground_m=None,
            )
            with self._lock:
                should_emit = (not self._stopping) and self._active_key == key
                if should_emit:
                    self._completed_key = key
                    self._failed_key = None
            if should_emit:
                self._emit_variant(
                    active_points,
                    mode=mode,
                    sea_points=sea_points,
                    dem_points=dem_points,
                    water_polygon_count=len(scope_cache.footprints),
                    source="Water: sea mask",
                )
        except DownloadCancelledError:
            return
        except Exception as exc:
            logger.warning("Water surface update failed: %s", exc)
            if cached_scope is not None:
                fallback_variant = self._select_cached_variant(
                    cached_scope,
                    use_dem_ground=bool(use_dem_ground),
                    observer_ground_m=float(observer_ground_m),
                )
                if fallback_variant is not None and self._active_key == key:
                    self._emit_variant(
                        fallback_variant["points"],
                        mode=fallback_variant["mode"],
                        sea_points=cached_scope.sea_points,
                        dem_points=cached_scope.dem_points,
                        water_polygon_count=len(cached_scope.footprints),
                        source="Water: cache-stale",
                    )
                with self._lock:
                    if self._active_key == key:
                        self._completed_key = key
                        self._failed_key = None
                return
            with self._lock:
                should_emit = (not self._stopping) and self._active_key == key
                if should_emit:
                    self._failed_key = key
            if should_emit:
                self.water_failed.emit({"banner": f"Water: {exc}"})
        finally:
            with self._lock:
                self._running = False
                pending_request = self._pending_request
                self._pending_request = None
            if pending_request is not None and not self._stopping:
                pending_viewer_data, pending_ground_m, pending_use_dem, pending_reason = pending_request
                self.update(
                    viewer_data=pending_viewer_data,
                    observer_ground_m=pending_ground_m,
                    use_dem_ground=pending_use_dem,
                    reason=pending_reason,
                )

    def _ensure_scope_cache(
        self,
        *,
        scope_key: str,
        lat_deg: float,
        lon_deg: float,
        scan_radius_km: float,
        cached_scope: _WaterOverlayScopeCache | None,
        now_utc: datetime,
    ) -> _WaterOverlayScopeCache:
        with self._lock:
            cached = self._scope_cache.get(scope_key)
            if cached is not None:
                return cached
        try:
            cache = _WaterOverlayScopeCache(
                footprints=(),
                fetched_at_utc=now_utc,
                sea_points=None,
                dem_points=None,
                dem_ground_m=None,
            )
            with self._lock:
                self._scope_cache[scope_key] = cache
            return cache
        except DownloadCancelledError:
            raise
        except Exception:
            raise

    def _build_requested_variants(
        self,
        scope_cache: _WaterOverlayScopeCache,
        *,
        observer_lat_deg: float,
        observer_lon_deg: float,
        observer_height_m: float,
        observer_ground_m: float,
        use_dem_ground: bool,
        scan_radius_km: float,
    ) -> tuple[tuple, tuple | None, tuple | None]:
        sea_points = scope_cache.sea_points
        if sea_points is None:
            sea_points = sample_water_surface_interface_points(
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                observer_height_m=float(observer_height_m) + float(observer_ground_m),
                max_distance_km=scan_radius_km,
            )
        active_points = sea_points
        dem_points = sea_points if use_dem_ground else None
        return active_points, sea_points, dem_points

    def _build_target_ground_sampler(
        self,
        *,
        observer_lat_deg: float,
        observer_lon_deg: float,
        scan_radius_km: float,
    ) -> Callable[[float, float], float] | None:
        try:
            download = fetch_copernicus_dem(
                observer_lat_deg=float(observer_lat_deg),
                observer_lon_deg=float(observer_lon_deg),
                max_distance_km=scan_radius_km,
                margin_km=10.0,
                cache_dir=self._dem_cache_dir,
                abort_event=self._download_abort_event,
            )
            dem = GeoTiffDem(download.paths, default_elevation_m=0.0)
        except DownloadCancelledError:
            raise
        except Exception:
            return None

        try:
            bbox = build_download_bbox(
                lat_deg=float(observer_lat_deg),
                lon_deg=float(observer_lon_deg),
                radius_km=scan_radius_km + 10.0,
            )
            dem_grid = dem.build_grid(bbox)
        except DownloadCancelledError:
            raise
        except Exception:
            dem.close()
            return None

        dem.close()

        def sampler(latitude_deg: float, longitude_deg: float) -> float:
            return sample_ground_elevation(
                dem_grid,
                latitude_deg=float(latitude_deg),
                longitude_deg=float(longitude_deg),
                dem_resampling="bilinear",
            )

        return sampler

    def _store_scope_cache(
        self,
        scope_key: str,
        scope_cache: _WaterOverlayScopeCache,
        *,
        sea_points: tuple | None,
        dem_points: tuple | None,
        dem_ground_m: float | None,
    ) -> None:
        with self._lock:
            if sea_points is not None:
                scope_cache.sea_points = sea_points
            if dem_points is not None:
                scope_cache.dem_points = dem_points
                scope_cache.dem_ground_m = dem_ground_m
            self._scope_cache[scope_key] = scope_cache

    def _select_cached_variant(
        self,
        scope_cache: _WaterOverlayScopeCache,
        *,
        use_dem_ground: bool,
        observer_ground_m: float,
    ) -> dict[str, object] | None:
        if use_dem_ground:
            if scope_cache.dem_points is not None and (
                scope_cache.dem_ground_m is None
                or abs(float(scope_cache.dem_ground_m) - float(observer_ground_m)) < 1e-6
            ):
                return {"mode": "dem", "points": scope_cache.dem_points}
            return None
        if scope_cache.sea_points is not None:
            return {"mode": "sea", "points": scope_cache.sea_points}
        return None

    def _emit_variant(
        self,
        points: tuple,
        *,
        mode: str,
        sea_points: tuple | None,
        dem_points: tuple | None,
        water_polygon_count: int,
        source: str,
    ) -> None:
        payload = {
            "points": list(points),
            "mode": mode,
            "source": source,
            "water_polygon_count": int(water_polygon_count),
            "water_point_count": len(points),
        }
        if sea_points is not None:
            payload["sea_points"] = list(sea_points)
        if dem_points is not None:
            payload["dem_points"] = list(dem_points)
        self.water_ready.emit(payload)
