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
from ..water_overlay import (
    DEFAULT_WATER_RADIUS_KM,
    DEFAULT_WATER_OVERPASS_ENDPOINT,
    DEFAULT_WATER_USER_AGENT,
    DEFAULT_WATER_SAMPLE_STEP_M,
    DEFAULT_WATER_AZIMUTH_STEP_DEG,
    bbox_from_point,
    extract_water_polygons,
    fetch_overpass_json,
    sample_water_overlay_points,
)
from ..terrain import GeoTiffDem
from ..terrain import build_download_bbox
from ..terrain import fetch_copernicus_dem
from ..terrain import sample_ground_elevation
from ..paths import CACHE_PATH
from .water_overlay_cache import (
    WaterOverlayCacheSnapshot,
    WATER_OVERLAY_CACHE_RETENTION_SECONDS,
    load_water_overlay_cache,
    save_water_overlay_cache,
    water_overlay_cache_is_recent,
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
        self._cache_retention_seconds = int(WATER_OVERLAY_CACHE_RETENTION_SECONDS)
        self._running = False
        self._stopping = False
        self._active_key: Optional[tuple[float, float, float, float, bool]] = None
        self._completed_key: Optional[tuple[float, float, float, float, bool]] = None
        self._failed_key: Optional[tuple[float, float, float, float, bool]] = None
        self._pending_request: tuple[ViewerData, float, bool, str] | None = None
        self._scope_cache: dict[str, _WaterOverlayScopeCache] = {}
        self._active_workers: set[threading.Thread] = set()
        self._lock = threading.Lock()

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        with self._lock:
            self._stopping = True
        self._wait_for_workers(wait_timeout_s)

    def update(
        self,
        *,
        viewer_data: ViewerData,
        observer_ground_m: float,
        use_dem_ground: bool,
        reason: str = "manual",
    ) -> bool:
        now = datetime.now(timezone.utc)
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
            radius_km=self._radius_km,
        )
        with self._lock:
            in_memory_scope = self._scope_cache.get(scope_key)
        cached_scope = in_memory_scope
        snapshot = self._load_scope_snapshot(scope_key, now=now)
        if snapshot is not None:
            snapshot_scope = self._scope_cache_from_snapshot(snapshot)
            if (
                cached_scope is None
                or cached_scope.fetched_at_utc is None
                or (
                    snapshot_scope.fetched_at_utc is not None
                    and cached_scope.fetched_at_utc < snapshot_scope.fetched_at_utc
                )
            ):
                cached_scope = snapshot_scope
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
        cached_scope: _WaterOverlayScopeCache | None,
    ) -> None:
        try:
            if reason == "initial":
                logger.info("Fetching initial water surface data...")
            else:
                logger.info("Fetching water surface data...")

            scope_cache = self._ensure_scope_cache(
                scope_key=scope_key,
                lat_deg=float(lat_deg),
                lon_deg=float(lon_deg),
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
            )
            mode = "dem" if use_dem_ground and dem_points is not None else "sea"
            self._store_scope_cache(
                scope_key,
                scope_cache,
                sea_points=sea_points,
                dem_points=dem_points,
                dem_ground_m=float(observer_ground_m) if use_dem_ground else None,
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
                    source="Water: Overpass",
                )
        except Exception as exc:
            logger.warning("Water surface update failed: %s", exc, exc_info=True)
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
                self.water_failed.emit({"banner": f"Water: failed ({exc})"})
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
        cached_scope: _WaterOverlayScopeCache | None,
        now_utc: datetime,
    ) -> _WaterOverlayScopeCache:
        with self._lock:
            cached = self._scope_cache.get(scope_key)
            if cached is not None:
                return cached
        snapshot = cached_scope if cached_scope is not None else self._load_scope_snapshot_any(scope_key)
        if snapshot is not None and not isinstance(snapshot, _WaterOverlayScopeCache):
            snapshot = self._scope_cache_from_snapshot(snapshot)
        if snapshot is not None and water_overlay_cache_is_recent(
            snapshot,
            now_utc=now_utc,
            max_age_seconds=self._cache_retention_seconds,
        ):
            cache = snapshot
            with self._lock:
                self._scope_cache[scope_key] = cache
            return cache

        if snapshot is None:
            snapshot = WaterOverlayCacheSnapshot(footprints=(), water_polygon_count=0, fetched_at_utc=None)

        try:
            bbox = bbox_from_point(float(lat_deg), float(lon_deg), self._radius_km)
            payload = fetch_overpass_json(
                bbox=bbox,
                endpoint=self._endpoint or DEFAULT_WATER_OVERPASS_ENDPOINT,
                user_agent=self._user_agent or DEFAULT_WATER_USER_AGENT,
                timeout_s=self._timeout_s,
            )
            elements = payload.get("elements")
            if not isinstance(elements, list):
                raise RuntimeError("Overpass payload missing elements")
            footprints = extract_water_polygons(elements, bbox=bbox)
            fresh_snapshot = WaterOverlayCacheSnapshot(
                footprints=footprints,
                water_polygon_count=len(footprints),
                fetched_at_utc=now_utc,
            )
            save_water_overlay_cache(scope_key, fresh_snapshot)
            cache = _WaterOverlayScopeCache(
                footprints=footprints,
                fetched_at_utc=now_utc,
            )
            with self._lock:
                self._scope_cache[scope_key] = cache
            return cache
        except Exception:
            if snapshot.footprints:
                cache = _WaterOverlayScopeCache(
                    footprints=snapshot.footprints,
                    fetched_at_utc=snapshot.fetched_at_utc,
                )
                with self._lock:
                    self._scope_cache[scope_key] = cache
                return cache
            raise

    def _load_scope_snapshot(
        self,
        scope_key: str,
        *,
        now: datetime,
    ) -> WaterOverlayCacheSnapshot | None:
        snapshot = load_water_overlay_cache(scope_key)
        if snapshot is None:
            return None
        if water_overlay_cache_is_recent(
            snapshot,
            now_utc=now,
            max_age_seconds=self._cache_retention_seconds,
        ):
            return snapshot
        return None

    def _load_scope_snapshot_any(self, scope_key: str) -> WaterOverlayCacheSnapshot | None:
        return load_water_overlay_cache(scope_key)

    def _scope_cache_from_snapshot(
        self,
        snapshot: WaterOverlayCacheSnapshot,
    ) -> _WaterOverlayScopeCache:
        return _WaterOverlayScopeCache(
            footprints=snapshot.footprints,
            fetched_at_utc=snapshot.fetched_at_utc,
        )

    def _build_requested_variants(
        self,
        scope_cache: _WaterOverlayScopeCache,
        *,
        observer_lat_deg: float,
        observer_lon_deg: float,
        observer_height_m: float,
        observer_ground_m: float,
        use_dem_ground: bool,
    ) -> tuple[tuple, tuple | None, tuple | None]:
        observer_absolute_height_m = float(observer_height_m) + float(observer_ground_m)
        sea_points = scope_cache.sea_points
        if sea_points is None:
            sea_points = sample_water_overlay_points(
                scope_cache.footprints,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                observer_height_m=observer_absolute_height_m,
                fallback_surface_height_m=observer_ground_m,
                max_distance_km=self._radius_km,
                sample_step_m=self._sample_step_m,
                azimuth_step_deg=self._azimuth_step_deg,
            )
        dem_points = scope_cache.dem_points
        if use_dem_ground and dem_points is None:
            target_ground_sampler = self._build_target_ground_sampler(
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
            )
            dem_points = sample_water_overlay_points(
                scope_cache.footprints,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                observer_height_m=observer_absolute_height_m,
                fallback_surface_height_m=observer_ground_m,
                target_ground_elevation_m_sampler=target_ground_sampler,
                max_distance_km=self._radius_km,
                sample_step_m=self._sample_step_m,
                azimuth_step_deg=self._azimuth_step_deg,
            )
        active_points = dem_points if use_dem_ground and dem_points is not None else sea_points
        return active_points, sea_points, dem_points

    def _build_target_ground_sampler(
        self,
        *,
        observer_lat_deg: float,
        observer_lon_deg: float,
    ) -> Callable[[float, float], float] | None:
        try:
            download = fetch_copernicus_dem(
                observer_lat_deg=float(observer_lat_deg),
                observer_lon_deg=float(observer_lon_deg),
                max_distance_km=self._radius_km,
                margin_km=10.0,
                cache_dir=self._dem_cache_dir,
            )
            dem = GeoTiffDem(download.paths, default_elevation_m=0.0)
        except Exception:
            return None

        try:
            bbox = build_download_bbox(
                lat_deg=float(observer_lat_deg),
                lon_deg=float(observer_lon_deg),
                radius_km=self._radius_km + 10.0,
            )
            dem_grid = dem.build_grid(bbox)
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
