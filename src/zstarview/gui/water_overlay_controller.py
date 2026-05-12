from __future__ import annotations

import logging
import threading
import time
from datetime import datetime, timezone
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
from .water_overlay_cache import (
    WaterOverlayCacheSnapshot,
    WATER_OVERLAY_CACHE_RETENTION_SECONDS,
    load_water_overlay_cache,
    save_water_overlay_cache,
    water_overlay_cache_is_recent,
    water_overlay_cache_scope_key,
)

logger = logging.getLogger(__name__)


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
        self._cache_retention_seconds = int(WATER_OVERLAY_CACHE_RETENTION_SECONDS)
        self._running = False
        self._stopping = False
        self._completed_key: Optional[tuple[float, float, float]] = None
        self._failed_key: Optional[tuple[float, float, float]] = None
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
        reason: str = "manual",
    ) -> bool:
        now = datetime.now(timezone.utc)
        key = (
            float(viewer_data.lat_deg),
            float(viewer_data.lon_deg),
            float(viewer_data.observer_height_m),
        )
        scope_key = water_overlay_cache_scope_key(
            observer_lat_deg=float(viewer_data.lat_deg),
            observer_lon_deg=float(viewer_data.lon_deg),
            observer_height_m=float(viewer_data.observer_height_m),
            radius_km=self._radius_km,
            sample_step_m=self._sample_step_m,
            azimuth_step_deg=self._azimuth_step_deg,
        )
        cached = load_water_overlay_cache(scope_key)
        if cached is not None:
            if water_overlay_cache_is_recent(
                cached,
                now_utc=now,
                max_age_seconds=self._cache_retention_seconds,
            ):
                self._emit_cached_result(cached)
                with self._lock:
                    self._completed_key = key
                    self._failed_key = None
                return False
        with self._lock:
            if self._stopping or self._running:
                return False
            if self._completed_key == key:
                return False
            if self._failed_key == key:
                return False
            self._running = True

        self.water_started.emit({"banner": "Water overlay: loading..."})
        try:
            self._spawn_worker(
                target=self._run_update,
                kwargs={
                    "lat_deg": float(viewer_data.lat_deg),
                    "lon_deg": float(viewer_data.lon_deg),
                    "observer_height_m": float(viewer_data.observer_height_m),
                    "reason": reason,
                    "key": key,
                    "scope_key": scope_key,
                    "cached": cached,
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
        reason: str,
        key: tuple[float, float, float],
        scope_key: str,
        cached: WaterOverlayCacheSnapshot | None,
    ) -> None:
        try:
            if reason == "initial":
                logger.info("Fetching initial water overlay data...")
            else:
                logger.info("Fetching water overlay data...")

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
            footprints = extract_water_polygons(elements)
            points = sample_water_overlay_points(
                footprints,
                observer_lat_deg=float(lat_deg),
                observer_lon_deg=float(lon_deg),
                observer_height_m=float(observer_height_m),
                max_distance_km=self._radius_km,
                sample_step_m=self._sample_step_m,
                azimuth_step_deg=self._azimuth_step_deg,
            )
            snapshot = WaterOverlayCacheSnapshot(
                points=points,
                water_polygon_count=len(footprints),
                water_point_count=len(points),
                fetched_at_utc=datetime.now(timezone.utc),
            )
            save_water_overlay_cache(scope_key, snapshot)
            with self._lock:
                if not self._stopping:
                    self._completed_key = key
                    self._failed_key = None
                should_emit = not self._stopping
            if should_emit:
                self.water_ready.emit(
                    {
                        "points": list(points),
                        "source": "Water: Overpass",
                        "water_polygon_count": len(footprints),
                        "water_point_count": len(points),
                }
            )
        except Exception as exc:
            logger.warning("Water overlay update failed: %s", exc, exc_info=True)
            if cached is not None:
                self._emit_cached_result(cached)
                with self._lock:
                    self._completed_key = key
                    self._failed_key = None
                return
            with self._lock:
                self._failed_key = key
                should_emit = not self._stopping
            if should_emit:
                self.water_failed.emit({"banner": f"Water overlay: failed ({exc})"})
        finally:
            with self._lock:
                self._running = False

    def _emit_cached_result(
        self,
        cached: WaterOverlayCacheSnapshot,
    ) -> None:
        self.water_ready.emit(
            {
                "points": list(cached.points),
                "source": "Water: cache",
                "water_polygon_count": cached.water_polygon_count,
                "water_point_count": cached.water_point_count,
            }
        )
