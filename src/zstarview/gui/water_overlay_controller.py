from __future__ import annotations

import logging
import threading
import time
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

from PySide6.QtCore import QObject, Signal

from ..types import ViewerData
from ..clouddisc.types import DownloadCancelledError
from ..water_overlay import (
    DEFAULT_WATER_OVERPASS_ENDPOINT,
    DEFAULT_WATER_USER_AGENT,
    DEFAULT_WATER_RADIUS_KM,
    DEFAULT_WATER_SAMPLE_STEP_M,
    DEFAULT_WATER_AZIMUTH_STEP_DEG,
    bbox_from_point,
    extract_water_polygons,
    fetch_overpass_json,
    WaterOverlayPoint,
    fetch_water_overlay_footprints,
    resolve_water_scan_radius_km,
    sample_water_overlay_points,
)
from ..water_mask_interface import (
    WaterSurfaceBandStats,
    sample_water_surface_interface_points_with_stats,
)
from ..terrain import GeoTiffDem
from ..terrain import build_download_bbox
from ..terrain import fetch_copernicus_dem
from ..terrain import sample_ground_elevation
from ..paths import CACHE_PATH
from .water_overlay_cache import (
    WaterOverlayCacheSnapshot,
    WATER_OVERLAY_CACHE_RETENTION_SECONDS,
    water_overlay_cache_scope_key,
    load_water_overlay_cache,
    save_water_overlay_cache,
    water_overlay_cache_is_recent,
)

logger = logging.getLogger(__name__)


def _water_overlay_band_counts(dots: tuple[WaterOverlayPoint, ...] | list[WaterOverlayPoint]) -> tuple[int, int, int]:
    counts = Counter(str(getattr(dot, "water_category", "")).strip().lower() for dot in dots)
    return (
        int(counts.get("sea-125", 0)),
        int(counts.get("sea-250", 0)),
        int(counts.get("sea-500", 0)),
    )


def _water_overlay_band_stats_text(stats: WaterSurfaceBandStats) -> str:
    return (
        f"{stats.band_name} tiles={int(stats.loaded_tile_count)} "
        f"raw={int(stats.raw_point_count)} "
        f"collapsed={int(stats.collapsed_point_count)} "
        f"visible={int(stats.visible_point_count)}"
    )


@dataclass(slots=True)
class _WaterOverlayScopeCache:
    footprints: tuple
    fetched_at_utc: datetime | None
    uses_dem_sampler: bool = False
    sea_mask_dots: tuple | None = None
    inland_dots: tuple | None = None
    sea_dots: tuple | None = None
    sea_band_stats: tuple[WaterSurfaceBandStats, ...] | None = None
    dem_dots: tuple | None = None
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
        terrain_horizon_profile_altaz: list[tuple[float, float]] | None = None,
        terrain_horizon_profile_distances_m: list[float] | None = None,
        terrain_horizon_secondary_profile_altaz_layers: list[list[tuple[float, float]]] | None = None,
        terrain_horizon_secondary_profile_distances_m_layers: list[list[float]] | None = None,
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
                cached_sea_dots = cached_scope.sea_mask_dots or cached_scope.sea_dots
                for band_stat in cached_scope.sea_band_stats or ():
                    logger.info("Water band stats: %s", _water_overlay_band_stats_text(band_stat))
                self._emit_variant(
                    cached_variant["dots"],
                    mode=cached_variant["mode"],
                    sea_dots=cached_sea_dots,
                    inland_dots=cached_scope.inland_dots,
                    dem_dots=cached_scope.dem_dots,
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
                    "terrain_horizon_profile_altaz": terrain_horizon_profile_altaz,
                    "terrain_horizon_profile_distances_m": terrain_horizon_profile_distances_m,
                    "terrain_horizon_secondary_profile_altaz_layers": terrain_horizon_secondary_profile_altaz_layers,
                    "terrain_horizon_secondary_profile_distances_m_layers": terrain_horizon_secondary_profile_distances_m_layers,
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
        terrain_horizon_profile_altaz: list[tuple[float, float]] | None = None,
        terrain_horizon_profile_distances_m: list[float] | None = None,
        terrain_horizon_secondary_profile_altaz_layers: list[list[tuple[float, float]]] | None = None,
        terrain_horizon_secondary_profile_distances_m_layers: list[list[float]] | None = None,
    ) -> None:
        now_utc = datetime.now(timezone.utc)
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
                now_utc=now_utc,
            )
            (
                active_dots,
                sea_mask_dots,
                sea_dots,
                inland_dots,
                dem_dots,
                band_stats,
                footprints,
            ) = self._build_requested_variants(
                scope_cache,
                observer_lat_deg=float(lat_deg),
                observer_lon_deg=float(lon_deg),
                observer_height_m=float(observer_height_m),
                observer_ground_m=float(observer_ground_m),
                use_dem_ground=bool(use_dem_ground),
                scan_radius_km=scan_radius_km,
                target_ground_sampler=(
                    self._build_target_ground_sampler(
                        observer_lat_deg=float(lat_deg),
                        observer_lon_deg=float(lon_deg),
                        scan_radius_km=scan_radius_km,
                    )
                    if use_dem_ground
                    else None
                ),
                key=key,
                scope_key=scope_key,
                terrain_horizon_profile_altaz=terrain_horizon_profile_altaz,
                terrain_horizon_profile_distances_m=terrain_horizon_profile_distances_m,
                terrain_horizon_secondary_profile_altaz_layers=terrain_horizon_secondary_profile_altaz_layers,
                terrain_horizon_secondary_profile_distances_m_layers=terrain_horizon_secondary_profile_distances_m_layers,
            )
            nearest_distance_km = min((float(dot.distance_km) for dot in active_dots), default=None)
            band_100_count, band_250_count, band_500_count = _water_overlay_band_counts(active_dots)
            for band_stat in band_stats:
                logger.info("Water band stats: %s", _water_overlay_band_stats_text(band_stat))
            if nearest_distance_km is None:
                logger.info(
                    "Water mask dots: 0 visible, nearest sea dot n/a, bands: 125m=%d 250m=%d 500m=%d",
                    band_100_count,
                    band_250_count,
                    band_500_count,
                )
            else:
                logger.info(
                    "Water mask dots: %d visible, nearest sea dot %.3f km, bands: 125m=%d 250m=%d 500m=%d",
                    len(active_dots),
                    nearest_distance_km,
                    band_100_count,
                    band_250_count,
                    band_500_count,
                )
                self._store_scope_cache(
                scope_key,
                scope_cache,
                footprints=footprints,
                uses_dem_sampler=scope_cache.uses_dem_sampler,
                sea_mask_dots=sea_mask_dots,
                sea_dots=sea_dots,
                inland_dots=inland_dots,
                sea_band_stats=band_stats,
                dem_dots=dem_dots,
                dem_ground_m=float(observer_ground_m) if dem_dots is not None else None,
            )
            with self._lock:
                should_emit = (not self._stopping) and self._active_key == key
                if should_emit:
                    self._completed_key = key
                    self._failed_key = None
            if should_emit:
                self._emit_variant(
                    active_dots,
                    mode="dem" if use_dem_ground and dem_dots is not None else "sea",
                    sea_dots=sea_mask_dots,
                    inland_dots=inland_dots,
                    dem_dots=dem_dots,
                    water_polygon_count=len(scope_cache.footprints),
                    source="Water: sea mask + OSM",
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
                    cached_sea_dots = cached_scope.sea_mask_dots or cached_scope.sea_dots
                    self._emit_variant(
                        fallback_variant["dots"],
                        mode=fallback_variant["mode"],
                        sea_dots=cached_sea_dots,
                        inland_dots=cached_scope.inland_dots,
                        dem_dots=cached_scope.dem_dots,
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
        snapshot = cached_scope if cached_scope is not None else self._load_scope_snapshot(
            scope_key,
            now=now_utc,
        )
        if snapshot is not None and water_overlay_cache_is_recent(
            snapshot,
            now_utc=now_utc,
            max_age_seconds=self._cache_retention_seconds,
        ):
            if isinstance(snapshot, _WaterOverlayScopeCache):
                cache = snapshot
            else:
                cache = self._scope_cache_from_snapshot(snapshot)
            with self._lock:
                self._scope_cache[scope_key] = cache
            return cache

        if snapshot is None:
            snapshot = WaterOverlayCacheSnapshot(
                footprints=(),
                water_polygon_count=0,
                fetched_at_utc=None,
            )

        try:
            bbox = bbox_from_point(float(lat_deg), float(lon_deg), float(scan_radius_km))
            payload = fetch_overpass_json(
                bbox=bbox,
                endpoint=self._endpoint or DEFAULT_WATER_OVERPASS_ENDPOINT,
                user_agent=self._user_agent or DEFAULT_WATER_USER_AGENT,
                timeout_s=self._timeout_s,
                abort_event=self._download_abort_event,
            )
            elements = payload.get("elements")
            if not isinstance(elements, list):
                raise RuntimeError("Overpass payload missing elements")
            footprints = extract_water_polygons(elements, bbox=bbox, abort_event=self._download_abort_event)
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
        except DownloadCancelledError:
            raise
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
        scan_radius_km: float,
        target_ground_sampler: Callable[[float, float], float] | None,
        key: tuple[float, float, float, float, bool],
        scope_key: str,
        terrain_horizon_profile_altaz: list[tuple[float, float]] | None = None,
        terrain_horizon_profile_distances_m: list[float] | None = None,
        terrain_horizon_secondary_profile_altaz_layers: list[list[tuple[float, float]]] | None = None,
        terrain_horizon_secondary_profile_distances_m_layers: list[list[float]] | None = None,
    ) -> tuple[
        tuple,
        tuple | None,
        tuple | None,
        tuple | None,
        tuple | None,
        tuple[WaterSurfaceBandStats, ...],
        tuple,
    ]:
        use_target_sampler = bool(use_dem_ground and target_ground_sampler is not None)
        sea_mask_dots = scope_cache.sea_mask_dots
        if sea_mask_dots is None or bool(scope_cache.uses_dem_sampler) != use_target_sampler:
            sea_mask_dots, band_stats = sample_water_surface_interface_points_with_stats(
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                observer_height_m=float(observer_height_m) + float(observer_ground_m),
                max_distance_km=scan_radius_km,
                target_ground_elevation_m_sampler=target_ground_sampler if use_target_sampler else None,
            )
            scope_cache.sea_mask_dots = sea_mask_dots
            scope_cache.uses_dem_sampler = use_target_sampler
        else:
            band_stats = ()

        partial_sea_dots = tuple(sea_mask_dots)
        with self._lock:
            should_emit_partial = (not self._stopping) and self._active_key == key
        if should_emit_partial:
            self._store_scope_cache(
                scope_key,
                scope_cache,
                footprints=None,
                uses_dem_sampler=scope_cache.uses_dem_sampler,
                sea_mask_dots=partial_sea_dots,
                sea_dots=None,
                inland_dots=None,
                sea_band_stats=band_stats,
                dem_dots=None,
                dem_ground_m=None,
            )
            self._emit_variant(
                partial_sea_dots,
                mode="sea",
                sea_dots=partial_sea_dots,
                inland_dots=None,
                dem_dots=None,
                water_polygon_count=len(scope_cache.footprints),
                source="Water: sea mask",
            )

        footprints = scope_cache.footprints
        if not footprints:
            footprints = fetch_water_overlay_footprints(
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                max_distance_km=scan_radius_km,
                endpoint=self._endpoint or DEFAULT_WATER_OVERPASS_ENDPOINT,
                user_agent=self._user_agent or DEFAULT_WATER_USER_AGENT,
                timeout_s=self._timeout_s,
                abort_event=self._download_abort_event,
            )
            scope_cache.footprints = footprints

        observer_absolute_height_m = float(observer_height_m) + float(observer_ground_m)
        inland_dots = scope_cache.inland_dots
        if inland_dots is None or bool(scope_cache.uses_dem_sampler) != use_target_sampler:
            inland_dots = sample_water_overlay_points(
                footprints,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                observer_height_m=observer_absolute_height_m,
                fallback_surface_height_m=float(observer_ground_m),
                target_ground_elevation_m_sampler=target_ground_sampler if use_target_sampler else None,
                max_distance_km=scan_radius_km,
                abort_event=self._download_abort_event,
            )
            scope_cache.inland_dots = inland_dots
            scope_cache.uses_dem_sampler = use_target_sampler
        combined_sea_dots = tuple(sea_mask_dots) + tuple(inland_dots)

        if use_target_sampler:
            inland_dem_dots = sample_water_overlay_points(
                footprints,
                observer_lat_deg=observer_lat_deg,
                observer_lon_deg=observer_lon_deg,
                observer_height_m=observer_absolute_height_m,
                fallback_surface_height_m=float(observer_ground_m),
                target_ground_elevation_m_sampler=target_ground_sampler,
                max_distance_km=scan_radius_km,
                abort_event=self._download_abort_event,
            )
            dem_dots: tuple | None = tuple(sea_mask_dots) + tuple(inland_dem_dots)
        elif use_dem_ground:
            dem_dots = combined_sea_dots
        else:
            dem_dots = None
        active_dots = dem_dots if use_dem_ground and dem_dots is not None else combined_sea_dots
        return (
            active_dots,
            tuple(sea_mask_dots),
            combined_sea_dots,
            tuple(inland_dots),
            dem_dots,
            band_stats,
            footprints,
        )

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
        footprints: tuple | None,
        uses_dem_sampler: bool | None,
        sea_mask_dots: tuple | None,
        inland_dots: tuple | None,
        sea_dots: tuple | None,
        sea_band_stats: tuple[WaterSurfaceBandStats, ...] | None,
        dem_dots: tuple | None,
        dem_ground_m: float | None,
    ) -> None:
        with self._lock:
            if footprints is not None:
                scope_cache.footprints = footprints
            if uses_dem_sampler is not None:
                scope_cache.uses_dem_sampler = bool(uses_dem_sampler)
            if sea_mask_dots is not None:
                scope_cache.sea_mask_dots = sea_mask_dots
            if inland_dots is not None:
                scope_cache.inland_dots = inland_dots
            if sea_dots is not None:
                scope_cache.sea_dots = sea_dots
            if sea_band_stats is not None:
                scope_cache.sea_band_stats = sea_band_stats
            if dem_dots is not None:
                scope_cache.dem_dots = dem_dots
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
            if scope_cache.dem_dots is not None and (
                scope_cache.dem_ground_m is None
                or abs(float(scope_cache.dem_ground_m) - float(observer_ground_m)) < 1e-6
            ):
                return {"mode": "dem", "dots": scope_cache.dem_dots}
            return None
        if scope_cache.sea_dots is not None:
            return {"mode": "sea", "dots": scope_cache.sea_dots}
        return None

    def _emit_variant(
        self,
        dots: tuple,
        *,
        mode: str,
        sea_dots: tuple | None,
        inland_dots: tuple | None,
        dem_dots: tuple | None,
        water_polygon_count: int,
        source: str,
    ) -> None:
        payload = {
            "dots": list(dots),
            "mode": mode,
            "source": source,
            "water_polygon_count": int(water_polygon_count),
            "water_dot_count": len(dots),
        }
        if sea_dots is not None:
            payload["sea_dots"] = list(sea_dots)
        if inland_dots is not None:
            payload["inland_dots"] = list(inland_dots)
        if dem_dots is not None:
            payload["dem_dots"] = list(dem_dots)
        self.water_ready.emit(payload)
