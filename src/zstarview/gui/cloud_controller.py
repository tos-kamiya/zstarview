"""
Cloud update controller for UI layer.

This module moves cloud fetching orchestration out of SkyWindow while keeping
all UI painting/state updates in the window class.
"""

from __future__ import annotations

import logging
import threading
import time
from collections.abc import Callable
from concurrent.futures import Future
from datetime import datetime, timezone

import numpy as np
from PySide6.QtCore import QObject, Signal

from ..clouddisc import (
    CloudDisc,
    CloudDiscError,
    DataNotFoundError,
    DownloadCancelledError,
    DownloadError,
    RenderError,
    TimeoutError,
    VisibilityError,
    cleanup_satellite_cache,
)
from ..clouddisc.altaz_grid import CloudAltAzGrid
from ..clouddisc.altaz_render import (
    render_altaz_grid_circles,
    render_altaz_missing_mask,
)
from ..clouddisc.providers.select import pick_satellite
from ..clouddisc.types import CloudSourceData, round_down_utc_to_slot
from ..clouddisc.workers.cloud_source import build_cloud_source_fetch_request
from ..clouddisc.workers.cloud_source_worker import run_cloud_source_worker_process
from .application_services import ApplicationServices, wait_for_gui_futures

logger = logging.getLogger(__name__)
DEFAULT_CLOUD_FOV_OVERSCAN_DEG = 2.0


class CloudController(QObject):
    """Orchestrates cloud source fetch and render work."""

    cloud_started = Signal(object)
    cloud_source_ready = Signal(object)
    cloud_ready = Signal(object)
    cloud_failed = Signal(object)

    def __init__(
        self,
        clouddisc: CloudDisc,
        services: ApplicationServices | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._owns_services = services is None
        self._services = services or ApplicationServices()
        self._clouddisc = clouddisc
        self._source_is_running = False
        self._render_is_running = False
        self._active_source_request_key: tuple[object, ...] | None = None
        self._pending_source_request_key: tuple[object, ...] | None = None
        self._pending_source_request: dict[str, object] | None = None
        self._active_render_request_key: tuple[object, ...] | None = None
        self._pending_render_request_key: tuple[object, ...] | None = None
        self._pending_render_request: dict[str, object] | None = None
        self._latest_source: CloudSourceData | None = None
        self._latest_source_request_id = 0
        self._latest_render_request_id = 0
        self._stopping = False
        self._cleanup_counter = 0
        self._cleanup_interval = 10
        self._active_workers: set[Future[None]] = set()
        self._lock = threading.Lock()
        self._download_abort_event = threading.Event()

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        with self._lock:
            self._stopping = True
            self._pending_source_request = None
            self._pending_source_request_key = None
            self._pending_render_request = None
            self._pending_render_request_key = None
        self._download_abort_event.set()
        self._wait_for_workers(wait_timeout_s)
        if self._owns_services:
            self._services.shutdown(wait=True)

    def invalidate_pending_render_results(self) -> None:
        """Mark in-flight render results as stale and drop queued render work."""
        with self._lock:
            if self._stopping:
                return
            self._latest_render_request_id += 1
            self._pending_render_request = None
            self._pending_render_request_key = None

    def has_in_flight_update(self) -> bool:
        """Return True while any cloud fetch or render work is still running."""
        with self._lock:
            return (
                self._source_is_running
                or self._render_is_running
                or self._pending_source_request is not None
                or self._pending_render_request is not None
                or bool(self._active_workers)
            )

    def has_source_data(self) -> bool:
        with self._lock:
            return self._latest_source is not None

    def update(
        self,
        *,
        lat: float,
        lon: float,
        reason: str = "manual",
    ) -> bool:
        """Fetch a new cloud source, without rendering it."""
        return self.update_source(lat=lat, lon=lon, reason=reason)

    def update_source(
        self,
        *,
        lat: float,
        lon: float,
        reason: str = "manual",
    ) -> bool:
        request_key = self._source_request_key(lat=lat, lon=lon)
        request = {
            "lat": float(lat),
            "lon": float(lon),
            "reason": str(reason),
        }
        run_cleanup = False
        with self._lock:
            if self._stopping:
                return False
            if self._source_is_running:
                if request_key in {
                    self._active_source_request_key,
                    self._pending_source_request_key,
                }:
                    return False
                self._pending_source_request = dict(request)
                self._pending_source_request_key = request_key
                return False
            self._latest_source_request_id += 1
            request_id = int(self._latest_source_request_id)
            request["request_id"] = request_id
            self._source_is_running = True
            self._active_source_request_key = request_key
            self._pending_source_request = None
            self._pending_source_request_key = None
            run_cleanup = self._tick_cleanup()

        if run_cleanup:
            self._spawn_worker(target=self._cleanup_cache, kwargs={})

        sat = self._predicted_satellite(request["lat"], request["lon"])
        self.cloud_started.emit({"satellite": sat, "banner": "Clouds: downloading..."})
        self._spawn_worker(
            target=self._run_source_update, kwargs=request
        )
        return True

    def update_render(
        self,
        *,
        alt: float,
        az: float,
        radius_px: int,
        content_fov_deg: float,
        reason: str = "manual",
        render_generation: int = 0,
    ) -> bool:
        with self._lock:
            source = self._latest_source
        if source is None:
            return False
        source_key = getattr(source, "source_key", None)
        request_key = self._render_request_key(
            source_key=source_key,
            alt=alt,
            az=az,
            radius_px=radius_px,
            content_fov_deg=content_fov_deg,
            render_generation=render_generation,
        )
        request = {
            "alt": float(alt),
            "az": float(az),
            "radius_px": int(radius_px),
            "content_fov_deg": float(content_fov_deg),
            "reason": str(reason),
            "render_generation": int(render_generation),
        }
        with self._lock:
            if self._stopping:
                return False
            if self._render_is_running:
                if request_key in {
                    self._active_render_request_key,
                    self._pending_render_request_key,
                }:
                    return False
                self._pending_render_request = dict(request)
                self._pending_render_request_key = request_key
                return False
            self._latest_render_request_id += 1
            request_id = int(self._latest_render_request_id)
            request["request_id"] = request_id
            request["source_id"] = id(source)
            request["source_key"] = source_key
            self._render_is_running = True
            self._active_render_request_key = request_key
            self._pending_render_request = None
            self._pending_render_request_key = None

        self._spawn_worker(
            target=self._run_render_update, kwargs=request
        )
        return True

    def _tick_cleanup(self) -> bool:
        run = (self._cleanup_counter % self._cleanup_interval) == 0
        self._cleanup_counter += 1
        return run

    def _predicted_satellite(self, lat: float, lon: float) -> str:
        return pick_satellite(lat, lon, ("AUTO",))

    def _source_request_key(self, *, lat: float, lon: float) -> tuple[object, ...]:
        return (
            "source",
            round(float(lat), 6),
            round(float(lon), 6),
            round_down_utc_to_slot(datetime.now(timezone.utc)).isoformat(),
        )

    def _render_request_key(
        self,
        *,
        source_key: object | None,
        alt: float,
        az: float,
        radius_px: int,
        content_fov_deg: float,
        render_generation: int,
    ) -> tuple[object, ...]:
        return (
            "render",
            source_key,
            round(float(alt), 6),
            round(float(az) % 360.0, 6),
            int(radius_px),
            round(float(content_fov_deg), 6),
            int(render_generation),
        )

    def _spawn_worker(
        self,
        *,
        target: Callable[..., None],
        kwargs: dict[str, object],
    ) -> None:
        def runner() -> None:
            target(**kwargs)

        worker = self._services.submit(runner)
        with self._lock:
            if self._stopping:
                return
            self._active_workers.add(worker)
            if worker.done():
                self._active_workers.discard(worker)
                return
        worker.add_done_callback(self._unregister_worker)

    def _unregister_worker(self, worker: Future[None]) -> None:
        with self._lock:
            self._active_workers.discard(worker)

    def _wait_for_workers(self, wait_timeout_s: float | None) -> None:
        deadline = (
            None
            if wait_timeout_s is None
            else time.monotonic() + max(0.0, float(wait_timeout_s))
        )
        while True:
            with self._lock:
                workers = tuple(self._active_workers)
            if not workers:
                return
            if deadline is None:
                wait_for_gui_futures(workers, None)
                continue
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                logger.warning(
                    "Timed out waiting for %d cloud worker task(s) to finish during shutdown",
                    len(workers),
                )
                return
            wait_for_gui_futures(workers, remaining)

    def _cleanup_cache(self) -> None:
        try:
            logger.info("Running satellite cache cleanup...")
            cleanup_satellite_cache(self._clouddisc.cfg.cache_root())
            logger.info("Done: Satellite cache cleanup.")
        except Exception:
            logger.exception("Error during cache cleanup")

    def _run_source_update(
        self,
        *,
        lat: float,
        lon: float,
        reason: str,
        request_id: int,
    ) -> None:
        next_req: dict[str, object] | None = None
        try:
            if reason == "initial":
                logger.info("Fetching initial cloud data (reason=%s)...", reason)
            else:
                logger.info("Fetching cloud source data (reason=%s)...", reason)

            try:
                source_request = build_cloud_source_fetch_request(lat=lat, lon=lon)
                source = run_cloud_source_worker_process(
                    self._clouddisc,
                    source_request,
                    request_id=request_id,
                    abort_event=self._download_abort_event,
                )
                with self._lock:
                    is_latest = (
                        not self._stopping
                        and request_id == self._latest_source_request_id
                    )
                    if is_latest:
                        self._latest_source = source
                if is_latest:
                    satellite = str(getattr(source, "satellite", "")).strip()
                    product = str(getattr(source, "product", "")).strip()
                    self.cloud_source_ready.emit(
                        {
                            "source": source,
                            "satellite": satellite,
                            "product": product,
                            "source_key": getattr(source, "source_key", None),
                            "refreshed_at_utc": datetime.now(timezone.utc),
                            "banner": "Clouds: projecting...",
                            "altaz_grid": getattr(source, "altaz_grid", None),
                        }
                    )
            except DownloadCancelledError:
                logger.info("Cloud source download cancelled")
            except VisibilityError as e:
                logger.error("Invalid params for cloud-disc image generation: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: unsupported region"})
            except DownloadError as e:
                logger.warning("Network/S3 download error: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: download failed"})
            except TimeoutError as e:
                logger.warning("Clouds fetch timed out: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: timed out"})
            except DataNotFoundError as e:
                logger.info("Satellite data not found in search window: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: unavailable"})
            except RenderError as e:
                logger.error("Failed to decode/render satellite data: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: decode failed"})
            except CloudDiscError as e:
                logger.error("Unexpected clouddisc error: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: unavailable"})
        except Exception:
            logger.exception("Cloud source update failed")
        finally:
            with self._lock:
                self._source_is_running = False
                self._active_source_request_key = None
                if not self._stopping and self._pending_source_request is not None:
                    next_req = dict(self._pending_source_request)
                    self._active_source_request_key = self._pending_source_request_key
                    self._pending_source_request = None
                    self._pending_source_request_key = None
                    self._source_is_running = True
            if next_req is not None:
                sat = self._predicted_satellite(next_req["lat"], next_req["lon"])
                self.cloud_started.emit(
                    {"satellite": sat, "banner": "Clouds: downloading..."}
                )
                self._spawn_worker(
                    target=self._run_source_update, kwargs=next_req
                )

    def _run_render_update(
        self,
        *,
        alt: float,
        az: float,
        radius_px: int,
        content_fov_deg: float,
        reason: str,
        request_id: int,
        source_id: int,
        source_key: object | None = None,
        render_generation: int = 0,
    ) -> None:
        next_req: dict[str, object] | None = None
        try:
            with self._lock:
                source = self._latest_source
            if source is None:
                return

            altaz_grid = getattr(source, "altaz_grid", None)
            if not isinstance(altaz_grid, CloudAltAzGrid):
                raise RuntimeError("cloud source is missing alt/az grid")
            with self._services.native_work_lock:
                cloud_rgba = render_altaz_grid_circles(
                    altaz_grid,
                    width=int(round(radius_px * 2 + 1)),
                    height=int(round(radius_px * 2 + 1)),
                    center_alt_deg=alt,
                    center_az_deg=az,
                    edge_fov_deg=content_fov_deg + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
                    mask_fov_deg=content_fov_deg + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
                )
                missing_mask = render_altaz_missing_mask(
                    altaz_grid,
                    width=int(round(radius_px * 2 + 1)),
                    height=int(round(radius_px * 2 + 1)),
                    center_alt_deg=alt,
                    center_az_deg=az,
                    edge_fov_deg=content_fov_deg + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
                    mask_fov_deg=content_fov_deg + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
                )
            meta = getattr(altaz_grid, "meta", None)
            if meta is None:
                from ..clouddisc.types import CloudMeta

                meta = CloudMeta(
                    satellite=altaz_grid.satellite,
                    product=altaz_grid.product,
                    time_utc=altaz_grid.time_utc,
                    src_paths=[],
                )
            coverage_ratio = float(altaz_grid.coverage_ratio)
            logger.info(
                "Cloud render ready (request_id=%s, reason=%s, sat=%s, product=%s, data_time=%s, coverage=%.1f%%)",
                request_id,
                reason,
                getattr(meta, "satellite", "?"),
                getattr(meta, "product", "?"),
                getattr(meta, "time_utc", "?"),
                float(coverage_ratio) * 100.0,
            )
            missing_alpha = np.where(missing_mask > 0, 255, 0).astype(np.uint8)
            cloud_amount_field = None
            finished_at_utc = datetime.now(timezone.utc)

            with self._lock:
                current_source = self._latest_source
                current_source_id = (
                    id(current_source) if current_source is not None else None
                )
                current_source_key = getattr(current_source, "source_key", None)
                is_latest = (
                    not self._stopping
                    and request_id == self._latest_render_request_id
                    and source_id == current_source_id
                    and (source_key is None or source_key == current_source_key)
                )
            if not is_latest:
                logger.debug(
                    "Discard stale cloud render result request_id=%s latest=%s source_id=%s current_source_id=%s",
                    request_id,
                    self._latest_render_request_id,
                    source_id,
                    current_source_id,
                )
                return

            self.cloud_ready.emit(
                {
                    "image": cloud_rgba,
                    "meta": meta,
                    "az": az,
                    "time_utc": finished_at_utc,
                    "finished_at_utc": finished_at_utc,
                    "cloud_amount_field": cloud_amount_field,
                    "altaz_grid": altaz_grid
                    if isinstance(altaz_grid, CloudAltAzGrid)
                    else None,
                    "missing_mask": missing_alpha,
                    "coverage_ratio": coverage_ratio,
                    "request_id": request_id,
                    "source_key": getattr(source, "source_key", None),
                    "render_generation": int(render_generation),
                    "source_expected_count": getattr(
                        source, "source_expected_count", None
                    ),
                    "source_available_count": getattr(
                        source, "source_available_count", None
                    ),
                    "source_completeness_ratio": getattr(
                        source, "source_completeness_ratio", None
                    ),
                }
            )
        except Exception:
            logger.exception("Cloud render update failed")
            self.cloud_failed.emit({"banner": "Clouds: render failed"})
        finally:
            with self._lock:
                self._render_is_running = False
                self._active_render_request_key = None
                if not self._stopping and self._pending_render_request is not None:
                    next_req = dict(self._pending_render_request)
                    self._active_render_request_key = self._pending_render_request_key
                    self._pending_render_request = None
                    self._pending_render_request_key = None
                    self._render_is_running = True
            if next_req is not None:
                self._spawn_worker(
                    target=self._run_render_update, kwargs=next_req
                )
