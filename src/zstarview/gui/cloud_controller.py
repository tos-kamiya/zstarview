# -*- coding: utf-8 -*-
"""
Cloud update controller for UI layer.

This module moves cloud fetching orchestration out of SkyWindow while keeping
all UI painting/state updates in the window class.
"""
from __future__ import annotations

import logging
import threading
import time
from concurrent.futures import Future
from datetime import datetime, timezone
from typing import Callable, Optional

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
from ..clouddisc.providers.select import pick_satellite
from ..clouddisc.types import CloudSourceData
from ..paths import CLOUD_SHELLS_KM
from .composite import build_cloud_amount_field_from_rgba
from .native_work_lock import HEAVY_NATIVE_WORK_LOCK
from .worker_pool import submit_gui_work, wait_for_gui_futures

logger = logging.getLogger(__name__)
DEFAULT_CLOUD_FOV_OVERSCAN_DEG = 2.0


class CloudController(QObject):
    """Orchestrates cloud source fetch and render work."""

    cloud_started = Signal(object)
    cloud_source_ready = Signal(object)
    cloud_ready = Signal(object)
    cloud_failed = Signal(object)

    def __init__(self, clouddisc: CloudDisc, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._clouddisc = clouddisc
        self._source_is_running = False
        self._render_is_running = False
        self._pending_source_request: Optional[dict[str, object]] = None
        self._pending_render_request: Optional[dict[str, object]] = None
        self._latest_source: Optional[CloudSourceData] = None
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
            self._pending_render_request = None
        self._download_abort_event.set()
        self._wait_for_workers(wait_timeout_s)

    def invalidate_pending_render_results(self) -> None:
        """Mark in-flight render results as stale and drop queued render work."""
        with self._lock:
            if self._stopping:
                return
            self._latest_render_request_id += 1
            self._pending_render_request = None

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
        request = {
            "lat": float(lat),
            "lon": float(lon),
            "reason": str(reason),
        }
        run_cleanup = False
        with self._lock:
            if self._stopping:
                return False
            self._latest_source_request_id += 1
            request_id = int(self._latest_source_request_id)
            request["request_id"] = request_id
            if self._source_is_running:
                self._pending_source_request = dict(request)
                return False
            self._source_is_running = True
            run_cleanup = self._tick_cleanup()

        if run_cleanup:
            self._spawn_worker(target=self._cleanup_cache, kwargs={}, label="cleanup")

        sat = self._predicted_satellite(request["lat"], request["lon"])
        self.cloud_started.emit({"satellite": sat, "banner": "Clouds: downloading..."})
        self._spawn_worker(target=self._run_source_update, kwargs=request, label="source")
        return True

    def update_render(
        self,
        *,
        lat: float,
        lon: float,
        alt: float,
        az: float,
        radius_px: int,
        content_fov_deg: float,
        reason: str = "manual",
        render_generation: int = 0,
    ) -> bool:
        request = {
            "lat": float(lat),
            "lon": float(lon),
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
            source = self._latest_source
            if source is None:
                return False
            self._latest_render_request_id += 1
            request_id = int(self._latest_render_request_id)
            request["request_id"] = request_id
            request["source_id"] = id(source)
            request["source_key"] = getattr(source, "source_key", None)
            if self._render_is_running:
                self._pending_render_request = dict(request)
                return False
            self._render_is_running = True

        self._spawn_worker(target=self._run_render_update, kwargs=request, label="render")
        return True

    def _tick_cleanup(self) -> bool:
        run = (self._cleanup_counter % self._cleanup_interval) == 0
        self._cleanup_counter += 1
        return run

    def _predicted_satellite(self, lat: float, lon: float) -> str:
        return pick_satellite(lat, lon, ("AUTO",))

    def _spawn_worker(
        self,
        *,
        target: Callable[..., None],
        kwargs: dict[str, object],
        label: str,
    ) -> None:
        def runner() -> None:
            target(**kwargs)

        worker = submit_gui_work(runner)
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
        deadline = None if wait_timeout_s is None else time.monotonic() + max(0.0, float(wait_timeout_s))
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
        except Exception as e:
            logger.error("Error during cache cleanup: %s", e, exc_info=True)

    def _run_source_update(
        self,
        *,
        lat: float,
        lon: float,
        reason: str,
        request_id: int,
    ) -> None:
        next_req: Optional[dict[str, object]] = None
        try:
            if reason == "initial":
                logger.info("Fetching initial cloud data (reason=%s)...", reason)
            else:
                logger.info("Fetching cloud source data (reason=%s)...", reason)

            try:
                with HEAVY_NATIVE_WORK_LOCK:
                    source = self._clouddisc.fetch_source(
                        lat=lat,
                        lon=lon,
                        abort_event=self._download_abort_event,
                    )
                with self._lock:
                    is_latest = not self._stopping and request_id == self._latest_source_request_id
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
        except Exception as e:
            logger.error("Cloud source update failed: %s", e, exc_info=True)
        finally:
            with self._lock:
                self._source_is_running = False
                if not self._stopping and self._pending_source_request is not None:
                    next_req = dict(self._pending_source_request)
                    self._pending_source_request = None
            if next_req is not None:
                sat = self._predicted_satellite(next_req["lat"], next_req["lon"])
                self.cloud_started.emit({"satellite": sat, "banner": "Clouds: downloading..."})
                with self._lock:
                    if not self._stopping:
                        self._source_is_running = True
                self._spawn_worker(target=self._run_source_update, kwargs=next_req, label="source")

    def _run_render_update(
        self,
        *,
        lat: float,
        lon: float,
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
        next_req: Optional[dict[str, object]] = None
        try:
            with self._lock:
                source = self._latest_source
            if source is None:
                return

            with HEAVY_NATIVE_WORK_LOCK:
                cloud_rgba, meta, missing_mask, coverage_ratio = self._clouddisc.render_from_source_with_coverage(
                    source=source,
                    lat=lat,
                    lon=lon,
                    alt=alt,
                    az=az,
                    radius_px=radius_px,
                    edge_fov_deg=content_fov_deg + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
                    mask_fov_deg=content_fov_deg + DEFAULT_CLOUD_FOV_OVERSCAN_DEG,
                    cloud_shells_km=CLOUD_SHELLS_KM,
                )
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
            cloud_amount_field = build_cloud_amount_field_from_rgba(cloud_rgba)
            finished_at_utc = datetime.now(timezone.utc)

            with self._lock:
                current_source = self._latest_source
                current_source_id = id(current_source) if current_source is not None else None
                current_source_key = getattr(current_source, "source_key", None)
                is_latest = (
                    not self._stopping
                    and request_id == self._latest_render_request_id
                    and source_id == current_source_id
                    and (
                        source_key is None
                        or source_key == current_source_key
                    )
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
                    "missing_mask": missing_alpha,
                    "coverage_ratio": coverage_ratio,
                    "request_id": request_id,
                    "source_key": getattr(source, "source_key", None),
                    "render_generation": int(render_generation),
                    "source_expected_count": getattr(source, "source_expected_count", None),
                    "source_available_count": getattr(source, "source_available_count", None),
                    "source_completeness_ratio": getattr(source, "source_completeness_ratio", None),
                }
            )
        except Exception as e:
            logger.error("Cloud render update failed: %s", e, exc_info=True)
        finally:
            with self._lock:
                self._render_is_running = False
                if (
                    not self._stopping
                    and self._pending_render_request is not None
                ):
                    next_req = dict(self._pending_render_request)
                    self._pending_render_request = None
                    self._render_is_running = True
            if next_req is not None:
                self._spawn_worker(target=self._run_render_update, kwargs=next_req, label="render")
