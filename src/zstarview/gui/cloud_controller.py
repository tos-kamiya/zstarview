# -*- coding: utf-8 -*-
"""
Cloud update controller for UI layer.

This module moves cloud fetching orchestration out of SkyWindow while keeping
all UI painting/state updates in the window class.
"""
from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from typing import Optional

import numpy as np
from PySide6.QtCore import QObject, Signal

from ..clouddisc import (
    CloudDisc,
    CloudDiscError,
    DataNotFoundError,
    DownloadError,
    RenderError,
    TimeoutError,
    VisibilityError,
    cleanup_satellite_cache,
)
from ..clouddisc.types import CloudSourceData
from ..clouddisc.providers.select import pick_satellite
from ..paths import CLOUD_SHELL_KM
from ..render.qt_image import np_rgba_to_qimage
from .composite import build_stripe_density_field

logger = logging.getLogger(__name__)


class CloudController(QObject):
    """Orchestrates cloud fetch/render and exposes state via Qt signals."""

    cloud_started = Signal(object)
    cloud_ready = Signal(object)
    cloud_failed = Signal(object)

    def __init__(self, clouddisc: CloudDisc, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._clouddisc = clouddisc
        self._source_is_running = False
        self._render_is_running = False
        self._pending_source_request: Optional[dict] = None
        self._pending_render_request: Optional[dict] = None
        self._last_render_request: Optional[dict] = None
        self._latest_source: Optional[CloudSourceData] = None
        self._latest_request_id = 0
        self._stopping = False
        self._cleanup_counter = 0
        self._cleanup_interval = 10
        self._lock = threading.Lock()

    def shutdown(self) -> None:
        with self._lock:
            self._stopping = True
            self._pending_source_request = None
            self._pending_render_request = None

    def update(
        self,
        *,
        lat: float,
        lon: float,
        alt: float,
        az: float,
        radius_px: int,
        content_fov_deg: float = 100.0,
        reason: str = "manual",
        render_generation: int = 0,
    ) -> None:
        render_req = {
            "lat": float(lat),
            "lon": float(lon),
            "alt": float(alt),
            "az": float(az),
            "radius_px": int(radius_px),
            "content_fov_deg": float(content_fov_deg),
            "reason": reason,
            "render_generation": int(render_generation),
        }
        source_req = {
            "lat": float(lat),
            "lon": float(lon),
            "reason": reason,
        }
        start_render_req: Optional[dict] = None
        start_source_req: Optional[dict] = None
        run_cleanup = False

        with self._lock:
            if self._stopping:
                return
            self._latest_request_id += 1
            request_id = int(self._latest_request_id)
            render_req["request_id"] = request_id
            self._last_render_request = dict(render_req)

            if self._render_is_running:
                self._pending_render_request = dict(render_req)
            else:
                start_render_req = dict(render_req)
                self._render_is_running = True

            need_source = self._latest_source is None or self._should_refresh_source(reason)
            if need_source:
                if self._source_is_running:
                    self._pending_source_request = dict(source_req)
                else:
                    self._source_is_running = True
                    start_source_req = dict(source_req)
                    run_cleanup = self._tick_cleanup()

        if run_cleanup:
            cleanup_thread = threading.Thread(target=self._cleanup_cache, daemon=True)
            cleanup_thread.start()

        if start_source_req is not None:
            sat = self._predicted_satellite(start_source_req["lat"], start_source_req["lon"])
            self.cloud_started.emit({"satellite": sat, "banner": "Clouds: downloading..."})
            source_worker = threading.Thread(target=self._run_source_update, kwargs=start_source_req, daemon=True)
            source_worker.start()

        if start_render_req is not None:
            render_worker = threading.Thread(target=self._run_render_update, kwargs=start_render_req, daemon=True)
            render_worker.start()

    def _tick_cleanup(self) -> bool:
        run = (self._cleanup_counter % self._cleanup_interval) == 0
        self._cleanup_counter += 1
        return run

    @staticmethod
    def _should_refresh_source(reason: str) -> bool:
        return reason in {"initial", "timer", "toggle-on", "manual"}

    def _predicted_satellite(self, lat: float, lon: float) -> str:
        return pick_satellite(lat, lon, ("AUTO",))

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
    ) -> None:
        next_req: Optional[dict] = None
        try:
            if reason == "initial":
                logger.info("Fetching initial cloud data...")
            else:
                logger.info("Fetching cloud source data...")

            try:
                source = self._clouddisc.fetch_source(
                    lat=lat,
                    lon=lon,
                )
                with self._lock:
                    if not self._stopping:
                        self._latest_source = source

                rerender_req = None
                with self._lock:
                    if not self._stopping:
                        if self._pending_render_request is not None:
                            rerender_req = dict(self._pending_render_request)
                            self._pending_render_request = None
                        elif self._last_render_request is not None:
                            rerender_req = dict(self._last_render_request)
                        if rerender_req is not None:
                            if not self._render_is_running:
                                self._render_is_running = True
                            else:
                                # Keep latest rerender request queued while current render is running.
                                self._pending_render_request = dict(rerender_req)
                                rerender_req = None
                if rerender_req is not None:
                    worker = threading.Thread(target=self._run_render_update, kwargs=rerender_req, daemon=True)
                    worker.start()
            except VisibilityError as e:
                logger.error("Invalid params for cloud-disc image generation: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: no supported satellite for this region"})
            except DownloadError as e:
                logger.warning("Network/S3 download error: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Network/S3 download error"})
            except TimeoutError as e:
                logger.warning("Clouds fetch timed out: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Clouds fetch timed out"})
            except DataNotFoundError as e:
                logger.info("Satellite data not found in search window: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Satellite data not found in search window"})
            except RenderError as e:
                logger.error("Failed to decode/render satellite data: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Failed to decode/render satellite data"})
            except CloudDiscError as e:
                logger.error("Unexpected clouddisc error: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Unexpected clouddisc error"})
        except Exception as e:
            logger.error("Cloud source update failed: %s", e, exc_info=True)
        finally:
            with self._lock:
                self._source_is_running = False
                if not self._stopping and self._pending_source_request is not None:
                    next_req = self._pending_source_request
                    self._pending_source_request = None
            if next_req is not None:
                sat = self._predicted_satellite(next_req["lat"], next_req["lon"])
                self.cloud_started.emit({"satellite": sat, "banner": "Clouds: downloading..."})
                worker = threading.Thread(target=self._run_source_update, kwargs=next_req, daemon=True)
                with self._lock:
                    if not self._stopping:
                        self._source_is_running = True
                worker.start()

    def _run_render_update(
        self,
        *,
        lat: float,
        lon: float,
        alt: float,
        az: float,
        radius_px: int,
        content_fov_deg: float = 100.0,
        reason: str,
        request_id: int,
        render_generation: int = 0,
    ) -> None:
        next_req: Optional[dict] = None
        try:
            with self._lock:
                source = self._latest_source
            if source is None:
                return

            cloud_rgba, meta, missing_mask, coverage_ratio = self._clouddisc.render_from_source_with_coverage(
                source=source,
                lat=lat,
                lon=lon,
                alt=alt,
                az=az,
                radius_px=radius_px,
                edge_fov_deg=90,
                mask_fov_deg=content_fov_deg,
                cloud_shell_km=CLOUD_SHELL_KM,
            )
            logger.info(
                "Cloud render ready (request_id=%s, sat=%s, product=%s, data_time=%s, coverage=%.1f%%)",
                request_id,
                getattr(meta, "satellite", "?"),
                getattr(meta, "product", "?"),
                getattr(meta, "time_utc", "?"),
                float(coverage_ratio) * 100.0,
            )
            qimg = np_rgba_to_qimage(cloud_rgba)
            missing_rgba = np.zeros((missing_mask.shape[0], missing_mask.shape[1], 4), dtype=np.uint8)
            missing_rgba[..., 3] = np.where(missing_mask > 0, 255, 0).astype(np.uint8)
            missing_qimg = np_rgba_to_qimage(missing_rgba)
            stripe_density = build_stripe_density_field(qimg)

            with self._lock:
                is_latest = (request_id == self._latest_request_id)
            if not is_latest:
                logger.debug("Discard stale cloud render result request_id=%s latest=%s", request_id, self._latest_request_id)
                return

            self.cloud_ready.emit(
                {
                    "image": qimg,
                    "meta": meta,
                    "az": az,
                    "time_utc": datetime.now(timezone.utc),
                    "stripe_density": stripe_density,
                    "missing_mask": missing_qimg,
                    "coverage_ratio": coverage_ratio,
                    "request_id": request_id,
                    "source_key": getattr(source, "source_key", None),
                    "render_generation": int(render_generation),
                }
            )
        except Exception as e:
            logger.error("Cloud render update failed: %s", e, exc_info=True)
        finally:
            with self._lock:
                self._render_is_running = False
                if not self._stopping and self._pending_render_request is not None:
                    next_req = self._pending_render_request
                    self._pending_render_request = None
                    self._render_is_running = True
            if next_req is not None:
                worker = threading.Thread(target=self._run_render_update, kwargs=next_req, daemon=True)
                worker.start()
