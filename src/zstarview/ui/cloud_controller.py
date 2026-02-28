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
from ..clouddisc.providers.select import pick_satellite
from ..paths import CLOUD_SHELL_KM
from ..utils.qt import pil_to_qimage
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
        self._is_running = False
        self._pending_request: Optional[dict] = None
        self._stopping = False
        self._cleanup_counter = 0
        self._cleanup_interval = 10
        self._lock = threading.Lock()

    def shutdown(self) -> None:
        with self._lock:
            self._stopping = True
            self._pending_request = None

    def update(
        self,
        *,
        lat: float,
        lon: float,
        alt: float,
        az: float,
        radius_px: int,
        reason: str = "manual",
    ) -> None:
        req = {
            "lat": float(lat),
            "lon": float(lon),
            "alt": float(alt),
            "az": float(az),
            "radius_px": int(radius_px),
            "reason": reason,
        }

        with self._lock:
            if self._stopping:
                return
            if self._is_running:
                self._pending_request = req
                return
            self._is_running = True

        if self._tick_cleanup():
            cleanup_thread = threading.Thread(target=self._cleanup_cache, daemon=True)
            cleanup_thread.start()

        sat = self._predicted_satellite(req["lat"], req["lon"])
        self.cloud_started.emit({"satellite": sat, "banner": "Clouds: downloading…"})

        worker = threading.Thread(target=self._run_update, kwargs=req, daemon=True)
        worker.start()

    def _tick_cleanup(self) -> bool:
        run = (self._cleanup_counter % self._cleanup_interval) == 0
        self._cleanup_counter += 1
        return run

    def _predicted_satellite(self, lat: float, lon: float) -> str:
        return pick_satellite(lat, lon, ("AUTO",))

    def _cleanup_cache(self) -> None:
        try:
            logger.info("Running satellite cache cleanup...")
            cleanup_satellite_cache(self._clouddisc.cfg.cache_root())
            logger.info("Done: Satellite cache cleanup.")
        except Exception as e:
            logger.error("Error during cache cleanup: %s", e, exc_info=True)

    def _run_update(
        self,
        *,
        lat: float,
        lon: float,
        alt: float,
        az: float,
        radius_px: int,
        reason: str,
    ) -> None:
        try:
            if reason == "initial":
                logger.info("Fetching initial cloud data...")
            else:
                logger.info("Updating cloud data...")

            try:
                pil_img, meta = self._clouddisc.render_now(
                    lat=lat,
                    lon=lon,
                    alt=alt,
                    az=az,
                    radius_px=radius_px,
                    edge_fov_deg=90,
                    mask_fov_deg=93,
                    cloud_shell_km=CLOUD_SHELL_KM,
                )
                qimg = pil_to_qimage(pil_img)
                stripe_density = build_stripe_density_field(qimg)
                self.cloud_ready.emit(
                    {
                        "image": qimg,
                        "meta": meta,
                        "az": az,
                        "time_utc": datetime.now(timezone.utc),
                        "stripe_density": stripe_density,
                    }
                )
            except VisibilityError as e:
                logger.error("Invalid params for cloud-disc image generation: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: no supported satellite for this region", "clear_image": True})
            except DownloadError as e:
                logger.warning("Network/S3 download error: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Network/S3 download error", "clear_image": True})
            except TimeoutError as e:
                logger.warning("Clouds fetch timed out: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Clouds fetch timed out", "clear_image": True})
            except DataNotFoundError as e:
                logger.info("Satellite data not found in search window: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Satellite data not found in search window", "clear_image": True})
            except RenderError as e:
                logger.error("Failed to decode/render satellite data: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Failed to decode/render satellite data", "clear_image": True})
            except CloudDiscError as e:
                logger.error("Unexpected clouddisc error: %s", e)
                self.cloud_failed.emit({"banner": "Clouds: Unexpected clouddisc error", "clear_image": True})
        except Exception as e:
            logger.error("Cloud update failed: %s", e, exc_info=True)
        finally:
            next_req: Optional[dict] = None
            with self._lock:
                self._is_running = False
                if not self._stopping and self._pending_request is not None:
                    next_req = self._pending_request
                    self._pending_request = None
            if next_req is not None:
                self.update(**next_req)
