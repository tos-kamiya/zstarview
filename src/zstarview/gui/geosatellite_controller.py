from __future__ import annotations

import logging
import threading
import time
from concurrent.futures import Future
from datetime import datetime, timezone
from types import SimpleNamespace
from typing import Callable, Optional

import numpy as np
from PySide6.QtCore import QObject, Signal

from ..geosatellite.pipeline import run_geo_satellite_pipeline
from ..geosatellite.projection import render_gray_image_to_cloud_rgba
from .composite import build_cloud_amount_field_from_rgba
from .native_work_lock import HEAVY_NATIVE_WORK_LOCK
from .worker_pool import submit_gui_work, wait_for_gui_futures

logger = logging.getLogger(__name__)


class GeoSatelliteController(QObject):
    geo_started = Signal(object)
    geo_source_ready = Signal(object)
    geo_ready = Signal(object)
    geo_failed = Signal(object)

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._running = False
        self._stopping = False
        self._pending_request: Optional[dict[str, object]] = None
        self._latest_request_id = 0
        self._active_workers: set[Future[None]] = set()
        self._lock = threading.Lock()

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        with self._lock:
            self._stopping = True
            self._pending_request = None
        self._wait_for_workers(wait_timeout_s)

    def has_in_flight_update(self) -> bool:
        with self._lock:
            return bool(self._running or self._pending_request is not None or self._active_workers)

    def update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        alt: float,
        az: float,
        fov_deg: float,
        render_generation: int = 0,
        reason: str = "manual",
    ) -> bool:
        request = {
            "observer_lat": float(observer_lat),
            "observer_lon": float(observer_lon),
            "alt": float(alt),
            "az": float(az),
            "fov_deg": float(fov_deg),
            "render_generation": int(render_generation),
            "reason": str(reason),
        }
        with self._lock:
            if self._stopping:
                return False
            self._latest_request_id += 1
            request["request_id"] = int(self._latest_request_id)
            if self._running:
                self._pending_request = dict(request)
                return False
            self._running = True

        self.geo_started.emit({"banner": "Geo-sat + Downloading"})
        self._spawn_worker(target=self._run_update, kwargs=request)
        return True

    def _spawn_worker(
        self,
        *,
        target: Callable[..., None],
        kwargs: dict[str, object],
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
                    "Timed out waiting for %d Geo-satellite worker task(s) to finish during shutdown",
                    len(workers),
                )
                return
            wait_for_gui_futures(workers, remaining)

    def _run_update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        alt: float,
        az: float,
        fov_deg: float,
        render_generation: int,
        reason: str,
        request_id: int,
    ) -> None:
        next_req: Optional[dict[str, object]] = None
        try:
            logger.info("Fetching Geo-satellite data (%s)...", reason)

            def _status_callback(stage: str, payload: dict[str, object]) -> None:
                with self._lock:
                    is_latest = not self._stopping and request_id == self._latest_request_id
                if not is_latest:
                    return
                if stage == "projecting":
                    download = payload.get("download")
                    captured_at_utc = payload.get("captured_at_utc")
                    fetched_at_utc = payload.get("fetched_at_utc")
                    self.geo_source_ready.emit(
                        {
                            "download": download,
                            "captured_at_utc": captured_at_utc,
                            "fetched_at_utc": fetched_at_utc,
                            "refreshed_at_utc": datetime.now(timezone.utc),
                            "banner": "Geo-sat + Projecting",
                        }
                    )

            with HEAVY_NATIVE_WORK_LOCK:
                result = run_geo_satellite_pipeline(
                    observer_lat=observer_lat,
                    observer_lon=observer_lon,
                    alt=alt,
                    az=az,
                    fov_deg=fov_deg,
                    status_callback=_status_callback,
                )

            cloud_rgba = render_gray_image_to_cloud_rgba(result.disc_gray)
            cloud_amount_field = build_cloud_amount_field_from_rgba(
                cloud_rgba,
                source_cache_key=int(result.intermediate.raw_digest[:16], 16),
            )
            download = result.download
            captured_at_utc = download.captured_at_utc or download.fetched_at_utc
            meta = SimpleNamespace(
                satellite="Geo-sat",
                product=str(download.kind),
                time_utc=captured_at_utc,
                src_paths=[],
            )
            with self._lock:
                is_latest = not self._stopping and request_id == self._latest_request_id
            if is_latest:
                self.geo_ready.emit(
                    {
                        "image": cloud_rgba,
                        "meta": meta,
                        "az": float(az),
                        "time_utc": captured_at_utc,
                        "cloud_amount_field": cloud_amount_field,
                        "missing_mask": None,
                        "coverage_ratio": float(np.count_nonzero(cloud_rgba[..., 3]) / max(1, cloud_rgba[..., 3].size)),
                        "source_key": None,
                        "request_id": request_id,
                        "render_generation": int(render_generation),
                        "captured_at_utc": captured_at_utc,
                        "fetched_at_utc": download.fetched_at_utc,
                        "refreshed_at_utc": datetime.now(timezone.utc),
                        "banner": "",
                    }
                )
        except Exception as exc:
            logger.warning("Geo-satellite update failed: %s", exc, exc_info=True)
            with self._lock:
                is_latest = not self._stopping and request_id == self._latest_request_id
            if is_latest:
                self.geo_failed.emit({"banner": "Geo-sat + error"})
        finally:
            with self._lock:
                self._running = False
                if not self._stopping and self._pending_request is not None:
                    next_req = dict(self._pending_request)
                    self._pending_request = None
            if next_req is not None:
                self.geo_started.emit({"banner": "Geo-sat + Downloading"})
                with self._lock:
                    if not self._stopping:
                        self._running = True
                self._spawn_worker(target=self._run_update, kwargs=next_req)
