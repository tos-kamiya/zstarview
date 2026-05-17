from __future__ import annotations

import logging
import threading
import time
from concurrent.futures import Future
from datetime import datetime, timezone
from typing import Callable, Optional
from urllib.error import URLError

from PySide6.QtCore import QObject, Signal

from ..search.jpl import extract_horizons_state_vector
from ..satellites import fetch_horizons_vector_csv
from ..search.models import SearchJumpTarget
from .worker_pool import submit_gui_work, wait_for_gui_futures

logger = logging.getLogger(__name__)


class JplSmallBodyController(QObject):
    jpl_started = Signal(object)
    jpl_ready = Signal(object)
    jpl_failed = Signal(object)

    def __init__(self, *, parent: QObject | None = None) -> None:
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

    def update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        observer_height_m: float,
        target: SearchJumpTarget,
        target_time_utc: datetime,
        reason: str = "manual",
    ) -> bool:
        request = {
            "observer_lat": float(observer_lat),
            "observer_lon": float(observer_lon),
            "observer_height_m": float(observer_height_m),
            "target": target,
            "target_time_utc": target_time_utc.astimezone(timezone.utc),
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

        self.jpl_started.emit({"banner": "JPL: fetching small-body ephemeris..."})
        self._spawn_worker(target=self._run_update, kwargs=request, label="jpl")
        return True

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
                    "Timed out waiting for %d JPL worker task(s) to finish during shutdown",
                    len(workers),
                )
                return
            wait_for_gui_futures(workers, remaining)

    def _run_update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        observer_height_m: float,
        target: SearchJumpTarget,
        target_time_utc: datetime,
        reason: str,
        request_id: int,
    ) -> None:
        next_request: Optional[dict[str, object]] = None
        try:
            target_label = str(getattr(target, "label", "")).strip() or "<unnamed>"
            command = str(target.command).strip()
            if not command:
                command = f"DES={target.object_key};" if target.object_key else ""
            if not command:
                raise RuntimeError("JPL small-body target has no usable command")
            logger.info(
                "Fetching JPL small-body state vector (%s): target=%s command=%s target_time_utc=%s",
                reason,
                target_label,
                command,
                target_time_utc.astimezone(timezone.utc).isoformat(),
            )
            vector_rows = fetch_horizons_vector_csv(
                command,
                target_time_utc=target_time_utc,
            )
            state_vector = extract_horizons_state_vector(vector_rows)
            if state_vector is None:
                raise RuntimeError("JPL vector table did not contain a state vector sample")
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                position_km, velocity_km_s = state_vector
                payload = {
                    "target": target,
                    "target_time_utc": target_time_utc.astimezone(timezone.utc),
                    "refreshed_at_utc": datetime.now(timezone.utc),
                    "rows": vector_rows,
                    "reason": reason,
                    "horizons_epoch_utc": target_time_utc.astimezone(timezone.utc),
                    "horizons_position_km": position_km,
                    "horizons_velocity_km_s": velocity_km_s,
                }
                self.jpl_ready.emit(
                    payload
                )
        except Exception as exc:
            logger.warning(
                "JPL small-body update failed (%s): target=%s command=%s target_time_utc=%s error=%s exception=%r",
                reason,
                str(getattr(target, "label", "")).strip() or "<unnamed>",
                str(getattr(target, "command", "")).strip() or "<missing>",
                target_time_utc.astimezone(timezone.utc).isoformat(),
                str(exc).strip() or exc.__class__.__name__,
                exc,
                exc_info=not isinstance(exc, URLError),
            )
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                self.jpl_failed.emit(
                    {
                        "target": target,
                        "target_time_utc": target_time_utc.astimezone(timezone.utc),
                        "refreshed_at_utc": datetime.now(timezone.utc),
                        "banner": f"JPL: {exc}",
                        "error": str(exc),
                        "reason": reason,
                    }
                )
        finally:
            with self._lock:
                self._running = False
                if not self._stopping and self._pending_request is not None:
                    next_request = dict(self._pending_request)
                    self._pending_request = None
                    self._running = True
            if next_request is not None:
                self.jpl_started.emit({"banner": "JPL: fetching small-body ephemeris..."})
                self._spawn_worker(target=self._run_update, kwargs=next_request, label="jpl")
