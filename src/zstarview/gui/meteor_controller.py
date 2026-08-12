from __future__ import annotations

import logging
import threading
from collections.abc import Callable
from concurrent.futures import Future
from datetime import datetime, timezone

from PySide6.QtCore import QObject, Signal

from ..meteors import load_celestial_meteor_trails
from ..meteors.types import MeteorWindowResult
from .worker_pool import submit_gui_work, wait_for_gui_futures

logger = logging.getLogger(__name__)
MeteorLoader = Callable[..., MeteorWindowResult]


class MeteorController(QObject):
    meteor_started = Signal(object)
    meteor_ready = Signal(object)
    meteor_failed = Signal(object)

    def __init__(self, *, loader: MeteorLoader | None = None, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._loader = loader or load_celestial_meteor_trails
        self._stopping = False
        self._running = False
        self._latest_request_id = 0
        self._active_workers: set[Future[None]] = set()
        self._lock = threading.Lock()

    def has_in_flight_update(self) -> bool:
        with self._lock:
            return bool(self._running or self._active_workers)

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        with self._lock:
            self._stopping = True
        wait_for_gui_futures(tuple(self._active_workers), wait_timeout_s)

    def update(self, *, display_time_utc: datetime, observer_lat: float,
               observer_lon: float, observer_height_m: float, reason: str = "manual") -> bool:
        with self._lock:
            if self._stopping or self._running:
                return False
            self._running = True
            self._latest_request_id += 1
            request_id = self._latest_request_id
        self.meteor_started.emit({"banner": "GMN meteors: loading..."})

        def runner() -> None:
            try:
                logger.info("Loading GMN meteor trails (%s)...", reason)
                result = self._loader(
                    display_time_utc,
                    observer_lat=observer_lat,
                    observer_lon=observer_lon,
                    observer_height_m=observer_height_m,
                    now_utc=datetime.now(timezone.utc),
                )
                with self._lock:
                    emit = not self._stopping and request_id == self._latest_request_id
                if emit:
                    self.meteor_ready.emit({"result": result})
            except Exception as exc:
                logger.warning("GMN meteor update failed: %s", exc)
                with self._lock:
                    emit = not self._stopping and request_id == self._latest_request_id
                if emit:
                    self.meteor_failed.emit({"banner": "GMN meteors: unavailable"})
            finally:
                with self._lock:
                    self._running = False

        future = submit_gui_work(runner)
        with self._lock:
            self._active_workers.add(future)
        future.add_done_callback(self._worker_done)
        return True

    def _worker_done(self, future: Future[None]) -> None:
        with self._lock:
            self._active_workers.discard(future)
