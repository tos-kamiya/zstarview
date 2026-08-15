"""Background loading and in-memory coordination for Sun hover images."""

from __future__ import annotations

import threading
from concurrent.futures import Future
from datetime import datetime, timezone

from PySide6.QtCore import QObject, Signal

from ..solar_hover import (
    SolarHoverImage,
    fetch_solar_hover_image,
    normalize_solar_hover_time,
)
from .worker_pool import submit_gui_work


class SolarHoverController(QObject):
    """Load AIA 193 images without blocking the GUI thread."""

    image_ready = Signal(object)

    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._lock = threading.Lock()
        self._images: dict[datetime, SolarHoverImage] = {}
        self._inflight: dict[datetime, Future] = {}
        self._failed_until: dict[datetime, float] = {}

    def request(self, target_time: datetime) -> SolarHoverImage | None:
        key = normalize_solar_hover_time(target_time)
        with self._lock:
            cached = self._images.get(key)
            if cached is not None:
                return cached
            if key in self._inflight:
                return None
            if self._failed_until.get(key, 0.0) > datetime.now(timezone.utc).timestamp():
                return None
            future = submit_gui_work(fetch_solar_hover_image, target_time)
            self._inflight[key] = future
        future.add_done_callback(lambda completed: self._complete(key, completed))
        return None

    def _complete(self, key: datetime, future: Future) -> None:
        try:
            image = future.result()
        except Exception:
            with self._lock:
                self._inflight.pop(key, None)
                self._failed_until[key] = datetime.now(timezone.utc).timestamp() + 30.0
            self.image_ready.emit((key, None))
            return
        with self._lock:
            self._inflight.pop(key, None)
            self._images[key] = image
        self.image_ready.emit((key, image))


__all__ = ["SolarHoverController"]
