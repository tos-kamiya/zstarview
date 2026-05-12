from __future__ import annotations

import logging
import threading
import time
from pathlib import Path
from typing import Callable, Optional

from PySide6.QtCore import QObject, Signal

from ..types import ViewerData
from ..water_overlay_layer import resolve_water_overlay_for_viewer

logger = logging.getLogger(__name__)


class WaterOverlayController(QObject):
    water_started = Signal(object)
    water_ready = Signal(object)
    water_failed = Signal(object)

    def __init__(
        self,
        *,
        dem_cache_dir: Path,
        radius_km: float = 2.5,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._dem_cache_dir = Path(dem_cache_dir)
        self._radius_km = float(radius_km)
        self._running = False
        self._stopping = False
        self._completed_key: Optional[str] = None
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
        dataset_name = f"water_r{self._radius_km:.1f}_{float(viewer_data.lat_deg):.5f}_{float(viewer_data.lon_deg):.5f}"
        with self._lock:
            if self._stopping or self._running:
                return False
            if self._completed_key == dataset_name:
                return False
            self._running = True

        if reason == "initial":
            self.water_started.emit({"banner": "Water overlay: downloading..."})
        else:
            self.water_started.emit({"banner": "Water overlay: updating..."})

        self._spawn_worker(
            target=self._run_update,
            kwargs={
                "viewer_data": viewer_data,
                "dataset_name": dataset_name,
                "reason": reason,
            },
            label="water",
        )
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
                    "Timed out waiting for %d water worker thread(s) to finish during shutdown",
                    len(workers),
                )
                return
            for worker in workers:
                worker.join(timeout=remaining)

    def _run_update(
        self,
        *,
        viewer_data: ViewerData,
        dataset_name: str,
        reason: str,
    ) -> None:
        try:
            logger.info("Resolving water overlay (%s)...", reason)
            surfaces = resolve_water_overlay_for_viewer(
                viewer_data,
                radius_km=self._radius_km,
                dem_cache_dir=self._dem_cache_dir,
            )
            if not surfaces:
                raise RuntimeError("Water overlay: no polygons found")
            with self._lock:
                if not self._stopping:
                    self._completed_key = dataset_name
            if not self._stopping:
                self.water_ready.emit(
                    {
                        "surfaces": surfaces,
                        "source": "Water: cache",
                    }
                )
        except Exception as exc:
            logger.warning("Water overlay update failed: %s", exc, exc_info=True)
            if not self._stopping:
                self.water_failed.emit({"banner": "Water overlay: unavailable"})
        finally:
            with self._lock:
                self._running = False
