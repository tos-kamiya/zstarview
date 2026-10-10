from __future__ import annotations

import logging
import threading
import time
from collections.abc import Callable
from concurrent.futures import Future
from datetime import datetime, timezone

from PySide6.QtCore import QObject, Signal

from ..aircraft import (
    AircraftBoundingBox,
    AircraftSnapshot,
    CachedAircraftSnapshotSet,
    build_observer_bbox,
    fetch_cached_opensky_states,
    load_aircraft_cache_for_time,
)
from .application_services import ApplicationServices, wait_for_gui_futures

logger = logging.getLogger(__name__)

AircraftFetcher = Callable[[AircraftBoundingBox], list[AircraftSnapshot] | CachedAircraftSnapshotSet]


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


class AircraftController(QObject):
    """Background fetch controller for the planned aircraft overlay."""

    aircraft_started = Signal(object)
    aircraft_ready = Signal(object)
    aircraft_failed = Signal(object)

    def __init__(
        self,
        *,
        fetcher: AircraftFetcher | None = None,
        projector: object | None = None,
        services: ApplicationServices | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._owns_services = services is None
        self._services = services or ApplicationServices()
        self._fetcher = fetcher or fetch_cached_opensky_states
        self._projector = projector
        self._running = False
        self._stopping = False
        self._pending_request: dict[str, object] | None = None
        self._latest_request_id = 0
        self._active_workers: set[Future[None]] = set()
        self._lock = threading.Lock()

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        with self._lock:
            self._stopping = True
            self._pending_request = None
        self._wait_for_workers(wait_timeout_s)
        if self._owns_services:
            self._services.shutdown(wait=True)

    def has_in_flight_update(self) -> bool:
        with self._lock:
            return bool(self._running or self._pending_request is not None or self._active_workers)

    def update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        reason: str = "manual",
        target_time_utc: datetime | None = None,
    ) -> bool:
        target = (
            _normalize_utc(target_time_utc)
            if target_time_utc is not None
            else None
        )
        historical = bool(
            target is not None
            and (datetime.now(timezone.utc) - target).total_seconds() > 300.0
        )
        request = {
            "observer_lat": float(observer_lat),
            "observer_lon": float(observer_lon),
            "reason": str(reason),
            "target_time_utc": target,
            "historical": historical,
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

        self.aircraft_started.emit(
            {
                "banner": "Aircraft: cache" if historical else "Aircraft: fetching OpenSky states...",
                "historical": historical,
            }
        )
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
                    "Timed out waiting for %d aircraft worker task(s) to finish during shutdown",
                    len(workers),
                )
                return
            wait_for_gui_futures(workers, remaining)

    def _run_update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        reason: str,
        request_id: int,
        target_time_utc: datetime | None = None,
        historical: bool = False,
    ) -> None:
        next_request: dict[str, object] | None = None
        try:
            if historical:
                bbox = build_observer_bbox(observer_lat, observer_lon)
                fetched = load_aircraft_cache_for_time(
                    bbox,
                    target_time_utc or datetime.now(timezone.utc),
                )
                if fetched is None:
                    raise RuntimeError("Aircraft: historical cache unavailable")
                snapshots = fetched.snapshots
                refreshed_at_utc = fetched.fetched_at_utc
                source = fetched.source
                banner = "Aircraft: cache"
            elif reason == "initial":
                logger.info("Fetching initial aircraft state...")
            else:
                logger.info("Fetching aircraft state...")
            if not historical:
                bbox = build_observer_bbox(observer_lat, observer_lon)
                fetched = self._fetcher(bbox)
                if isinstance(fetched, list):
                    snapshots = fetched
                    refreshed_at_utc = datetime.now(timezone.utc)
                    source = "opensky"
                    is_stale = False
                else:
                    snapshots = fetched.snapshots
                    refreshed_at_utc = fetched.fetched_at_utc
                    source = fetched.source
                    is_stale = bool(fetched.is_stale)
                if source == "rate-limited-skip":
                    banner = "Aircraft: deferred"
                elif is_stale:
                    banner = "Aircraft: cache"
                else:
                    banner = ""
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                self.aircraft_ready.emit(
                    {
                        "snapshots": snapshots,
                        "bbox": bbox,
                        "refreshed_at_utc": refreshed_at_utc,
                        "banner": banner,
                        "source": source,
                    }
                )
        except Exception as exc:
            logger.warning("Aircraft update failed: %s", exc)
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                self.aircraft_failed.emit(
                    {"banner": "Aircraft: unavailable", "historical": historical}
                )
        finally:
            with self._lock:
                self._running = False
                if not self._stopping and self._pending_request is not None:
                    next_request = dict(self._pending_request)
                    self._pending_request = None
                    self._running = True
            if next_request is not None:
                historical = bool(next_request.get("historical", False))
                self.aircraft_started.emit(
                    {
                        "banner": "Aircraft: cache"
                        if historical
                        else "Aircraft: fetching OpenSky states...",
                        "historical": historical,
                    }
                )
                self._spawn_worker(target=self._run_update, kwargs=next_request)
