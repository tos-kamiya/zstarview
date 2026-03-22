from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from typing import Callable, Optional

import astropy.time
from PySide6.QtCore import QObject, Signal

from ..aircraft import (
    AircraftBoundingBox,
    CachedAircraftSnapshotSet,
    AircraftOverlayPoint,
    AircraftSnapshot,
    build_observer_bbox,
    fetch_cached_opensky_states,
    project_aircraft_snapshots,
)

logger = logging.getLogger(__name__)

AircraftFetcher = Callable[[AircraftBoundingBox], list[AircraftSnapshot] | CachedAircraftSnapshotSet]
AircraftProjector = Callable[..., list[AircraftOverlayPoint]]


class AircraftController(QObject):
    """Background fetch controller for the planned aircraft overlay."""

    aircraft_started = Signal(object)
    aircraft_ready = Signal(object)
    aircraft_failed = Signal(object)

    def __init__(
        self,
        *,
        fetcher: AircraftFetcher | None = None,
        projector: AircraftProjector | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._fetcher = fetcher or fetch_cached_opensky_states
        self._projector = projector or project_aircraft_snapshots
        self._running = False
        self._stopping = False
        self._pending_request: Optional[dict[str, object]] = None
        self._latest_request_id = 0
        self._lock = threading.Lock()

    def shutdown(self) -> None:
        with self._lock:
            self._stopping = True
            self._pending_request = None

    def update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        observer_height_m: float,
        time_obj: astropy.time.Time,
        reason: str = "manual",
    ) -> bool:
        request = {
            "observer_lat": float(observer_lat),
            "observer_lon": float(observer_lon),
            "observer_height_m": float(observer_height_m),
            "time_obj": time_obj,
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

        self.aircraft_started.emit({"banner": "Aircraft: fetching OpenSky states..."})
        worker = threading.Thread(target=self._run_update, kwargs=request, daemon=True)
        worker.start()
        return True

    def _run_update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        observer_height_m: float,
        time_obj: astropy.time.Time,
        reason: str,
        request_id: int,
    ) -> None:
        next_request: Optional[dict[str, object]] = None
        try:
            if reason == "initial":
                logger.info("Fetching initial aircraft state...")
            else:
                logger.info("Fetching aircraft state...")
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
            overlay_points = self._projector(
                snapshots,
                observer_lat=observer_lat,
                observer_lon=observer_lon,
                observer_height_m=observer_height_m,
                time_obj=time_obj,
            )
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                self.aircraft_ready.emit(
                    {
                        "snapshots": snapshots,
                        "overlay_points": overlay_points,
                        "bbox": bbox,
                        "refreshed_at_utc": refreshed_at_utc,
                        "banner": (
                            "Aircraft: using stale cached OpenSky states"
                            if is_stale
                            else ""
                        ),
                        "source": source,
                    }
                )
        except Exception as exc:
            logger.warning("Aircraft update failed: %s", exc, exc_info=True)
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                self.aircraft_failed.emit({"banner": f"Aircraft: {exc}"})
        finally:
            with self._lock:
                self._running = False
                if not self._stopping and self._pending_request is not None:
                    next_request = dict(self._pending_request)
                    self._pending_request = None
                    self._running = True
            if next_request is not None:
                self.aircraft_started.emit({"banner": "Aircraft: fetching OpenSky states..."})
                worker = threading.Thread(target=self._run_update, kwargs=next_request, daemon=True)
                worker.start()
