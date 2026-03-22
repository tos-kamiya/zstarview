from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from typing import Callable, Optional

import astropy.time
from PySide6.QtCore import QObject, Signal

from ..satellites import (
    CachedSatelliteElementSet,
    fetch_cached_satellite_elements,
    project_satellite_records,
)
from ..satellites.types import SatelliteOmmRecord, SatelliteOverlayPoint

logger = logging.getLogger(__name__)

SatelliteFetcher = Callable[[str], CachedSatelliteElementSet]
SatelliteProjector = Callable[..., list[SatelliteOverlayPoint]]


class SatelliteController(QObject):
    satellite_started = Signal(object)
    satellite_ready = Signal(object)
    satellite_failed = Signal(object)

    def __init__(
        self,
        *,
        fetcher: Callable[..., CachedSatelliteElementSet] | None = None,
        projector: SatelliteProjector | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._fetcher = fetcher or fetch_cached_satellite_elements
        self._projector = projector or project_satellite_records
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
        enabled_groups: tuple[str, ...],
        reason: str = "manual",
    ) -> bool:
        request = {
            "observer_lat": float(observer_lat),
            "observer_lon": float(observer_lon),
            "observer_height_m": float(observer_height_m),
            "time_obj": time_obj,
            "enabled_groups": tuple(enabled_groups),
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

        self.satellite_started.emit({"banner": "Satellites: fetching orbital elements..."})
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
        enabled_groups: tuple[str, ...],
        reason: str,
        request_id: int,
    ) -> None:
        next_request: Optional[dict[str, object]] = None
        try:
            logger.info("Fetching satellite element sets (%s)...", reason)
            records_by_group: dict[str, list[SatelliteOmmRecord]] = {}
            refreshed_at_utc = datetime.now(timezone.utc)
            failed_groups: list[str] = []
            for group_key in enabled_groups:
                try:
                    fetched = self._fetcher(group_key)
                except Exception:
                    logger.warning("Satellite element fetch failed for %s", group_key, exc_info=True)
                    failed_groups.append(group_key)
                    continue
                records_by_group[str(group_key)] = list(fetched.records)
                if fetched.fetched_at_utc > refreshed_at_utc:
                    refreshed_at_utc = fetched.fetched_at_utc
            if not records_by_group:
                raise RuntimeError(
                    "Satellites: failed to fetch orbital elements"
                    if failed_groups
                    else "Satellites: no enabled groups"
                )
            overlay_points = self._projector(
                records_by_group,
                observer_lat=observer_lat,
                observer_lon=observer_lon,
                observer_height_m=observer_height_m,
                time_obj=time_obj,
            )
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                banner = ""
                if failed_groups:
                    banner = f"Satellites: partial ({', '.join(sorted(failed_groups))} unavailable)"
                self.satellite_ready.emit(
                    {
                        "records_by_group": records_by_group,
                        "overlay_points": overlay_points,
                        "refreshed_at_utc": refreshed_at_utc,
                        "banner": banner,
                    }
                )
        except Exception as exc:
            logger.warning("Satellite update failed: %s", exc, exc_info=True)
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                self.satellite_failed.emit({"banner": f"Satellites: {exc}"})
        finally:
            with self._lock:
                self._running = False
                if not self._stopping and self._pending_request is not None:
                    next_request = dict(self._pending_request)
                    self._pending_request = None
                    self._running = True
            if next_request is not None:
                self.satellite_started.emit({"banner": "Satellites: fetching orbital elements..."})
                worker = threading.Thread(target=self._run_update, kwargs=next_request, daemon=True)
                worker.start()
