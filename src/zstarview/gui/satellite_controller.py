from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from typing import Callable, Optional

import astropy.time
from PySide6.QtCore import QObject, Signal

from ..overlay_time import TimeMode, classify_target_time
from ..satellites import (
    CachedSatelliteElementSet,
    resolve_satellite_elements_for_time,
    project_satellite_records,
)
from ..satellites.types import SatelliteOmmRecord, SatelliteOverlayPoint

logger = logging.getLogger(__name__)

SatelliteFetcher = Callable[..., CachedSatelliteElementSet]
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
        self._fetcher = fetcher or resolve_satellite_elements_for_time
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
            element_epoch_utc: datetime | None = None
            failed_groups: list[str] = []
            target_time_utc = time_obj.to_datetime(timezone=timezone.utc)
            time_mode: TimeMode = classify_target_time(target_time_utc)
            for group_key in enabled_groups:
                try:
                    fetched = self._fetcher(
                        group_key,
                        target_time_utc=target_time_utc,
                        time_mode=time_mode,
                    )
                except Exception:
                    logger.warning("Satellite element fetch failed for %s", group_key, exc_info=True)
                    failed_groups.append(group_key)
                    continue
                records_by_group[str(group_key)] = list(fetched.records)
                if element_epoch_utc is None or fetched.element_epoch_utc > element_epoch_utc:
                    element_epoch_utc = fetched.element_epoch_utc
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
                        "element_epoch_utc": element_epoch_utc or target_time_utc,
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
