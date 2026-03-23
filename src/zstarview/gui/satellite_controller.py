from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from typing import Callable, Optional
from urllib.error import URLError

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
EXPECTED_FETCH_FAILURE_MESSAGES = {
    "Satellites: time-shifted view is not supported",
}


def _is_expected_fetch_failure(exc: Exception) -> bool:
    return str(exc) in EXPECTED_FETCH_FAILURE_MESSAGES


def _is_timeout_url_error(exc: Exception) -> bool:
    if isinstance(exc, URLError):
        reason = getattr(exc, "reason", None)
        if isinstance(reason, TimeoutError):
            return True
    return "timed out" in str(exc).lower()


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
            failure_messages: list[str] = []
            target_time_utc = time_obj.to_datetime(timezone=timezone.utc)
            time_mode: TimeMode = classify_target_time(target_time_utc)
            for group_key in enabled_groups:
                try:
                    fetched = self._fetcher(
                        group_key,
                        target_time_utc=target_time_utc,
                        time_mode=time_mode,
                    )
                except Exception as exc:
                    if _is_expected_fetch_failure(exc):
                        logger.info("Satellite element fetch unavailable for %s: %s", group_key, exc)
                    elif _is_timeout_url_error(exc):
                        logger.warning("Satellite element fetch timed out for %s: %s", group_key, exc)
                        failed_groups.append(group_key)
                        failure_messages.append(str(exc))
                        break
                    else:
                        logger.warning("Satellite element fetch failed for %s", group_key, exc_info=True)
                    failed_groups.append(group_key)
                    failure_messages.append(str(exc))
                    continue
                records_by_group[str(group_key)] = list(fetched.records)
                if element_epoch_utc is None or fetched.element_epoch_utc > element_epoch_utc:
                    element_epoch_utc = fetched.element_epoch_utc
            if not records_by_group:
                if failed_groups:
                    unique_messages = []
                    for message in failure_messages:
                        if message not in unique_messages:
                            unique_messages.append(message)
                    if len(unique_messages) == 1:
                        raise RuntimeError(unique_messages[0])
                    raise RuntimeError("Satellites: failed to fetch orbital elements")
                raise RuntimeError("Satellites: no enabled groups")
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
                        "refreshed_at_utc": datetime.now(timezone.utc),
                        "banner": banner,
                    }
                )
        except Exception as exc:
            if _is_timeout_url_error(exc):
                logger.warning("Satellite update failed: %s", exc)
            else:
                logger.warning("Satellite update failed: %s", exc, exc_info=True)
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                message = str(exc).strip()
                banner = message if message.startswith("Satellites:") else f"Satellites: {message}"
                self.satellite_failed.emit({"banner": banner})
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
