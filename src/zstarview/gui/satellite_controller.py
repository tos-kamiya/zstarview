from __future__ import annotations

import inspect
import logging
import threading
import time
from collections.abc import Callable
from concurrent.futures import Future
from datetime import datetime, timezone
from urllib.error import URLError

import astropy.time
from PySide6.QtCore import QObject, Signal

from ..overlay_time import TimeMode, classify_target_time
from ..satellites import (
    CachedSatelliteElementSet,
    resolve_satellite_elements_for_time,
)
from ..satellites.fetch import SatelliteFetchCancelled
from ..satellites.types import SatelliteOmmRecord
from .worker_pool import submit_gui_work, wait_for_gui_futures

logger = logging.getLogger(__name__)

SatelliteFetcher = Callable[..., CachedSatelliteElementSet]
EXPECTED_FETCH_FAILURE_MESSAGES = {
    "Satellites: time-shifted view is not supported",
    "Horizons fetch returned no spacecraft records",
}


def _is_expected_fetch_failure(exc: Exception) -> bool:
    return str(exc) in EXPECTED_FETCH_FAILURE_MESSAGES


def _is_timeout_url_error(exc: Exception) -> bool:
    if isinstance(exc, URLError):
        reason = getattr(exc, "reason", None)
        if isinstance(reason, TimeoutError):
            return True
    return "timed out" in str(exc).lower()


def _is_expected_network_failure(exc: Exception) -> bool:
    if isinstance(exc, (ConnectionError, TimeoutError, URLError)):
        return True
    return isinstance(exc, OSError) and not isinstance(exc, FileNotFoundError)


def _call_fetcher_with_supported_kwargs(
    fetcher: SatelliteFetcher,
    group_key: str,
    **kwargs: object,
) -> CachedSatelliteElementSet:
    signature = inspect.signature(fetcher)
    if any(param.kind == inspect.Parameter.VAR_KEYWORD for param in signature.parameters.values()):
        return fetcher(group_key, **kwargs)
    supported_kwargs = {
        name: value
        for name, value in kwargs.items()
        if name in signature.parameters
    }
    return fetcher(group_key, **supported_kwargs)


class SatelliteController(QObject):
    satellite_started = Signal(object)
    satellite_ready = Signal(object)
    satellite_failed = Signal(object)

    def __init__(
        self,
        *,
        fetcher: Callable[..., CachedSatelliteElementSet] | None = None,
        projector: object | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._fetcher = fetcher or resolve_satellite_elements_for_time
        self._projector = projector
        self._running = False
        self._stopping = False
        self._pending_request: dict[str, object] | None = None
        self._latest_request_id = 0
        self._active_workers: set[Future[None]] = set()
        self._lock = threading.Lock()
        self._abort_event = threading.Event()

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        with self._lock:
            self._stopping = True
            self._pending_request = None
        self._abort_event.set()
        self._wait_for_workers(wait_timeout_s)

    def has_in_flight_update(self) -> bool:
        with self._lock:
            return bool(self._running or self._pending_request is not None or self._active_workers)

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
        self._spawn_worker(target=self._run_update, kwargs=request, label="satellite")
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
                    "Timed out waiting for %d satellite worker task(s) to finish during shutdown",
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
        time_obj: astropy.time.Time,
        enabled_groups: tuple[str, ...],
        reason: str,
        request_id: int,
    ) -> None:
        next_request: dict[str, object] | None = None
        try:
            logger.info("Fetching satellite element sets (%s)...", reason)
            records_by_group: dict[str, list[SatelliteOmmRecord]] = {}
            element_epoch_utc: datetime | None = None
            failed_groups: list[str] = []
            failure_messages: list[str] = []
            used_stale_cache = False
            target_time_utc = time_obj.to_datetime(timezone=timezone.utc)
            time_mode: TimeMode = classify_target_time(target_time_utc)

            def emit_partial_horizons_record(record: SatelliteOmmRecord) -> None:
                horizons_group_key = "horizons"
                horizons_records = list(records_by_group.get(horizons_group_key, []))
                record_name = str(record.get("OBJECT_NAME", "")).strip()
                if record_name:
                    horizons_records = [
                        existing
                        for existing in horizons_records
                        if str(existing.get("OBJECT_NAME", "")).strip() != record_name
                    ]
                horizons_records.append(dict(record))
                records_by_group[horizons_group_key] = horizons_records
                with self._lock:
                    should_emit_partial = (
                        not self._stopping and request_id == self._latest_request_id
                    )
                if should_emit_partial:
                    self.satellite_ready.emit(
                        {
                            "records_by_group": {
                                key: list(value)
                        for key, value in records_by_group.items()
                            },
                            "element_epoch_utc": target_time_utc,
                            "refreshed_at_utc": datetime.now(timezone.utc),
                            "banner": "Satellites: partial",
                        }
                    )

            for group_key in enabled_groups:
                try:
                    fetched = _call_fetcher_with_supported_kwargs(
                        self._fetcher,
                        group_key,
                        target_time_utc=target_time_utc,
                        time_mode=time_mode,
                        observer_lat=observer_lat,
                        observer_lon=observer_lon,
                        observer_height_m=observer_height_m,
                        horizons_record_callback=emit_partial_horizons_record
                        if group_key == "horizons"
                        else None,
                        abort_event=self._abort_event,
                    )
                except Exception as exc:
                    if isinstance(exc, SatelliteFetchCancelled):
                        with self._lock:
                            shutting_down = self._stopping
                        if shutting_down:
                            logger.debug("Satellite element fetch cancelled during shutdown")
                            return
                    if _is_expected_fetch_failure(exc):
                        logger.info("Satellite element fetch unavailable for %s: %s", group_key, exc)
                        if str(exc) == "Satellites: time-shifted view is not supported":
                            failure_messages.append(str(exc))
                    elif _is_timeout_url_error(exc):
                        logger.warning("Satellite element fetch timed out for %s: %s", group_key, exc)
                        failed_groups.append(group_key)
                        failure_messages.append("Satellites: timed out")
                        break
                    elif _is_expected_network_failure(exc):
                        logger.warning("Satellite element fetch unavailable for %s: %s", group_key, exc)
                    else:
                        logger.warning("Satellite element fetch failed for %s", group_key, exc_info=True)
                    failed_groups.append(group_key)
                    continue
                records_by_group[str(group_key)] = list(fetched.records)
                used_stale_cache = used_stale_cache or fetched.source in {"cache-stale", "cache-backoff"}
                if element_epoch_utc is None or fetched.element_epoch_utc > element_epoch_utc:
                    element_epoch_utc = fetched.element_epoch_utc
            if not records_by_group:
                if failed_groups:
                    banner = failure_messages[0] if failure_messages else "unavailable"
                    if not banner.startswith("Satellites:"):
                        banner = f"Satellites: {banner}"
                    with self._lock:
                        should_emit = not self._stopping and request_id == self._latest_request_id
                    if should_emit:
                        self.satellite_failed.emit({"banner": banner})
                    return
                raise RuntimeError("Satellites: no enabled groups")
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                banner = ""
                if failed_groups:
                    banner = "Satellites: partial"
                elif used_stale_cache:
                    banner = "Satellites: cache"
                self.satellite_ready.emit(
                    {
                        "records_by_group": records_by_group,
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
                banner = "Satellites: timed out" if _is_timeout_url_error(exc) else "Satellites: unavailable"
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
                self._spawn_worker(target=self._run_update, kwargs=next_request, label="satellite")
