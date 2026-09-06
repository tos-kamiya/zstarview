from __future__ import annotations

import logging
import json
import sys
import threading
import time
from collections.abc import Callable
from concurrent.futures import Future
from datetime import datetime, timedelta, timezone
from pathlib import Path

from PySide6.QtCore import QObject, QTimer, Signal

from ..paths import TROPICAL_CYCLONE_CACHE_DIR
from ..tropical_cyclones.cache import (
    TROPICAL_CYCLONE_CACHE_TTL_SECONDS,
    TROPICAL_CYCLONE_CACHE_VERSION,
    TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS,
    TropicalCycloneCacheEntry,
    is_tropical_cyclone_cache_current,
    is_tropical_cyclone_cache_stale,
    load_tropical_cyclone_cache,
    save_tropical_cyclone_cache,
)
from ..tropical_cyclones.client import (
    DEFAULT_SERVICE_URL,
    DEFAULT_TIMEOUT_S,
    DEFAULT_USER_AGENT,
    TropicalCycloneFetchError,
    fetch_active_hurricanes_snapshot,
    fetch_latest_observed_feature,
)
from ..tropical_cyclones.models import TropicalCycloneSnapshotCollection
from ..processes import JobRequest, ProcessJobSupervisor
from .application_services import ApplicationServices, wait_for_gui_futures

logger = logging.getLogger(__name__)
_EMPTY_OBSERVED_POSITION_MESSAGE = "No observed position features returned"


class TropicalCycloneController(QObject):
    cyclone_started = Signal(object)
    cyclone_ready = Signal(object)
    cyclone_failed = Signal(object)

    def __init__(
        self,
        *,
        service_url: str = DEFAULT_SERVICE_URL,
        cache_root: Path | str = TROPICAL_CYCLONE_CACHE_DIR,
        timeout_s: float = DEFAULT_TIMEOUT_S,
        user_agent: str = DEFAULT_USER_AGENT,
        services: ApplicationServices | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._owns_services = services is None
        self._services = services or ApplicationServices()
        self._service_url = str(service_url)
        self._cache_root = Path(cache_root)
        self._timeout_s = float(timeout_s)
        self._user_agent = str(user_agent)
        self._running = False
        self._stopping = False
        self._pending_request: dict[str, object] | None = None
        self._latest_request_id = 0
        self._active_workers: set[Future[None]] = set()
        self._lock = threading.Lock()
        self._process_supervisor = ProcessJobSupervisor(
            self._cache_root / ".process_jobs",
            max_restarts=2,
        )
        self._process_poll_timer = QTimer(self)
        self._process_poll_timer.setInterval(50)
        self._process_poll_timer.timeout.connect(self._poll_process_worker)

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        with self._lock:
            self._stopping = True
            self._pending_request = None
        self._process_poll_timer.stop()
        self._process_supervisor.close()
        self._wait_for_workers(wait_timeout_s)
        if self._owns_services:
            self._services.shutdown(wait=True)

    def has_in_flight_update(self) -> bool:
        with self._lock:
            return bool(self._running or self._pending_request is not None or self._active_workers)

    def update(self, *, reason: str = "manual") -> bool:
        request = {"reason": str(reason)}
        with self._lock:
            if self._stopping:
                return False
            self._latest_request_id += 1
            request["request_id"] = int(self._latest_request_id)
            if self._running:
                self._pending_request = dict(request)
                started = False
            else:
                self._running = True
                started = True

        process_request = JobRequest(
            session_id=self._process_supervisor.session_id,
            worker_epoch=1,
            request_id=int(request["request_id"]),
            job_kind="tropical-cyclone-fetch",
            layer_generation=1,
            view_generation=0,
            input_revision=str(TROPICAL_CYCLONE_CACHE_VERSION),
            payload={
                "reason": str(reason),
                "service_url": self._service_url,
                "cache_root": str(self._cache_root),
                "timeout_s": self._timeout_s,
                "user_agent": self._user_agent,
            },
        )
        try:
            self._process_supervisor.submit(
                process_request,
                (sys.executable, "-m", "zstarview.tropical_cyclones.worker"),
                timeout_s=max(1.0, self._timeout_s + 10.0),
            )
        except Exception as exc:
            logger.warning("Failed to launch tropical cyclone worker: %s", exc)
            with self._lock:
                self._running = False
            self._emit_failed("Typhoon: unavailable", request_id=int(request["request_id"]))
            return False
        self._process_poll_timer.start()
        if started:
            self.cyclone_started.emit({"banner": "Typhoon: checking..."})
        return started

    def _poll_process_worker(self) -> None:
        events = self._process_supervisor.poll()
        for event in events:
            try:
                request_id = event.request.request_id
                if event.failure is not None:
                    logger.warning(
                        "Tropical cyclone process failed (%s): %s",
                        event.failure.kind,
                        event.failure,
                    )
                    self._emit_failed("Typhoon: unavailable", request_id=request_id)
                elif event.result is None or event.result.status != "ok":
                    message = "worker returned an unsuccessful result"
                    if event.result is not None and event.result.error_message:
                        message = event.result.error_message
                    logger.warning("Tropical cyclone process failed: %s", message)
                    self._emit_failed("Typhoon: unavailable", request_id=request_id)
                else:
                    try:
                        payload = self._load_process_payload(event)
                    except Exception as exc:
                        logger.warning("Invalid tropical cyclone worker payload: %s", exc)
                        self._emit_failed("Typhoon: unavailable", request_id=request_id)
                    else:
                        self._emit_ready(payload, request_id=request_id)
            finally:
                self._process_supervisor.release(event)
        if not self._process_supervisor.has_pending_work:
            with self._lock:
                self._running = False
                self._pending_request = None
            self._process_poll_timer.stop()
        else:
            with self._lock:
                had_pending = self._pending_request is not None
                self._pending_request = None
            if events and had_pending:
                self.cyclone_started.emit({"banner": "Typhoon: checking..."})

    def _load_process_payload(self, event) -> dict[str, object]:
        if event.result is None or len(event.result.artifacts) != 1:
            raise ValueError("worker payload artifact is missing")
        artifact = event.result.artifacts[0]
        if event.directory is None:
            raise ValueError("worker job directory is missing")
        job_dir = event.directory
        payload_path = (job_dir / artifact.relative_path).resolve()
        raw = json.loads(payload_path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ValueError("worker payload must be an object")
        collection_raw = raw.get("snapshot_collection")
        if not isinstance(collection_raw, dict):
            raise ValueError("worker payload is missing snapshot collection")
        collection = TropicalCycloneSnapshotCollection.from_dict(collection_raw)
        if collection is None:
            raise ValueError("worker snapshot collection is invalid")

        def parse_datetime(value: object) -> datetime:
            if not isinstance(value, str):
                raise ValueError("worker payload datetime is missing")
            text = value[:-1] + "+00:00" if value.endswith("Z") else value
            parsed = datetime.fromisoformat(text)
            if parsed.tzinfo is None:
                parsed = parsed.replace(tzinfo=timezone.utc)
            return parsed.astimezone(timezone.utc)

        return {
            "snapshot_collection": collection.to_dict(),
            "cached_at_utc": parse_datetime(raw["cached_at_utc"]),
            "last_checked_utc": parse_datetime(raw["last_checked_utc"]),
            "next_check_utc": parse_datetime(raw["next_check_utc"]),
            "next_refresh_utc": parse_datetime(raw["next_refresh_utc"]),
            "banner": str(raw.get("banner", "")),
            "service_url": str(raw.get("service_url", self._service_url)),
        }

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
                    "Timed out waiting for %d tropical-cyclone worker task(s) to finish during shutdown",
                    len(workers),
                )
                return
            wait_for_gui_futures(workers, remaining)

    def _cache_payload(
        self,
        snapshot_collection: TropicalCycloneSnapshotCollection,
        *,
        cached_at_utc: datetime,
        last_checked_utc: datetime,
        next_check_utc: datetime,
        next_refresh_utc: datetime,
        banner: str = "",
    ) -> dict[str, object]:
        return {
            "snapshot_collection": snapshot_collection.to_dict(),
            "cached_at_utc": cached_at_utc,
            "last_checked_utc": last_checked_utc,
            "next_check_utc": next_check_utc,
            "next_refresh_utc": next_refresh_utc,
            "banner": banner,
            "service_url": self._service_url,
        }

    def _payload_from_cache_entry(
        self,
        entry: TropicalCycloneCacheEntry,
        *,
        now_utc: datetime,
        next_check_utc: datetime | None = None,
        banner: str = "",
    ) -> dict[str, object]:
        cached_at = entry.cached_at_utc
        next_check = (
            next_check_utc
            if next_check_utc is not None
            else max(
                now_utc,
                cached_at + timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
            )
        )
        next_refresh = cached_at + timedelta(seconds=TROPICAL_CYCLONE_CACHE_TTL_SECONDS)
        return self._cache_payload(
            entry.snapshot_collection,
            cached_at_utc=cached_at,
            last_checked_utc=now_utc,
            next_check_utc=next_check,
            next_refresh_utc=next_refresh,
            banner=banner,
        )

    def _emit_ready(
        self,
        payload: dict[str, object],
        *,
        request_id: int,
    ) -> None:
        with self._lock:
            should_emit = not self._stopping and request_id == self._latest_request_id
        if should_emit:
            self.cyclone_ready.emit(payload)

    def _emit_failed(self, banner: str, *, request_id: int) -> None:
        with self._lock:
            should_emit = not self._stopping and request_id == self._latest_request_id
        if should_emit:
            self.cyclone_failed.emit({"banner": banner})

    def _emit_empty_overlay(
        self,
        *,
        request_id: int,
        now: datetime,
        banner: str = "Typhoon: none",
    ) -> None:
        empty_collection = TropicalCycloneSnapshotCollection(
            snapshots=(),
            source_url=self._service_url,
            service_name="",
            refreshed_at_utc=now,
        )
        cached_at = datetime.now(timezone.utc)
        entry = TropicalCycloneCacheEntry(
            snapshot_collection=empty_collection,
            cached_at_utc=cached_at,
            cache_version=TROPICAL_CYCLONE_CACHE_VERSION,
        )
        try:
            save_tropical_cyclone_cache(entry, cache_root=self._cache_root)
        except Exception:
            logger.warning("Failed to write tropical cyclone cache", exc_info=True)
        payload = self._cache_payload(
            empty_collection,
            cached_at_utc=cached_at,
            last_checked_utc=cached_at,
            next_check_utc=cached_at + timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
            next_refresh_utc=cached_at + timedelta(seconds=TROPICAL_CYCLONE_CACHE_TTL_SECONDS),
            banner=banner,
        )
        self._emit_ready(payload, request_id=request_id)

    def _run_update(self, *, reason: str, request_id: int) -> None:
        next_request: dict[str, object] | None = None
        try:
            logger.info("Updating tropical cyclone overlay (%s)...", reason)
            now = datetime.now(timezone.utc)
            cached_entry = load_tropical_cyclone_cache(self._cache_root)
            cached_is_stale = False
            if cached_entry is not None:
                cached_is_stale = is_tropical_cyclone_cache_stale(cached_entry, now_utc=now)
            if (
                cached_entry is not None
                and not cached_is_stale
                and is_tropical_cyclone_cache_current(cached_entry)
                and now < cached_entry.cached_at_utc + timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS)
            ):
                self._emit_ready(
                    self._payload_from_cache_entry(
                        cached_entry,
                        now_utc=now,
                        next_check_utc=cached_entry.cached_at_utc
                        + timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
                    ),
                    request_id=request_id,
                )
                return

            latest_feature = fetch_latest_observed_feature(
                service_url=self._service_url,
                timeout_s=self._timeout_s,
                user_agent=self._user_agent,
            )
            if latest_feature is None:
                logger.warning("No observed tropical cyclone positions returned; treating overlay as empty.")
                self._emit_empty_overlay(request_id=request_id, now=now)
                return
            latest_attrs = latest_feature.get("attributes")
            if not isinstance(latest_attrs, dict):
                raise TropicalCycloneFetchError("Observed position payload missing attributes")
            latest_storm_name = latest_attrs.get("STORMNAME")
            latest_basin = latest_attrs.get("BASIN")
            latest_advdate = latest_attrs.get("ADVDATE")
            latest_advdate_int = int(latest_advdate) if isinstance(latest_advdate, (int, float)) else None

            cached_collection = cached_entry.snapshot_collection if cached_entry is not None else None
            cached_snapshot = (
                cached_collection.snapshots[0]
                if cached_collection is not None and cached_collection.snapshots
                else None
            )
            if (
                cached_snapshot is not None
                and not cached_is_stale
                and cached_entry is not None
                and is_tropical_cyclone_cache_current(cached_entry)
                and cached_snapshot.has_projectable_timeline()
                and isinstance(latest_storm_name, str)
                and latest_storm_name == cached_snapshot.storm_name
                and (
                    (latest_basin is None and cached_snapshot.basin is None)
                    or (
                        isinstance(latest_basin, str)
                        and latest_basin == cached_snapshot.basin
                    )
                )
                and latest_advdate_int is not None
                and cached_snapshot.advdate_utc is not None
            ):
                cached_advdate_ms = int(cached_snapshot.advdate_utc.timestamp() * 1000.0)
                if cached_advdate_ms == latest_advdate_int:
                    self._emit_ready(
                        self._payload_from_cache_entry(
                            cached_entry,
                            now_utc=now,
                            next_check_utc=now
                            + timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
                        ),
                        request_id=request_id,
                    )
                    return

            snapshot = fetch_active_hurricanes_snapshot(
                service_url=self._service_url,
                timeout_s=self._timeout_s,
                user_agent=self._user_agent,
            )
            cached_at = datetime.now(timezone.utc)
            entry = TropicalCycloneCacheEntry(
                snapshot_collection=snapshot,
                cached_at_utc=cached_at,
                cache_version=TROPICAL_CYCLONE_CACHE_VERSION,
            )
            try:
                save_tropical_cyclone_cache(entry, cache_root=self._cache_root)
            except Exception:
                logger.warning("Failed to write tropical cyclone cache", exc_info=True)
            payload = self._cache_payload(
                snapshot,
                cached_at_utc=cached_at,
                last_checked_utc=cached_at,
                next_check_utc=cached_at + timedelta(seconds=TROPICAL_CYCLONE_CHECK_INTERVAL_SECONDS),
                next_refresh_utc=cached_at + timedelta(seconds=TROPICAL_CYCLONE_CACHE_TTL_SECONDS),
            )
            self._emit_ready(payload, request_id=request_id)
        except Exception as exc:
            if isinstance(exc, TropicalCycloneFetchError) and str(exc) == _EMPTY_OBSERVED_POSITION_MESSAGE:
                logger.warning("No observed tropical cyclone positions returned; treating overlay as empty.")
                self._emit_empty_overlay(request_id=request_id, now=datetime.now(timezone.utc))
                return
            logger.warning("Tropical cyclone update failed: %s", exc)
            self._emit_failed("Typhoon: unavailable", request_id=request_id)
        finally:
            with self._lock:
                self._running = False
                if not self._stopping and self._pending_request is not None:
                    next_request = dict(self._pending_request)
                    self._pending_request = None
                    self._running = True
            if next_request is not None:
                self.cyclone_started.emit({"banner": "Typhoon: checking..."})
                self._spawn_worker(target=self._run_update, kwargs=next_request)
