from __future__ import annotations

import logging
import json
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path

from PySide6.QtCore import QObject, QTimer, Signal

from ..paths import TROPICAL_CYCLONE_CACHE_DIR
from ..tropical_cyclones.cache import TROPICAL_CYCLONE_CACHE_VERSION
from ..tropical_cyclones.client import DEFAULT_SERVICE_URL, DEFAULT_TIMEOUT_S, DEFAULT_USER_AGENT
from ..tropical_cyclones.models import TropicalCycloneSnapshotCollection
from ..processes import JobRequest, ProcessJobSupervisor
from .application_services import ApplicationServices

logger = logging.getLogger(__name__)
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
        # Keep the argument for caller compatibility; this controller no longer
        # submits Python work to the shared GUI executor.
        del services
        self._service_url = str(service_url)
        self._cache_root = Path(cache_root)
        self._timeout_s = float(timeout_s)
        self._user_agent = str(user_agent)
        self._running = False
        self._stopping = False
        self._pending_request: dict[str, object] | None = None
        self._latest_request_id = 0
        self._lock = threading.Lock()
        self._process_supervisor = ProcessJobSupervisor(
            self._cache_root / ".process_jobs",
            max_restarts=2,
        )
        self._process_poll_timer = QTimer(self)
        self._process_poll_timer.setInterval(50)
        self._process_poll_timer.timeout.connect(self._poll_process_worker)

    def shutdown(self, *, wait_timeout_s: float | None = None) -> None:
        del wait_timeout_s
        with self._lock:
            self._stopping = True
            self._pending_request = None
        self._process_poll_timer.stop()
        self._process_supervisor.close()

    def has_in_flight_update(self) -> bool:
        with self._lock:
            return bool(self._running or self._pending_request is not None)

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
