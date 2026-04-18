from __future__ import annotations

import logging
import threading
from datetime import datetime, timezone
from typing import Optional
from urllib.error import URLError

from PySide6.QtCore import QObject, Signal

from ..satellites import fetch_horizons_observer_csv
from .famous_star_shortcuts import SearchJumpTarget

logger = logging.getLogger(__name__)


class JplSmallBodyController(QObject):
    jpl_started = Signal(object)
    jpl_ready = Signal(object)
    jpl_failed = Signal(object)

    def __init__(self, *, parent: QObject | None = None) -> None:
        super().__init__(parent)
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
        target: SearchJumpTarget,
        target_time_utc: datetime,
        reason: str = "manual",
    ) -> bool:
        request = {
            "observer_lat": float(observer_lat),
            "observer_lon": float(observer_lon),
            "observer_height_m": float(observer_height_m),
            "target": target,
            "target_time_utc": target_time_utc.astimezone(timezone.utc),
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

        self.jpl_started.emit({"banner": "JPL: fetching small-body ephemeris..."})
        worker = threading.Thread(target=self._run_update, kwargs=request, daemon=True)
        worker.start()
        return True

    def _run_update(
        self,
        *,
        observer_lat: float,
        observer_lon: float,
        observer_height_m: float,
        target: SearchJumpTarget,
        target_time_utc: datetime,
        reason: str,
        request_id: int,
    ) -> None:
        next_request: Optional[dict[str, object]] = None
        try:
            command = str(target.command).strip()
            if not command:
                command = f"DES={target.object_key};" if target.object_key else ""
            if not command:
                raise RuntimeError("JPL small-body target has no usable command")
            logger.info("Fetching JPL small-body ephemeris (%s)...", reason)
            rows = fetch_horizons_observer_csv(
                command,
                target_time_utc=target_time_utc,
                observer_lat=observer_lat,
                observer_lon=observer_lon,
                observer_height_m=observer_height_m,
            )
            alt_az = self._extract_altaz(rows)
            if alt_az is None:
                raise RuntimeError("JPL observer table did not contain an alt/az sample")
            alt_deg, az_deg = alt_az
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                self.jpl_ready.emit(
                    {
                        "target": target,
                        "target_time_utc": target_time_utc.astimezone(timezone.utc),
                        "refreshed_at_utc": datetime.now(timezone.utc),
                        "alt_deg": float(alt_deg),
                        "az_deg": float(az_deg) % 360.0,
                        "rows": rows,
                        "reason": reason,
                    }
                )
        except Exception as exc:
            if isinstance(exc, URLError):
                logger.warning("JPL small-body update failed: %s", exc)
            else:
                logger.warning("JPL small-body update failed: %s", exc, exc_info=True)
            with self._lock:
                should_emit = not self._stopping and request_id == self._latest_request_id
            if should_emit:
                self.jpl_failed.emit(
                    {
                        "target": target,
                        "target_time_utc": target_time_utc.astimezone(timezone.utc),
                        "refreshed_at_utc": datetime.now(timezone.utc),
                        "banner": f"JPL: {exc}",
                        "error": str(exc),
                        "reason": reason,
                    }
                )
        finally:
            with self._lock:
                self._running = False
                if not self._stopping and self._pending_request is not None:
                    next_request = dict(self._pending_request)
                    self._pending_request = None
                    self._running = True
            if next_request is not None:
                self.jpl_started.emit({"banner": "JPL: fetching small-body ephemeris..."})
                worker = threading.Thread(target=self._run_update, kwargs=next_request, daemon=True)
                worker.start()

    def _extract_altaz(self, rows: list[list[str]]) -> tuple[float, float] | None:
        for row in rows:
            numeric_values: list[float] = []
            for value in row:
                try:
                    numeric_values.append(float(str(value).strip()))
                except (TypeError, ValueError):
                    continue
            if len(numeric_values) >= 2:
                return numeric_values[-1], numeric_values[-2]
        return None
