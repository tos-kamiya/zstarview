from __future__ import annotations

import threading
from concurrent.futures import Future
from datetime import datetime, timezone

from PySide6.QtCore import QObject, Signal

from ..precipitation import (
    PrecipitationSnapshot,
    fetch_open_meteo_precipitation,
    generate_precipitation_request_samples,
    precipitation_cache_key,
    precipitation_snapshot_is_fresh,
    project_precipitation_columns,
)
from ..types import ViewerData


class PrecipitationController(QObject):
    precipitation_started = Signal(object)
    precipitation_ready = Signal(object)
    precipitation_failed = Signal(object)

    def __init__(self, *, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._future: Future[None] | None = None
        self._stopping = False
        self._cache_key: tuple[object, ...] | None = None
        self._snapshot: PrecipitationSnapshot | None = None

    def shutdown(self) -> None:
        self._stopping = True

    def has_in_flight_update(self) -> bool:
        return self._future is not None and not self._future.done()

    def update(self, *, viewer_data: ViewerData, reason: str = "manual") -> bool:
        if self._stopping or self.has_in_flight_update():
            return False
        samples = generate_precipitation_request_samples(
            viewer_data.lat_deg, viewer_data.lon_deg
        )
        key = precipitation_cache_key(samples)
        now = datetime.now(timezone.utc)
        if (
            key == self._cache_key
            and self._snapshot is not None
            and precipitation_snapshot_is_fresh(self._snapshot, now_utc=now)
        ):
            return False
        self.precipitation_started.emit({"reason": reason})
        future: Future[None] = Future()
        self._future = future
        thread = threading.Thread(
            target=self._run,
            args=(future, viewer_data, samples, key),
            name="precipitation",
            daemon=True,
        )
        thread.start()
        future.add_done_callback(self._finished)
        return True

    def _run(self, future, viewer_data, samples, key) -> None:
        try:
            snapshot = fetch_open_meteo_precipitation(samples)
            columns = project_precipitation_columns(snapshot, viewer_data)
            self._cache_key = key
            self._snapshot = snapshot
            first_value = snapshot.values[0]
            if not self._stopping:
                self.precipitation_ready.emit(
                    {
                        "columns": list(columns),
                        "forecast_time_utc": first_value.forecast_time_utc,
                        "interval_seconds": first_value.interval_seconds,
                        "fetched_at_utc": snapshot.fetched_at_utc,
                    }
                )
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            future.set_exception(exc)
        else:
            future.set_result(None)

    def _finished(self, future: Future[None]) -> None:
        self._future = None
        if self._stopping:
            return
        try:
            future.result()
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            self._snapshot = None
            self.precipitation_failed.emit({"error": str(exc)})
