"""Keep cyclic Qt/Python object finalization on the GUI's owning thread."""

from __future__ import annotations

import gc
import threading
import time
from types import TracebackType

from PySide6.QtCore import QCoreApplication, QThread, QTimer


class GuiGarbageCollector:
    """Replace automatic cyclic GC with bounded-frequency GUI-thread callbacks.

    Reference counting stays active. A full collection periodically visits old
    cycles too. Native-work contention skips a tick rather than blocking Qt.
    This controls GC execution, not explicit finalizers or third-party collects.
    """

    def __init__(self, native_work_lock: threading.Lock) -> None:
        self._native_work_lock = native_work_lock
        self._owner_ident = threading.get_ident()
        self._was_enabled = False
        self._active = False
        self._last_full_collection = time.monotonic()
        self._timer: QTimer | None = None

    def __enter__(self) -> GuiGarbageCollector:
        self._require_owner_thread()
        if self._active:
            raise RuntimeError("GUI garbage collector is already active")
        self._was_enabled = gc.isenabled()
        self._active = True
        gc.disable()
        return self

    def start_timer(self, app: QCoreApplication) -> None:
        self._require_owner_thread()
        if app.thread() != QThread.currentThread():
            raise RuntimeError("GUI garbage collection must attach on the application thread")
        if not self._active or not self._was_enabled or self._timer is not None:
            return
        self._timer = QTimer(app)
        self._timer.setInterval(1000)
        self._timer.timeout.connect(self.collect_if_due)
        self._timer.start()

    def _require_owner_thread(self) -> None:
        if threading.get_ident() != self._owner_ident:
            raise RuntimeError("GUI garbage collection must run on its owner thread")

    def collect_if_due(self) -> None:
        self._require_owner_thread()
        if not self._active or not self._was_enabled:
            return
        now = time.monotonic()
        full_due = now - self._last_full_collection >= 30.0
        if not full_due and gc.get_count()[0] < gc.get_threshold()[0]:
            return
        if not self._native_work_lock.acquire(blocking=False):
            return
        try:
            gc.collect(2 if full_due else 0)
            if full_due:
                self._last_full_collection = time.monotonic()
        finally:
            self._native_work_lock.release()

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self._require_owner_thread()
        self._active = False
        try:
            if self._timer is not None:
                self._timer.stop()
                self._timer.timeout.disconnect(self.collect_if_due)
                self._timer = None
        finally:
            if self._was_enabled:
                gc.enable()
            else:
                gc.disable()
