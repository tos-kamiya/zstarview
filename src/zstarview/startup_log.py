from __future__ import annotations

import logging
import threading
from typing import Callable


class BufferedStartupLogHandler(logging.Handler):
    """Buffer startup log lines until a UI consumer is attached."""

    def __init__(self) -> None:
        super().__init__()
        self._lock = threading.Lock()
        self._consumer: Callable[[str, int], None] | None = None
        self._pending: list[tuple[str, int]] = []
        self.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))

    def set_consumer(self, consumer: Callable[[str, int], None] | None) -> None:
        pending: list[tuple[str, int]] = []
        with self._lock:
            self._consumer = consumer
            if consumer is not None and self._pending:
                pending = self._pending
                self._pending = []
        for line, levelno in pending:
            consumer(line, levelno)

    def emit(self, record: logging.LogRecord) -> None:
        try:
            line = self.format(record)
            with self._lock:
                consumer = self._consumer
                if consumer is None:
                    self._pending.append((line, int(record.levelno)))
                    return
            consumer(line, int(record.levelno))
        except Exception:
            self.handleError(record)
