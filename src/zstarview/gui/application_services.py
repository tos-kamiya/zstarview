from __future__ import annotations

import threading
from collections.abc import Callable
from concurrent.futures import Future, ThreadPoolExecutor
from typing import Any

from ..astro import load_ephemeris_uncached


class EphemerisProvider:
    """Application-owned lazy ephemeris provider."""

    def __init__(self, loader: Callable[[], Any] = load_ephemeris_uncached) -> None:
        self._loader = loader
        self._lock = threading.Lock()
        self._value: Any | None = None

    def load(self) -> Any:
        if self._value is not None:
            return self._value
        with self._lock:
            if self._value is None:
                self._value = self._loader()
        return self._value


class ApplicationServices:
    """Own shared background resources for one GUI application lifetime."""

    def __init__(
        self,
        *,
        max_workers: int = 2,
        ephemeris_provider: EphemerisProvider | None = None,
    ) -> None:
        self._worker_pool = ThreadPoolExecutor(
            max_workers=max(1, int(max_workers)),
            thread_name_prefix="zstarview-gui",
        )
        self.native_work_lock = threading.Lock()
        self.ephemeris = ephemeris_provider or EphemerisProvider()
        self._shutdown = False

    def submit(self, target, /, *args, **kwargs) -> Future:
        if self._shutdown:
            raise RuntimeError("application services are shut down")
        return self._worker_pool.submit(target, *args, **kwargs)

    def shutdown(self, *, wait: bool = True) -> None:
        if self._shutdown:
            return
        self._shutdown = True
        self._worker_pool.shutdown(wait=wait, cancel_futures=False)
