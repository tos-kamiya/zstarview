from __future__ import annotations

import os
import sys
from datetime import datetime, timezone
from threading import Lock

_DEBUG_TIMING_ENV = "ZSTARVIEW_DEBUG_EARTH_GUIDE_TIMING"
_TIMING_LOCK = Lock()
_TIMING_COUNTER = 0


def earth_guide_timing_enabled() -> bool:
    raw = os.getenv(_DEBUG_TIMING_ENV, "").strip().casefold()
    return raw in {"1", "true", "yes", "on"}


def log_earth_guide_timing(source: str, message: str) -> None:
    if not earth_guide_timing_enabled():
        return
    global _TIMING_COUNTER
    with _TIMING_LOCK:
        _TIMING_COUNTER += 1
        seq = _TIMING_COUNTER
    timestamp = datetime.now(timezone.utc).strftime("%H:%M:%S.%f")[:-3]
    print(
        f"[earth-guide-timing {seq:04d} {timestamp}] {source}: {message}",
        file=sys.stderr,
        flush=True,
    )
