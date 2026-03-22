from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Any

SatelliteOmmRecord = dict[str, Any]


@dataclass(frozen=True)
class CachedSatelliteElementSet:
    group_key: str
    fetched_at_utc: datetime
    source: str
    records: list[SatelliteOmmRecord]
