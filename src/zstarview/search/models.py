from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime


@dataclass(frozen=True)
class SearchJumpTarget:
    label: str
    kind: str
    sort_key: tuple[float, str]
    ra_hours: float = 0.0
    dec_deg: float = 0.0
    subtitle: str = ""
    object_key: str = ""
    latitude_deg: float | None = None
    longitude_deg: float | None = None
    command: str = ""
    alt_deg: float | None = None
    az_deg: float | None = None
    target_time_utc: datetime | None = None
    jpl_group: str = ""
    persistent_keep_marker: bool = False


@dataclass(frozen=True)
class SearchRequest:
    query: str
    list_only: bool = False
    keep_marker: bool = False


@dataclass(frozen=True)
class SearchResolution:
    query: str
    candidates: tuple[SearchJumpTarget, ...]
    selected_target: SearchJumpTarget | None = None
    status: str = ""

