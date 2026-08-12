"""End-to-end GMN loading through fixed celestial trail generation."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Protocol

from .constants import GMN_WINDOW
from .projection import project_meteor_observations_to_altaz
from .repository import GmnMeteorRepository
from .types import GmnLoadResult, MeteorWindowResult


class MeteorRepository(Protocol):
    def load_latest_window(
        self,
        display_time_utc: datetime,
        *,
        now_utc: datetime | None = None,
    ) -> GmnLoadResult: ...


def load_celestial_meteor_trails(
    display_time_utc: datetime,
    *,
    observer_lat: float,
    observer_lon: float,
    observer_height_m: float,
    repository: MeteorRepository | None = None,
    now_utc: datetime | None = None,
) -> MeteorWindowResult:
    display_time = _normalize_utc(display_time_utc)
    repo = repository or GmnMeteorRepository()
    loaded = repo.load_latest_window(display_time, now_utc=now_utc)
    window_end = loaded.window_end_utc
    if window_end is None:
        raise ValueError("GMN contains no observations in the latest search window")
    window_start = window_end - GMN_WINDOW
    trails = project_meteor_observations_to_altaz(
        loaded.observations,
        observer_lat=observer_lat,
        observer_lon=observer_lon,
        observer_height_m=observer_height_m,
    )
    return MeteorWindowResult(
        trails=trails,
        display_time_utc=display_time,
        window_start_utc=window_start,
        window_end_utc=window_end,
        source_files=loaded.source_files,
        unavailable_files=loaded.unavailable_files,
        used_stale_index=loaded.used_stale_index,
        used_stale_files=loaded.used_stale_files,
    )


def _normalize_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)
