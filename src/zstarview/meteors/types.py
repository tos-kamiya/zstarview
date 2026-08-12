"""Data structures for GMN meteor observations."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime


@dataclass(frozen=True, slots=True)
class MeteorObservation:
    trajectory_id: str
    beginning_utc: datetime
    begin_lat_deg: float
    begin_lon_deg: float
    begin_height_km: float
    end_lat_deg: float
    end_lon_deg: float
    end_height_km: float
    duration_s: float | None = None
    peak_abs_magnitude: float | None = None
    initial_speed_km_s: float | None = None
    shower_code: str | None = None


@dataclass(frozen=True, slots=True)
class MeteorTrail:
    trajectory_id: str
    beginning_utc: datetime
    begin_alt_deg: float
    begin_az_deg: float
    end_alt_deg: float
    end_az_deg: float
    duration_s: float | None = None
    peak_abs_magnitude: float | None = None
    shower_code: str | None = None


CelestialMeteorTrail = MeteorTrail


@dataclass(frozen=True, slots=True)
class GmnLoadResult:
    observations: tuple[MeteorObservation, ...]
    source_files: tuple[str, ...]
    unavailable_files: tuple[str, ...]
    used_stale_index: bool = False
    used_stale_files: bool = False
    window_end_utc: datetime | None = None


@dataclass(frozen=True, slots=True)
class MeteorWindowResult:
    trails: tuple[MeteorTrail, ...]
    display_time_utc: datetime
    window_start_utc: datetime
    window_end_utc: datetime
    source_files: tuple[str, ...]
    unavailable_files: tuple[str, ...]
    used_stale_index: bool = False
    used_stale_files: bool = False


@dataclass(frozen=True, slots=True)
class GmnDailyFile:
    filename: str
    nominal_date: date
