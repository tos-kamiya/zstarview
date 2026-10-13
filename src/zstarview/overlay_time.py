from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Literal

from .config import cache_retention_seconds

TimeMode = Literal["past", "present", "future"]

OVERLAY_PRESENT_TOLERANCE_SECONDS = 5 * 60
# Default/compatibility value; runtime checks use ``overlay_history_seconds``.
OVERLAY_HISTORY_SECONDS = 24 * 60 * 60


def overlay_history_seconds() -> int:
    """Return the configured past-time window for live overlay data."""
    return cache_retention_seconds()


def _within_past_history(delta_t: timedelta, history_seconds: int) -> bool:
    age_seconds = -delta_t.total_seconds()
    return 0.0 < age_seconds <= history_seconds


@dataclass(frozen=True)
class OverlayAvailability:
    cloud: bool
    aircraft: bool
    satellite: bool
    tropical_cyclone: bool
    precipitation: bool


def current_utc_time(*, now_utc: datetime | None = None) -> datetime:
    value = now_utc or datetime.now(timezone.utc)
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def target_time_utc_from_delta(
    delta_t: timedelta,
    *,
    now_utc: datetime | None = None,
) -> datetime:
    return current_utc_time(now_utc=now_utc) + delta_t


def overlay_source_time_utc_from_delta(
    delta_t: timedelta,
    *,
    now_utc: datetime | None = None,
    history_seconds: int | None = None,
) -> datetime:
    """Use current source data while displaying an eligible past sky time."""
    now = current_utc_time(now_utc=now_utc)
    target = now + delta_t
    if delta_t.total_seconds() < 0.0:
        retention = overlay_history_seconds() if history_seconds is None else history_seconds
        if _within_past_history(delta_t, retention):
            return now
    return target


def classify_target_time(
    target_time_utc: datetime,
    *,
    now_utc: datetime | None = None,
    present_tolerance_seconds: int = OVERLAY_PRESENT_TOLERANCE_SECONDS,
) -> TimeMode:
    now = current_utc_time(now_utc=now_utc)
    target = current_utc_time(now_utc=target_time_utc)
    delta_seconds = (target - now).total_seconds()
    if abs(delta_seconds) <= max(0, int(present_tolerance_seconds)):
        return "present"
    return "future" if delta_seconds > 0.0 else "past"


def classify_delta_t(
    delta_t: timedelta,
    *,
    now_utc: datetime | None = None,
    present_tolerance_seconds: int = OVERLAY_PRESENT_TOLERANCE_SECONDS,
) -> TimeMode:
    return classify_target_time(
        target_time_utc=target_time_utc_from_delta(delta_t, now_utc=now_utc),
        now_utc=now_utc,
        present_tolerance_seconds=present_tolerance_seconds,
    )


def overlay_availability_for_time_mode(time_mode: TimeMode) -> OverlayAvailability:
    if time_mode == "present":
        return OverlayAvailability(
            cloud=True,
            aircraft=True,
            satellite=True,
            tropical_cyclone=True,
            precipitation=True,
        )
    return OverlayAvailability(
        cloud=False,
        aircraft=False,
        satellite=False,
        tropical_cyclone=False,
        precipitation=False,
    )


def overlay_availability_for_delta(
    delta_t: timedelta,
    *,
    now_utc: datetime | None = None,
    present_tolerance_seconds: int = OVERLAY_PRESENT_TOLERANCE_SECONDS,
    allow_cached_history: bool = False,
    history_seconds: int | None = None,
) -> OverlayAvailability:
    now = current_utc_time(now_utc=now_utc)
    time_mode = classify_delta_t(
        delta_t,
        now_utc=now,
        present_tolerance_seconds=present_tolerance_seconds,
    )
    if time_mode == "present":
        return overlay_availability_for_time_mode(time_mode)
    if allow_cached_history and time_mode == "past":
        retention = overlay_history_seconds() if history_seconds is None else history_seconds
        if _within_past_history(delta_t, retention):
            return OverlayAvailability(
                cloud=True,
                aircraft=True,
                satellite=True,
                tropical_cyclone=False,
                precipitation=False,
            )
    return overlay_availability_for_time_mode(time_mode)
