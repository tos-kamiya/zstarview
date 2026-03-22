from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Literal

TimeMode = Literal["past", "present", "future"]

OVERLAY_PRESENT_TOLERANCE_SECONDS = 5 * 60


@dataclass(frozen=True)
class OverlayAvailability:
    cloud: bool
    aircraft: bool
    satellite: bool


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
        return OverlayAvailability(cloud=True, aircraft=True, satellite=True)
    return OverlayAvailability(cloud=False, aircraft=False, satellite=False)


def overlay_availability_for_delta(
    delta_t: timedelta,
    *,
    now_utc: datetime | None = None,
    present_tolerance_seconds: int = OVERLAY_PRESENT_TOLERANCE_SECONDS,
) -> OverlayAvailability:
    return overlay_availability_for_time_mode(
        classify_delta_t(
            delta_t,
            now_utc=now_utc,
            present_tolerance_seconds=present_tolerance_seconds,
        )
    )
