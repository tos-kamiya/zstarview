"""Parse launch-time inputs without handling location resolution."""

from __future__ import annotations

import logging
import re
from datetime import datetime, timedelta, timezone
from typing import Optional, Tuple

from .utils.timezone_parser import parse_tz_string

logger = logging.getLogger(__name__)


class LaunchSetupError(Exception):
    """Abort the startup sequence (handled by callers to show concise UI errors)."""


def _parse_flexible_time(time_str: str) -> Tuple[int, int, int]:
    """Parse time string that may omit minutes/seconds."""
    match = re.fullmatch(r"\s*(\d{1,2})(?::(\d{1,2}))?(?::(\d{1,2}))?\s*", time_str)
    if not match:
        raise ValueError(f"Invalid time: {time_str!r}. Use HH, HH:MM, or HH:MM:SS.")
    hour = int(match.group(1))
    minute = int(match.group(2)) if match.group(2) is not None else 0
    second = int(match.group(3)) if match.group(3) is not None else 0

    if not (0 <= hour <= 23):
        raise ValueError(f"Hour out of range (0-23): {hour}")
    if not (0 <= minute <= 59):
        raise ValueError(f"Minute out of range (0-59): {minute}")
    if not (0 <= second <= 59):
        raise ValueError(f"Second out of range (0-59): {second}")
    return hour, minute, second


def parse_launch_time_arguments(
    args_datetime: Optional[str],
    args_days: float,
    args_hours: float,
    *,
    timezone_name: str = "UTC",
    timezone_override: str | None = None,
) -> timedelta:
    """
    Parse time-related arguments and return a timedelta from now.

    Supports --datetime in 'YYYY-MM-DD HH[:MM[:SS]] [TZ]' (TZ optional).
    """
    if not args_datetime:
        return timedelta(days=args_days, hours=args_hours)

    if args_hours != 0 or args_days != 0:
        logger.error("Invalid option: --datetime cannot be used with --hours or --days.")
        raise LaunchSetupError()

    try:
        parts = args_datetime.split()
        if len(parts) < 2:
            raise ValueError("Missing time part. Use 'YYYY-MM-DD HH[:MM[:SS]] [TZ]'.")

        date_str = parts[0]
        time_str = parts[1]
        tz_str = parts[2] if len(parts) >= 3 else None

        try:
            date_only = datetime.strptime(date_str, "%Y-%m-%d")
        except ValueError as exc:
            raise ValueError(f"Invalid date: {date_str!r}. Use YYYY-MM-DD.") from exc

        hour, minute, second = _parse_flexible_time(time_str)
        dt_naive = date_only.replace(hour=hour, minute=minute, second=second)

        effective_tz_str = timezone_override if timezone_override is not None else tz_str

        if effective_tz_str:
            try:
                tz = parse_tz_string(effective_tz_str)
                dt_local = dt_naive.replace(tzinfo=tz)
                target_time_utc = dt_local.astimezone(timezone.utc)
            except Exception as exc:
                logger.error("Invalid timezone '%s'. %s", effective_tz_str, exc)
                raise LaunchSetupError() from exc
        else:
            try:
                tz = parse_tz_string(timezone_name)
            except Exception as exc:
                logger.error("Invalid timezone '%s'. %s", timezone_name, exc)
                raise LaunchSetupError() from exc
            target_time_utc = dt_naive.replace(tzinfo=tz).astimezone(timezone.utc)

        now_utc = datetime.now(timezone.utc)
        return target_time_utc - now_utc
    except ValueError as exc:
        logger.error("%s Input was: %r", exc, args_datetime)
        raise LaunchSetupError() from exc
