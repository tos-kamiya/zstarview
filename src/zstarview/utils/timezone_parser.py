"""
Time zone string parsing utility.

This module provides a robust function for converting various time zone string
formats into standard `tzinfo` objects from Python's `datetime` or `zoneinfo`
libraries.
"""

import re
from datetime import timedelta, timezone
from zoneinfo import ZoneInfo

# A mapping of common, unambiguous time zone abbreviations to their corresponding
# IANA database names. This provides a convenient shortcut for users.
# DST (Daylight Saving Time) is handled automatically by the IANA zone.
TZ_ABBREV_MAP: dict[str, str] = {
    # Core
    "UTC": "UTC",
    "GMT": "Etc/GMT",
    "JST": "Asia/Tokyo",
    # East Asia
    "KST": "Asia/Seoul",
    "HKT": "Asia/Hong_Kong",
    # Pacific / Oceania
    "HST": "Pacific/Honolulu",  # UTC-10, no DST
    "AWST": "Australia/Perth",  # UTC+8, no DST
    "ACST": "Australia/Adelaide",  # UTC+9:30, auto-switches to ACDT
    "AEST": "Australia/Sydney",  # UTC+10, auto-switches to AEDT
    "NZST": "Pacific/Auckland",  # Standard time
    "NZDT": "Pacific/Auckland",  # Daylight time handled by same zone
    # Europe / Africa
    "MSK": "Europe/Moscow",  # UTC+3, no DST
    "EAT": "Africa/Nairobi",  # UTC+3, no DST
}

# A regular expression to parse UTC offset formats like "UTC+9", "UTC-07", "UTC+09:30".
UTC_OFFSET_RE = re.compile(
    r"^UTC(?P<sign>[+-])(?P<h>\d{1,2})(?::?(?P<m>\d{2}))?$",
    re.IGNORECASE,
)


def parse_tz_string(tz_str: str) -> timezone | ZoneInfo:
    """
    Parses a time zone string and returns a corresponding tzinfo object.

    This function attempts to parse the string in the following order:
    1. A whitelisted, unambiguous abbreviation (e.g., "JST", "PST").
    2. A UTC offset format (e.g., "UTC+9", "UTC-07:00").
    3. A full IANA time zone identifier (e.g., "Asia/Tokyo", "America/New_York").

    Args:
        tz_str: The time zone string to parse.

    Returns:
        A `datetime.timezone` object for fixed offsets or a `zoneinfo.ZoneInfo`
        object for IANA zones.

    Raises:
        ValueError: If the time zone string is unknown or in an invalid format.
    """
    s = tz_str.strip()
    s_up = s.upper()

    # 1. Check against the map of whitelisted abbreviations.
    if s_up in TZ_ABBREV_MAP:
        return ZoneInfo(TZ_ABBREV_MAP[s_up])

    # 2. Try to match a UTC offset format (e.g., "UTC+9:30").
    m = UTC_OFFSET_RE.match(s_up)
    if m:
        sign = -1 if m.group("sign") == "-" else 1
        h = int(m.group("h"))
        m_str = m.group("m")
        mm = int(m_str) if m_str is not None else 0

        # Validate the hour and minute ranges.
        if not (0 <= h <= 14):  # Practical bounds of the IANA database
            raise ValueError(f"UTC offset hour out of range: {h}")
        if not (0 <= mm < 60):
            raise ValueError(f"UTC offset minute out of range: {mm}")

        delta = timedelta(hours=sign * h, minutes=sign * mm)
        return timezone(delta)

    # 3. As a fallback, try to interpret the string as a direct IANA ID.
    try:
        return ZoneInfo(s)
    except Exception:
        raise ValueError(
            f"Unknown time zone '{tz_str}'. "
            f"Use an IANA ID (e.g., 'Asia/Tokyo'), a supported abbreviation "
            f"({', '.join(sorted(TZ_ABBREV_MAP.keys()))}), or a UTC offset "
            f"(e.g., 'UTC+9', 'UTC+09:30')."
        )
