from datetime import timezone, timedelta
from zoneinfo import ZoneInfo
import re

# Mapping of safe, unambiguous abbreviations to IANA zone names.
# DST handling (if any) is delegated to the IANA zone itself.
TZ_ABBREV_MAP: dict[str, str] = {
    # Core
    "UTC": "UTC",
    "GMT": "Etc/GMT",
    "JST": "Asia/Tokyo",

    # East Asia
    "KST": "Asia/Seoul",
    "HKT": "Asia/Hong_Kong",

    # Pacific / Oceania
    "HST": "Pacific/Honolulu",     # UTC-10, no DST
    "AWST": "Australia/Perth",     # UTC+8, no DST
    "ACST": "Australia/Adelaide",  # UTC+9:30, auto-switches to ACDT
    "AEST": "Australia/Sydney",    # UTC+10, auto-switches to AEDT
    "NZST": "Pacific/Auckland",    # Standard time
    "NZDT": "Pacific/Auckland",    # Daylight time handled by same zone

    # Europe / Africa
    "MSK": "Europe/Moscow",        # UTC+3, no DST
    "EAT": "Africa/Nairobi",       # UTC+3, no DST
}

# Regex to support UTC±offset formats like "UTC+9", "UTC-07", "UTC+09:30"
UTC_OFFSET_RE = re.compile(
    r"^UTC(?P<sign>[+-])(?P<h>\d{1,2})(?::?(?P<m>\d{2}))?$",
    re.IGNORECASE,
)

def parse_tz_string(tz_str: str):
    """
    Parse a time zone string into tzinfo.

    Supports:
    - Whitelisted abbreviations (see TZ_ABBREV_MAP)
    - UTC±offset (e.g., UTC+9, UTC-07, UTC+09:30)
    - Full IANA IDs (e.g., Asia/Tokyo)
    """
    s = tz_str.strip()
    s_up = s.upper()

    # 1) Whitelisted abbreviations
    if s_up in TZ_ABBREV_MAP:
        return ZoneInfo(TZ_ABBREV_MAP[s_up])

    # 2) UTC±offset
    m = UTC_OFFSET_RE.match(s_up)
    if m:
        sign = -1 if m.group("sign") == "-" else 1
        h = int(m.group("h"))
        m_str = m.group("m")
        mm = int(m_str) if m_str is not None else 0

        if not (0 <= h <= 14):  # Practical bound of IANA DB
            raise ValueError(f"UTC offset hour out of range: {h}")
        if not (0 <= mm < 60):
            raise ValueError(f"UTC offset minute out of range: {mm}")

        delta = timedelta(hours=sign * h, minutes=sign * mm)
        return timezone(delta)

    # 3) Fallback: try as IANA time zone ID
    try:
        return ZoneInfo(s)
    except Exception:
        raise ValueError(
            f"Unknown time zone '{tz_str}'. "
            f"Use IANA ID (e.g., 'Asia/Tokyo'), a supported abbreviation "
            f"({', '.join(sorted(TZ_ABBREV_MAP.keys()))}), or UTC±offset "
            f"(e.g., 'UTC+9', 'UTC+09:30')."
        )
