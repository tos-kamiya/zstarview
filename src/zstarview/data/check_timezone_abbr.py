"""Build a mapping from timezone abbreviations (e.g., JST, PDT) to IANA zones.

The script samples multiple dates across the year to capture DST/STD variants,
then keeps only abbreviations that uniquely map to exactly one IANA timezone.
Finally, a few well-known exceptions are explicitly whitelisted/overridden.
"""

from datetime import datetime
from zoneinfo import ZoneInfo, available_timezones

# If you want a fallback, you can switch to pytz as a last resort.
# import pytz

# Representative dates to probe both DST and standard time across regions.
SAMPLE_DATES = [
    datetime(2025, 1, 1, 12, 0, 0),
    datetime(2025, 3, 1, 12, 0, 0),
    datetime(2025, 7, 1, 12, 0, 0),
    datetime(2025, 10, 1, 12, 0, 0),
]


def build_abbrev_map() -> dict[str, str]:
    """Construct a dict mapping TZ abbreviations to unique IANA zones.

    Strategy:
        1) Iterate all available IANA timezones.
        2) For each zone, derive its abbreviation on several dates in the year.
        3) Collect a reverse map {abbr -> set(zones)}.
        4) Keep only abbreviations that resolve to a single zone.
        5) Apply explicit overrides for common abbreviations (e.g., JST, UTC, GMT).

    Returns:
        Dict[str, str]: Abbreviation to IANA zone for those that are unambiguous.
    """
    zones = available_timezones()

    abbrev_to_zones: dict[str, set[str]] = {}
    for z in zones:
        tz = ZoneInfo(z)
        for dt in SAMPLE_DATES:
            abbr = dt.replace(tzinfo=tz).strftime("%Z")
            if not abbr:
                continue
            abbrev_to_zones.setdefault(abbr, set()).add(z)

    # Keep only abbreviations that uniquely identify a single zone
    unique: dict[str, str] = {abbr: next(iter(zs)) for abbr, zs in abbrev_to_zones.items() if len(zs) == 1}

    # Explicit overrides for widely used abbreviations (policy choice)
    unique.update(
        {
            "JST": "Asia/Tokyo",
            "UTC": "UTC",
            "GMT": "Etc/GMT",
        }
    )
    return unique


# Build at import time to mirror the original behavior (prints below).
TZ_ABBREV_MAP = build_abbrev_map()

# Keep the original behavior: print the resulting mapping when the script runs.
print(TZ_ABBREV_MAP)
