from datetime import datetime
from zoneinfo import ZoneInfo, available_timezones

# フォールバックしたい場合は pytz を利用
# import pytz

SAMPLE_DATES = [
    datetime(2025, 1, 1, 12, 0, 0),
    datetime(2025, 3, 1, 12, 0, 0),
    datetime(2025, 7, 1, 12, 0, 0),
    datetime(2025, 10, 1, 12, 0, 0),
]


def build_abbrev_map():
    try:
        zones = available_timezones()
    except Exception:
        # zones = set(pytz.all_timezones)  # フォールバック案
        raise

    abbrev_to_zones = {}
    for z in zones:
        tz = ZoneInfo(z)
        for dt in SAMPLE_DATES:
            abbr = dt.replace(tzinfo=tz).strftime("%Z")
            if not abbr:
                continue
            abbrev_to_zones.setdefault(abbr, set()).add(z)

    # 一意に決まるものだけ採用
    unique = {abbr: next(iter(zs)) for abbr, zs in abbrev_to_zones.items() if len(zs) == 1}

    # 例外的に許したいものをここで上書き
    unique.update(
        {
            "JST": "Asia/Tokyo",
            "UTC": "UTC",
            "GMT": "Etc/GMT",
        }
    )
    return unique


TZ_ABBREV_MAP = build_abbrev_map()

print(TZ_ABBREV_MAP)
