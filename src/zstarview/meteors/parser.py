"""Parser for GMN trajectory summary text files and directory indexes."""

from __future__ import annotations

import math
import re
from datetime import datetime, timezone
from html.parser import HTMLParser

from .types import GmnDailyFile, MeteorObservation

_DAILY_FILENAME_RE = re.compile(
    r"^traj_summary_(\d{8})_solrange_[0-9.]+-[0-9.]+\.txt$"
)

# Column positions in the documented GMN trajectory summary schema. Sigma
# columns are interleaved with physical values, so positional parsing is less
# ambiguous than the two-row display header.
_UTC_TIME = 2
_SHOWER_CODE = 4
_INITIAL_SPEED = 59
_BEGIN_LAT = 63
_BEGIN_LON = 65
_BEGIN_HEIGHT = 67
_END_LAT = 69
_END_LON = 71
_END_HEIGHT = 73
_DURATION = 75
_PEAK_ABS_MAG = 76
_MIN_COLUMNS = 77


class _DailyIndexParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.filenames: set[str] = set()

    def handle_starttag(
        self,
        tag: str,
        attrs: list[tuple[str, str | None]],
    ) -> None:
        if tag.casefold() != "a":
            return
        href = dict(attrs).get("href")
        if href and _DAILY_FILENAME_RE.fullmatch(href):
            self.filenames.add(href)


def parse_gmn_daily_index(text: str) -> tuple[GmnDailyFile, ...]:
    parser = _DailyIndexParser()
    parser.feed(text)
    files: list[GmnDailyFile] = []
    for filename in parser.filenames:
        match = _DAILY_FILENAME_RE.fullmatch(filename)
        if match is None:
            continue
        nominal_date = datetime.strptime(match.group(1), "%Y%m%d").date()
        files.append(GmnDailyFile(filename=filename, nominal_date=nominal_date))
    return tuple(sorted(files, key=lambda item: (item.nominal_date, item.filename)))


def parse_gmn_trajectory_summary(text: str) -> tuple[MeteorObservation, ...]:
    observations: list[MeteorObservation] = []
    for raw_line in text.splitlines():
        line = raw_line.strip().lstrip("\ufeff")
        if not line or line.startswith("#"):
            continue
        fields = [field.strip() for field in line.split(";")]
        observation = _parse_row(fields)
        if observation is not None:
            observations.append(observation)
    return tuple(observations)


def _parse_row(fields: list[str]) -> MeteorObservation | None:
    if len(fields) < _MIN_COLUMNS:
        return None
    try:
        trajectory_id = fields[0]
        beginning_utc = _parse_utc(fields[_UTC_TIME])
        begin_lat = _required_float(fields[_BEGIN_LAT])
        begin_lon = _required_float(fields[_BEGIN_LON])
        begin_height = _required_float(fields[_BEGIN_HEIGHT])
        end_lat = _required_float(fields[_END_LAT])
        end_lon = _required_float(fields[_END_LON])
        end_height = _required_float(fields[_END_HEIGHT])
    except (ValueError, OverflowError):
        return None
    if not trajectory_id or not _valid_geodetic(begin_lat, begin_lon, begin_height):
        return None
    if not _valid_geodetic(end_lat, end_lon, end_height):
        return None
    shower_code = fields[_SHOWER_CODE] or None
    if shower_code == "...":
        shower_code = None
    return MeteorObservation(
        trajectory_id=trajectory_id,
        beginning_utc=beginning_utc,
        begin_lat_deg=begin_lat,
        begin_lon_deg=begin_lon,
        begin_height_km=begin_height,
        end_lat_deg=end_lat,
        end_lon_deg=end_lon,
        end_height_km=end_height,
        duration_s=_optional_float(fields[_DURATION]),
        peak_abs_magnitude=_optional_float(fields[_PEAK_ABS_MAG]),
        initial_speed_km_s=_optional_float(fields[_INITIAL_SPEED]),
        shower_code=shower_code,
    )


def _parse_utc(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _required_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError("non-finite GMN value")
    return parsed


def _optional_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _valid_geodetic(lat_deg: float, lon_deg: float, height_km: float) -> bool:
    return -90.0 <= lat_deg <= 90.0 and -180.0 <= lon_deg <= 180.0 and height_km >= 0.0
