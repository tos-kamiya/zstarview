"""Indexed history selection shared by coverage and future image export."""

from __future__ import annotations

from bisect import bisect_right
from datetime import datetime, timedelta

from .models import MINUTE, HistoryRecord, ceil_minute


class HistoryIndex:
    def __init__(self, records: list[HistoryRecord]) -> None:
        self.records = sorted(records, key=lambda r: (r.starts_at, str(r.path)))
        self.starts = [r.starts_at for r in self.records]
        self.max_age = max(
            (r.expires_at - r.starts_at for r in self.records), default=timedelta(0)
        )

    def select(self, target: datetime) -> HistoryRecord | None:
        stop = bisect_right(self.starts, target)
        if stop and self.records[stop - 1].layer in {
            "cloud-native",
            "cloud-geo",
            "satellite-iss",
            "satellite-horizons",
        }:
            latest = self.records[stop - 1]
            if latest.empty and latest.layer.startswith("satellite-"):
                return None
            # Successful acquisitions replace the displayed data. Keep the
            # previous source until that replacement, independent of source
            # age. Only the final record needs a bounded trailing interval.
            return (
                latest if stop < len(self.records) or latest.available(target) else None
            )
        if stop and self.records[stop - 1].layer == "aircraft":
            latest = self.records[stop - 1]
            return latest if latest.available(target) else None
        begin = bisect_right(
            self.starts, target - self.max_age - timedelta(microseconds=1)
        )
        candidates = self.records[begin:stop]
        if not candidates:
            return None
        return max(
            (record for record in candidates if record.available(target)),
            key=lambda r: (r.source_time, r.fetched_at, str(r.path)),
            default=None,
        )


def coverage(
    index: HistoryIndex,
    start: datetime,
    end: datetime,
) -> list[tuple[datetime, datetime]]:
    ranges: list[tuple[datetime, datetime]] = []
    target = ceil_minute(start)
    while target < end:
        if index.select(target) is not None:
            if ranges and ranges[-1][1] == target:
                ranges[-1] = (ranges[-1][0], target + MINUTE)
            else:
                ranges.append((target, target + MINUTE))
        target += MINUTE
    return ranges


def intersect(
    left: list[tuple[datetime, datetime]],
    right: list[tuple[datetime, datetime]],
) -> list[tuple[datetime, datetime]]:
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        begin = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        if begin < end:
            result.append((begin, end))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return result
